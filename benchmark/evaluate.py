from argparse import ArgumentParser
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from functools import reduce, total_ordering
import gc
from itertools import permutations
import logging
import os
from pathlib import Path
from textwrap import dedent
from typing import Any, Literal
import numpy as np
import numpy.typing as npt
import polars as pl
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# -------------------------------------------------------------------------------------------------
# Basic static variables and constants.
# -------------------------------------------------------------------------------------------------
LABELS: list[str] = ['D5', 'D6', 'F7', 'M8', 'T1', 'T2', 'T3', 'T4', 'HF', 'BC', 'BL']
REPS: list[str] = ['1', '2']
SAMPLES: list[str] = [f'{label}_{rep}' for label in LABELS for rep in REPS]


# -------------------------------------------------------------------------------------------------
# Directories contining annotations.
# -------------------------------------------------------------------------------------------------
TRUSET_CPG_DIR = Path('/mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/c2cpg')
TRUSET_C_DIR = Path('/mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/calibrated')
TRUSET_DMC_DIR = Path('/mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/dmc')
UNION_CPG_DIR = Path('/mnt/eqa/zhangyuanfeng/methylation/retrieve/UNION')
UNION_C_DIR = Path('/mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE/UNION')


# -------------------------------------------------------------------------------------------------
# The following dataclass is used to store a sample dataset.
# -------------------------------------------------------------------------------------------------
@total_ordering
@dataclass(frozen=True)
class SampleData:
    lab: str
    label: Literal['D5', 'D6', 'F7', 'M8', 'T1', 'T2', 'T3', 'T4', 'HF', 'BC', 'BL']
    rep: Literal['1', '2']
    f_path: Path
    sample: str = field(init=False)
    name: str = field(init=False)

    def __post_init__(self):
        object.__setattr__(self, 'name', '_'.join([self.lab, self.label, self.rep]))
        object.__setattr__(self, 'sample', '_'.join([self.label, self.rep]))

    def __lt__(self, other: Any) -> bool:
        if not isinstance(other, SampleData):
            return NotImplemented
        return self.name < other.name

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, SampleData):
            return NotImplemented
        return self.name == other.name


# -------------------------------------------------------------------------------------------------
# The following classes are used to create a lab profile and validate the input paths.
# -------------------------------------------------------------------------------------------------
@total_ordering
class LabDataSet:
    def __init__(self, path_list: list[Path]):
        self.lab: str | None = None
        self.samples: dict[str, SampleData] = {}
        self._parse_paths(path_list)

    def _parse_paths(self, paths: list[Path]):
        if not any(p.exists() for p in paths):
            raise FileNotFoundError('None of input paths exist.')

        labs: set[str] = {p.name.split('.')[0].split('_')[0] for p in paths}

        if len(labs) > 1:
            raise ValueError('only one lab is allowed in the input paths.')

        sample_dict: dict[str, Path] = {
            '_'.join(p.name.split('.')[0].split('_')[1:3]): p
            for p in paths}

        invalid_samples: set[str] = set(sample_dict.keys()) - set(SAMPLES)

        if invalid_samples:
            raise ValueError('Invalid sample identifiers found: '
                             f'{", ".join(sorted(invalid_samples))}')

        self.lab = list(labs)[0]

        for name, path in sample_dict.items():
            if path.exists():
                label, rep = name.split('_')
                self.samples[name] = SampleData(lab=self.lab,
                                                label=label,  # type: ignore
                                                rep=rep,  # type: ignore
                                                f_path=path)

    def sample(self, name: str) -> SampleData:
        if name in self.samples:
            return self.samples[name]
        else:
            raise AttributeError(f'No sample "{name}" in {self.__class__.__name__}')

    def __getattr__(self, name: str) -> SampleData:
        if name in self.samples:
            return self.samples[name]
        else:
            raise AttributeError(f'No sample "{name}" in {self.__class__.__name__}')

    def __iter__(self) -> Iterator[SampleData]:
        yield from self.samples.values()

    def __repr__(self) -> str:
        valid = sum(1 if s.f_path.exists() else 0 for s in self.samples.values())
        return f"<LabDataSet lab={self.lab} total={len(self.samples)} valid={valid}>"

    def __lt__(self, other: Any) -> bool:
        if not isinstance(other, LabDataSet):
            return NotImplemented
        # levae None at the end of queue
        if self.lab is None:
            return False
        if other.lab is None:
            return True
        return self.lab < other.lab

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, LabDataSet):
            return NotImplemented
        return self.lab == other.lab


@total_ordering
@dataclass
class LabProfile:
    Name: str
    Level: Literal['c', 'cpg']
    DataSet: LabDataSet

    def __lt__(self, other: Any) -> bool:
        if not isinstance(other, LabProfile):
            return NotImplemented
        return self.Name < other.Name

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, LabProfile):
            return NotImplemented
        return self.Name == other.Name


class LabProfileFactory:
    @staticmethod
    def create_profile(name: str,
                       level: Literal['c', 'cpg'],
                       paths: set[Path]) -> LabProfile:
        dataset = LabDataSet(paths)
        return LabProfile(name, level, dataset)


# -------------------------------------------------------------------------------------------------
# The following classes are used to read and format the truset.
# -------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class CutParams:
    before_cut: str
    after_cut: str
    breaks: Sequence[int | float]
    labels: Sequence[str]


@dataclass(frozen=True)
class QCutParams:
    before_cut: str
    after_cut: str
    quantiles: Sequence[float]
    labels: Sequence[str]


@dataclass(frozen=True)
class StrSliceParams:
    before_slice: str
    after_slice: str
    offset: int
    length: int


def parse_comma_separated(s: str | None) -> list[str] | None:
    if s is None or s.strip() == '':
        return None
    return [item.strip() for item in s.split(',')]

def get_feature_cols(choice: str, pyro: bool = False) -> list[str]:
    if pyro:
        beta_col = 'beta_pyro'
    else:
        beta_col = 'beta'
    feature_cols: dict[str, list[str]] = {
        'seq_cols': ['seq_2'],
        'complexity_cols': ['gc_bin', 'gc_skew_bin', 'cpg_prop_bin', 'se_bin', 'bwt_bin'],
        'genomic_cols': [
            'promoter', 'enhancer', 'lncrna', 'cpg', 'location', 'BadPromoters', 'AllHomopolymers',
            'AllTandemRepeats', 'Satellites', 'DinucleotideTandemRepeats', 'diTR_10to49', 'diTR_50to149',
            'diTR_ge150', 'triTR_14to49', 'triTR_50to149', 'triTR_ge150', 'quadTR_19to49', 'quadTR_50to149',
            'quadTR_ge150', 'Homopolymer_4to6_AT', 'Homopolymer_4to6_GC', 'Homopolymer_4to6',
            'Homopolymer_7to11_AT', 'Homopolymer_7to11_GC', 'Homopolymer_7to11', 'Homopolymer_ge12_AT',
            'Homopolymer_ge12_GC', 'Homopolymer_ge12', 'Homopolymer_ge21_AT', 'Homopolymer_ge21_GC',
            'Homopolymer_ge21', 'LowMappability', 'Not_LowMappability', 'LowMappability_150bp',
            'Not_LowMappability_150bp', 'Contigs_lt500', 'Gaps_15kb', 'MHC', 'VDJ', 'Human_LINE_1']
    }

    feature_pairs: dict[str, list[str]] = (
        {
            f'{fixed_col}:{col}': [fixed_col, col]
            for fixed_col in [f'{beta_col}_bin', 'depth_bin']
            for col in (
                feature_cols['complexity_cols'] + ['cpg'])
        }
    )
    feature_pairs[f'{beta_col}_bin:depth_bin'] = [f'{beta_col}_bin', 'depth_bin']

    feature_triples: dict[str, list[str]] = (
        {
            f'{beta_col}_bin:depth_bin:{col}': [f'{beta_col}_bin', 'depth_bin', col]
            for col in (feature_cols['complexity_cols'] + ['cpg'])
        }
    )

    feature_cols['all'] = (
        ['stratum', f'{beta_col}_bin', 'depth_bin'] +  # noqa:W504
        feature_cols['complexity_cols'] + feature_cols['genomic_cols'] +  # noqa:W504
        list(feature_pairs.keys()) + list(feature_triples.keys()))

    feature_cols['two_sample'] = ['beta1_bin', 'beta2_bin', 'depth1_bin', 'depth2_bin',
                                  'beta1_bin:beta2_bin', 'depth1_bin:depth2_bin',
                                  'beta1_bin:depth1_bin', 'beta2_bin:depth2_bin', 'global']

    feature_cols['two_sample_features'] = [
        f'depth1_bin:depth2_bin:{col}'
        for col in ('seq_2', 'seq_3', 'gc_bin', 'gc_skew_bin',
                    'cpg_prop_bin', 'se_bin', 'bwt_bin',
                    'b1', 'b2', 'b3', 'a2', 'a3',
                    'promoter', 'enhancer', 'cpg', 'location')
    ]

    feature_cols['compact_union'] = [
        'compact_beta1_bin', 'compact_beta2_bin', 'compact_depth1_bin', 'compact_depth2_bin',
        'compact_beta1_bin:compact_beta2_bin', 'compact_depth1_bin:compact_depth2_bin',
        'compact_beta1_bin:compact_depth1_bin', 'compact_beta2_bin:compact_depth2_bin', 'global']

    feature_cols['c_curated'] = [
        'beta_pyro_bin', 'curated_depth_bin', 'beta_pyro_bin:curated_depth_bin',
        'beta_pyro_bin:curated_depth_bin:cpg', 'beta_pyro_bin:curated_depth_bin:promoter',
        'beta_pyro_bin:curated_depth_bin:enhancer', 'beta_pyro_bin:curated_depth_bin:location',
        'beta_pyro_bin:curated_depth_bin:seq_3', 'beta_pyro_bin:curated_depth_bin:b1',
        'beta_pyro_bin:curated_depth_bin:b2', 'beta_pyro_bin:curated_depth_bin:a2',
        'beta_pyro_bin:curated_depth_bin:b3', 'beta_pyro_bin:curated_depth_bin:a3',
        'beta_pyro_bin:curated_depth_bin:se_bin', 'beta_pyro_bin:curated_depth_bin:bwt_bin',
        'beta_pyro_bin:curated_depth_bin:gc_bin', 'global']

    feature_cols['compact_curated'] = [
        'global', 'compact_depth_bin', 'compact_beta_bin', 'compact_depth_bin:compact_beta_bin',
        'compact_beta_bin:compact_depth_bin:cpg', 'compact_beta_bin:compact_depth_bin:promoter',
        'compact_beta_bin:compact_depth_bin:enhancer', 'compact_beta_bin:compact_depth_bin:location',
        'compact_beta_bin:compact_depth_bin:seq_3', 'compact_beta_bin:compact_depth_bin:b1',
        'compact_beta_bin:compact_depth_bin:b2', 'compact_beta_bin:compact_depth_bin:a2',
        'compact_beta_bin:compact_depth_bin:b3', 'compact_beta_bin:compact_depth_bin:a3',
        'compact_beta_bin:compact_depth_bin:b4', 'compact_beta_bin:compact_depth_bin:a4',
        'compact_beta_bin:compact_depth_bin:b5', 'compact_beta_bin:compact_depth_bin:a5',
        'compact_beta_bin:compact_depth_bin:se_bin', 'compact_beta_bin:compact_depth_bin:bwt_bin',
        'compact_beta_bin:compact_depth_bin:gc_bin']

    match choice:
        # case 'combined_features_ref':
        #     return {
        #         k: v for k, v in (
        #             {**feature_pairs,
        #              **feature_triples}
        #         ).items() if 'depth' not in k
        #     }
        # case 'combined_features_lab':
        #     return {
        #         k: v for k, v in (
        #             {**feature_pairs,
        #              **feature_triples}
        #         ).items() if 'depth' in k
        #     }
        case 'c_ref':
            return [i for i in feature_cols['all'] if 'depth' not in i]
        # case 'cpg_ref':
        #     return [i for i in feature_cols['all']
        #             if not any(b in i for b in ('depth', 'stratum'))]

        case 'c_lab':
            return [i for i in feature_cols['all'] if 'strand' not in i] + ['global']
        case 'cpg_ref':
            return [i for i in feature_cols['all']
                    if not any(b in i for b in ('strand', 'stratum'))] + ['global']

        case 'union':
            # return [*feature_cols['two_sample'],
            #         *feature_cols['seq_cols'],
            #         *feature_cols['complexity_cols'],
            #         *feature_cols['genomic_cols']]
            return feature_cols['two_sample_features'] + feature_cols['two_sample']
        case 'union_without_ref':
            return feature_cols['two_sample']

        case 'compact_union':
            return feature_cols['compact_union']

        case 'dmc':
            return feature_cols['genomic_cols']

        case 'c_curated':
            return feature_cols['c_curated']

        case 'compact_curated':
            return feature_cols['compact_curated']

        case _:
            raise ValueError('Invalid choice')


def resolve_feature(df: pl.DataFrame, feature_col: str) -> pl.LazyFrame:
    cut_params: dict[str, CutParams] = {
        'depth_bin': CutParams(before_cut='depth', after_cut='depth_bin',
                               breaks=list(range(2, 100, 2)),
                               labels=[f'{i - 2}–{i}x' for i in (range(2, 102, 2))] + ['≥100x']),
        'curated_depth_bin': CutParams(before_cut='depth', after_cut='curated_depth_bin',
                                       breaks=[2, 3, 4, 5, 10, 20, 30, 40, 50, 60],
                                       labels=['1x', '2x', '3x', '4x', '5–10x', '10–20x', '20–30x',
                                               '30–40x', '40–50x', '50–60x', '>60x']),
        'compact_depth_bin': CutParams(before_cut='depth', after_cut='compact_depth_bin',
                                       breaks=[2, 3, 4, 5, 10],
                                       labels=['1x', '2x', '3x', '4x', '5–10x', '≥10x']),
        'grouped_depth_bin': CutParams(before_cut='depth', after_cut='grouped_depth_bin',
                                       breaks=[5, 10],
                                       labels=['1–4x', '5–10x', '≥10x']),
        'depth1_bin': CutParams(before_cut='depth1', after_cut='depth1_bin',
                                breaks=list(range(2, 102, 2)),
                                labels=[f'{i - 2}–{i}x' for i in (range(2, 102, 2))] + ['≥100x']),
        'depth2_bin': CutParams(before_cut='depth2', after_cut='depth2_bin',
                                breaks=list(range(2, 102, 2)),
                                labels=[f'{i - 2}–{i}x' for i in (range(2, 102, 2))] + ['≥100x']),
        'compact_depth1_bin': CutParams(before_cut='depth1', after_cut='compact_depth1_bin',
                                        breaks=[2, 3, 4, 5, 10],
                                        labels=['1x', '2x', '3x', '4x', '5–10x', '≥10x']),
        'compact_depth2_bin': CutParams(before_cut='depth2', after_cut='compact_depth2_bin',
                                        breaks=[2, 3, 4, 5, 10],
                                        labels=['1x', '2x', '3x', '4x', '5–10x', '≥10x']),
        's1_depth_bin': CutParams(before_cut='s1_depth', after_cut='s1_depth_bin',
                                  breaks=list(range(2, 102, 2)),
                                  labels=[f'{i - 2}–{i}x' for i in (range(2, 102, 2))] + ['≥100x']),
        's2_depth_bin': CutParams(before_cut='s2_depth', after_cut='s2_depth_bin',
                                  breaks=list(range(2, 102, 2)),
                                  labels=[f'{i - 2}–{i}x' for i in (range(2, 102, 2))] + ['≥100x']),
        'beta_bin': CutParams(before_cut='beta', after_cut='beta_bin',
                              breaks=list(range(5, 100, 5)),
                              labels=[f'{i}_{i + 5}' for i in range(0, 100, 5)]),
        'beta_pyro_bin': CutParams(before_cut='beta_pyro', after_cut='beta_pyro_bin',
                                   breaks=list(range(5, 100, 5)),
                                   labels=[f'{i}_{i + 5}' for i in range(0, 100, 5)]),
        'compact_beta_bin': CutParams(before_cut='beta_pyro', after_cut='compact_depth_bin',
                                      breaks=list(range(10, 100, 10)),
                                      labels=[f'{i}–{i + 10}' for i in range(0, 100, 10)]),
        'beta1_bin': CutParams(before_cut='beta1', after_cut='beta1_bin',
                               breaks=list(range(5, 100, 5)),
                               labels=[f'{i}–{i + 5}' for i in range(0, 100, 5)]),
        'beta2_bin': CutParams(before_cut='beta2', after_cut='beta2_bin',
                               breaks=list(range(5, 100, 5)),
                               labels=[f'{i}–{i + 5}' for i in range(0, 100, 5)]),
        'compact_beta1_bin': CutParams(before_cut='beta1', after_cut='compact_beta1_bin',
                                       breaks=list(range(10, 100, 10)),
                                       labels=[f'{i}–{i + 10}' for i in range(0, 100, 10)]),
        'compact_beta2_bin': CutParams(before_cut='beta2', after_cut='compact_beta2_bin',
                                       breaks=list(range(10, 100, 10)),
                                       labels=[f'{i}–{i + 10}' for i in range(0, 100, 10)]),
        # 'beta_bin': CutParams(before_cut='beta', after_cut='beta_bin',
        #                       breaks=list(range(1, 100)),
        #                       labels=[f'{i - 1}_{i}' for i in (range(1, 100))] + ['99_100']),
        # 'beta_pyro_bin': CutParams(before_cut='beta_pyro', after_cut='beta_pyro_bin',
        #                            breaks=list(range(1, 100)),
        #                            labels=[f'{i - 1}_{i}' for i in (range(1, 100))] + ['99_100']),
        # 'beta1_bin': CutParams(before_cut='beta1', after_cut='beta1_bin',
        #                        breaks=list(range(1, 100)),
        #                        labels=[f'{i - 1}_{i}' for i in (range(1, 100))] + ['99_100']),
        # 'beta2_bin': CutParams(before_cut='beta2', after_cut='beta2_bin',
        #                        breaks=list(range(1, 100)),
        #                        labels=[f'{i - 1}_{i}' for i in (range(1, 100))] + ['99_100']),
        'gc_bin': CutParams(before_cut='GC%_70', after_cut='gc_bin',
                            breaks=[i / 100 for i in range(10, 100, 10)],
                            labels=[f'{i}–{i + 10}' for i in range(0, 100, 10)]),
        'gc_skew_bin': CutParams(before_cut='GC_skew_70', after_cut='gc_skew_bin',
                                 breaks=[i / 10 for i in range(-9, 10, 1)],
                                 labels=[f'{i / 10} – {(i + 1) / 10}' for i in range(-10, 10, 1)]),
        'cpg_prop_bin': CutParams(before_cut='CpG_GC_ratio_70', after_cut='cpg_prop_bin',
                                  breaks=[i / 100 for i in range(2, 42, 2)],
                                  labels=[f'{i / 100}–{(i + 2) / 100}' for i in range(0, 42, 2)]),
        'SNV_dist_bin': CutParams(before_cut='nearest_SNV_dist', after_cut='SNV_dist_bin',
                                  breaks=list(range(0, 150, 10)),
                                  labels=['None'] + [f'{i}–{i + 10}' for i in range(0, 150, 10)]),
        'sINDEL_dist_bin': CutParams(before_cut='nearest_sINDEL_dist', after_cut='sINDEL_dist_bin',
                                     breaks=list(range(0, 150, 10)),
                                     labels=['None'] + [f'{i}–{i + 10}' for i in range(0, 150, 10)]),
        '1000bp_SVs_bin': CutParams(before_cut='1000bp_SVs', after_cut='1000bp_SVs_bin',
                                    breaks=[1, 5, 10, 15, 20],
                                    labels=['None', '1–5', '5–10', '10–15', '15–20', '>20']),
        'SV_dist_bin': CutParams(before_cut='nearest_SV_dist', after_cut='SV_dist_bin',
                                 breaks=list(range(0, 1000, 100)),
                                 labels=['None'] + [f'{i}–{i + 100}' for i in range(0, 1000, 100)]),
        'SV_size_bin': CutParams(before_cut='nearest_SV_size', after_cut='SV_size_bin',
                                 breaks=[1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7],
                                 labels=['None', '1–10', '10–100', '100–1e3', '1e3–1e4',
                                         '1e4–1e5', '1e5–1e6', '1e6–1e7', '>1e7'])}

    qcut_params: dict[str, QCutParams] = {
        'se_bin': QCutParams(before_cut='ShannonEntropy_70', after_cut='se_bin',
                             quantiles=[0.25, 0.5, 0.75],
                             labels=['low', 'medium', 'high', 'very_high']),
        'bwt_bin': QCutParams(before_cut='BWT_ratio_70', after_cut='bwt_bin',
                              quantiles=[0.25, 0.5, 0.75],
                              labels=['low', 'medium', 'high', 'very_high']),
    }

    str_slice_params: dict[str, StrSliceParams] = {
        'seq_2': StrSliceParams(before_slice='seq_5', after_slice='seq_2',
                                offset=3, length=5),
        'seq_3': StrSliceParams(before_slice='seq_5', after_slice='seq_3',
                                offset=2, length=7),
        'b5': StrSliceParams(before_slice='seq_5', after_slice='b1',
                             offset=0, length=1),
        'b4': StrSliceParams(before_slice='seq_5', after_slice='b2',
                             offset=1, length=1),
        'b3': StrSliceParams(before_slice='seq_5', after_slice='b3',
                             offset=2, length=1),
        'b2': StrSliceParams(before_slice='seq_5', after_slice='b4',
                             offset=3, length=1),
        'b1': StrSliceParams(before_slice='seq_5', after_slice='b5',
                             offset=4, length=1),
        'a1': StrSliceParams(before_slice='seq_5', after_slice='a1',
                             offset=6, length=1),
        'a2': StrSliceParams(before_slice='seq_5', after_slice='a2',
                             offset=7, length=1),
        'a3': StrSliceParams(before_slice='seq_5', after_slice='a3',
                             offset=8, length=1),
        'a4': StrSliceParams(before_slice='seq_5', after_slice='a4',
                             offset=9, length=1),
        'a5': StrSliceParams(before_slice='seq_5', after_slice='a5',
                             offset=10, length=1)
    }

    grouped_beta_expr: pl.Expr = (pl.when(pl.col('beta_pyro') < 10)
                                    .then(pl.lit('a'))
                                    .when(pl.col('beta_pyro') >= 90)
                                    .then(pl.lit('a'))
                                    .when(pl.col('beta_pyro').is_between(lower_bound=10,
                                                                         upper_bound=20,
                                                                         closed='left'))
                                    .then(pl.lit('b'))
                                    .when(pl.col('beta_pyro').is_between(lower_bound=20,
                                                                         upper_bound=30,
                                                                         closed='left'))
                                    .then(pl.lit('c'))
                                    .when(pl.col('beta_pyro').is_between(lower_bound=80,
                                                                         upper_bound=90,
                                                                         closed='left'))
                                    .then(pl.lit('c'))
                                    .when(pl.col('beta_pyro').is_between(lower_bound=30,
                                                                         upper_bound=40,
                                                                         closed='left'))
                                    .then(pl.lit('d'))
                                    .otherwise(pl.lit('e')))

    existing_cols = ['global', 'stratum', 'seq_5', 'promoter', 'enhancer', 'lncrna', 'cpg',
                     'location', 'BadPromoters', 'AllHomopolymers', 'AllTandemRepeats',
                     'Satellites', 'DinucleotideTandemRepeats', 'diTR_10to49', 'diTR_50to149',
                     'diTR_ge150', 'triTR_14to49', 'triTR_50to149', 'triTR_ge150', 'quadTR_19to49',
                     'quadTR_50to149', 'quadTR_ge150', 'Homopolymer_4to6_AT', 'Homopolymer_4to6_GC',
                     'Homopolymer_4to6', 'Homopolymer_7to11_AT', 'Homopolymer_7to11_GC',
                     'Homopolymer_7to11', 'Homopolymer_ge12_AT', 'Homopolymer_ge12_GC', 'Homopolymer_ge12',
                     'Homopolymer_ge21_AT', 'Homopolymer_ge21_GC', 'Homopolymer_ge21', 'LowMappability',
                     'Not_LowMappability', 'LowMappability_150bp', 'Not_LowMappability_150bp',
                     'Contigs_lt500', 'Gaps_15kb', 'MHC', 'VDJ', 'Human_LINE_1',
                     '150bp_SNVs', 'nearest_SNV_ref', 'nearest_SNV_alt',
                     '150bp_sINDELs', 'nearest_sINDEL_type', 'nearest_sINDEL_size', '1000bp_SVs',
                     'nearest_SV_type', '1000bp_SVs',
                     'D5', 'D6', 'F7', 'M8', 'T1', 'T2', 'T3', 'T4', 'BC', 'BL', 'HF']
    if feature_col in existing_cols:
        return df.lazy()
    elif feature_col in cut_params.keys():
        return (df.lazy()
                  .with_columns(pl.col(cut_params[feature_col].before_cut)
                                  .cut(breaks=cut_params[feature_col].breaks,
                                       labels=cut_params[feature_col].labels,
                                       left_closed=True).alias(feature_col)))
    elif feature_col in qcut_params.keys():
        return (df.lazy()
                  .with_columns(pl.col(qcut_params[feature_col].before_cut)
                                  .qcut(quantiles=qcut_params[feature_col].quantiles,
                                        labels=qcut_params[feature_col].labels,
                                        left_closed=True).alias(feature_col)))
    elif feature_col in str_slice_params.keys():
        return (df.lazy()
                  .with_columns(pl.col('seq_5')
                                  .str
                                  .slice(offset=str_slice_params[feature_col].offset,
                                         length=str_slice_params[feature_col].length)
                                  .alias(feature_col)))
    elif feature_col == 'grouped_beta_bin':
        return (df.lazy()
                  .with_columns(grouped_beta_expr.alias('grouped_beta_bin')))
    else:
        if ':' in feature_col:
            df = df.lazy()  # type: ignore
            f_parts = feature_col.split(':')
            for f_part in f_parts:
                if f_part in existing_cols:
                    pass
                elif f_part in cut_params.keys():
                    df = (df.with_columns(pl.col(cut_params[f_part].before_cut)
                                            .cut(breaks=cut_params[f_part].breaks,
                                                 labels=cut_params[f_part].labels,
                                                 left_closed=True).alias(f_part)))
                elif f_part in qcut_params.keys():
                    df = (df.with_columns(pl.col(qcut_params[f_part].before_cut)
                                            .qcut(quantiles=qcut_params[f_part].quantiles,
                                                  labels=qcut_params[f_part].labels,
                                                  left_closed=True).alias(f_part)))
                elif f_part in str_slice_params.keys():
                    df = (df.with_columns(pl.col('seq_5')
                                            .str
                                            .slice(offset=str_slice_params[f_part].offset,
                                                   length=str_slice_params[f_part].length)
                                            .alias(f_part)))
                elif f_part == 'grouped_beta_bin':
                    df = (df.with_columns(grouped_beta_expr.alias('grouped_beta_bin')))
                else:
                    raise ValueError(f'how to prepare it is not specified: {f_part}')
            df = (df.with_columns(pl.concat_str(f_parts, separator=':').alias(feature_col))
                    .drop(f_parts))
            return df  # type: ignore
        else:
            raise ValueError(f'how to prepare it is not specified: {feature_col}')


class ReferenceResolver:
    def __init__(self,
                 cpg_dir: Path, c_dir: Path, uc_dir: Path, ucpg_dir: Path, dmc_dir: Path,
                 cpg_suffix: str = '.parquet.lz4', c_suffix: str = '.parquet.lz4',
                 uc_suffix: str = '.parquet.lz4', ucpg_suffix: str = '.parquet.zst',
                 dmc_suffix: str = '.parquet.lz4'):
        self.cpg_dir: Path = cpg_dir
        self.c_dir: Path = c_dir
        self.uc_dir: Path = uc_dir
        self.ucpg_dir: Path = ucpg_dir
        self.dmc_dir: Path = dmc_dir

        self.type_mapping: dict[str, tuple[Path, str]] = ({
            'cpg': (cpg_dir, cpg_suffix),
            'c': (c_dir, c_suffix),
            'uc': (uc_dir, uc_suffix),
            'ucpg': (ucpg_dir, ucpg_suffix),
            'dmc': (dmc_dir, dmc_suffix)
        })

        self.mixup_props: dict[str, tuple[float, float]] = ({
            'T1': (0.8, 0.2), 'T2': (0.6, 0.4), 'T3': (0.4, 0.6), 'T4': (0.2, 0.8)
        })

        self.mixup_cols: list[str] = [
            'chrom', 'start', 'type', 'end', 'strand', 'in_hcr', 'seq_2', 'seq_3',
            'GC%_75', 'GC_skew_75', 'CpG_GC_ratio_75', 'ShannonEntropy_75', 'BWT_ratio_75',
            'gene_symbols', 'promoter', 'enhancer', 'lncrna', 'cpg', 'location',
            'BadPromoters', 'AllHomopolymers', 'AllTandemRepeats', 'Satellites',
            'DinucleotideTandemRepeats', 'diTR_10to49', 'diTR_50to149', 'diTR_ge150',
            'triTR_14to49', 'triTR_50to149', 'triTR_ge150', 'quadTR_19to49', 'quadTR_50to149',
            'quadTR_ge150', 'Homopolymer_4to6_AT', 'Homopolymer_4to6_GC', 'Homopolymer_4to6',
            'Homopolymer_7to11_AT', 'Homopolymer_7to11_GC', 'Homopolymer_7to11', 'Homopolymer_ge12_AT',
            'Homopolymer_ge12_GC', 'Homopolymer_ge12', 'Homopolymer_ge21_AT', 'Homopolymer_ge21_GC',
            'Homopolymer_ge21', 'LowMappability', 'Not_LowMappability', 'LowMappability_150bp',
            'Not_LowMappability_150bp', 'Contigs_lt500', 'Gaps_15kb', 'MHC', 'VDJ', 'Human_LINE_1']

        for path in (cpg_dir, c_dir, dmc_dir):
            if not path.exists() or not path.is_dir():
                raise FileNotFoundError(f'No such dir: {path}')

    def get(self,
            ref_type: Literal['c', 'cpg', 'uc', 'ucpg', 'dmc'],
            labels: list[str]) -> pl.LazyFrame:
        ref_dir: Path = self.type_mapping[ref_type][0]
        ref_suffix: str = self.type_mapping[ref_type][1]
        label_str: str = '_'.join(labels)

        if 'parquet' in ref_suffix:
            def load_func(path: Path):
                return pl.scan_parquet(path)
        elif 'bedgraph' in ref_suffix:
            def load_func(path: Path):
                return pl.scan_csv(path, separator='\t', has_header=False)
        elif 'csv' in ref_suffix:
            def load_func(path: Path):
                return pl.scan_csv(path)
        else:
            raise ValueError(f'Unsupported file suffix: {ref_suffix}')

        to_return: pl.LazyFrame
        if ref_type in ('c', 'cpg'):
            if ref_type == 'c':
                filter_condition = pl.col('in_hcr')
            else:
                filter_condition = (pl.col('in_hcr-') & pl.col('in_hcr+'))

            if label_str in ['D5', 'D6', 'F7', 'M8', 'BC', 'BL', 'HF']:
                to_return = (load_func(ref_dir / f'{label_str}{ref_suffix}')
                             .filter(filter_condition))
            elif label_str in ['T1', 'T2', 'T3', 'T4']:
                if (ref_dir / f'mixup{ref_suffix}').exists():
                    _MIXUP_COLS: set[str] = {'D5', 'D5_pyro', 'D6', 'D6_pyro',
                                             'T1', 'T1_pyro', 'T2', 'T2_pyro',
                                             'T3', 'T3_pyro', 'T4', 'T4_pyro'}
                    to_return = (load_func(ref_dir / f'mixup{ref_suffix}')
                                 .filter(filter_condition)
                                 .drop({c for c in _MIXUP_COLS if label_str not in c})
                                 .rename({label_str: 'beta', f'{label_str}_pyro': 'beta_pyro'}))
                else:
                    _a, _b = self.mixup_props[label_str]
                    (_d5_tmp,
                     _d6_tmp) = [(load_func(ref_dir / f'{_l}{ref_suffix}')
                                  .rename({'beta': _l, 'beta_pyro': f'{_l}_pyro'})
                                  .filter(filter_condition)) for _l in ('D5', 'D6')]
                    to_return = (_d5_tmp.join(other=_d6_tmp, on=['chrom', 'start'])
                                 .drop_nulls(subset=['D5', 'D5_pyro', 'D6', 'D6_pyro'])
                                 .with_columns([((pl.col(f'D5{_s}') * _a) + (pl.col(f'D6{_s}') * _b)).alias(f'beta{_s}')
                                                for _s in ('', '_pyro')])
                                 .drop(['D5', 'D6', 'D5_pyro', 'D6_pyro'])) #  .rename({label_str: 'beta', f'{label_str}_pyro': 'beta_pyro'}))
            else:
                raise ValueError(f'how to prepare this label is not specified: {label_str}')
        elif ref_type == 'uc':
            to_return = load_func(ref_dir / f'{label_str}{ref_suffix}')
        else:
            raise ValueError(f'how to prepare this ref_type is not specified: {ref_type}')

        return (to_return.with_columns(pl.when(pl.col('seq_5')
                                                 .str
                                                 .slice(offset=5, length=1) == 'G')
                                         .then(pl.col('seq_5')
                                                 .str
                                                 .reverse()
                                                 .replace_strict(old=['A', 'C', 'G', 'T'],
                                                                 new=['T', 'G', 'C', 'A'],
                                                                 default='N',
                                                                 return_dtype=pl.String))
                                         .otherwise(pl.col('seq_5'))))

    @staticmethod
    def help():
        print(dedent("""\
        c ref bedgraph contains columns:
        (1) basic columns:
            chrom, start, end, strand, type, in_hcr
        (2) beta value columns:
            beta, beta_pyro, 95ci, sd, mad, cv, credibility
        (3) sequence pattern columns:
            seq_2, seq_3, seq_5, GC%_75, GC_skew_75, CpG_GC_ratio_75, ShannonEntropy_75, BWT_ratio_75
        (4) genomic feature columns:
            gene_symbols, promoter, enhancer, lncrna, cpg, location
        (5) GIAB stratification columns:
            BadPromoters, AllHomopolymers, AllTandemRepeats, Satellites, DinucleotideTandemRepeats,
            diTR_10to49, diTR_50to149, diTR_ge150, triTR_14to49, triTR_50to149, triTR_ge150, quadTR_19to49,
            quadTR_50to149, quadTR_ge150, Homopolymer_4to6_AT, Homopolymer_4to6_GC, Homopolymer_4to6,
            Homopolymer_7to11_AT, Homopolymer_7to11_GC, Homopolymer_7to11, Homopolymer_ge12_AT,
            Homopolymer_ge12_GC, Homopolym er_ge12, Homopolymer_ge21_AT, Homopolymer_ge21_GCHomopolymer_ge21,
            LowMappability, Not_LowMappability, LowMappability_150bp, Not_LowMappability_150bp, Contigs_lt500,
            Gaps_15kb, MHC, VDJ, Human_LINE_1

        cpg ref bedgraph contains columns:
        (1) basic columns:
            chrom, start, end, strand, type, in_hcr+, in_hcr-
        (2) beta value columns:
            beta, beta+, beta-, beta_pyro, beta_pyro+, beta_pyro-, 95ci+, 95ci-,
            sd+, sd-, mad+, mad-, cv+, cv-, credibility+, credibility-
        (3) sequence pattern columns:
            seq_2, seq_3, seq_5, GC%_75, GC_skew_75, CpG_GC_ratio_75, ShannonEntropy_75, BWT_ratio_75
        (4) genomic feature columns:
            gene_symbols, promoter, enhancer, lncrna, cpg, location
        (5) GIAB stratification columns:
            BadPromoters, AllHomopolymers, AllTandemRepeats, Satellites, DinucleotideTandemRepeats,
            diTR_10to49, diTR_50to149, diTR_ge150, triTR_14to49, triTR_50to149, triTR_ge150, quadTR_19to49,
            quadTR_50to149, quadTR_ge150, Homopolymer_4to6_AT, Homopolymer_4to6_GC, Homopolymer_4to6,
            Homopolymer_7to11_AT, Homopolymer_7to11_GC, Homopolymer_7to11, Homopolymer_ge12_AT,
            Homopolymer_ge12_GC, Homopolymer_ge12, Homopolymer_ge21_AT, Homopolymer_ge21_GCHomopolymer_ge21,
            LowMappability, Not_LowMappability, LowMappability_150bp, Not_LowMappability_150bp, Contigs_lt500,
            Gaps_15kb, MHC, VDJ, Human_LINE_1

        union cpg and union c bedgraph contains columns:
        (1) basic columns:
            chrom, start, end, strand
        (2) sequence pattern columns:
            seq_2, seq_3, seq_5, GC%_75, GC_skew_75, CpG_GC_ratio_75, ShannonEntropy_75, BWT_ratio_75
        (3) genomic feature columns:
            gene_symbols, promoter, enhancer, lncrna, cpg, location
        (4) GIAB stratification columns:
            BadPromoters, AllHomopolymers, AllTandemRepeats, Satellites, DinucleotideTandemRepeats,
            diTR_10to49, diTR_50to149, diTR_ge150, triTR_14to49, triTR_50to149, triTR_ge150, quadTR_19to49,
            quadTR_50to149, quadTR_ge150, Homopolymer_4to6_AT, Homopolymer_4to6_GC, Homopolymer_4to6,
            Homopolymer_7to11_AT, Homopolymer_7to11_GC, Homopolymer_7to11, Homopolymer_ge12_AT,
            Homopolymer_ge12_GC, Homopolymer_ge12, Homopolymer_ge21_AT, Homopolymer_ge21_GCHomopolymer_ge21,
            LowMappability, Not_LowMappability, LowMappability_150bp, Not_LowMappability_150bp, Contigs_lt500,
            Gaps_15kb, MHC, VDJ, Human_LINE_1

        dmc ref bedgraph contains columns:
        (1) basic columns:
            chrom, start, end, strand, type
        (2) beta value columns:
            beta_D5, beta_D6, beta_F7, beta_M8, beta_T1, beta_T2, beta_T3, beta_T4, beta_BC, beta_BL
        (3) dmc columns:
            D5_D6, D5_F7, D5_M8, D6_F7, D6_M8, F7_M8, BC_BL,
            D5_T1, D5_T2, D5_T3, D5_T4, D6_T1, D6_T2, D6_T3, D6_T4,
            T1_T3, T1_T4, T2_T3, T2_T4, T3_T4
        (4) genomic feature columns:
            gene_symbols, promoter, enhancer, lncrna, cpg, location
        (5) GIAB stratification columns:
            BadPromoters, AllHomopolymers, AllTandemRepeats, Satellites, DinucleotideTandemRepeats,
            diTR_10to49, diTR_50to149, diTR_ge150, triTR_14to49, triTR_50to149, triTR_ge150, quadTR_19to49,
            quadTR_50to149, quadTR_ge150, Homopolymer_4to6_AT, Homopolymer_4to6_GC, Homopolymer_4to6,
            Homopolymer_7to11_AT, Homopolymer_7to11_GC, Homopolymer_7to11, Homopolymer_ge12_AT,
            Homopolymer_ge12_GC, Homopolym er_ge12, Homopolymer_ge21_AT, Homopolymer_ge21_GCHomopolymer_ge21,
            LowMappability, Not_LowMappability, LowMappability_150bp, Not_LowMappability_150bp, Contigs_lt500,
            Gaps_15kb, MHC, VDJ, Human_LINE_1
        """))


# -------------------------------------------------------------------------------------------------
# The following classes are used to calculate the methylation metrics.

# ! To make it fast, naive polars api is used
# ! instead of scikit-learn.metrics.root_mean_squared_error and scipy.stat.spearmanr.
# -------------------------------------------------------------------------------------------------
class RootMeanSquareError:
    def __init__(self, sample: SampleData, ref: pl.DataFrame, output_dir: Path,
                 level: Literal['cpg', 'c'], depth: int, lab_beta: str = 'beta',
                 pyro: bool = False, feature: bool = True, curated: bool = True,
                 compact: bool = True, tmp_features: list[str] | None = None,
                 force: bool = False):
        self.sample: SampleData = sample
        self.ref: pl.DataFrame = ref
        self.out_f: Path = Path(f'{output_dir.as_posix()}/{self.sample.name}.csv')
        self.level: Literal['cpg', 'c'] = level
        self.lab_beta: str = lab_beta
        self.pyro: bool = pyro
        self.feature: bool = feature
        self.curated: bool = curated
        self.compact: bool = compact
        self.tmp_features: list[str] | None = tmp_features
        self.depth: int = depth
        self.force: bool = force

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        logging.debug('🛠️ RMSE object initialized')
        logging.debug(str(self))

    def __repr__(self) -> str:
        return (f'RootMeanSquareError('
                f'sample={self.sample!r}, '
                f'out_f={self.out_f!r}, '
                f'pyro={self.pyro!r}, '
                f'feature={self.feature!r}, '
                f'depth={self.depth})')

    def __str__(self) -> str:
        return (f'RootMeanSquareError\n'
                f'-------------------\n'
                f'Sample Name: {self.sample.name}\n'
                f'Output Directory: {self.out_f}\n'
                f'Beta Col to Evalute: {self.lab_beta}\n'
                f'Pyro Mode: {'Enabled' if self.pyro else 'Disabled'}\n'
                f'Feature Mode: {'Enabled' if self.feature or self.tmp_features else 'Disabled'}\n'
                f'Temporary Features: {'Enabled' if self.tmp_features else 'Disabled'}\n'
                f'Depth Filter: depth >= {self.depth}')

    def _read_parquet(self) -> pl.LazyFrame:
        to_select: list[str]
        lzdf: pl.LazyFrame
        match (self.level, self.sample.lab[: 2]):
            case ('cpg', 'MA'):
                to_select = ['chrom', 'start', 'beta_pred']
                lzdf = pl.scan_parquet(self.sample.f_path)
            case ('c', 'MA'):
                to_select = ['chrom', 'start', 'beta_pred']
                lzdf = (pl.concat([pl.scan_parquet(self.sample.f_path)
                                     .with_columns((pl.col('end') - 1)),
                                   pl.scan_parquet(self.sample.f_path)
                                     .with_columns((pl.col('start') + 1))]))
            case _:
                to_select = ['chrom', 'start', 'beta_pred', 'depth']
                lzdf = (pl.scan_parquet(self.sample.f_path)
                          .filter(pl.col('depth') >= self.depth))

        return (lzdf.with_columns(pl.col('beta')
                                    .cast(pl.Float64)
                                    .alias('beta_pred'))
                    .select(to_select))

    @staticmethod
    def _process_feature(df: pl.DataFrame, feature: str, ref_beta: str) -> pl.LazyFrame:
        return (resolve_feature(df=df, feature_col=feature)
                .group_by(feature)
                .agg(pl.len().alias('count'),
                     (pl.col(ref_beta) - pl.col('beta_pred')).pow(2)
                     .mean()
                     .sqrt()
                     .alias('rmse'))
                .rename({feature: 'fgroup'})
                .with_columns(pl.lit(feature).alias('feature')))

    def _compute(self):
        # to_select = ['chrom', 'start', 'beta_pred', 'depth']
        # if self.sample.lab[:2] == 'MA':
        #     to_select = to_select[:-1]
        #     if self.level == 'cpg':
        #         lab_tmp: pl.LazyFrame = (pl.scan_parquet(self.sample.f_path)
        #                                    .rename({self.lab_beta: 'beta_pred'})
        #                                    .with_columns(pl.col('beta_pred').cast(pl.Float64)))

        # lab_tmp: pl.LazyFrame = (pl.scan_parquet(self.sample.f_path)
        #                            .rename({self.lab_beta: 'beta_pred'})
        #                            .with_columns(pl.col('beta_pred').cast(pl.Float64)))
        # if self.depth > 1 and self.sample.lab[: 2] != 'MA':
        #     lab_tmp = lab_tmp.filter(pl.col('depth') >= self.depth)
        lab_tmp: pl.LazyFrame = self._read_parquet()

        to_compute = (self.ref
                          .lazy()
                          .join(other=lab_tmp,
                                on=['chrom', 'start'], how='inner')
                          .with_columns(pl.lit('global').alias('global'))
                          .collect())

        if to_compute.is_empty():
            logging.warning(f'⚠️ No data to compute RMSE for {self.sample}')
            (pl.DataFrame({'lab': [self.sample.lab], 'sample': [self.sample.sample],
                           'feature': ['global'], 'fgroup': ['global'], 'count': [0], 'rmse': [-1]})
               .write_csv(self.out_f))
            return
        else:
            logging.debug('🛠️ sample and ref dataframe joined, '
                          f'containing {to_compute.shape[0]} rows')

        feature_cols: set[str]
        if self.tmp_features:
            feature_cols = set(self.tmp_features)
        elif self.compact:
            feature_cols = set(get_feature_cols(choice='compact_curated', pyro=self.pyro))
        else:
            if self.feature:
                if self.level == 'c':
                    if self.curated:
                        feature_cols = set(get_feature_cols(choice='c_curated', pyro=self.pyro))
                    else:
                        feature_cols = {*get_feature_cols(choice='c_lab', pyro=self.pyro),
                                        *get_feature_cols(choice='c_ref', pyro=self.pyro)}  # type: ignore
                    # if self.sample.label == 'BC':
                    #     _variant_cols: set[str] = ({
                    #         '150bp_SNVs', 'SNV_dist_bin', 'nearest_SNV_ref', 'nearest_SNV_alt',
                    #         'nearest_SNV_ref:nearest_SNV_alt', '150bp_sINDELs', 'nearest_sINDEL_type',
                    #         'sINDEL_dist_bin', 'nearest_sINDEL_size', '1000bp_SVs_bin', 'nearest_SV_type',
                    #         'SV_dist_bin', 'SV_size_bin'})
                    #     feature_cols.update(_variant_cols)
                    #     feature_cols.update({f'beta_bin:depth_bin:{_v_col}' for _v_col in _variant_cols})
                else:
                    feature_cols = {*get_feature_cols(choice='cpg_lab', pyro=self.pyro),
                                    *get_feature_cols(choice='cpg_ref', pyro=self.pyro)}  # type: ignore
            else:
                feature_cols = {'global'}

        if self.sample.lab[:2] == 'MA':
            feature_cols = {i for i in feature_cols if 'depth' not in i}

        feature_dfs: list[pl.DataFrame] = []
        logging.debug(f'🛠️ collecting {len(feature_cols)} feature results')
        for feature in feature_cols:
            logging.debug(f'🛠️ processing feature: {feature}')
            feature_dfs.append(self._process_feature(df=to_compute,
                                                     feature=feature,
                                                     ref_beta='beta_pyro' if self.pyro else 'beta')
                                   .with_columns(pl.col('fgroup').cast(pl.String))
                                   .collect())
            logging.debug(f'🛠️ {feature} processed')
        logging.debug('🛠️ feature dfs collected')

        (pl.concat(feature_dfs)
           .with_columns([pl.lit(self.sample.lab).alias('lab'),
                          pl.lit(self.sample.sample).alias('sample')])
           .select(['lab', 'sample', 'feature', 'fgroup', 'count', 'rmse'])
           .write_csv(self.out_f))
        logging.debug(f'🛠️ dataframes concatenated and saved at {self.out_f}')

    def run(self):
        if self.out_f.exists() and not self.force:
            logging.info(f'ℹ️ RMSE already calculated for {self.sample.name}')
            return
        logging.info(f'ℹ️ Calculating RMSE for {self.sample.name}')
        self._compute()
        logging.info(f'✅ RMSE calculated for {self.sample.name}')


class SpearmanR:
    def __init__(self, rep1: SampleData, rep2: SampleData, union_ref: pl.DataFrame | None,
                 output_dir: Path, depth: int, beta_col: str, feature: bool = True,
                 compact: bool = True, tmp_features: list[str] | None = None,
                 force: bool = False):
        self.rep1: SampleData = rep1
        self.rep2: SampleData = rep2
        self.ref: pl.DataFrame | None = union_ref
        self.feature: bool = feature
        self.out_f: Path = Path(f'{output_dir.as_posix()}/{self.rep1.lab}_{self.rep1.label}.csv')
        self.feature: bool = feature
        self.compact: bool = compact
        self.tmp_features: list[str] | None = tmp_features
        self.beta_col: str = beta_col
        self.depth: int = depth
        self.force: bool = force

        output_dir.mkdir(parents=True, exist_ok=True)

        logging.debug('🛠️ SpearmanR object initialized')
        logging.debug(str(self))

    def __repr__(self) -> str:
        return (f'SpearmanR('
                f'rep1={self.rep1!r}, '
                f'rep2={self.rep2!r}, '
                f'out_f={self.out_f!r}, '
                f'ref={self.ref!r}, '
                f'beta_col={self.beta_col!r}, '
                f'feature={self.feature!r}, '
                f'depth={self.depth})')

    def __str__(self) -> str:
        return (f'SpearmanR\n'
                f'---------\n'
                f'Replicate 1 Name: {self.rep1.name}\n'
                f'Replicate 2 Name: {self.rep2.name}\n'
                f'Output File: {self.out_f}\n'
                f'Beta Col to Evaluate: {self.beta_col}\n'
                f'Feature Mode: {'Enabled' if self.feature or self.tmp_features else 'Disabled'}\n'
                f'Temporary Features: {'Enabled' if self.tmp_features else 'Disabled'}\n'
                f'Depth Filter: depth >= {self.depth}')

    @staticmethod
    def _read_parquet(f_path: Path, suffix: str, lab: str,
                      depth: int, beta: str) -> pl.LazyFrame:
        tmp_df: pl.LazyFrame
        if lab[: 2] == 'MA':
            tmp_df = (pl.scan_parquet(f_path)
                      .with_columns(pl.col(beta).cast(pl.Float64))
                      .rename({beta: f'beta{suffix}'})
                      .select('chrom', 'start', f'beta{suffix}'))
        else:
            tmp_df = (pl.scan_parquet(f_path)
                        .with_columns(pl.col(beta).cast(pl.Float64))
                        .rename({beta: f'beta{suffix}',
                                 'depth': f'depth{suffix}'}))
            if depth > 1:
                tmp_df = tmp_df.filter(pl.col(f'depth{suffix}') >= depth)

            tmp_df = tmp_df.select('chrom', 'start', f'beta{suffix}', f'depth{suffix}')

        return tmp_df

    @staticmethod
    def _process_feature(df: pl.DataFrame, feature: str) -> pl.LazyFrame:
        return (resolve_feature(df=df, feature_col=feature)
                .group_by(feature)
                .agg(pl.len().alias('count'),
                     pl.corr(a=pl.col('beta1'),
                             b=pl.col('beta2'),
                             method='spearman')
                     .alias('spearmanr'))
                .rename({feature: 'fgroup'})
                .with_columns(pl.lit(feature).alias('feature')))

    def _compute(self):
        logging.info(f'ℹ️ Calculating Spearman correlation for {self.rep1.name} vs {self.rep2.name}')
        to_compute: pl.DataFrame
        _df = (self._read_parquet(f_path=self.rep1.f_path, lab=self.rep1.lab,
                                  suffix='1', depth=self.depth, beta=self.beta_col)
                   .join(other=self._read_parquet(f_path=self.rep2.f_path, lab=self.rep1.lab,
                                                  suffix='2', depth=self.depth, beta=self.beta_col),
                         on=['chrom', 'start'], how='inner'))
        if self.ref is not None:
            _df = (_df.join(other=self.ref.lazy(), on=['chrom', 'start'], how='left'))

        to_compute = (_df.drop_nulls(subset=['beta1', 'beta2'])
                         .with_columns(pl.lit('global').alias('global'))
                         .collect())
        if to_compute.is_empty():
            logging.warning(f'⚠️ No data to compute SpearmanR for {self.rep1.name} vs {self.rep2.name}')
            (pl.DataFrame({'lab': [self.rep1.lab], 'label': [self.rep1.label],
                           'feature': ['global'], 'fgroup': ['global'], 'count': [0], 'spearmanr': [-1]})
               .write_csv(self.out_f))
            return
        else:
            logging.debug(f'🛠️ Replicates and reference dataframe joined, containing {to_compute.shape[0]} rows')

        feature_cols: set[str]
        if self.tmp_features:
            feature_cols = set(self.tmp_features)
        elif self.compact:
            feature_cols = set(get_feature_cols(choice='compact_union'))
        else:
            match (self.feature, self.ref):
                case (True, pl.DataFrame()):
                    feature_cols = set(get_feature_cols(choice='union'))  # type: ignore
                case (True, None):
                    feature_cols = set(get_feature_cols(choice='union_without_ref'))
                case _:
                    feature_cols = {'global'}

        if self.rep1.lab[:2] == 'MA':
            feature_cols = {i for i in feature_cols if 'depth' not in i}

        feature_dfs: list[pl.DataFrame] = []

        logging.debug(f'🛠️ Collecting {len(feature_cols)} feature results')
        for feature_col in feature_cols:
            logging.debug(f'🛠️ Processing feature: {feature_col}')
            feature_dfs.append(self._process_feature(df=to_compute,
                                                     feature=feature_col)
                                   .with_columns(pl.col('fgroup').cast(pl.String))
                                   .collect())
            logging.debug(f'🛠️ {feature_col} processed')
        (pl.concat(feature_dfs)
           .with_columns([pl.lit(self.rep1.lab).alias('lab'),
                          pl.lit(self.rep1.label).alias('label')])
           .select(['lab', 'label', 'feature', 'fgroup', 'count', 'spearmanr'])
           .write_csv(self.out_f))
        logging.debug(f'🛠️ dataframes concatenated and saved at {self.out_f}')

    def run(self):
        if self.out_f.exists() and not self.force:
            logging.info(f'ℹ️ SpearmanR alreadly calcualted for {self.rep1.lab} {self.rep1.label}')
        logging.info(f'ℹ️ Calculating SpearmanR for {self.rep1.lab} {self.rep1.label}')
        self._compute()
        logging.info(f'✅ SpearmanR calculated for {self.rep1.lab} {self.rep1.label}')


class SignalNoiseRatio:
    def __init__(self, lab: str, samples: dict[str, SampleData],
                 output_dir: Path, depth: int, force: bool = False):
        self.lab: str = lab
        self.samples: dict[str, SampleData] = samples
        self.out_f: Path = Path(output_dir / f'{self.lab}.csv')
        self.depth: int = depth
        self.force: bool = force

        logging.debug('🛠️ SignalNoiseRatio object initialized')
        logging.debug(str(self))

    def __repr__(self) -> str:
        return (f'SignalNoiseRatio('
                f'samples={", ".join(f"{k!r}: {v!r}" for k, v in self.samples.items())},\n'
                f'depth={self.depth})')

    def __str__(self) -> str:
        return (f'SignalNoiseRatio\n'
                f'----------------\n'
                f'Sample Count: {len(self.samples)}\n'
                f'Samples:\n'
                f'{chr(10).join(f"  - {k}: {v.name}" for k, v in self.samples.items())}\n'
                'Depth Filter: {self.depth}')

    def _preprocess(self):
        logging.debug('🛠️ Preprocessing samples')
        if self.depth > 1 and next(iter(self.samples.values())).lab[: 2] != 'MA':
            return {
                sample: (pl.scan_parquet(sample_data.f_path)
                           .with_columns(pl.col('beta').cast(pl.Float64))
                           .rename({'beta': sample})
                           .filter(pl.col('depth') >= self.depth)
                           .select('chrom', 'start', sample))
                for sample, sample_data in self.samples.items()
            }
        else:
            return {
                sample: (pl.scan_parquet(sample_data.f_path)
                           .with_columns(pl.col('beta').cast(pl.Float64))
                           .rename({'beta': sample})
                           .select('chrom', 'start', sample))
                for sample, sample_data in self.samples.items()
            }

    def _compute(self):
        sample_count = len(self.samples)

        if sample_count < 2:
            raise ValueError('at least two samples are required to calculate SNR.')
        elif 2 <= sample_count < 8:
            logging.warning('SNR is not reliable with less than 8 samples:\n'
                            'e.g.: D5_1 D5_2 D6_1 D6_2 F7_1 F7_2 M8_1 M8_2')
        elif 8 <= sample_count < 16:
            logging.warning('it is recommended to use 16 samples:\n.'
                            'e.g.: D5_1 D5_2 D6_1 D6_2 F7_1 F7_2 M8_1 M8_2\n'
                            '      T1_1 T1_2 T2_1 T2_2 T3_1 T3_2 T4_1 T4_2')
        elif sample_count > 16:
            logging.warning('including non-quartet samples will increase the observed SNR')
        else:
            pass

        logging.debug('🛠️ Merging sample data')
        merged: pl.DataFrame
        merged = (reduce(lambda acc, x: acc.join(x,
                                                 on=['chrom', 'start'],
                                                 how='inner'),
                         [df for df in self._preprocess().values()])
                  .drop_nulls()
                  .drop('chrom', 'start')
                  .with_columns([pl.col(c).cast(pl.Int64)  # cast to int64 to reduce memory and CPU usage.
                                 for c in self.samples.keys()])
                  .filter(pl.concat_list(pl.all()).list.var() != 0)
                  .collect(engine='gpu'))  # ! GPU

        if merged.is_empty():
            logging.warning('No valid data after merging samples.')
            return (pl.DataFrame({'lab': [list(self.samples.values())[0].lab],
                                  'samples': [len(self.samples)],
                                  'count': [0], 'snr': [None]})
                      .lazy())

        logging.debug('🛠️ Sample data merged')
        logging.debug('🛠️ Performing PCA')
        pca = PCA(n_components=min(merged.shape[0], 8),
                  svd_solver='full')  # use 8 as alignment to previous Quartet studies.

        pcs_array: npt.NDArray[np.float64]
        pcs: dict[str, dict[str, float]]

        pcs_array = pca.fit_transform(StandardScaler().fit_transform(merged.transpose()  # type:ignore
                                                                           .to_pandas()))
        logging.debug('🛠️ PCA fitted')

        merged_length: int = merged.shape[0]
        merged_columns: list[str] = merged.columns

        del merged
        gc.collect()

        logging.debug('🛠️ Merged DataFrame deleted from RAM')

        pcs = {
            merged_columns[index]: {'PC1': pc[0], 'PC2': pc[1]}
            for index, pc in enumerate(pcs_array)
        }

        weights: npt.NDArray[np.float64] = pca.explained_variance_ratio_  # type:ignore
        weight1: float = float(weights[0])
        weight2: float = float(weights[1])

        dt_dist: pl.DataFrame
        dt_dist = (pl.DataFrame(data=[[sample1, sample1.split('_')[0],
                                       sample2, sample2.split('_')[0],
                                       'intra'
                                       if sample1.split('_')[0] == sample2.split('_')[0]
                                       else 'inter',
                                       weight1 * (pcs[sample1]['PC1'] - pcs[sample2]['PC1']) ** 2 +  # noqa:W504
                                       weight2 * (pcs[sample1]['PC2'] - pcs[sample2]['PC2']) ** 2]
                                      for sample1, sample2 in permutations(merged_columns, 2)],
                                orient='row',
                                schema={'sample1': pl.String, 'label1': pl.String,
                                        'sample2': pl.String, 'label2': pl.String,
                                        'relation': pl.String, 'dist': pl.Float64})
                   .group_by('relation')
                   .agg(pl.col('dist').mean().alias('mean_dist')))

        logging.debug('🛠️ Distance matrix calculated')

        snr: float = 10 * np.log10(dt_dist.select(pl.col('mean_dist')
                                                  .filter(pl.col('relation') == 'inter') /  # noqa:W504
                                                  pl.col('mean_dist')
                                                  .filter(pl.col('relation') == 'intra'))
                                   .item())
        logging.info(f'🛠️ SNR calculated: {snr}')
        logging.info(f'✅ SNR calculated for {self.lab}: {snr}')
        (pl.DataFrame({'lab': [list(self.samples.values())[0].lab],
                       'samples': [len(self.samples)],
                       'count': [merged_length],
                       'snr': [snr]})
            .write_csv(self.out_f))
        return None

    def run(self):
        if self.out_f.exists() and not self.force:
            logging.info(f'ℹ️ SNR already calculated for {self.lab}')
            return
        self._compute()


class DifferentialMethylation:
    pass


# -------------------------------------------------------------------------------------------------
# The following classes are used to generate and run the tasks
# -------------------------------------------------------------------------------------------------
class TaskManager:
    def __init__(self, output_dir: Path, bedgraph_files: list[Path], ref: ReferenceResolver,
                 level: Literal['cpg', 'c'], depth: int, beta_col: str,
                 pyro: bool = True, feature: bool = True, curated: bool = False,
                 compact: bool = False, tmp_rmse_features: list[str] | None = None,
                 tmp_spearmanr_features: list[str] | None = None,
                 rmse: bool = True, spearmanr: bool = True, spearmanr_with_ref: bool = True,
                 snr: bool = True, dmc: bool = False, gather: bool = False,
                 force: bool = True):
        self.rmse_dir: Path = output_dir / f'{level}_rmse'
        self.spearmanr_dir: Path = output_dir / f'{level}_spearmanr'
        self.snr_dir: Path = output_dir / f'{level}_snr'
        self.dmc_dir: Path = output_dir / f'{level}_dmc'
        self.ref: ReferenceResolver = ref
        self.level: Literal['cpg', 'c'] = level
        self.pyro: bool = pyro
        self.feature: bool = feature
        self.curated: bool = curated
        self.compact: bool = compact
        self.tmp_rmse_features: list[str] | None = tmp_rmse_features
        self.tmp_spearmanr_features: list[str] | None = tmp_spearmanr_features
        self.spearmanr_with_ref: bool = spearmanr_with_ref
        self.gather: bool = gather
        self.depth: int = depth
        self.beta_col: str = beta_col
        self.force: bool = force

        task_categories: set[str] = set()
        if rmse:
            task_categories.add('rmse')
        if spearmanr:
            task_categories.add('spearmanr')
        if snr:
            logging.warning('It is recommended to only run SNR for one lab one time')
            task_categories.add('snr')
        if dmc:
            task_categories.add('dmc')
        self.task_categories = task_categories

        for d in (self.rmse_dir, self.spearmanr_dir, self.snr_dir, self.dmc_dir):
            d.mkdir(parents=True, exist_ok=True)

        self.labs: list[LabProfile] = self._parse_bedgraphs(bedgraph_files)

    def _parse_bedgraphs(self, bedgraph_files: list[Path]):
        _legal_bedgraphs: dict[str, set[Path]] = {
            lab: {f for f in bedgraph_files if f.name.split('_')[0] == lab}
            for lab in {'BS1', 'BS2', 'BS3', 'BS4', 'EM1', 'EM2', 'EM3', 'EM4',
                        'PS1', 'PS2', 'PS3', 'NP1', 'RM1', 'RR1', 'TB1', 'MA1', 'MA2', 'MA3'}
            if len({f for f in bedgraph_files if f.name.split('_')[0] == lab}) > 0
        }

        return [
            LabProfileFactory.create_profile(name=lab, level=self.level, paths=lab_files)
            for lab, lab_files in _legal_bedgraphs.items()
        ]

    def _gather_results(self, gather_job: Literal['rmse', 'spearmanr', 'snr', 'dmc']):
        match gather_job:
            case 'rmse':
                logging.debug('🛠️ All RMSE calculated, gathering...')
                (pl.concat([pl.scan_csv(f, has_header=True,
                                        schema={'lab': pl.String, 'sample': pl.String,
                                                'feature': pl.String, 'fgroup': pl.String,
                                                'count': pl.Int64, 'rmse': pl.Float64})
                            for f in self.rmse_dir.glob('*.csv')])
                   .sort(by=['lab', 'sample', 'feature', 'fgroup'])
                   .sink_csv(self.rmse_dir / 'rmse.csv'))
                logging.debug(f'🛠️ All RMSE gathered and saved at {self.rmse_dir / "rmse.csv"}')

            case 'spearmanr':
                logging.debug('🛠️ All SpearmanR calculated, gathering...')
                (pl.concat([pl.scan_csv(f, has_header=True,
                                        schema={'lab': pl.String, 'label': pl.String,
                                                'feature': pl.String, 'fgroup': pl.String,
                                                'count': pl.Int64, 'spearmanr': pl.Float64})
                            for f in self.spearmanr_dir.glob('*.csv')])
                   .sort(by=['lab', 'label', 'feature', 'fgroup'])
                   .sink_csv(self.spearmanr_dir / 'spearmanr.csv'))
                logging.debug(f'🛠️ All SpearmanR gathered and saved at {self.spearmanr_dir / "spearmanr.csv"}')

            case _:
                raise NotImplementedError(f'Gathering for {gather_job} is not implemented')

    def run_tasks(self):
        available_labels: set[str] = {
            sample_data.label
            for lab in self.labs
            for sample_data in lab.DataSet.samples.values()
        }

        logging.debug(f'🛠️ Available labels: {available_labels}')

        if 'rmse' in self.task_categories:
            for label in sorted([_l for _l in available_labels if _l != 'HF']):
                label_ref = self.ref.get(ref_type=self.level, labels=[label]).collect()

                rmse_task_collection: set[SampleData] = {
                    _sd
                    for _l in self.labs
                    for _sd in sorted(_l.DataSet.samples.values())
                    if _sd.label == label
                }
                for sample_data in rmse_task_collection:
                    RootMeanSquareError(sample=sample_data, ref=label_ref,
                                        output_dir=self.rmse_dir, level=self.level,
                                        pyro=self.pyro, feature=self.feature,
                                        depth=self.depth, curated=self.curated,
                                        compact=self.compact,
                                        tmp_features=self.tmp_rmse_features,
                                        lab_beta=self.beta_col,
                                        force=self.force).run()
            if self.gather:
                self._gather_results(gather_job='rmse')

        if 'spearmanr' in self.task_categories:
            for label in available_labels:
                if self.spearmanr_with_ref:
                    union_ref = self.ref.get(ref_type='uc' if self.level == 'c' else 'ucpg',
                                             labels=[label]).collect()
                else:
                    union_ref = None

                spearmanr_task_collections: set[tuple[SampleData, SampleData]] = set()
                for lab in sorted(self.labs):
                    rep1, rep2 = None, None
                    for sample in sorted(lab.DataSet.samples.values()):
                        if sample.label == label:
                            if sample.rep == '1':
                                rep1 = sample
                            elif sample.rep == '2':
                                rep2 = sample
                            else:
                                pass
                            if rep1 is not None and rep2 is not None:
                                break
                    if rep1 and rep2:
                        spearmanr_task_collections.add((rep1, rep2))
                    else:
                        logging.warning(f'⚠️ No valid replicate pair found for {lab.Name} {label}')

                for rep1, rep2 in spearmanr_task_collections:
                    SpearmanR(rep1=rep1, rep2=rep2,
                              union_ref=union_ref,
                              output_dir=self.spearmanr_dir,
                              feature=self.feature,
                              compact=self.compact,
                              tmp_features=self.tmp_spearmanr_features,
                              depth=self.depth,
                              beta_col=self.beta_col,
                              force=self.force).run()
            if self.gather:
                self._gather_results(gather_job='spearmanr')

        if 'snr' in self.task_categories:
            for lab in self.labs:
                SignalNoiseRatio(lab=lab.Name,
                                 samples={k: v for k, v in lab.DataSet.samples.items()
                                          if k in ('D5_1', 'D5_2', 'D6_1', 'D6_2',
                                                   'F7_1', 'F7_2', 'M8_1', 'M8_2',
                                                   'T1_1', 'T1_2', 'T2_1', 'T2_2',
                                                   'T3_1', 'T3_2', 'T4_1', 'T4_2')},
                                 output_dir=self.snr_dir,
                                 depth=self.depth,
                                 force=self.force).run()
            if self.gather:
                self._gather_results(gather_job='snr')


def main():
    arg_parser = ArgumentParser()
    arg_parser.add_argument('-i', '--input-dir', dest='input_dir', type=Path,
                            required=True, help='input directory for the bedgraph files')
    arg_parser.add_argument('-o', '--output-dir', dest='output_dir', type=Path,
                            required=True, help='output directory for the results')
    arg_parser.add_argument('-l', '--level', dest='level', type=str,
                            choices=['cpg', 'c'], default='c',
                            help='level of the reference files')
    arg_parser.add_argument('--pyro', action='store_true', dest='pyro',
                            help='use pyro beta values')
    arg_parser.add_argument('--beta-col', type=str, dest='beta_col', default='beta',
                            help='The beta col in lab parquet lz4 file to evaluate')
    arg_parser.add_argument('--feature', action='store_true', dest='feature',
                            help='use feature columns')
    arg_parser.add_argument('--curated', action='store_true',
                            help='use curated features')
    arg_parser.add_argument('--compact', action='store_true',
                            help='use compact features')
    arg_parser.add_argument('-trf', '--tmp-rmse-features', type=parse_comma_separated,
                            dest='tmp_rmse_features', default=None,
                            help='temporary features for RMSE evaluation to use')
    arg_parser.add_argument('-tsf', '--tmp-spearmanr-features', type=parse_comma_separated,
                            dest='tmp_spearmanr_features', default=None,
                            help='temporary features for SpearmanR evaluation to use')
    arg_parser.add_argument('-d', '--depth', dest='depth', type=int, default=1,
                            help=('only depth greater than this value will be used,'
                                  'default is 1 (depth filter disabled)'))
    arg_parser.add_argument('--rmse', action='store_true', dest='rmse',
                            help='calculate RMSE')
    arg_parser.add_argument('--spearmanr', action='store_true', dest='spearmanr',
                            help='calculate SpearmanR')
    arg_parser.add_argument('--snr', action='store_true', dest='snr',
                            help='calculate SNR. !!!PLEASE DO NOT TRY THIS HERE UNLESS '
                                 'YOU HAVE RAM ≥1TB FOR THIS TASK!!!')
    arg_parser.add_argument('--gather', action='store_true', dest='gather',
                            help='gather all results')
    arg_parser.add_argument('--sf', action='store_true', dest='sf',
                            help='calculate spearman with ref or not')
    arg_parser.add_argument('--cpg-ref-dir', dest='truset_cpg_dir',
                            type=Path, default=TRUSET_CPG_DIR,
                            help='the directory of truSet cpg reference files')
    arg_parser.add_argument('--c-ref-dir', dest='truset_c_dir',
                            type=Path, default=TRUSET_C_DIR,
                            help='the directory of truSet c reference files')
    arg_parser.add_argument('--ucpg-ref-dir', dest='union_cpg_dir',
                            type=Path, default=UNION_CPG_DIR,
                            help='the directory of union cpg reference files')
    arg_parser.add_argument('--uc-ref-dir', dest='union_c_dir',
                            type=Path, default=UNION_C_DIR,
                            help='the directory of union c reference files')
    arg_parser.add_argument('--dmc-dir', dest='truset_dmc_dir',
                            type=Path, default=TRUSET_DMC_DIR,
                            help='the directory of truSet dmc reference files')
    arg_parser.add_argument('-s', '--suffix', type=str, dest='suffix', default='parquet.lz4',
                            help='the suffix to locate bedgraph files')
    arg_parser.add_argument('-w', '--whitelist', dest='whitelist', type=str,
                            default='', help='pattern of files that should only be included')
    arg_parser.add_argument('-b', '--blacklist', dest='blacklist', type=str,
                            default='', help='pattern of files that should not be included')
    arg_parser.add_argument('-v', '--log-level', dest='log_level',
                            choices=['debug', 'info', 'warning', 'error', 'critical'],
                            default='info')
    arg_parser.add_argument('-t', '--threads', dest='threads', type=int, default=8)
    arg_parser.add_argument('-f', '--force', dest='force', action='store_true',
                            help='ignore the existing output files and run.')

    args = arg_parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    os.environ['NUMEXPR_MAX_THREADS'] = str(args.threads)

    bedgraph_files = list(Path(args.input_dir).glob(f'*.{args.suffix}'))

    whitelist: list[str] = args.whitelist.split(',')
    blacklist: list[str] = args.blacklist.split(',')

    if whitelist:
        bedgraph_files: list[Path] = [i for i in bedgraph_files
                                      if any(white_word in i.name for white_word in whitelist)]
    else:
        if blacklist:
            bedgraph_files: list[Path] = [i for i in bedgraph_files
                                          if not any(black_word in i.name for black_word in blacklist)]
    if len(bedgraph_files) == 0:
        print('no bedgraph files found')
        return

    ref_collection = ReferenceResolver(cpg_dir=args.truset_cpg_dir, c_dir=args.truset_c_dir,
                                       ucpg_dir=args.union_cpg_dir, uc_dir=args.union_c_dir,
                                       dmc_dir=args.truset_dmc_dir)

    manager = TaskManager(output_dir=Path(args.output_dir),
                          bedgraph_files=bedgraph_files, ref=ref_collection,
                          level=args.level, pyro=args.pyro, feature=args.feature,
                          rmse=args.rmse, spearmanr=args.spearmanr, snr=args.snr,
                          gather=args.gather, spearmanr_with_ref=args.sf,
                          depth=args.depth, beta_col=args.beta_col,
                          curated=args.curated, compact=args.compact,
                          tmp_rmse_features=args.tmp_rmse_features,
                          tmp_spearmanr_features=args.tmp_spearmanr_features,
                          force=args.force)
    manager.run_tasks()


if __name__ == '__main__':
    main()
