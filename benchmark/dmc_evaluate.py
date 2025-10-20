from argparse import ArgumentParser, Namespace
from dataclasses import dataclass
import logging
from pathlib import Path
from typing import Literal
import pandas as pd
import polars as pl
import polars_ds as pds
from evaluate import resolve_feature


SAMPLE_PAIRS: list[str] = ['D6_vs_D5', 'D6_vs_F7', 'D6_vs_M8', 'BC_vs_BL']

tn: pl.Expr = pl.col('tn')
tp: pl.Expr = pl.col('tp')
fn: pl.Expr = pl.col('fn')
fp: pl.Expr = pl.col('fp')
recall: pl.Expr = tp / (tp + fn)
precision: pl.Expr = tp / (tp + fp)
f1_score: pl.Expr = 2 * tp / (2 * tp + fp + fn)
ppv: pl.Expr = tp / (tp + fp)
npv: pl.Expr = tn / (tn + fn)


@dataclass
class DmCallerInfo:
    dmc_output_pattern: str
    dmc_ref_suffix: dict[str, str]
    dmc_evaluated_fname: str


methylkit_info: DmCallerInfo
methylsig_beta_binomial_info: DmCallerInfo
methylsig_binomial_info: DmCallerInfo
methylsig_dss_info: DmCallerInfo
cpgtools_info: DmCallerInfo

methylkit_info = DmCallerInfo(dmc_output_pattern='methylkit',
                              dmc_ref_suffix={'p': 'methylkit.p.parquet.lz4',
                                              'q': 'methylkit.q.parquet.lz4'},
                              dmc_evaluated_fname='methylkit')
methylsig_beta_binomial_info = DmCallerInfo(dmc_output_pattern='beta_binomial',
                                            dmc_ref_suffix={'p': 'beta_binomial.p.parquet.lz4',
                                                            'q': 'beta_binomial.q.parquet.lz4'},
                                            dmc_evaluated_fname='methylsig.beta_binomial')
methylsig_binomial_info = DmCallerInfo(dmc_output_pattern='binomial',
                                       dmc_ref_suffix={'p': 'binomial.p.parquet.lz4',
                                                       'q': 'binomial.q.parquet.lz4'},
                                       dmc_evaluated_fname='methylsig.binomial')
methylsig_dss_info = DmCallerInfo(dmc_output_pattern='dss',
                                  dmc_ref_suffix={'p': 'dss.p.parquet.lz4',
                                                  'q': 'dss.q.parquet.lz4'},
                                  dmc_evaluated_fname='methylsig.dss')
cpgtools_info = DmCallerInfo(dmc_output_pattern='cpgtools',
                             dmc_ref_suffix={'p': 'cpgtools.p.parquet.lz4',
                                             'q': 'cpgtools.q.parquet.lz4'},
                             dmc_evaluated_fname='cpgtools')

TOOL_MAPPING: dict[str, DmCallerInfo] = {
    'methylkit': methylkit_info,
    'methylsig-beta-binomial': methylsig_beta_binomial_info,
    'methylsig-binomial': methylsig_binomial_info,
    'methylsig-DSS': methylsig_dss_info,
    'cpgtools': cpgtools_info
}


@dataclass
class LabInfo:
    lab: str
    pair: str
    s1: str
    s2: str
    lab_df: pl.LazyFrame


def mcc_pd(row: pd.Series):
    tp = row['tp']
    tn = row['tn']
    fp = row['fp']
    fn = row['fn']

    if tp + fp == 0 or tp + fn == 0 or tn + fp == 0 or tn + fn == 0:
        return None

    mcc_fz = tp * tn - fp * fn
    mcc_fm1 = (tp + fp) ** 0.5
    mcc_fm2 = (tp + fn) ** 0.5
    mcc_fm3 = (tn + fp) ** 0.5
    mcc_fm4 = (tn + fn) ** 0.5

    return mcc_fz / mcc_fm1 / mcc_fm2 / mcc_fm3 / mcc_fm4


def parse_input(file: Path, cutoff: float, p_type: Literal['p', 'q']) -> LabInfo:
    lab_df = pl.scan_csv(file, separator='\t', infer_schema_length=10000)

    P_TYPE_DICT: dict[str, dict[str, str]] = {
        'methylkit': {'p': 'pvalue', 'q': 'qvalue'},
        'methylsig': {'p': 'pvalue', 'q': 'fdr'},
        'cpgtools': {'p': 'pval', 'q': 'adj.pval'}
    }

    s1 = ''
    s2 = ''
    for pair in SAMPLE_PAIRS:
        if pair in file.name:
            s1 = pair.split('_vs_')[0]
            s2 = pair.split('_vs_')[1]
    if not s1 or not s2:
        raise ValueError('No sample pair found')

    lab_df: pl.LazyFrame

    if 'methylkit' in file.name:
        lab_df = lab_df.rename({P_TYPE_DICT['methylkit'][p_type]: 'p'})
    elif 'methylsig' in file.name:
        lab_df = lab_df.rename({P_TYPE_DICT['methylsig'][p_type]: 'p'})
    else:
        lab_df = (lab_df.rename({P_TYPE_DICT['cpgtools'][p_type]: 'p'})
                        .with_columns(pl.col('cgID')
                                        .str.split_exact('_', 1)
                                        .struct.rename_fields(['chrom', 'start'])
                                        .alias('chrom_start'))
                        .unnest('chrom_start'))
    lab_df = (lab_df.with_columns(pl.when(pl.col('p') < cutoff)
                                    .then(pl.lit(True))
                                    .otherwise(pl.lit(False))
                                    .alias('lab_DMC'))
                    .with_columns(pl.col('start').cast(pl.Int64))
                    .select('chrom', 'start', 'lab_DMC'))
    return LabInfo(lab=file.name[: 3], pair=f'{s1}_vs_{s2}', s1=s1, s2=s2, lab_df=lab_df)


def get_depth_info(cytosine_dir: Path, lab: str, s1_str: str, s2_str: str) -> pl.DataFrame:
    (s1r1,
     s1r2,
     s2r1,
     s2r2) = ([pl.scan_parquet(cytosine_dir / f'{lab}_{_s}_{_r}.parquet.lz4')
                 .with_columns(pl.col('start').cast(pl.Int64))
                 .select('chrom', 'start', 'depth')
                 .rename({'depth': f'{_s}_{_r}'})
               for _s, _r in [(s1_str, 1), (s1_str, 2), (s2_str, 1), (s2_str, 2)]])

    (s1,
     s2) = ([r1.join(r2, on=['chrom', 'start'], how='inner')
               .with_columns(pl.mean_horizontal(f'{_s}_1', f'{_s}_2').alias(f'{_str}_depth'))
               .select('chrom', 'start', f'{_str}_depth')
             for r1, r2, _s, _str in [(s1r1, s1r2, s1_str, 's1'), (s2r1, s2r2, s2_str, 's2')]])
    return s1.join(s2, on=['chrom', 'start'], how='inner').collect()


def get_feature_cols(bc: bool = False):
    depth_cols: list[str] = ['s1_depth_bin', 's2_depth_bin']
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
    if bc:
        feature_cols['genomic_cols'] += ([
            '150bp_SNVs', 'SNV_dist_bin', 'nearest_SNV_ref', 'nearest_SNV_alt',
            'nearest_SNV_ref:nearest_SNV_alt', '150bp_sINDELs', 'nearest_sINDEL_type',
            'sINDEL_dist_bin', 'nearest_sINDEL_size', '1000bp_SVs_bin', 'nearest_SV_type',
            'SV_dist_bin', 'SV_size_bin'])
    feature_pairs: list[str] = (
        [f's1_depth_bin:s2_depth_bin:{_col}'
         for _col in feature_cols['seq_cols'] + feature_cols['complexity_cols'] + feature_cols['genomic_cols']]
    )

    return set(depth_cols + ['s1_depth_bin:s2_depth_bin', 'global'] +  # noqa:W504
               feature_cols['seq_cols'] + feature_cols['complexity_cols'] +  # noqa:W504
               feature_cols['genomic_cols'] + feature_pairs)


def process_feature(df: pl.DataFrame, feature: str) -> pl.DataFrame:
    logging.debug(f'calculating for feature: {feature}')
    return (resolve_feature(df=df,
                            feature_col=feature)
            .lazy()
            .group_by(feature)
            .agg(pl.len().alias('count'),
                 (pl.col('DMC') & pl.col('lab_DMC')).sum().alias('tp'),
                 (pl.col('DMC') & ~pl.col('lab_DMC')).sum().alias('fn'),
                 (~pl.col('DMC') & pl.col('lab_DMC')).sum().alias('fp'),
                 (~pl.col('DMC') & ~pl.col('lab_DMC')).sum().alias('tn'),
                 pds.query_roc_auc('DMC_numeric', 'lab_DMC_numeric').alias('roc_auc'))
            .rename({feature: 'fgroup'})
            .with_columns(pl.lit(feature).alias('feature'))
            .collect())


def evaluate_file(ref_df: pl.DataFrame,
                  lab_file: LabInfo,
                  cytosine_dir: Path,
                  global_only: bool = False) -> pl.DataFrame:
    feature_cols: set[str]
    to_evaluate: pl.DataFrame
    if global_only:
        feature_cols = {'global'}
        to_evaluate = (ref_df.lazy()
                             .join(lab_file.lab_df, on=['chrom', 'start'], how='left')
                             .drop_nulls(subset=['lab_DMC', 'DMC'])
                             .with_columns([pl.col(c).cast(pl.Int64).alias(f'{c}_numeric')
                                            for c in ('DMC', 'lab_DMC')])
                             .with_columns(pl.lit('global').alias('global'))
                             .collect())
    else:
        logging.debug(f'retrieving depth info of {lab_file.lab} {lab_file.s1} {lab_file.s2} '
                      f'from {cytosine_dir.as_posix()}')
        depth_df: pl.DataFrame = get_depth_info(cytosine_dir=cytosine_dir, lab=lab_file.lab,
                                                s1_str=lab_file.s1, s2_str=lab_file.s2)
        to_evaluate = (ref_df.lazy()
                             .join(lab_file.lab_df, on=['chrom', 'start'], how='left')
                             .join(depth_df.lazy(), on=['chrom', 'start'], how='left')
                             .drop_nulls(subset=['lab_DMC', 'DMC'])
                             .with_columns([pl.col(c).cast(pl.Int64).alias(f'{c}_numeric')
                                            for c in ('DMC', 'lab_DMC')])
                             .with_columns(pl.lit('global').alias('global'))
                             .collect())
        logging.debug(f'dataframe to evaluate prepared, containing columns:\n'
                      f'{to_evaluate.columns}\n')
        if lab_file.s1 == 'BC':
            feature_cols = get_feature_cols(bc=True)
        else:
            feature_cols = get_feature_cols(bc=False)
    logging.debug(f'feature cols: {feature_cols}')

    feature_df: pl.DataFrame

    feature_df = (pl.concat([process_feature(df=to_evaluate,
                                             feature=feature)
                             .with_columns(pl.col('fgroup').cast(pl.String))
                             for feature in feature_cols])
                    .with_columns([pl.lit(lab_file.lab).alias('lab'),
                                   recall.alias('recall'),
                                   precision.alias('precision'),
                                   ppv.alias('ppv'),
                                   npv.alias('npv')])
                    .with_columns(f1_score.alias('f1_score')))

    return feature_df


def evaluate_pair(files: list[LabInfo],
                  cytosine_dir: Path,
                  ref_f: Path,
                  global_only: bool) -> pl.DataFrame:
    return pl.concat([evaluate_file(lab_file=_f, cytosine_dir=cytosine_dir,
                                    ref_df=pl.read_parquet(ref_f),
                                    global_only=global_only)
                      for _f in files])


def main():
    arg_parser: ArgumentParser = ArgumentParser()
    arg_parser.add_argument('-i', '--input-dir', dest='input_dir',
                            type=str, required=True,
                            help='DMC result file')
    arg_parser.add_argument('-r', '--ref-dir', dest='ref_dir',
                            type=str, required=True,
                            help='directory of DMC reference files')
    arg_parser.add_argument('-c', '--cytosine-dir', dest='cytosine_dir',
                            type=str, required=True,
                            help='directory of formatted cytosine files')
    arg_parser.add_argument('-o', '--output-dir', dest='output_dir',
                            type=str, required=True,
                            help='directory of evaluated results')
    arg_parser.add_argument('-m', '--method', type=str, required=True,
                            choices=['methylkit', 'methylsig-beta-binomial',
                                     'methylsig-binomial', 'methylsig-DSS', 'cpgtools'],
                            help='The DMC method to use')
    arg_parser.add_argument('-g', '--global-only', dest='global_only', action='store_true',
                            help='Only evaluate the global feature, not the pairwise features')
    arg_parser.add_argument('--p-type', type=str, dest='p_type', choices={'p', 'q'}, default='p',
                            help='whether to use p or q value for evaluation')
    arg_parser.add_argument('--cutoff', type=float, dest='cutoff', default=0.05,
                            help='p value cutoff for DMC detection')
    arg_parser.add_argument('-v', '--log-level', dest='log_level',
                            choices=['debug', 'info', 'warning', 'error', 'critical'],
                            default='info')

    args: Namespace = arg_parser.parse_args()
    input_dir: Path = Path(args.input_dir)
    ref_dir: Path = Path(args.ref_dir)
    cytosine_dir: Path = Path(args.cytosine_dir)
    output_dir: Path = Path(args.output_dir)
    method: str = args.method
    cutoff: float = args.cutoff
    p_type: Literal['p', 'q'] = args.p_type
    global_only: bool = args.global_only

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    output_dir.mkdir(parents=True, exist_ok=True)

    dmcaller_info: DmCallerInfo = TOOL_MAPPING[method]

    logging.info('runtime config:\n'
                 f'input_dir:\t{input_dir.as_posix()}\n'
                 f'ref_dir:\t{ref_dir.as_posix()}\n'
                 f'cytosine_dir:\t{cytosine_dir.as_posix()}\n'
                 f'output_dir:\t{output_dir.as_posix()}\n'
                 f'method:\t{method}\n'
                 f'cutoff:\t{cutoff}\n'
                 f'p_type:\t{p_type}\n'
                 f'granularity:\t{"global" if global_only else "features"}\n')

    input_list: list[LabInfo]
    if method == 'methylsig-binomial':
        input_list = [(parse_input(file=i, cutoff=cutoff, p_type=p_type))
                      for i in input_dir.glob(f'*{dmcaller_info.dmc_output_pattern}*')
                      if not any(f in i.name for f in ('BS0', 'EM0')) and 'beta_' not in i.name]
    else:
        input_list = [(parse_input(file=i, cutoff=cutoff, p_type=p_type))
                      for i in input_dir.glob(f'*{dmcaller_info.dmc_output_pattern}*')
                      if not any(f in i.name for f in ('BS0', 'EM0'))]

    sample_pair_results: list[pl.DataFrame] = []
    for pair in SAMPLE_PAIRS:
        pair_files = [i for i in input_list if pair == i.pair]
        if pair_files:
            logging.info(f'Evaluating {len(pair_files)} files of {pair}')
            sample_pair_results.append((evaluate_pair(cytosine_dir=cytosine_dir,
                                                      files=pair_files,
                                                      ref_f=ref_dir / f'{pair}.{dmcaller_info.dmc_ref_suffix[p_type]}',
                                                      global_only=global_only))
                                       .with_columns(pl.lit(pair).alias('sample_pair')))
        else:
            logging.warning(f'No files for {pair}, skipping')

    df = (pl.concat(sample_pair_results)
            .with_columns(pl.lit(p_type).alias('p_type'))
            .select('sample_pair', 'lab', 'p_type', 'feature', 'fgroup', 'count',
                    'tp', 'tn', 'fn', 'fp', 'recall', 'precision', 'f1_score',
                    'ppv', 'npv', 'roc_auc')
            .to_pandas())
    df['mcc'] = df.apply(mcc_pd, axis=1)

    df[['sample_pair', 'lab', 'p_type', 'feature', 'fgroup', 'count',
        'tp', 'tn', 'fn', 'fp', 'recall', 'precision', 'f1_score', 'mcc',
        'ppv', 'npv', 'roc_auc']].to_csv(output_dir / f'{dmcaller_info.dmc_evaluated_fname}.{p_type}.csv',
                                         index=False)


if __name__ == '__main__':
    main()
