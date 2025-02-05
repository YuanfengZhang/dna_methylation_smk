from pathlib import Path
from typing import List, Set
import numpy as np
import pandas as pd
# ? We have to use pandas since it is the default choice in snakemake.
# ? And the sample sheet is too tiny to use the polars library.

"""
To explain...
# 0. The reason why the sample sheet is parsed here instead of main smk file is
to 1) make the main smk file clean;
   2) decouple, since the snakemake API is not stable;
   3) scrutiny the sample sheet in detail;
   4) provide a clear error message when the sample sheet is wrong.
# 1. FNAME should be the stem of R1.fq.gz and R2.fq.gz files. E.g.,
for BS_HCC1395_2.R1.fq.gz and BS_HCC1395_2.R2.fq.gz, FNAME should be BS_HCC1395_2.
The files will be searched in the input dir, which you should specify in the config.yaml,
using the pathlib.Path().glob() method. If the Snakemake complains about missing files
or more than 2 files found, please debug yourself.

# 2. TRIMMER is the tool to trim the adapters as well as low quality reads.
# 3. QC_REPORTER is the statistic tool to check out the sequencing quality.
# 4. ALIGNER is the alignment tool to use for this sample.
# 5. DEDUPER is the tool to remove duplicated reads which are labeled by the aligner.
# 6. RECALIBRATE is the choice to apply BQSR on the aligned bam file.
# 7. COUNTER is the tool to get per-CpG methylation level.
# 8. STAS is the choice to get the mapping quality and methylation level statistics.

#! 9. The names of tools are in the restrict vocabulary below.
# 10. It's possible to run multiple tools for a single sample. Let's say we want to compare
the results of bismark_bowtie2 and bwa-meth. We can add two rows in the sample sheet,
with same other columns but different ALIGNER values.
#! In other words, one row = one run. Multiple tools in one cell is not allowed,
because some tools are not compatible and all combinations once is a horrible choice.
"""

REQUIRED_COLS: List[str] = ['FNAME', 'TRIMMER', 'QC_REPORTER',
                            'ALIGNER', 'DEDUPER', 'RECALIBRATE',
                            'COUNTER', 'STATS']

# ! To get rid of spelling issues,
# ! only lower case names and hyphen-minus (-) are allowed.
AVAILABLE_TRIMMERS: Set[str] = {'fastp', 'trim-galore'}
AVAILABLE_QC_REPORTERS: Set[str] = {'fastqc', 'falco'}  # no fastp here
AVAILABLE_ALIGNERS: Set[str] = {'bismark-bowtie2', 'bismark_-hisat2',
                                'bwa-meth', 'bwa-meme', 'dnmtools',
                                'biscuit', 'bsmapz'}
AVAILABLE_DEDUPERS: Set[str] = {'bismark', 'gatk', 'samtools', 'sambamba'}
RECALIBRATE: Set[bool] = {True, False}
COUNTER: Set[str] = {'bismark', 'biscuit', 'methyldackel',
                     'dnmtools', 'bs_seeker2', 'astair'}
STATS: Set[bool] = {True, False}


class OneRun:
    def __init__(self, row: pd.Series):
        self.row: pd.Series = row
        self.fname: str = row['FNAME']
        self.trimmer: str | None = row['TRIMMER']
        self.qc_reporter: str | None = row['QC_REPORTER']
        self.aligner: str | None = row['ALIGNER']
        self.deduper: str | None = row['DEDUPER']
        self.recalibrate: bool = row['RECALIBRATE']
        self.counter: str | None = row['COUNTER']
        self.stats: bool = row['STATS']

    def __repr__(self):
        return (f'{self.fname}: {self.trimmer}, {self.qc_reporter},'
                f'{self.aligner}, {self.deduper},'
                f'{self.recalibrate}, {self.counter}, {self.stats}')

    def __str__(self):
        return (f'{self.fname}: {self.trimmer}, {self.qc_reporter},'
                f'{self.aligner}, {self.deduper},'
                f'{self.recalibrate}, {self.counter}, {self.stats}')

    def _validate_tool(self,
                       tool: str | None,
                       available_tools: Set[str],
                       tool_name: str) -> bool:
        if tool in available_tools:
            return True
        else:
            raise ValueError(f'Not implemented yet: {tool}.'
                             f'Available {tool_name}: {available_tools}')
    """
    # TODO: The output file paths are absolutely wrong and need fix.
    There should be a single output dir for each fname, e,g., result/{self.fname}
    And the output files should be in nested dirs,
    e.g., result/{self.fname}/fastp/bismark-bowtie2/bismark/bqsr/bismark/{self.fname}.bedgraph.gz
                    ^        ^          ^          ^       ^     ^
                    fname    trimmer    aligner    deduper recal counter
    or result/{self.fname}/fastp/bismark-bowtie2/bismark/no_bqsr/bismark/{self.fname}.bedgraph.gz
                    ^        ^          ^          ^       ^     ^
                    fname    trimmer    aligner    deduper recal counter
    """
    def _generate_trimmed_file_ls(self) -> List[str]:
        trim_files: List[str]
        # check the trimmer first
        if self._validate_tool(self.trimmer,
                               available_tools=AVAILABLE_TRIMMERS,
                               tool_name='trimmer'):
            match self.trimmer:
                case 'fastp':
                    trim_files = [
                        f'result/{self.fname}/{self.trimmer}/{self.fname}.R1.fq.gz',
                        f'result/{self.fname}/{self.trimmer}/{self.fname}.R2.fq.gz',
                        f'result/{self.fname}/{self.trimmer}/{self.fname}.fastp.json',
                        f'result/{self.fname}/{self.trimmer}/{self.fname}.fastp.html'
                    ]
                case 'trim-galore':
                    trim_files = [
                        f'result/{self.fname}/{self.trimmer}/{self.fname}.R1.fq.gz',
                        f'result/{self.fname}/{self.trimmer}/{self.fname}.R2.fq.gz'
                    ]
                case _:
                    trim_files = []
        else:
            trim_files = []

        return trim_files

    def _generate_qc_report_file_ls(self) -> List[str]:
        qc_files: List[str]
        if self._validate_tool(self.qc_reporter,
                               available_tools=AVAILABLE_QC_REPORTERS,
                               tool_name='qc_reporter'):
            match self.qc_reporter:
                case 'fastqc':
                    qc_files = [
                        f'result/{self.fname}/{self.trimmer}/{self.fname}.R1_fastqc.html',
                        f'result/{self.fname}/{self.trimmer}/{self.fname}.R2_fastqc.html'
                    ]
                case 'falco':
                    qc_files = [
                        f'result/{self.fname}/{self.trimmer}/{self.fname}.R1_falco.html',
                        f'result/{self.fname}/{self.trimmer}/{self.fname}.R2_falco.html'
                    ]
                case _:
                    qc_files = []
        else:
            qc_files = []

        return qc_files

    def _generate_aligned_bam(self) -> List[str]:
        aligned_files: List[str]
        if self._validate_tool(self.aligner,
                               available_tools=AVAILABLE_ALIGNERS,
                               tool_name='aligner'):
            aligned_files = [
                (f'result/{self.fname}/{self.trimmer}/'
                 f'{self.aligner}/{self.fname}.bam')
            ]
        else:
            aligned_files = []

        return aligned_files

    def _generate_deduped_bam(self) -> List[str]:
        dedup_files: List[str]
        if self._validate_tool(self.deduper,
                               available_tools=AVAILABLE_DEDUPERS,
                               tool_name='deduper'):
            dedup_files = [
                (f'result/{self.fname}/{self.trimmer}/'
                 f'{self.aligner}/{self.deduper}/'
                 f'{self.fname}.bam')
            ]
        else:
            dedup_files = []

        return dedup_files

    def _generate_recalibrated_bam(self) -> List[str]:
        recalibrated_files: List[str]
        if self.recalibrate:
            recalibrated_files = [
                (f'result/{self.fname}/{self.trimmer}/'
                 f'{self.aligner}/{self.deduper}/'
                 f'{self.fname}.bqsr.bam')
            ]
        else:
            recalibrated_files = []

        return recalibrated_files

    def _generate_qualimap_files(self) -> List[str]:
        qualimap_files: List[str]
        if self.recalibrate:
            qualimap_files = [
                (f'result/{self.fname}/{self.trimmer}/'
                 f'{self.aligner}/{self.deduper}/'
                 f'{self.fname}.bqsr.qualimap.pdf')
            ]
        else:
            qualimap_files = [
                (f'result/{self.fname}/{self.trimmer}/'
                 f'{self.aligner}/{self.deduper}/'
                 f'{self.fname}.qualimap.pdf')
            ]

        return qualimap_files

    def _generate_methylation(self) -> List[str]:
        methylation_files: List[str]
        if self._validate_tool(self.counter,
                               available_tools=COUNTER,
                               tool_name='counter'):
            if self.recalibrate:
                match self.counter:
                    case 'bismark':
                        methylation_files = [
                            (f'result/{self.fname}/{self.trimmer}/'
                             f'{self.aligner}/{self.deduper}/'
                             f'{self.counter}/bqsr/'
                             f'{self.fname}.bedgraph.gz'),
                            (f'result/{self.fname}/{self.trimmer}'
                             f'/{self.aligner}/{self.deduper}'
                             f'/{self.counter}/bqsr/'
                             f'{self.fname}.ucsc.bedgraph.gz')
                        ]
                    case _:
                        methylation_files = []
            else:
                match self.counter:
                    case 'bismark':
                        methylation_files = [
                            (f'result/{self.fname}/{self.trimmer}/'
                             f'{self.aligner}/{self.deduper}/'
                             f'{self.counter}/no_bqsr/'
                             f'{self.fname}.bedgraph.gz'),
                            (f'result/{self.fname}/{self.trimmer}'
                             f'/{self.aligner}/{self.deduper}'
                             f'/{self.counter}/no_bqsr/'
                             f'{self.fname}.ucsc.bedgraph.gz')
                        ]
                    case _:
                        methylation_files = []
        else:
            methylation_files = []

        return methylation_files

    def _generate_stats(self) -> List[str]:
        stats_files: List[str]

        match (self.recalibrate, self.stats):
            case (_, False):
                stats_files = []
            case (False, True):
                stats_files = [
                    (f'result/{self.fname}/{self.trimmer}'
                     f'/{self.aligner}/{self.deduper}'
                     f'/{self.counter}/no_bqsr/'
                     f'{self.fname}.ucsc.bedgraph.stats')
                ]
            case (True, True):
                stats_files = [
                    (f'result/{self.fname}/{self.trimmer}'
                     f'/{self.aligner}/{self.deduper}'
                     f'/{self.counter}/bqsr/'
                     f'{self.fname}.ucsc.bedgraph.stats')
                ]
            case _:
                stats_files = []

        return stats_files

    def generate_output_files_ls(self) -> List[str]:
        output_ls = []

        # add qualimap files for aligned files.
        if any(x for x in (self.aligner, self.deduper,
                           self.counter, self.stats)):
            output_ls.extend(self._generate_qualimap_files())

        # early return since the intermediate files can be solved automatically.
        if self.counter:
            output_ls.extend(self._generate_methylation())
            if self.stats:
                output_ls.extend(self._generate_stats())
            return output_ls

        # handle the situation when the bedgraph is not needed.
        if self.recalibrate:
            output_ls.extend(self._generate_recalibrated_bam())
            return output_ls

        if self.deduper:
            output_ls.extend(self._generate_deduped_bam())
            return output_ls

        if self.aligner:
            output_ls.extend(self._generate_aligned_bam())
            return output_ls

        # handle the situation when the bam is not needed.
        if self.trimmer:
            output_ls.extend(self._generate_trimmed_file_ls())
            if self.qc_reporter:
                output_ls.extend(self._generate_qc_report_file_ls())
            return output_ls

        raise ValueError('No valid tools found in the sample sheet.')


def clear_row(row: pd.Series) -> pd.Series:
    for col_idx in range(len(row)):
        if pd.isna(row.iloc[col_idx]):
            row.iloc[col_idx + 1:] = np.nan
            break
    return row


def read_sample_sheet(tsv_path: Path) -> List[str]:
    df: pd.DataFrame
    # There should be no missing values in the FNAME col
    df = (pd.read_table(tsv_path)
            .apply(clear_row, axis=1)
            .dropna(subset=['FNAME'])
            .replace(to_replace=np.nan, value=None))
    # check if the sample sheet has the required columns
    if df.columns.tolist() != REQUIRED_COLS:
        raise ValueError('Wrong columns in sample sheet.'
                         f'Expected: {REQUIRED_COLS}'
                         f'Found: {df.columns.tolist()}')
    if df.empty:
        raise ValueError('Missing values in FNAME col is not allowed.')

    print(df)

    all_files = []

    for idx, row in df.iterrows():
        try:
            run = OneRun(row)
            all_files.extend(run.generate_output_files_ls())
        except ValueError as e:
            print(f'at {idx} there is an error:\n{e}')
            continue

    if not all_files:
        raise ValueError('No valid runs found in the sample sheet.')

    return all_files


def hook_test():
    print('hook test succeed')


if __name__ == '__main__':
    for i in sorted(read_sample_sheet(Path('../test/example.tsv'))):
        print(i)
