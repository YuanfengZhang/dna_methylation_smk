from pathlib import Path
from typing import List, Set
import warnings
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

REQUIRED_COLS: Set[str] = {'FNAME', 'TRIMMER', 'QC_REPORTER',
                           'ALIGNER', 'DEDUPER', 'RECALIBRATE',
                           'COUNTER', 'STATS'}

# ! To get rid of spelling issues,
# ! only lower case names and underlines are allowed.
AVAILABLE_TRIMMERS: Set[str] = {'fastp', 'trim_galore'}
AVAILABLE_QC_REPORTERS: Set[str] = {'fastqc', 'falco'}  # no fastp here
AVAILABLE_ALIGNERS: Set[str] = {'bismark_bowtie2', 'bismark_hisat2',
                                'bwa_meth', 'bwa_meme', 'dnmtools',
                                'biscuit', 'bsmapz'}
AVAILABLE_DEDUPERS: Set[str] = {'bismark', 'gatk', 'samtools', 'sambamba'}
RECALIBRATE: Set[bool] = {True, False}
COUNTER: Set[str] = {'bismark', 'biscuit', 'methyldackel',
                     'dnmtools', 'bs_seeker2', 'astair'}
STATS: Set[bool] = {True, False}


class OneRun:
    def __init__(self, row: pd.Series):
        self.row = row
        self.fname = row['FNAME']
        self.trimmer = row['TRIMMER']
        self.qc_reporter = row['QC_REPORTER']
        self.aligner = row['ALIGNER']
        self.deduper = row['DEDUPER']
        self.recalibrate = row['RECALIBRATE']
        self.counter = row['COUNTER']
        self.stats = row['STATS']

    def __repr__(self):
        return (f'{self.fname}: {self.trimmer}, {self.qc_reporter},'
                f'{self.aligner}, {self.deduper},'
                f'{self.recalibrate}, {self.counter}, {self.stats}')

    def __str__(self):
        return (f'{self.fname}: {self.trimmer}, {self.qc_reporter},'
                f'{self.aligner}, {self.deduper},'
                f'{self.recalibrate}, {self.counter}, {self.stats}')

    def _generate_trimmed_file_ls(self) -> List[str]:
        # check the trimmer first
        if pd.isna(self.trimmer):
            warnings.warn(message=('No trimmer specified!'
                                   'Skipping this row.'),
                          category=UserWarning)
            return []
        if self.trimmer not in AVAILABLE_TRIMMERS:
            raise ValueError(f'Unknown trimmer: {self.trimmer}.'
                             f'Available trimmers: {AVAILABLE_TRIMMERS}')
        match self.trimmer:
            case 'fastp':
                return [f'result/trimmed/fastp/{self.fname}.R1.fq.gz',
                        f'result/trimmed/fastp/{self.fname}.R2.fq.gz',
                        f'result/trimmed/fastp/{self.fname}.fastp.json',
                        f'result/trimmed/fastp/{self.fname}.fastp.html']
            case 'trim_galore':
                return [f'result/trimmed/trim_galore/{self.fname}.R1.fq.gz',
                        f'result/trimmed/trim_galore/{self.fname}.R2.fq.gz',
                        f'result/trimmed/trim_galore/{self.fname}_trimming_report.txt']
            case _:
                raise ValueError(f'Not implemented yet: {self.trimmer}')

    def _generate_qc_report_file_ls(self) -> List[str]:
        if any(pd.isna(_t) for _t in (self.trimmer, self.qc_reporter)):
            warnings.warn(message=('No trimmer or qc_reporter specified!'
                                   'Skipping this row.'),
                          category=UserWarning)
            return []
        if self.qc_reporter not in AVAILABLE_QC_REPORTERS:
            raise ValueError(f'Unknown qc_reporter: {self.qc_reporter}.'
                             f'Available qc_reporters: {AVAILABLE_QC_REPORTERS}')
        match self.qc_reporter:
            case 'fastqc':
                return [f'result/qc/fastqc/{self.fname}.R1_fastqc.html',
                        f'result/qc/fastqc/{self.fname}.R2_fastqc.html']
            case 'falco':
                return [f'result/qc/falco/{self.fname}.R1_falco.html',
                        f'result/qc/falco/{self.fname}.R2_falco.html']
            case _:
                raise ValueError(f'Not implemented yet: {self.qc_reporter}')

    def _generate_aligned_bam(self) -> List[str]:
        if any(pd.isna(_t) for _t in (self.trimmer,
                                      self.qc_reporter,
                                      self.aligner)):
            warnings.warn(message=('No trimmer or qc_reporter or aligner specified!'
                                   'Skipping this row.'),
                          category=UserWarning)
            return []
        if self.aligner not in AVAILABLE_ALIGNERS:
            raise ValueError(f'Unknown aligner: {self.aligner}.'
                             f'Available aligners: {AVAILABLE_ALIGNERS}')
        match self.aligner:
            case 'bismark_bowtie2':
                return [f'result/aligned/bismark_bowtie2/{self.fname}.bam']
            case 'bismark_hisat2':
                return [f'result/aligned/bismark_hisat2/{self.fname}.bam']
            case 'bwa_meth':
                return [f'result/aligned/bwa_meth/{self.fname}.bam']
            case 'bwa_meme':
                return [f'result/aligned/bwa_meme/{self.fname}.bam']
            case 'biscuit':
                return [f'result/aligned/biscuit/{self.fname}.bam']
            case 'bsmapz':
                return [f'result/aligned/bsmapz/{self.fname}.bam']
            case 'dnmtools':
                return [f'result/aligned/dnmtools/{self.fname}.bam']
            case _:
                raise ValueError(f'Not implemented yet: {self.aligner}')

    def generate_output_files_ls(self):
        return (self._generate_trimmed_file_ls() +
                self._generate_qc_report_file_ls() +
                self._generate_aligned_bam())


def read_sample_sheet(tsv_path: Path) -> pd.DataFrame:
    df = pd.read_table(tsv_path)
    # check if the sample sheet has the required columns
    if df.columns.tolist() != REQUIRED_COLS:
        raise ValueError('Wrong columns in sample sheet.'
                         f'Expected: {REQUIRED_COLS}'
                         f'Found: {df.columns.tolist()}')

    # There should be no missing values in the FNAME col.
    if df['FNAME'].isnull().any():
        raise ValueError('Missing values in FNAME col is not allowed.')

    return df
