from argparse import ArgumentParser, Namespace
from pathlib import Path
from typing import Dict, List, Set
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
# 8. STATS is the choice to get the mapping quality and methylation level statistics.

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

    def _generate_trimmed_file_ls(self) -> List[str]:
        trim_files: List[str]
        # check the trimmer first
        if self._validate_tool(self.trimmer,
                               available_tools=AVAILABLE_TRIMMERS,
                               tool_name='trimmer'):
            match self.trimmer:
                case 'fastp':
                    trim_files = [
                        f'result/{self.fname}/fastp/'
                        f'{self.fname}{ext}'
                        for ext in ('.R1.fq.gz',
                                    '.R2.fq.gz',
                                    '.fastp.json',
                                    '.fastp.html')]
                case 'trim-galore':
                    trim_files = [
                        f'result/{self.fname}/trim-galore/'
                        f'{self.fname}{ext}'
                        for ext in ('.R1.fq.gz',
                                    '.R2.fq.gz',
                                    '.R1.report',
                                    '.R2.report')]
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
                        f'result/{self.fname}/{self.trimmer}/'
                        f'{self.fname}{ext}'
                        for ext in ('.R1.fastqc.html',
                                    '.R2.fastqc.html')]
                case 'falco':
                    qc_files = [
                        f'result/{self.fname}/{self.trimmer}/'
                        f'{self.fname}{ext}'
                        for ext in ('.R1.falco.html',
                                    '.R2.falco.html')]
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
            match self.deduper:
                case 'bismark':
                    dedup_files = [
                        f'result/{self.fname}/{self.trimmer}/'
                        f'{self.aligner}/bismark/'
                        f'{self.fname}{ext}'
                        for ext in ('.bam', '.bam.bai')
                    ]
                case 'gatk':
                    dedup_files = [
                        f'result/{self.fname}/{self.trimmer}/'
                        f'{self.aligner}/gatk/{self.fname}{ext}'
                        for ext in ('.bam',
                                    '.bam.bai',
                                    '.dedup.metrics')
                    ]
                case 'samtools':
                    dedup_files = [
                        f'result/{self.fname}/{self.trimmer}/'
                        f'{self.aligner}/samtools/{self.fname}{ext}'
                        for ext in ('.bam',
                                    '.bam.bai',
                                    '.dup.stats')
                    ]
                case 'no_dedup':
                    dedup_files = [
                        f'result/{self.fname}/{self.trimmer}/'
                        f'{self.aligner}/no_dedup/{self.fname}{ext}'
                        for ext in ('.bam',
                                    '.bam.bai')
                    ]
                case _:
                    dedup_files = []
        else:
            dedup_files = []

        return dedup_files

    def _generate_recalibrated_bam(self) -> List[str]:
        recalibrated_files: List[str]
        if self.recalibrate:
            recalibrated_files = [
                f'result/{self.fname}/{self.trimmer}/'
                f'{self.aligner}/{self.deduper}/{self.fname}{ext}'
                for ext in ('.bqsr.bam',
                            '.bqsr.bam.bai',
                            '.before.table',
                            '.after.table',
                            '.AnalyzeCovariates.pdf')
            ]
        else:
            recalibrated_files = []

        return recalibrated_files

    def _generate_qualimap_files(self) -> List[str]:
        qualimap_files: List[str]

        if self.recalibrate:
            qualimap_files = [
                f'result/{self.fname}/{self.trimmer}/'
                f'{self.aligner}/{self.deduper}/qualimap_bqsr/{ext}'
                for ext in (f'{self.fname}.bqsr.qualimap.pdf',
                            'qualimapReport.html')
            ]
        else:
            qualimap_files = [
                f'result/{self.fname}/{self.trimmer}/'
                f'{self.aligner}/{self.deduper}/qualimap/{ext}'
                for ext in (f'{self.fname}.qualimap.pdf',
                            'qualimapReport.html')
            ]

        return qualimap_files

    def _generate_methylation(self) -> List[str]:
        pattern: str
        methylation_files: List[str]
        if self._validate_tool(self.counter,
                               available_tools=COUNTER,
                               tool_name='counter'):
            if self.recalibrate:
                pattern = '.bqsr'
            else:
                pattern = ''

            match self.counter:
                case 'bismark':
                    methylation_files = [
                        f'result/{self.fname}/{self.trimmer}/'
                        f'{self.aligner}/{self.deduper}/'
                        f'bismark/{self.fname}{pattern}{ext}'
                        for ext in ('.bedgraph.gz', '.ucsc_bedgraph.gz',
                                    '.CHG_OB.txt.gz', '.CHG_OT.txt.gz',
                                    '.CHH_OB.txt.gz', '.CHH_OT.txt.gz',
                                    '.CpG_OB.txt.gz', '.CpG_OT.txt.gz',
                                    '.bismark.cov.gz', '.splitting_report',
                                    '.c2c.cov.gz', '.c2c.report.gz',
                                    '.c2c.summary')]
                case 'astair':
                    methylation_files = [
                        f'result/{self.fname}/{self.trimmer}/'
                        f'{self.aligner}/{self.deduper}/'
                        f'astair/{self.fname}{pattern}{ext}'
                        for ext in ('.astair.mods.gz',
                                    '.astair.stats')]
                case 'methyldackel':
                    methylation_files = [
                        f'result/{self.fname}/{self.trimmer}/'
                        f'{self.aligner}/{self.deduper}/'
                        f'methyldackel/{self.fname}{pattern}{ext}'
                        for ext in ('.bedgraph.gz',
                                    '.merged.bedgraph.gz')]
                case _:
                    methylation_files = []

        else:
            methylation_files = []

        return methylation_files

    def _generate_stats(self) -> List[str]:
        pattern: str
        stats_files: List[str]

        if self.stats:
            if self.recalibrate:
                pattern = '.bqsr'
            else:
                pattern = ''

            if self.aligner == 'bismark':
                stats_files = [
                    f'result/{self.fname}/{self.trimmer}/'
                    f'{self.aligner}/{self.deduper}/'
                    f'{self.fname}{pattern}{ext}'
                    for ext in ('.samtools.stats.txt',
                                '.samtools.flagstats.txt',
                                '.methydackel_mbias_OT.svg',
                                '.methydackel_mbias_OB.svg'
                                '.bam2nuc.txt')
                ]
            else:
                stats_files = [
                    f'result/{self.fname}/{self.trimmer}/'
                    f'{self.aligner}/{self.deduper}/'
                    f'{self.fname}{pattern}{ext}'
                    for ext in ('.samtools.stats.txt',
                                '.samtools.flagstats.txt',
                                '.methydackel_mbias_OT.svg',
                                '.methydackel_mbias_OB.svg')]
                # well, it's because bam2nuc would recognize the bam file
                # from other aligners as single-end and throw an error.

                # also, astair mbias and IDbias files are excluded because it will takes years
                # to run on an actual ~30x pair-end WGBS sample.
                # '.astair.IDbias_abundance_10bp_mod_site.pdf',
                # '.astair.IDbias_abundance.pdf',
                # '.astair.IDbias_indel_rate.pdf',
                # '.astair.IDbias_mod_co-localize.pdf',
                # '.astair.IDbias.stats',
                # '.astair.Mbias.stats',
                # '.astair.Mbias.pdf'

        else:
            stats_files = []

        return stats_files

    def generate_output_files_ls(self) -> List[str]:
        output_ls = []

        # early return since the intermediate files can be solved automatically.
        if self.counter:
            output_ls.extend(self._generate_methylation())
            if self.stats:
                output_ls.extend(self._generate_stats())
                output_ls.extend(self._generate_qualimap_files())
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


def read_sample_sheet(csv_path: Path) -> pd.DataFrame:
    df: pd.DataFrame
    # There should be no missing values in the FNAME col
    df = (pd.read_csv(csv_path, comment='#',
                      dtype={'FNAME': 'str', 'TRIMMER': 'str',
                             'QC_REPORTER': 'str', 'ALIGNER': 'str',
                             'DEDUPER': 'str', 'RECALIBRATE': 'boolean',
                             'COUNTER': 'str', 'STATS': 'boolean'})
            .apply(clear_row, axis=1)
            .dropna(subset=['FNAME'])
            .replace(to_replace=np.nan, value=None))
    # We have to skip the deduplication for RRBS
    df['DEDUPER'] = (
        df.apply(lambda row: 'no_dedup'
                 if row['FNAME'].startswith('RR') and row['DEDUPER']
                 else row['DEDUPER'],
                 axis=1))
    # Check if the sample sheet has the required columns
    if df.columns.tolist() != REQUIRED_COLS:
        raise ValueError('Wrong columns in sample sheet.'
                         f'Expected: {REQUIRED_COLS}'
                         f'Found: {df.columns.tolist()}')
    if df.empty:
        raise ValueError('Missing values in FNAME col is not allowed.')

    return df


def generate_file_ls(csv_path: Path) -> List[str]:
    df = read_sample_sheet(csv_path)

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


def generate_tool_ls(csv_path: Path) -> Dict[str, List[str]]:
    df = read_sample_sheet(csv_path)
    tool_dict = {
        'TRIMMER': [i for i in df['TRIMMER'].unique() if i],
        'QC_REPORTER': [i for i in df['QC_REPORTER'].unique() if i],
        'ALIGNER': [i for i in df['ALIGNER'].unique() if i],
        'DEDUPER': [i for i in df['DEDUPER'].unique() if i],
        'COUNTER': [i for i in df['COUNTER'].unique() if i],
        'RECALIBRATE': df['RECALIBRATE'].unique().tolist(),
        'STATS': df['STATS'].unique().tolist()}
    return tool_dict


def hook_test():
    print('hook found')


if __name__ == '__main__':
    arg_parser: ArgumentParser = ArgumentParser()
    arg_parser.add_argument('-i', '--input', type=str, required=True,
                            help='csv file to investigate')
    args: Namespace = arg_parser.parse_args()
    for i in sorted(generate_file_ls(Path(args.input))):
        print(i)
    print('\ntool list:')

    for k, v in generate_tool_ls(Path(args.input)).items():
        print(f'{k}: {v}')
