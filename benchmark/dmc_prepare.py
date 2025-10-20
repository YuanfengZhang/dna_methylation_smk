from argparse import ArgumentParser, Namespace
import logging
from pathlib import Path
import polars as pl


def get_ref_locs(ref_dir: Path, label1: str, label2: str, chrom: str) -> set[int]:
    logging.info(f'reading ref parquets of {label1} and {label2}')
    try:
        ref1_f: Path = list(ref_dir.glob(f'{label1}*parquet*'))[0]
        ref2_f: Path = list(ref_dir.glob(f'{label2}*parquet*'))[0]
    except IndexError:
        raise FileNotFoundError(f'No parquet file found for {label1} and {label2} in {ref_dir}.')

    (ref1_locs,
     ref2_locs) = [set(pl.read_parquet(_f).filter(pl.col('chrom') == chrom,
                                                  pl.col('in_hcr'))['start'].to_list())
                   for _f in (ref1_f, ref2_f)]

    return ref1_locs.intersection(ref2_locs)


def read_for_methylkit(input_dir: Path, beta_col: str,
                       lab: str, label: str, rep: str,
                       cutoff: int, chrom: str = 'chr17') -> pl.LazyFrame:
    try:
        sample_f: Path = list(input_dir.glob(f'{lab}_{label}_{rep}*parquet*'))[0]
    except IndexError:
        raise FileNotFoundError(f'No parquet file found for {lab}_{label}_{rep}.')

    logging.debug(f'reading {lab}_{label}_{rep}: {sample_f}')

    return (pl.scan_parquet(sample_f)
              .filter(pl.col('depth') >= cutoff,
                      pl.col('chrom') == chrom)
              .unique(subset=['start'])
              .rename({'m': f'numCs_{label}_{rep}',
                       'u': f'numTs_{label}_{rep}',
                       beta_col: f'freqC_{label}_{rep}',
                       'depth': f'coverage_{label}_{rep}'})
              .with_columns((100 - pl.col(f'freqC_{label}_{rep}')).alias(f'freqT_{label}_{rep}'))
              .select('start', f'coverage_{label}_{rep}',
                      f'freqC_{label}_{rep}', f'freqT_{label}_{rep}'))


def read_for_methylsig(input_dir: Path,
                       lab: str, label: str, rep: str,
                       cutoff: int, chrom: str = 'chr17') -> pl.LazyFrame:
    try:
        sample_f: Path = list(input_dir.glob(f'{lab}_{label}_{rep}*parquet*'))[0]
    except IndexError:
        raise FileNotFoundError(f'No parquet file found for {lab}_{label}_{rep}.')

    logging.debug(f'reading {lab}_{label}_{rep}: {sample_f}')

    return (pl.scan_parquet(sample_f)
              .filter(pl.col('depth') >= cutoff,
                      pl.col('chrom') == chrom)
              .unique(subset=['start'])
              .rename({'m': f'M_{label}_{rep}',
                       'depth': f'Cov_{label}_{rep}'})
              .select('start', f'M_{label}_{rep}', f'Cov_{label}_{rep}'))


def read_for_cpgtools(input_dir: Path,
                      lab: str, label: str, rep: str,
                      cutoff: int, chrom: str = 'chr17') -> pl.LazyFrame:
    try:
        sample_f: Path = list(input_dir.glob(f'{lab}_{label}_{rep}*parquet*'))[0]
    except IndexError:
        raise FileNotFoundError(f'No parquet file found for {lab}_{label}_{rep}.')

    logging.debug(f'reading {lab}_{label}_{rep}: {sample_f}')

    return (pl.scan_parquet(sample_f)
              .filter(pl.col('depth') >= cutoff,
                      pl.col('chrom') == chrom)
              .unique(subset=['start'])
              .with_columns(pl.concat_str('m', pl.lit(','), 'depth').alias(f'{label}_{rep}'))
              .select('start', f'{label}_{rep}'))


def prepare_methylkit(input_dir: Path, output_dir: Path, lab: str, labels: list[str],
                      cutoff: int, chrom: str, ref_locs: set[int], beta_col: str) -> None:
    (l1r1,
     l1r2,
     l2r1,
     l2r2) = [read_for_methylkit(input_dir=input_dir, lab=lab, label=label, rep=rep,
                                 cutoff=cutoff, chrom=chrom, beta_col=beta_col)
              for label in labels for rep in ('1', '2')]

    logging.debug(f'{lab} {labels} parquets merging')

    merged: pl.DataFrame
    if ref_locs:
        merged = (l1r1.join(l1r2, on=['start'], how='inner')
                      .join(l2r1, on=['start'], how='inner')
                      .join(l2r2, on=['start'], how='inner')
                      .filter(pl.col('start').is_in(ref_locs))
                      .with_columns([pl.lit(chrom).alias('chrom'),
                                     pl.lit('+').alias('strand')])
                      .with_columns(pl.concat_str('chrom', pl.lit('.'), 'start').alias('chrBase'))
                      .collect())
    else:
        merged = (l1r1.join(l1r2, on=['start'], how='inner')
                      .join(l2r1, on=['start'], how='inner')
                      .join(l2r2, on=['start'], how='inner')
                      .with_columns([pl.lit(chrom).alias('chrom'),
                                     pl.lit('+').alias('strand')])
                      .with_columns(pl.concat_str('chrom', pl.lit('.'), 'start').alias('chrBase'))
                      .collect())
    logging.debug(f'{lab} {labels} parquets merged')

    if merged.shape[0] == 0:
        raise ValueError(f'No overlapping locations in HCR found for {lab} and {labels}.')

    for label in labels:
        for rep in ('1', '2'):
            (merged.select('chrBase', 'chrom', 'start', 'strand',
                           f'coverage_{label}_{rep}', f'freqC_{label}_{rep}', f'freqT_{label}_{rep}')
                   .sort('chrom', 'start')
                   .rename({'chrom': 'chr', 'start': 'base',
                            f'coverage_{label}_{rep}': 'coverage',
                            f'freqC_{label}_{rep}': 'freqC',
                            f'freqT_{label}_{rep}': 'freqT'})
                   .write_csv(output_dir / f'{lab}_{labels[0]}_vs_{labels[1]}_{label}_{rep}.bed',
                              separator='\t', include_header=True))
            (output_dir / f'{lab}_{labels[0]}_vs_{labels[1]}_{label}_{rep}.bed').chmod(0o777)
            logging.info(f'methylkit input file {lab}_{label}_{rep}.bed at {output_dir.as_posix()} created')
    logging.info(f'{lab} {labels} methylkit finished\n------------------------------------------')


def prepare_methylsig(input_dir: Path, output_dir: Path, lab: str, labels: list[str],
                      cutoff: int, chrom: str, ref_locs: set[int]) -> None:
    (l1r1,
     l1r2,
     l2r1,
     l2r2) = [read_for_methylsig(input_dir=input_dir, lab=lab, label=label, rep=rep,
                                 cutoff=cutoff, chrom=chrom)
              for label in labels for rep in ('1', '2')]

    logging.debug(f'{lab} {labels} parquets merging')

    merged: pl.DataFrame

    if ref_locs:
        merged = (l1r1.join(l1r2, on=['start'], how='inner')
                      .join(l2r1, on=['start'], how='inner')
                      .join(l2r2, on=['start'], how='inner')
                      .filter(pl.col('start').is_in(ref_locs))
                      .with_columns([pl.lit(chrom).alias('chrom'),
                                     pl.lit('+').alias('strand')])
                      .collect())
    else:
        merged = (l1r1.join(l1r2, on=['start'], how='inner')
                      .join(l2r1, on=['start'], how='inner')
                      .join(l2r2, on=['start'], how='inner')
                      .with_columns([pl.lit(chrom).alias('chrom'),
                                     pl.lit('+').alias('strand')])
                      .collect())
    logging.debug(f'{lab} {labels} parquets merged')

    if merged.shape[0] == 0:
        raise ValueError(f'No overlapping locations in HCR found for {lab} and {labels}.')

    for label in labels:
        for rep in ('1', '2'):
            (merged.rename({f'M_{label}_{rep}': 'M', f'Cov_{label}_{rep}': 'Cov'})
                   .select('chrom', 'start', 'strand', 'M', 'Cov')
                   .sort('chrom', 'start')
                   .write_csv(output_dir / f'{lab}_{labels[0]}_vs_{labels[1]}_{label}_{rep}.bed',
                              separator='\t', include_header=True))
            (output_dir / f'{lab}_{labels[0]}_vs_{labels[1]}_{label}_{rep}.bed').chmod(0o777)
            logging.info(f'methylsig input file {lab}_{label}_{rep}.bed at {output_dir.as_posix()} created')
    logging.info(f'{lab} {labels} methylsig finished\n------------------------------------------')


def prepare_cpgtools(input_dir: Path, output_dir: Path, lab: str, labels: list[str],
                     cutoff: int, chrom: str, ref_locs: set[int]) -> None:
    (l1r1,
     l1r2,
     l2r1,
     l2r2) = [read_for_cpgtools(input_dir=input_dir, lab=lab, label=label, rep=rep,
                                cutoff=cutoff, chrom=chrom)
              for label in labels for rep in ('1', '2')]

    logging.debug(f'{lab} {labels} parquets merging')

    merged: pl.DataFrame

    if ref_locs:
        merged = (l1r1.join(l1r2, on=['start'], how='inner')
                      .join(l2r1, on=['start'], how='inner')
                      .join(l2r2, on=['start'], how='inner')
                      .filter(pl.col('start').is_in(ref_locs))
                      .with_columns(pl.concat_str(pl.lit(chrom),
                                                  pl.lit('_'),
                                                  'start').alias('cgID'))
                      .drop('start')
                      .collect())
    else:
        merged = (l1r1.join(l1r2, on=['start'], how='inner')
                      .join(l2r1, on=['start'], how='inner')
                      .join(l2r2, on=['start'], how='inner')
                      .with_columns(pl.concat_str(pl.lit(chrom),
                                                  pl.lit('_'),
                                                  'start').alias('cgID'))
                      .drop('start')
                      .collect())
    logging.debug(f'{lab} {labels} parquets merged')

    if merged.shape[0] == 0:
        raise ValueError(f'No overlapping locations in HCR found for {lab} and {labels}.')

    (merged.select('cgID', f'{labels[0]}_1', f'{labels[0]}_2', f'{labels[1]}_1', f'{labels[1]}_2')
           .write_csv(output_dir / f'{lab}_{labels[0]}_vs_{labels[1]}.bed',
                      separator='\t', include_header=True))
    (pl.DataFrame({'sampleID': [f'{labels[0]}_1', f'{labels[0]}_2',
                                f'{labels[1]}_1', f'{labels[1]}_2'],
                   'group': [labels[0], labels[0], labels[1], labels[1]]})
       .write_csv(output_dir / f'{lab}_{labels[0]}_vs_{labels[1]}.group',
                  include_header=True))
    (output_dir / f'{lab}_{labels[0]}_vs_{labels[1]}.bed').chmod(0o777)
    (output_dir / f'{lab}_{labels[0]}_vs_{labels[1]}.group').chmod(0o777)
    logging.info(f'cpgtools input file {lab}_{labels[0]}_vs_{labels[1]}.bed '
                 f'and group csv at {output_dir.as_posix()} created')
    logging.info(f'{lab} {labels} cpgtools finished\n------------------------------------------')


def main():
    arg_parser: ArgumentParser = ArgumentParser(usage='This script is made for preparing '
                                                      'the tsv files for DMC analysis')
    arg_parser.add_argument('-i', '--input-dir', dest='input_dir', type=str, required=True,
                            help='Input directory containing formatted parquet files.')
    arg_parser.add_argument('-r', '--ref-dir', dest='ref_dir', type=str, required=True,
                            help='Directory containing reference parquet files.')
    arg_parser.add_argument('--ref-only', dest='ref_only', action='store_true',
                            help='Whether to only generate input within the truset locs.')
    arg_parser.add_argument('-o', '--output-dir', dest='output_dir', type=str, required=True,
                            help='Output directory to save the prepared tsv files.')
    arg_parser.add_argument('--beta-col', dest='beta_col', type=str, default='beta',
                            help='Name of the beta column in the input lz4 files.')
    arg_parser.add_argument('-l', '--labs', dest='labs', type=str, required=True,
                            help='comma-separated labs to calculate for.')
    arg_parser.add_argument('-s', '--labels', dest='labels', type=str, required=True,
                            help='two comma-separated labels to calculate.')
    arg_parser.add_argument('-m', '--method', dest='method', type=str, required=True,
                            choices=('methylkit', 'methylsig', 'cpgtools'),
                            help='destination.')
    arg_parser.add_argument('-d', '--depth', dest='depth', type=int, default=5,
                            help='minimal depth to include.')
    arg_parser.add_argument('-c', '--chrom', dest='chrom', type=str, default='chr17',
                            choices=('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
                                     'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
                                     'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                                     'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM',
                                     'lambda', 'pUC19'),
                            help='one chromosome to include.')
    arg_parser.add_argument('-v', '--log-level', dest='log_level',
                            choices=['debug', 'info', 'warning', 'error', 'critical'],
                            default='info')
    arg_parser.add_argument('-a', '--action', dest='action', type=str, required=True,
                            choices=['prepare', 'remove'],
                            help='prepare or remove.')

    args: Namespace = arg_parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    input_dir: Path = Path(args.input_dir)
    ref_dir: Path = Path(args.ref_dir)
    output_dir: Path = Path(args.output_dir)
    beta_col: str = args.beta_col
    labs: list[str] = args.labs.split(',')
    method: str = args.method
    cutoff: int = args.depth
    chrom: str = args.chrom
    action: str = args.action
    ref_only: bool = args.ref_only

    if not input_dir.is_dir():
        raise FileNotFoundError(f'Input directory {input_dir} is not a valid directory.')

    output_dir.mkdir(parents=True, exist_ok=True)

    if args.labels.count(',') != 1:
        raise ValueError('Please provide two comma-separated labels.')

    order_mapping: dict[str, int] = {'D6': 0, 'D5': 1, 'T1': 2, 'T2': 3, 'T3': 4, 'T4': 5,
                                     'F7': 6, 'M8': 7, 'BC': 8, 'BL': 9, 'HF': 10}

    labels: list[str] = sorted(args.labels.split(','),
                               key=lambda i: order_mapping.get(i, 11))

    if action == 'remove':
        match method:
            case 'methylkit' | 'methylsig':
                for lab in labs:
                    logging.info(f'removing methylsig input files for {lab} {labels}')
                    for _l in labels:
                        for _r in ('1', '2'):
                            (output_dir / f'{lab}_{labels[0]}_vs_{labels[1]}_{_l}_{_r}.bed').unlink(missing_ok=True)
            case 'cpgtools':
                for lab in labs:
                    logging.info(f'removing cpgtools input files for {lab} {labels}')
                    (output_dir / f'{lab}_{labels[0]}_vs_{labels[1]}.bed').unlink(missing_ok=True)
            case _:
                raise ValueError(f'Unknown  dmc tool: {method}.')
    else:
        ref_locs: set[int] = set()
        if ref_only:
            ref_locs: set[int] = get_ref_locs(ref_dir=ref_dir, label1=labels[0], label2=labels[1], chrom=chrom)

        match method:
            case 'methylkit':
                logging.info('methylkit chosen')
                for lab in labs:
                    logging.info(f'preparing methlykit input files for {lab} {labels}')
                    prepare_methylkit(input_dir=input_dir, output_dir=output_dir,
                                      lab=lab, labels=labels, cutoff=cutoff, chrom=chrom,
                                      ref_locs=ref_locs, beta_col=beta_col)
            case 'methylsig':
                logging.info('methylsig chosen')
                for lab in labs:
                    logging.info(f'preparing methlykit input files for {lab} {labels}')
                    prepare_methylsig(input_dir=input_dir, output_dir=output_dir,
                                      lab=lab, labels=labels, cutoff=cutoff, chrom=chrom,
                                      ref_locs=ref_locs)

            case 'cpgtools':
                logging.info('cpgtools chosen')
                for lab in labs:
                    logging.info(f'preparing cpgtools input files for {lab} {labels}')
                    prepare_cpgtools(input_dir=input_dir, output_dir=output_dir,
                                     lab=lab, labels=labels, cutoff=cutoff, chrom=chrom,
                                     ref_locs=ref_locs)
            case _:
                raise ValueError(f'Unknown dmc tool: {method}.')


if __name__ == '__main__':
    main()
