from argparse import ArgumentParser, Namespace
from pathlib import Path
import polars as pl
from tqdm import tqdm


def conversion_rate_per_file(methylation_file: Path,
                             chromosomes: list[str],
                             depth_cutoff: int) -> dict[str, str | float]:
    f_name = methylation_file.stem.split('.')[0]
    lab: str = f_name.split('_')[0]
    sample: str = '_'.join(f_name.split('_')[1:])

    to_return: dict[str, str | float] = {'lab': lab, 'sample': sample}

    df: pl.DataFrame
    if 'parquet' in methylation_file.name:
        df = pl.read_parquet(methylation_file)
    elif 'csv' in methylation_file.name:
        df = pl.read_csv(methylation_file)
    else:
        df = pl.read_csv(methylation_file, separator='\t')

    available_chromosomes: list[str] = df['chrom'].unique().to_list()

    for chromosome in chromosomes:
        if chromosome in available_chromosomes:
            to_return[f'{chromosome}_avg_beta'] = (df.lazy()
                                                     .filter(pl.col('chrom') == chromosome,
                                                             pl.col('depth') >= depth_cutoff)
                                                     .select(pl.col('beta').mean())
                                                     .collect()
                                                     .item())
        else:
            to_return[f'{chromosome}_avg_beta'] = -1
    return to_return


def main():
    parser = ArgumentParser(description='calculate the conversion rate of BS/EM/PS methods')
    parser.add_argument('-i', '--input-dir', dest='input_dir', required=True, type=str,
                        help='Path to the folder containing formatted methylation files')
    parser.add_argument('-ro', '--raw-output', dest='raw_output', required=True, type=str,
                        help='Path to the output file containing conversion rates')
    parser.add_argument('-so', '--stat-output', dest='stat_output', required=True, type=str,
                        help='Path to the output file containing statistics of conversion rates')
    parser.add_argument('-s', '--suffix', dest='input_suffix', required=True, type=str,
                        help='suffix of the input files')
    parser.add_argument('-l', '--lab', dest='lab', type=str,
                        help='lab of the input files. If not provided, all files will be used')
    parser.add_argument('-c', '--chromosomes', dest='chromosomes', type=str,
                        default='chrM,lambda,pUC19', help='Chromosomes to use')
    parser.add_argument('-d', '--depth-cutoff', dest='depth_cutoff', type=int,
                        default=1, help='the lowest depth to include')

    args: Namespace = parser.parse_args()
    input_dir: Path = Path(args.input_dir)
    raw_output: Path = Path(args.raw_output)
    stat_output: Path = Path(args.stat_output)
    chromosomes: list[str] = args.chromosomes.split(',')
    depth_cutoff: int = args.depth_cutoff if args.depth_cutoff > 0 else 1

    if not (input_dir.exists() and input_dir.is_dir()):
        raise FileNotFoundError(f'Input directory {input_dir} does not exist or is not a directory.')

    methylation_files: list[Path] = [i for i in input_dir.glob(f'*{args.input_suffix}')]
    if args.lab:
        methylation_files: list[Path] = [i for i in methylation_files if i.stem.startswith(args.lab)]

    if not methylation_files:
        raise FileNotFoundError(f'No files with suffix {args.input_suffix} found in {input_dir}.')

    results: list[dict[str, str | float]] = [conversion_rate_per_file(methylation_file=m_f,
                                                                      chromosomes=chromosomes,
                                                                      depth_cutoff=depth_cutoff)
                                             for m_f in tqdm(methylation_files)]

    raw_cr: pl.DataFrame
    raw_cr = (pl.from_dicts(results)
                .with_columns([pl.when(pl.col('lab').is_in(['PS1', 'PS2', 'PS3']))
                                 .then(pl.col(f'{chromosome}_avg_beta'))
                                 .otherwise(100 - pl.col(f'{chromosome}_avg_beta'))
                                 .alias(chromosome)
                               for chromosome in chromosomes])
                .fill_null(value=-1.0)
                .with_columns(pl.when(pl.col('lab').str.contains('PS'))
                                .then(pl.col('lambda'))
                                .otherwise(100 - pl.col('lambda'))
                                .alias('lambda_gap'),
                              pl.when(pl.col('lab').str.contains('PS'))
                                .then(100 - pl.col('pUC19'))
                                .otherwise(pl.col('pUC19'))
                                .alias('pUC19_gap'))
                .with_columns(pl.sum_horizontal('lambda_gap',
                                                'pUC19_gap')
                                .alias('conversion_bias')))
    cr_stat: pl.DataFrame
    cr_stat = (raw_cr.filter(~pl.col('sample').str.contains('HF'))
                     .group_by('lab')
                     .agg(pl.col('chrM').min().alias('chrM_min'),
                          pl.col('chrM').max().alias('chrM_max'),
                          pl.col('chrM').mean().alias('chrM_mean'),
                          pl.col('chrM').median().alias('chrM_median'),
                          pl.col('lambda').min().alias('lambda_min'),
                          pl.col('lambda').max().alias('lambda_max'),
                          pl.col('lambda').mean().alias('lambda_mean'),
                          pl.col('lambda').median().alias('lambda_median'),
                          pl.col('pUC19').min().alias('pUC19_min'),
                          pl.col('pUC19').max().alias('pUC19_max'),
                          pl.col('pUC19').mean().alias('pUC19_mean'),
                          pl.col('pUC19').median().alias('pUC19_median'))
                     .sort(by=pl.col('lab')
                                .cast(pl.Enum(['BS1', 'BS2', 'BS3', 'BS4',
                                               'EM1', 'EM2', 'EM3', 'EM4',
                                               'RR1', 'PS1', 'PS2', 'PS3',
                                               'MA1', 'MA2', 'MA3']))))
    raw_cr.write_csv(raw_output, include_header=True)
    cr_stat.write_csv(stat_output, include_header=True)


if __name__ == '__main__':
    main()
