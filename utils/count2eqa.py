from argparse import ArgumentParser, Namespace
from pathlib import Path
import polars as pl


def main():
    parser = ArgumentParser(description='Unify format of methylation results')
    parser.add_argument('-i', '--input', dest='input', type=str, required=True,
                        help='Path to the input methylation file.')
    parser.add_argument('-o', '--outdir', dest='outdir', type=str,
                        help=('Path to the folder containing output Parquet file. '
                              'Default to the same dir as input.'))
    parser.add_argument('-d', '--depth', dest='depth', type=int,
                        default=5, help='Minimum depth of every cytosine.')
    parser.add_argument('-c', '--counter', dest='counter', type=str,
                        choices=['bismark_c2c', 'bismark_cov',
                                 'biscuit', 'methyldackel',
                                 'dnmtools', 'bs_seeker2', 'astair',
                                 'bsgenova', 'fame',
                                 'msuite2-bowtie2', 'msuite2-hisat2'],
                        help='Counter used to generate the input file.')
    parser.add_argument('-f', '--force', dest='force', action='store_true',
                        help='Force overwrite the output file if it exists.')

    args: Namespace = parser.parse_args()

    depth: int = args.depth if args.depth > 0 else 5

    input_f: Path = Path(args.input)
    if not input_f.exists():
        raise FileNotFoundError(f'Input file {input_f} does not exist.')

    fname: str = input_f.name.split('.')[0]

    output_f: Path
    if args.outdir:
        output_f = Path(args.outdir) / f'{fname}.parquet.lz4'
    else:
        output_f = input_f.parent / f'{fname}.parquet.lz4'

    if output_f.exists() and not args.force:
        return None

    match args.counter:
        case 'astair':  # mod.gz
            (pl.scan_csv(input_f, separator='\t', has_header=True)
               .filter(pl.col('TOTAL_DEPTH') >= depth)
               .with_columns(pl.sum_horizontal('UNMOD', 'MOD').alias('depth'))
               .with_columns((100 * pl.col('MOD') / pl.col('depth')).alias('beta'))
               .select('#CHROM', 'START', 'END', 'UNMOD', 'MOD', 'beta', 'depth')
               .rename({'#CHROM': 'chrom', 'START': 'start', 'END': 'end',
                        'UNMOD': 'u', 'MOD': 'm'})
               .sink_parquet(output_f, compression='lz4', sync_on_close='all'))
        case 'rastair':  # mod.gz
            (pl.scan_csv(input_f, separator='\t', has_header=True)
               .with_columns(pl.sum_horizontal('unmod', 'mod').alias('depth'))
               .with_columns((100 * pl.col('mod') / pl.col('depth')).alias('beta'))
               .filter(pl.col('depth') >= depth)
               .select('#chr', 'start', 'end', 'unmod', 'mod', 'strand', 'beta', 'depth')
               .rename({'#chr': 'chrom', 'unmod': 'u', 'mod': 'm'})
               .sink_parquet(output_f, compression='lz4', sync_on_close='all'))
        case 'bismark_c2c':  # c2c.cov.gz
            (pl.scan_csv(input_f, separator='\t', has_header=False,
                         schema={'chrom': pl.String, 'start': pl.Int64, 'strand': pl.String,
                                 'm': pl.Int64, 'u': pl.Int64, 'd': pl.String, 't': pl.String,
                                 'q': pl.String, 'p': pl.String})
               .with_columns([(pl.col('start') + 1).alias('end'),
                              pl.sum_horizontal('m', 'u').alias('depth')])
               .filter(pl.sum_horizontal('depth') >= depth)
               .with_columns((100 * pl.col('m') / pl.col('depth')).alias('beta'))
               .select('chrom', 'start', 'end', 'u', 'm', 'strand', 'beta', 'depth')
               .sink_parquet(output_f, compression='lz4', sync_on_close='all'))
        case 'bismark_cov':  # bismark.zero.cov
            (pl.scan_csv(input_f, separator='\t', has_header=False,
                         schema={'chrom': pl.String, 'start': pl.Int64, 'end': pl.String,
                                 'l': pl.Float64, 'u': pl.Int64, 'm': pl.Int64})
               .with_columns([pl.sum_horizontal('m', 'u').alias('depth'),
                              pl.lit('+').alias('strand')])
               .filter(pl.sum_horizontal('depth') >= depth)
               .with_columns((100 * pl.col('m') / pl.col('depth')).alias('beta'))
               .select('chrom', 'start', 'end', 'u', 'm', 'strand', 'beta', 'depth')
               .sink_parquet(output_f, compression='lz4', sync_on_close='all'))
        case 'dnmtools':  # dnmtools.gz
            (pl.scan_csv(input_f, separator='\t', has_header=False,
                         schema={'chrom': pl.String, 'start': pl.Int64, 'strand': pl.String,
                                 'type': pl.String, 'level': pl.Float64, 'depth': pl.Int64})
               .with_columns((pl.col('level') * pl.col('depth'))
                             .round(0)
                             .cast(pl.Int64)
                             .alias('m'))
               .with_columns([(pl.col('depth') - pl.col('m')).alias('u'),
                              (100 * pl.col('m') / pl.col('depth')).alias('beta'),
                              (pl.col('start') + 1).alias('end')])
               .select('chrom', 'start', 'end', 'u', 'm', 'strand', 'beta', 'depth')
               .sink_parquet(output_f, compression='lz4', sync_on_close='all'))
        case 'methyldackel':  # bedgraph.gz
            (pl.scan_csv(input_f, separator='\t', has_header=False,
                         schema={'chrom': pl.String, 'start': pl.Int64, 'end': pl.Int64,
                                 'level': pl.Float64, 'm': pl.Int64, 'u': pl.Int64})
               .with_columns(pl.sum_horizontal('m', 'u').alias('depth'))
               .with_columns((100 * pl.col('m') / pl.col('depth')).alias('beta'))
               .filter(pl.col('depth') >= depth)
               .select('chrom', 'start', 'end', 'u', 'm', 'beta', 'depth')
               .sink_parquet(output_f, compression='lz4', sync_on_close='all'))
        case _:
            raise NotImplementedError(f'Counter {args.counter} is not implemented yet.')


if __name__ == '__main__':
    main()
