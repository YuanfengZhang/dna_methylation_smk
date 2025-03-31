from argparse import ArgumentParser
import warnings
import polars as pl
from zstandard import ZstdCompressor


warnings.filterwarnings("ignore", message="Polars found a filename")


def main():
    parser = ArgumentParser(description='Format the FAME output')
    parser.add_argument('-i', '--input', help='The input FAME tsv')
    parser.add_argument('-o', '--output', help='The output FAME bedgraph file')
    parser.add_argument('-t', '--threads', help='The number of threads to use',
                        default=1, type=int)
    args = parser.parse_args()
    df: pl.DataFrame = (pl.scan_csv(args.input, separator='\t', has_header=False,
                                    schema={'chrom': pl.String, 'start': pl.Int64,
                                            'm+': pl.Int64, 'u+': pl.Int64, 'm-': pl.Int64, 'u-': pl.Int64})
                        .filter(pl.sum_horizontal('m+', 'u+', 'm-', 'u-') > 0)
                        .with_columns([(pl.col('start') + 2).alias('end'),
                                       (100 * pl.sum_horizontal('m+',
                                                                'm-') / pl.sum_horizontal('m+', 'u+',
                                                                                          'm-', 'u-'))
                                       .alias('beta')])
                        .select(['chrom', 'start', 'end', 'beta', 'm+', 'u+', 'm-', 'u-'])
                        .collect())

    zstd_compressor = ZstdCompressor(level=19, threads=args.threads)
    with zstd_compressor.stream_writer(writer=open(file=args.output,
                                                   mode='wb+')) as _writer:
        df.write_csv(_writer, separator='\t')


if __name__ == '__main__':
    main()
