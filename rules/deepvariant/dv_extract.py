from argparse import ArgumentParser, Namespace
from pathlib import Path
from warnings import filterwarnings
import polars as pl
from zstandard import ZstdCompressor

# ! This python script will costs a lot of RAM (> 50GB maybe),
# ! so use the awk approach instead, unless you have enough RAM and wanna it quicker.
filterwarnings("ignore", message="Polars found a filename")


def main():
    arg_parser: ArgumentParser = ArgumentParser()
    arg_parser.add_argument('-i', '--input', type=str, required=True,
                            help='Input DeepVariant VCF file')
    arg_parser.add_argument('-o', '--output', type=str, required=True,
                            help='Output Bedgraph file, compressed by zstd')
    arg_parser.add_argument('-m', '--method', type=str, required=True,
                            choices={'BS', 'EM', 'PS', 'RR'})
    arg_parser.add_argument('-p', '--parallel', type=int, default=1,
                            help='Number of parallel processes to use')

    args: Namespace = arg_parser.parse_args()
    input_vcf: Path = Path(args.input)

    if not input_vcf.exists():
        raise FileNotFoundError(f"Input file {input_vcf} does not exist.")

    output_bedgraph: Path = Path(args.output)

    if not output_bedgraph.parent.exists():
        raise FileNotFoundError(f"Output directory {output_bedgraph.parent} does not exist.")

    m_index: int
    u_index: int
    if args.method == 'PS':
        u_index = 0
        m_index = 1
    else:
        m_index = 0
        u_index = 1

    df: pl.DataFrame
    df = (pl.scan_csv(input_vcf, has_header=False, comment_prefix='#', separator='\t',
                      schema={'chrom': pl.String, 'end': pl.Int64, 'id': pl.String,
                              'ref': pl.String, 'alt': pl.String,
                              'qual': pl.Float64, 'filter': pl.String, 'info': pl.String,
                              'format': pl.String, 'sample': pl.String})
            .filter(
                ((
                    pl.col('ref') == 'C') & (pl.col('alt') == 'T')
                 ) | (
                    (pl.col('ref') == 'G') & (pl.col('alt') == 'A')
                ))
            .with_columns([(pl.col('end') - 1).alias('start'),
                           pl.col('sample').str.split(':').alias('fields')])
            .filter(pl.col('filter').is_in({'RefCall', 'PASS'}))
            .select('chrom', 'start', 'end',
                    pl.col('fields').list.get(2)
                      .cast(pl.Int64).alias('depth'),
                    pl.col('fields').list.get(3)
                      .str.split(',').list.get(m_index)
                      .cast(pl.Int32).alias('m'),
                    pl.col('fields').list.get(3)
                      .str.split(',').list.get(u_index)
                      .cast(pl.Int32).alias('u'))
            .collect())

    zstd_compressor = ZstdCompressor(level=19, threads=args.parallel)
    with open(output_bedgraph, mode='wb+') as f_out:
        with zstd_compressor.stream_writer(f_out) as _writer:
            df.write_csv(_writer, include_header=True, separator='\t', line_terminator='\n')


if __name__ == '__main__':
    main()
