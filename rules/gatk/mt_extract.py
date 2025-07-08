from argparse import ArgumentParser, Namespace
from pathlib import Path
import polars as pl


# ! This python script will costs a lot of RAM (> 20GB maybe),
# ! so use the awk approach instead, unless you have enough RAM and wanna it quicker.

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

    m_origin: str
    u_origin: str
    if args.method == 'PS':
        u_origin = 'UNCONVERTED_BASE_COV'
        m_origin = 'CONVERTED_BASE_COV'
    else:
        m_origin = 'UNCONVERTED_BASE_COV'
        u_origin = 'CONVERTED_BASE_COV'

    (pl.scan_csv(input_vcf, has_header=False, comment_prefix='#', separator='\t',
                 schema={'chrom': pl.String, 'end': pl.Int64, 'id': pl.String,
                         'ref': pl.String, 'alt': pl.String,
                         'qual': pl.String, 'filter': pl.String, 'info': pl.String,
                         'format': pl.String, 'sample': pl.String})
       .with_columns((pl.col('end') - 1).alias('start'))
       .with_columns([pl.col('info')
                        .str
                        .extract(regex, group_index=1)
                        .cast(pl.Int64 if field != "REFERENCE_CONTEXT" else pl.String)
                        .alias(field)
                      for field, regex in {'CONVERTED_BASE_COV': r'CONVERTED_BASE_COV=(\d+)',
                                           'REFERENCE_CONTEXT': r'REFERENCE_CONTEXT=([A-Z]+)',
                                           'UNCONVERTED_BASE_COV': r'UNCONVERTED_BASE_COV=(\d+)',
                                           'DP': r'DP=(\d+)'}.items()])
       .with_columns([pl.col(u_origin).alias('u'),
                      pl.col(m_origin).alias('m')])
       .with_columns(pl.sum_horizontal('u', 'm').alias('depth'))
       .filter(pl.col('depth') >= 5)
       .select('chrom', 'start', 'end', 'depth', 'u', 'm')
       .sink_csv(output_bedgraph, include_header=True, separator='\t', line_terminator='\n'))


if __name__ == '__main__':
    main()
