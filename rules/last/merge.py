from argparse import ArgumentParser
from pathlib import Path
import polars as pl


def main():
    parser = ArgumentParser(description=('Merge and unique BISF, BISR and '
                                         'NEAR sam files to output.sam file.'))
    parser.add_argument('--bisf', help='Path to the BISF sam file.')
    parser.add_argument('--bisr', help='Path to the BISR sam file.')
    parser.add_argument('--near', help='Path to the NEAR sam file.')
    parser.add_argument('--sam_tmp', help='Path to the tmp sam file.')
    args = parser.parse_args()

    for _f in (args.bisf, args.bisr, args.near):
        if not Path(_f).exists():
            raise FileNotFoundError(f'File {_f} does not exist.')

    Path(args.sam_tmp).parent.mkdir(exist_ok=True)

    # Read the CSV file
    reads_df: pl.LazyFrame = (
        pl.concat(items=[
            pl.scan_csv(args.bisf, separator='\t', has_header=False,
                        schema={'QNAME': pl.String, 'FLAG': pl.UInt16, 'RNAME': pl.String,
                                'POS': pl.Int64, 'MAPQ': pl.Int64, 'CIGAR': pl.String,
                                'RNEXT': pl.String, 'PNEXT': pl.Int64, 'TLEN': pl.Int64,
                                'SEQ': pl.String, 'QUAL': pl.String}),
            pl.scan_csv(args.bisr, separator='\t', has_header=False,
                        schema={'QNAME': pl.String, 'FLAG': pl.UInt16, 'RNAME': pl.String,
                                'POS': pl.Int64, 'MAPQ': pl.Int64, 'CIGAR': pl.String,
                                'RNEXT': pl.String, 'PNEXT': pl.Int64, 'TLEN': pl.Int64,
                                'SEQ': pl.String, 'QUAL': pl.String}),
            pl.scan_csv(args.near, separator='\t', has_header=False,
                        schema={'QNAME': pl.String, 'FLAG': pl.UInt16, 'RNAME': pl.String,
                                'POS': pl.Int64, 'MAPQ': pl.Int64, 'CIGAR': pl.String,
                                'RNEXT': pl.String, 'PNEXT': pl.Int64, 'TLEN': pl.Int64,
                                'SEQ': pl.String, 'QUAL': pl.String})],
                  how='vertical')
          .sort(by=['QNAME', 'FLAG', 'MAPQ'], descending=[False, False, True])
          .unique(subset=['QNAME', 'FLAG'], keep='first'))

    reads_df.sink_csv(path=args.sam_tmp, separator='\t', include_header=False)


if __name__ == '__main__':
    main()
