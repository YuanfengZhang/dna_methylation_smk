from argparse import ArgumentParser
from pathlib import Path
from tempfile import NamedTemporaryFile
import warnings
import polars as pl


warnings.filterwarnings("ignore", message="Polars found a filename")


def main():
    parser = ArgumentParser(description='Format the FAME output')
    parser.add_argument('-i', '--input', help='The input FAME tsv')
    parser.add_argument('-o', '--output', help='The output FAME bedgraph file')
    parser.add_argument('-t', '--tmp-dir', dest='tmp_dir', help='temporary directory for tmp bed files')
    parser.add_argument('-c', '--cores', help='The number of threads to use',
                        default=1, type=int)
    args = parser.parse_args()

    tmp_dir: Path = Path(args.tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    df: pl.DataFrame = (pl.scan_csv(args.input, separator='\t', has_header=False,
                                    schema={'chrom': pl.String, 'start': pl.Int64,
                                            'm+': pl.Int64, 'u+': pl.Int64, 'm-': pl.Int64, 'u-': pl.Int64})
                        .filter(pl.sum_horizontal('m+', 'u+', 'm-', 'u-') > 0)
                        .with_columns([(pl.col('start') + 2).alias('end')])
                        .select(['chrom', 'start', 'end', 'm+', 'u+', 'm-', 'u-'])
                        .collect())
    to_slice: pl.DataFrame = df.select(['chrom', 'start', 'end'])
    tmp_file_list: list[str] = []

    for sliced_df in to_slice.iter_slices(n_rows=1000000):
        with NamedTemporaryFile(delete=False, delete_on_close=False,
                                dir=tmp_dir, suffix='.bed') as tmp_file:
            sliced_df.write_csv(tmp_file.name, include_header=False,
                                separator='\t', line_terminator='\n')
            tmp_file_list.append(tmp_file.name)
            tmp_file.close()


if __name__ == '__main__':
    # main()
    print('this python script has been deprecated')
