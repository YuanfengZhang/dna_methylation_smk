from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
from typing import List


def process_fame_row(row: List):
    if (not row) or any(col == '0' for col in row[2:6]):
        return None
    try:
        new_row = (
            row[:2] + [           # 0-based start
                int(row[1]) + 2,  # 0-based end
                100 * (int(row[2]) + int(row[4])) / (sum(int(i) for i in row[2: 6]))  # beta value
            ] + row[2:])          # m+, u+, m-, u-
        return new_row
    except ValueError:
        return None


def main():
    parser = ArgumentParser(description='Format the FAME output')
    parser.add_argument('-i', '--input', help='The input FAME tsv')
    parser.add_argument('-o', '--output', help='The output FAME bedgraph file')
    parser.add_argument('-t', '--threads', help='The number of threads to use',
                        default=1, type=int)
    args = parser.parse_args()

    rows = [['chrom', 'start', 'end', 'beta', 'm+', 'u+', 'm-', 'u-']]

    with open(args.input, mode='r', newline='') as fame_f:
        reader = csv.reader(fame_f, delimiter='\t')
        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            rows += list(executor.map(process_fame_row, reader))

    with gzip.open(args.output, mode='wt', compresslevel=9) as bedgraph_f:
        writer = csv.writer(bedgraph_f, delimiter='\t')  # type: ignore
        writer.writerows(filter(None, rows))


if __name__ == '__main__':
    main()
