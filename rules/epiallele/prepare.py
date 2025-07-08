from argparse import ArgumentParser, Namespace
from pathlib import Path
from subprocess import run


def main():
    arg_parser: ArgumentParser = ArgumentParser(
        description='Prepare the input bam file for epiallele.R')
    arg_parser.add_argument(
        '-i', '--input-dir', dest='input_dir',
        type=str, required=True,
        help='Input directory containing the bam files')
    arg_parser.add_argument(
        '-b', '--base-name', dest='base_name',
        type=str, required=True,
        help='Basename of the bam files')
    arg_parser.add_argument(
        '-r', '--read', dest='read_string',
        type=str, required=True,
        help='The first read of original bam file')
    arg_parser.add_argument(
        '-t', '--threads', dest='threads',
        type=int, default=1,
        help='Number of threads to use for samtools')

    args: Namespace = arg_parser.parse_args()
    input_dir: Path = Path(args.input_dir)
    base_name: str = args.base_name
    read_string: str = args.read_string
    threads: int = args.threads

    if not input_dir.is_dir():
        raise ValueError(f'Input directory {input_dir} does not exist or is not a directory')

    output_dir: Path = input_dir / 'epiallele'

    original_bam: Path = input_dir / f'{base_name}.bam'
    copied_bam: Path = output_dir / f'{base_name}.bam'
    preprocessed_bam: Path = output_dir / f'{base_name}.preprocessed.bam'
    sortn_bam: Path = output_dir / f'{base_name}.sortn.bam'

    if sortn_bam.exists():
        print(f'Sortn bam file {sortn_bam.name} already exists, skipping sortn step')
        return

    if 'XM' in read_string:
        print('XM tag found in read string, prepare the sortn bam file now')
        run(f'samtools sort -n -@ {threads} -o {sortn_bam} {original_bam}',
            shell=True)
        run(f'touch {copied_bam} {preprocessed_bam}', shell=True)

    elif any(tag in read_string for tag in ['XG', 'YD', 'ZS']):
        print('XG, XH, or XP tag found in read string, prepare the copied bam file now')
        if copied_bam.exists():
            print(f'Copied bam file {copied_bam.name} already exists, skipping copy step')
            return
        else:
            run(f'cp {original_bam} {copied_bam}', shell=True)

    else:
        raise ValueError(
            'Read string must contain XM, XG, YD, or ZS tag. '
            f'Provided read string: {read_string}'
        )


if __name__ == '__main__':
    main()
