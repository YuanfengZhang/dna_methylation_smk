from argparse import ArgumentParser, Namespace
import logging
from pathlib import Path
from pickle import dump as pickle_dump
from pickle import load as pickle_load
from matplotlib.font_manager import fontManager, FontProperties
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
from scipy.stats import entropy, gaussian_kde, skew
from scipy.signal import find_peaks
import seaborn as sns


logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)


def draw_distribution(df: pd.DataFrame,
                      info_col: str,
                      fname: str,
                      output_img_path: Path,  # be aware that no suffix is added
                      img_format: list[str],
                      font_path: Path):
    info_mapping: dict[str, str] = {
        'beta': 'Beta Value (%)',
        'beta_pyro': 'Beta Value (%)',
        'depth': 'Sequencing Depth (x)'
    }

    info_str: str = info_mapping.get(info_col, info_col)

    # set font
    fontManager.addfont(font_path)
    prop = FontProperties(fname=font_path)
    sns.set_theme(style='whitegrid', font=prop.get_name())

    # init figure and gridspec
    fig: Figure = plt.figure(figsize=(15, 8))
    gs: GridSpec = fig.add_gridspec(2, 3)

    fig.suptitle(f'{info_str} Distribution',
                 fontsize=16, fontweight='bold')

    # draw chr1 ... chr22, chrX, chrY
    chromosome_df = df[(~df.chrom.isin(['chrM', 'pUC19', 'lambda'])) & (df.depth >= 10)].copy()
    ax0: Axes = fig.add_subplot(gs[0, :])
    if 'beta' in info_col:
        sns.kdeplot(data=chromosome_df,
                    x=info_col, ax=ax0, fill=True, cut=0,
                    color='#68b65c', alpha=.8, bw_adjust=1.0)
        ax0.xaxis.set_ticks(range(0, 105, 5))
    else:
        sns.kdeplot(data=chromosome_df,
                    x=info_col, ax=ax0, fill=True, cut=0,
                    color='#68b65c', alpha=.8)
    ax0.set_title('')
    ax0.set_yticks([])
    ax0.set_xlabel('Chromosomes 1-22 + XY')

    # add description for ax0
    desc = chromosome_df[info_col].describe()
    describe_str: str = (f'count: {round(desc["count"])}   mean: {desc["mean"]:.4f}\n'
                         f'std: {desc["std"]:.4f}    min: {desc["min"]:.4f}\n'
                         f'25%: {desc["25%"]:.4f}    50%: {desc["50%"]:.4f}\n'
                         f'75%: {desc["75%"]:.4f}    max: {desc["max"]:.4f}')
    ax0.text(x=0.4, y=0.95,
             s=(f'{info_str} Statistics:\n' + describe_str),
             transform=ax0.transAxes, fontsize=14, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # draw chrM, lambda, pUC19
    control_contigs: list[tuple[str, str]] = [('chrM', '#d57ca3'),
                                              ('lambda', '#776c96'),
                                              ('pUC19', '#d96436')]

    for index, control_contig in enumerate(control_contigs):
        ax: Axes = fig.add_subplot(gs[1, index])
        if control_contig[0] in df.chrom.unique() and all(n not in fname for n in ('RR1', 'RM1')):
            sns.kdeplot(data=df[(df['chrom'] == control_contig[0]) & df['depth'] >= 10], x=info_col, ax=ax,
                        fill=True, cut=0, color=control_contig[1], alpha=.8)
            ax.set_title('')
            ax.set_yticks([])
            ax.set_xlabel(control_contig[0])
        else:
            ax.text(x=0.5, y=0.5, s=f'{control_contig[0]} not found in {fname}',
                    ha='center', va='center', fontsize=16, color='red')

    # save figure
    output_img_path.parent.mkdir(parents=True, exist_ok=True)
    for fmt in img_format:
        if fmt == 'tiff':
            fig.savefig(output_img_path.with_suffix(f'.{fmt}'),
                        format=fmt, dpi=300, bbox_inches='tight',
                        pil_kwargs={'compression': 'tiff_lzw'})
        else:
            fig.savefig(output_img_path.with_suffix(f'.{fmt}'),
                        format=fmt, dpi=300, bbox_inches='tight')
    plt.close(fig)


def get_quantiles(series: pd.Series) -> list[float]:
    # get 1 to 99 percentiles
    return np.percentile(series, q=np.arange(1, 100, 1)).tolist()


def get_mad(series: pd.Series,
            median: float) -> float:
    mad_value = np.median(np.abs(series - median))
    return float(mad_value)


def get_description(series: pd.Series) -> pd.Series:
    return series.describe()


def get_peaks(series: pd.Series,
              origin_col: str,
              n: int) -> tuple[int, list[float | None]]:
    # init KDE
    if 'beta' in origin_col:
        kde = gaussian_kde(series.values, bw_method=1)
    else:
        kde = gaussian_kde(series.values)

    # create linspace and calculate density
    x_vals = np.linspace(series.min(), series.max(), 1000)
    spaced_density = kde(x_vals)

    # find peaks
    if 'beta' in origin_col:
        peaks_idxs, _ = find_peaks(spaced_density, distance=10, prominence=.01)
    elif 'depth' in origin_col:
        peaks_idxs, _ = find_peaks(spaced_density, prominence=.1)
    else:
        peaks_idxs, _ = find_peaks(spaced_density, prominence=0)

    peaks = x_vals[peaks_idxs]
    densities = spaced_density[peaks_idxs]

    sorted_peaks = peaks[np.argsort(-densities)]

    if len(sorted_peaks) >= n:
        selected_peaks = (sorted_peaks[: n]).tolist()
    else:
        selected_peaks = sorted_peaks.tolist() + [None] * (n - len(sorted_peaks))

    return len(peaks), selected_peaks


def run_stat(args: Namespace):
    input_path: Path = Path(args.input)
    logger.info(f'Parsed arguments: {args}')

    if input_path.exists() and input_path.is_file():
        pass
    else:
        raise FileNotFoundError(f'Input file {input_path} does not exist.')

    fname: str = input_path.name.split('.')[0]

    output_dir: Path = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    font_path: Path = Path(args.font_path)
    if not font_path.exists() or not font_path.is_file():
        raise FileNotFoundError(f'Font file {font_path} does not exist.')

    img_format_list: list[str] = []
    if args.svg:
        img_format_list.append('svg')
    if args.png:
        img_format_list.append('png')
    if args.pdf:
        img_format_list.append('pdf')
    if args.eps:
        img_format_list.append('eps')
    if args.tiff:
        try:
            import PIL
            img_format_list.append('tiff')
        except ImportError:
            logger.warning('Pillow library is not installed, TIFF format will not be available.')

    if not img_format_list:
        img_format_list = ['svg']

    info_cols: list[str] = args.info_col.split(',')

    logger.info(f'[{fname}] Run Config:\n'
                f'--input path:\t{input_path.as_posix()}\n'
                f'--output directory:\t{output_dir.as_posix()}\n'
                f'--font file path:\t{font_path.as_posix()}\n'
                f'--info columns:\t\t{info_cols}\n'
                f'--image formats:\t{img_format_list}\n'
                '=====================================================================')

    df: pd.DataFrame = pd.read_parquet(input_path)
    logger.info(f'[{fname}] Loaded data from {input_path} with shape {df.shape}')

    for info_col in info_cols:
        if info_col not in df.columns:
            raise ValueError(f'Column {info_col} not found in input file {input_path}.')

        logger.info(f'[{fname}] Processing column: {info_col}')

        # draw distribution
        output_img_path: Path = output_dir / f'{fname}-{info_col.strip()}'
        draw_distribution(df, info_col, fname, output_img_path, img_format_list, font_path)
        logger.info(f'[{fname}]\t{info_col}: Distribution image saved to {output_img_path.as_posix()}')

        # get description
        target_series: pd.Series = df[df['depth'] >= 10][info_col] if 'depth' in info_col else df[info_col]
        description: pd.Series = get_description(target_series)
        logger.info(f'[{fname}]\t{info_col}:\ndescription: {description}')

        # get mad
        mad = get_mad(target_series, description['50%'])

        # get cv
        cv = description['std'] / description['mean'] if description['mean'] != 0 else np.nan

        # get skew and kurtosis
        skew_value: float = float(skew(target_series.to_numpy(), nan_policy='omit'))
        kurtosis_value: float = df[info_col].kurt()

        # get quantiles
        quantiles: list[float] = get_quantiles(target_series)
        logger.info(f'[{fname}]\t{info_col}: quantiles calculated: {quantiles}')

        # get IQR
        iqr = quantiles[74] - quantiles[24]

        # get entropy
        hist, _ = np.histogram(target_series, bins=20, density=True)
        entropy_value: float = entropy(hist, base=2) if np.any(hist) else 0.0

        # get peaks
        peak_num: int
        peaks: list[float]
        peak_num, peaks = get_peaks(target_series, info_col, args.num_peaks)
        logger.info(f'[{fname}]\t{info_col}: {peak_num} peaks found, best 2 at {peaks}')

        # save metadata
        with open(output_img_path.with_suffix('.pkl'), 'wb') as f:
            pickle_dump({
                'quantiles': quantiles,
                'description': [
                    description.tolist()[i] for i in (0, 1, 2, 3, 7)
                ] + [
                    mad, cv, skew_value, kurtosis_value, iqr, entropy_value, peak_num
                ],
                'peaks': peaks}, f)
        logger.info(f'[{fname}]\t{info_col}: quantiles and peaks saved to '
                    f'{output_img_path.with_suffix(".pkl").as_posix()}')


def to_float(a_list) -> list[float]:
    return [float(x) if x is not None else None for x in a_list]


def run_merge(args: Namespace):
    output_dir: Path = Path(args.output_dir)
    if not output_dir.exists() or not output_dir.is_dir():
        raise FileNotFoundError(f'Output directory {output_dir} does not exist.')

    pkl_files: list[Path] = list(output_dir.glob('*.pkl'))
    if not pkl_files:
        raise FileNotFoundError(f'No .pkl files found in {output_dir.as_posix()}')

    merged_data: list[list[str | float | int]] = []
    for pkl_file in pkl_files:
        fname: str = pkl_file.parent.name
        info_col: str = pkl_file.stem
        with open(pkl_file, 'rb') as f:
            data = pickle_load(f)
            peaks: list[float | None] = data['peaks']

            if len(peaks) < args.num_peaks:
                padded_peaks: list[float | None] = peaks + [None] * (args.num_peaks - len(peaks))
                merged_data.append([fname, info_col] + to_float(data['quantiles']) + to_float(data['description']) + to_float(padded_peaks))
            else:
                merged_data.append([fname, info_col] + to_float(data['quantiles']) + to_float(data['description']) + to_float(data['peaks']))

    df: pd.DataFrame = pd.DataFrame(
        merged_data,
        columns=[
            'fname', 'info_col'
        ] + [
            f'q{i}' for i in range(1, 100)
        ] + [
            'count', 'mean', 'std', 'min', 'max', 'mad', 'cv', 'skew', 'kurtosis', 'iqr', 'entropy', 'peak_num'
        ] + [
            f'peak{i}' for i in range(1, args.num_peaks + 1)
        ])
    output_path: Path = output_dir / 'merged_stats'
    df.to_csv(output_path.with_suffix('.csv'), index=False)
    df.to_parquet(output_path.with_suffix('.parquet'))
    logger.info(f'Merged data saved to {output_path.as_posix()}.csv and .parquet')


def main():
    arg_parser: ArgumentParser = ArgumentParser(description='Plot and save the distribution '
                                                            'of info cols in formatted bedgraph')
    subparsers = arg_parser.add_subparsers(dest='subcommand', required=True)
    stat_parser = subparsers.add_parser('stat', help='Perform statistical analysis and plotting')
    stat_parser.add_argument('-i', '--input', dest='input',
                             type=str, required=True,
                             help='Input file path (formatted bedgraph)')
    stat_parser.add_argument('-o', '--output-dir', dest='output_dir',
                             type=str, required=True,
                             help='Root Output directory for statistics and images, '
                                  'there will be subdirs for every samples')
    stat_parser.add_argument('-f', '--font-path', dest='font_path',
                             type=str, required=True,
                             help='Font file path')
    stat_parser.add_argument('-c', '--info-col', dest='info_col',
                             type=str, default='beta,depth',
                             help='Info column to plot (e.g., beta,beta_pyro,depth)')
    stat_parser.add_argument('-n', '--num-peaks', dest='num_peaks',
                             type=int, default=2, help='Number of peaks to find')
    stat_parser.add_argument('--svg', dest='svg', action='store_true',
                             help='Generate SVG image')
    stat_parser.add_argument('--tiff', dest='tiff', action='store_true',
                             help='Generate TIFF image  # pillow lib needed')
    stat_parser.add_argument('--png', dest='png', action='store_true',
                             help='Generate PNG image')
    stat_parser.add_argument('--pdf', dest='pdf', action='store_true',
                             help='Generate PDF image')
    stat_parser.add_argument('--eps', dest='eps', action='store_true',
                             help='Generate EPS image')
    stat_parser.set_defaults(func=run_stat)

    merge_parser = subparsers.add_parser('merge', help='Merge all .pkl files into a single CSV')
    merge_parser.add_argument('-o', '--output-dir', dest='output_dir', type=str, required=True,
                              help='Root output directory to search for .pkl files')
    merge_parser.add_argument('-n', '--num-peaks', dest='num_peaks',
                              type=int, default=2, help='Number of peaks to find')
    merge_parser.set_defaults(func=run_merge)

    args = arg_parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
