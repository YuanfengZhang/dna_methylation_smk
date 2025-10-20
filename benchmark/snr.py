from argparse import ArgumentParser, Namespace, RawTextHelpFormatter
from functools import reduce
import gc
from pathlib import Path
from textwrap import dedent
from typing import Literal
# import cudf.pandas  # noqa:E402
# cudf.pandas.install()  # noqa:E402
# from cuml.accel.core import install as cuml_install   # noqa:E402
# cuml_install(disable_uvm=False)
import pandas as pd  # noqa:E402
import polars as pl  # noqa:E402
import numpy as np  # noqa:E402
from sklearn.decomposition import PCA  # noqa:E402
from sklearn.preprocessing import StandardScaler  # noqa:E402


DESCRIPTION = dedent("""
LABS=(
  "BS1" "BS2" "BS3" "BS4"
  "EM1" "EM2" "EM3" "EM4"
  "PS1" "PS2" "PS3"
  "RR1" "RM1")

for lab in "${LABS[@]}"; do
time python snr_old.py \
  --input-dir /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE \
  --tmp-dir /mnt/eqa/zhangyuanfeng/methylation/evaluated/c/c_snr/tmp \
  --output-dir /mnt/eqa/zhangyuanfeng/methylation/evaluated/c/c_snr \
  --type c --lab BS1 --job prepare
done

for lab in "${LABS[@]}"; do
time python snr_old.py \
  --input-dir /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE \
  --tmp-dir /mnt/eqa/zhangyuanfeng/methylation/evaluated/c/c_snr/tmp \
  --output-dir /mnt/eqa/zhangyuanfeng/methylation/evaluated/c/c_snr \
  --type c --lab BS1 --job calculate
done

time python snr_old.py \
  --input-dir /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE \
  --tmp-dir /mnt/eqa/zhangyuanfeng/methylation/evaluated/c/c_snr/tmp \
  --output-dir /mnt/eqa/zhangyuanfeng/methylation/evaluated/c/c_snr \
  --type c --job merge
""")


SAMPLE_LIST: list[str] = ['D5_1', 'D5_2', 'D6_1', 'D6_2',
                          'F7_1', 'F7_2', 'M8_1', 'M8_2',
                          'T1_1', 'T1_2', 'T2_1', 'T2_2',
                          'T3_1', 'T3_2', 'T4_1', 'T4_2']


def read_lab_file(f_path: Path, lab: str, level: Literal['c', 'cpg'],
                  min_depth: int, max_depth: int,
                  beta_col: str, float2int: bool = True) -> pl.LazyFrame:
    sample = '_'.join(f_path.name.split('.')[0].split('_')[1: 3])

    df: pl.LazyFrame = pl.scan_parquet(f_path)

    if lab == 'RM1' and level == 'cpg':
        df = df.filter(pl.col('ctype') == '5mC')

    df = (df.drop_nulls(subset=[beta_col])
            .unique(subset=['chrom', 'start', 'end'])
            .rename({beta_col: sample}))
    if 'MA' not in lab:
        if min_depth > 1:
            df = df.filter(pl.col('depth') >= min_depth)
        if max_depth > max(1, min_depth):
            df = df.filter(pl.col('depth') <= max_depth)
    if float2int:
        df = df.with_columns(pl.col(sample).cast(pl.Float64).cast(pl.Int64))
    return df.select('chrom', 'start', sample)


def lab_prepare(lab_files: list[Path],
                level: Literal['c', 'cpg'],
                tmp_dir: Path,
                min_depth: int,
                max_depth: int,
                beta_col: str):
    lab = lab_files[0].stem.split('_')[0]

    (reduce(lambda x, acc: x.join(acc,
                                  on=['chrom', 'start'],
                                  how='inner'),
            [read_lab_file(f_path=p, lab=lab,
                           level=level,
                           min_depth=min_depth,
                           max_depth=max_depth,
                           beta_col=beta_col) for p in lab_files])
     .drop('chrom', 'start')
     .drop_nulls()
     .collect()  # .collect(engine='gpu')
     .write_csv(tmp_dir / f'{lab}.csv', include_header=True))
    print(f'lab {lab} merge done.')


def lab_snr(f_path: Path, output_dir: Path):
    lab = f_path.stem
    matrix: pd.DataFrame = pd.read_csv(f_path, index_col=None)

    # Remove features which variance is zero to ensure the PCA
    matrix = matrix.loc[:, matrix.var() != 0]

    lens = matrix.__len__()
    print(lens)

    if lens == 0:
        (pd.DataFrame({'lab': [lab], 'PC1': [np.nan],
                       'PC2': [np.nan], 'SNR': [np.nan],
                       'length': [0]})
           .to_csv(output_dir / f'{lab}.snr.csv',
                   index=False, header=True))
        print(f'lab {lab} no available rows to calculate, prune now.')
        return

    samples = matrix.columns.tolist()
    group = pd.DataFrame({'code': samples})
    group['code'] = group['code'].str.split('_').str[0]

    samples_group_mat: pd.DataFrame = group
    samples_group_mat['samples'] = samples

    # 对数据集进行特征缩放
    scaler = StandardScaler()
    pca = PCA(n_components=8, svd_solver='full')
    print('正在应用PCA...')
    # 应用PCA到缩放后的数据
    pcs = pca.fit_transform(scaler.fit_transform(matrix.T))

    print('PCA应用完成.')
    del matrix
    gc.collect()

    # Create DataFrame with PC values and sample samples
    pcs_df = pd.DataFrame(pcs[:, :2], columns=['PC1', 'PC2'])
    pcs_df['Sample_id'] = samples
    group['Sample_id'] = samples
    pcs = pd.DataFrame(pcs)
    pcs['Sample_id'] = samples
    pcs.set_index(pcs.iloc[:, 8], inplace=True)

    # save PC1, PC2 and PC3 values for all samples
    (pcs_df.copy()
           .rename(columns={'Sample_id': 'sample'})[['sample', 'PC1', 'PC2']]
           .to_csv(output_dir / f'{lab}.pcs.csv', index=False, header=True))

    # Merge with group information
    pcs_df = pd.merge(pcs_df, group, on=['Sample_id'], how='inner')

    importance = pca.explained_variance_ratio_

    dt_perc_pcs = pd.DataFrame({
        'PCX': np.arange(1, len(importance) + 1),
        'Percent': importance,
        'AccumPercent': np.cumsum(importance)
    })

    # dt.dist
    dt_dist = pd.DataFrame({
        'ID.A': np.repeat(samples, len(samples)),
        'ID.B': np.tile(samples, len(samples))
    })

    # dt.dist$group.A
    dt_dist['group.A'] = dt_dist.iloc[:, 0].str.split('_').str[0]
    dt_dist['group.B'] = dt_dist.iloc[:, 1].str.split('_').str[0]

    # dt.dist$Type
    dt_dist['Type'] = np.where(dt_dist['ID.A'] == dt_dist['ID.B'], 'Same',
                               np.where(dt_dist['group.A'] == dt_dist['group.B'],
                                        'Intra',
                                        'Inter'))

    dt_dist['Dist'] = (
        dt_perc_pcs.loc[0, 'Percent'] * (  # type: ignore
            pcs.loc[dt_dist['ID.A'], 0].values - pcs.loc[dt_dist['ID.B'], 0].values) ** 2 +  # type: ignore
        dt_perc_pcs.loc[1, 'Percent'] * (  # type: ignore
            pcs.loc[dt_dist['ID.A'], 1].values - pcs.loc[dt_dist['ID.B'], 1].values) ** 2  # type: ignore
    )

    # dt.dist.stats
    dt_dist_stats: pd.DataFrame = (dt_dist.groupby('Type')  # type: ignore
                                          .agg(AvgDist=pd.NamedAgg(column='Dist',
                                                                   aggfunc='mean'))
                                          .reset_index())
    print(dt_dist_stats)
    dt_dist_stats.columns = ['Type', 'AvgDist']
    dt_dist_stats.set_index('Type', inplace=True)

    # signoise_db
    signoise: float = (10 * np.log10(  # type: ignore
        dt_dist_stats.loc['Inter', 'AvgDist'] / (dt_dist_stats.loc['Intra', 'AvgDist'])  # type: ignore
    ))

    print('SNR:', signoise)  # type: ignore
    pca_stat = pd.DataFrame({
        'lab': [lab],
        'PC1': [importance[0]],
        'PC2': [importance[1]],
        'SNR': [signoise],
        'length': [lens]
    })

    pca_stat.to_csv(output_dir / f'{lab}.snr.csv', index=False, header=True)


def main():
    arg_parser: ArgumentParser = ArgumentParser(usage=DESCRIPTION,
                                                formatter_class=RawTextHelpFormatter)
    arg_parser.add_argument('-i', '--input-dir', dest='input_dir', type=str, required=True,
                            help='Input directory containing the bedgraph files.')
    arg_parser.add_argument('-t', '--tmp-dir', dest='tmp_dir', type=str, required=True,
                            help='Tmp directory containing the merged bedgraph files.')
    arg_parser.add_argument('-o', '--output-dir', dest='output_dir', type=str, required=True,
                            help='Output directory to save the snr results.')
    arg_parser.add_argument('-j', '--job', dest='job', type=str,
                            choices=['prepare', 'calculate', 'merge'],
                            help='either merge or calculate.')
    arg_parser.add_argument('-T', '--type', dest='type', type=str,
                            choices=['c', 'cpg'], default='c', help='c or cpg to perform.')
    arg_parser.add_argument('-l', '--lab', dest='lab', type=str,
                            help='lab to calculate for.')
    arg_parser.add_argument('-ld', '--min-depth', dest='min_depth', type=int, default=-1,
                            help='minimal depth to filter.')
    arg_parser.add_argument('-hd', '--max-depth', dest='max_depth', type=int, default=0,
                            help='maximal depth to filter.')
    arg_parser.add_argument('-c', '--column', dest='beta_col', type=str, default='beta',
                            help='beta column to calculate, default as beta.')
    arg_parser.add_argument('-rm', '--remove', dest='remove', action='store_true',
                            help='remove temporary files after processing.')
    args: Namespace = arg_parser.parse_args()

    input_dir: Path = Path(args.input_dir)
    tmp_dir: Path = Path(args.tmp_dir)
    output_dir: Path = Path(args.output_dir)
    min_depth: int = args.min_depth
    max_depth: int = args.max_depth

    lab: str = args.lab

    if args.job == 'prepare':
        if min_depth > max_depth > 0:
            raise ValueError('min_depth must be less than or equal to max_depth')
        print(f'Merging bedgraph files for {lab}...')

        input_dir.mkdir(parents=True, exist_ok=True)
        tmp_dir.mkdir(parents=True, exist_ok=True)

        lab_files: list[Path] = [p for p in list(input_dir.glob(f'{lab}*.parquet.lz4'))
                                 if any(s in p.stem for s in SAMPLE_LIST)]

        lab_prepare(lab_files=lab_files, level=args.type, tmp_dir=tmp_dir,
                    min_depth=min_depth, max_depth=max_depth, beta_col=args.beta_col)

    elif args.job == 'calculate':
        output_dir.mkdir(parents=True, exist_ok=True)
        lab_snr(f_path=tmp_dir / f'{lab}.csv', output_dir=output_dir)
        print('SNR results calculated')

    else:
        snr_results: pd.DataFrame = pd.concat([pd.read_csv(f) for f in output_dir.glob('*snv.csv')],
                                              ignore_index=True)
        snr_results.to_csv((output_dir / 'snr.csv'), header=True, index=False,
                           lineterminator='\n', encoding='utf-8')
        print('SNR results saved')

        if args.remove:
            tmp_files = list(tmp_dir.glob('*.csv'))
            for tmp_file in tmp_files:
                tmp_file.unlink()

            print('Temporary files deleted')


if __name__ == '__main__':
    main()
