from argparse import ArgumentParser, Namespace
from pathlib import Path
import polars as pl
from tqdm import tqdm


contigs: list[str] = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM', 'lambda', 'pUC19']
types: list[str] = ['cytosine', 'reads']


def main():
    arg_parser: ArgumentParser = ArgumentParser(description='Count the cytosines and reads on chromosomes')
    arg_parser.add_argument('-i', '--input-dir', type=Path, help='Input directory')
    arg_parser.add_argument('-o', '--output-dir', type=Path, help='Output directory')
    args: Namespace = arg_parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    files = [i for i in args.input_dir.glob('*.parquet*')]

    chrom_cytosine_dfs: list[pl.DataFrame] = []
    chrom_reads_dfs: list[pl.DataFrame] = []
    chrom_depth_dfs: list[pl.DataFrame] = []
    chrom_depth_reads_dfs: list[pl.DataFrame] = []
    chrom_cytosine_counts: pl.DataFrame
    chrom_reads_counts: pl.DataFrame

    for f in tqdm(files):
        fname = f.name.split('.')[0]
        df = pl.read_parquet(f)
        chrom_groupby = df.group_by('chrom')
        (chrom_cytosine_dfs.append(chrom_groupby
                                   .agg(pl.len().cast(pl.Int64).alias('count'))
                                   .with_columns(pl.lit(fname)
                                                   .alias('fname'))))
        (chrom_reads_dfs.append(chrom_groupby
                                .agg(pl.col('depth').sum().cast(pl.Int64).alias('count'))
                                .with_columns(pl.lit(fname)
                                                .alias('fname'))))

        chrom_depth_groupby = (df.with_columns(pl.col('depth')
                                               .cut(breaks=[2, 3, 4, 5, 10, 20, 30, 40, 50, 60],
                                                    labels=['1x', '2x', '3x', '4x',
                                                            '5–10x', '10–20x', '20–30x', '30–40x',
                                                            '40–50x', '50–60x', '≥60x'])
                                               .alias('depth_bin')).group_by('chrom', 'depth_bin'))
        (chrom_depth_dfs.append(chrom_depth_groupby
                                .agg(pl.len().cast(pl.Int64).alias('count'))
                                .with_columns(pl.lit(fname).alias('fname'))))
        (chrom_depth_reads_dfs.append(chrom_depth_groupby
                                      .agg(pl.col('depth').sum().cast(pl.Int64).alias('count'))
                                      .with_columns(pl.lit(fname).alias('fname'))))

    chrom_cytosine_tmp: pl.DataFrame
    chrom_reads_tmp: pl.DataFrame
    (chrom_cytosine_tmp,
     chrom_reads_tmp) = ([pl.concat(dfs)
                            .with_columns(pl.col('count').fill_null(0))
                            .pivot(index='fname', on='chrom', values='count')
                            .with_columns(pl.sum_horizontal(pl.all().exclude('fname')).alias('total'))
                            .with_columns(pl.col('fname').str.head(n=3).alias('lab'),
                                          pl.col('fname').str.slice(offset=4, length=2).alias('label'),
                                          pl.col('fname').str.slice(offset=7, length=1).alias('rep'),
                                          pl.lit(types[idx]).alias('count_type'))
                          for idx, dfs in enumerate([chrom_cytosine_dfs, chrom_reads_dfs])])
    cytosine_missing_contigs: list[str] = [i for i in contigs if i not in chrom_cytosine_tmp.columns]
    reads_missing_contigs: list[str] = [i for i in contigs if i not in chrom_reads_tmp.columns]

    (chrom_cytosine_counts,
     chrom_reads_counts) = [(df.with_columns(pl.sum_horizontal(pl.all()
                                                                 .exclude('fname', 'lab',
                                                                          'label', 'rep',
                                                                          'count_type'))
                                               .alias('total'))
                               .with_columns([pl.lit(0).alias(contig)
                                              for contig in missing_contigs])
                               .with_columns(pl.when(pl.col('label').is_in({'D5', 'D6', 'M8', 'BC', 'BL',
                                                                            'T1', 'T2', 'T3', 'T4'}))
                                               .then(pl.col('total') - pl.col('chrY'))
                                               .otherwise(pl.col('total'))
                                               .alias('total'))
                               .with_columns(pl.when(pl.col('lab').is_in({'RR1', 'RM1'}))
                                               .then(pl.col('total') - pl.col('pUC19') - pl.col('lambda'))
                                               .otherwise(pl.col('total'))
                                               .alias('total'))
                               .select('lab', 'label', 'rep', 'count_type', 'fname', 'total', *contigs)
                               .sort(by=['count_type', 'lab', 'label', 'rep']))
                            for missing_contigs, df in zip([cytosine_missing_contigs, reads_missing_contigs],
                                                           [chrom_cytosine_tmp, chrom_reads_tmp])]

    (pl.concat([chrom_cytosine_counts, chrom_reads_counts])
       .write_parquet(args.output_dir / 'chrom_stats.parquet.lz4'))

    (pl.concat([i.with_columns(pl.lit('cytosine').alias('count_type')) for i in chrom_depth_dfs
                ] + [j.with_columns(pl.lit('reads').alias('count_type')) for j in chrom_depth_reads_dfs
                     ])
       .with_columns(pl.col('fname').str.head(n=3).alias('lab'),
                     pl.col('fname').str.slice(offset=4, length=2).alias('label'),
                     pl.col('fname').str.slice(offset=7, length=1).alias('rep'))
       .select('lab', 'label', 'rep', 'fname', 'chrom', 'depth_bin', 'count_type', 'count')
       .sort(by=['count_type', 'lab', 'label', 'rep', 'chrom', 'depth_bin', 'count'])
       .write_parquet(args.output_dir / 'chrom_depth_stats.parquet.lz4'))


if __name__ == '__main__':
    main()
