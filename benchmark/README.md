# Benchmark Scripts Documentation

This directory contains scripts for benchmarking DNA methylation analysis methods. The scripts evaluate conversion rates, prepare data for differential methylation calling (DMC), evaluate DMC results, and calculate quality metrics.

The scripts heavily rely on `polars` for efficient data processing and support various input/output formats including parquet, CSV, and BED files. To be able to run the scripts, please create a conda env:
```bash
mamba create -n data -c conda-forge pandas numpy scipy statsmodel matplotlib seaborn xlsxwriter openpyxl ipykernel jupyter autopep8 python=3.13
mamba activate data
pip install 'polars[numpy,pandas,pyarrow,openpyxl,xlsxwriter,xlsx2csv]'
```

## Overview

The benchmark workflow consists of several stages:

1. **Conversion Rate Analysis** (`cr.py`) - Calculate bisulfite/enzymatic conversion rates
2. **Data Preparation** (`dmc_prepare.py`) - Format methylation data for DMC tools
3. **DMC Analysis** (`methylkit.R`, `methylsig.R`) - Run differential methylation calling
4. **DMC Evaluation** (`dmc_evaluate.py`) - Evaluate DMC results against reference
5. **Quality Metrics** (`evaluate.py`, `snr.py`) - Calculate RMSE, Spearman correlation, and SNR

---

## Scripts

### 1. `cr.py` - Conversion Rate Calculator

Calculates bisulfite/enzymatic conversion rates from methylation data using control sequences (lambda, pUC19, chrM).

#### Usage

```bash
python cr.py \
  -i <input_dir> \
  -ro <raw_output_file> \
  -so <stat_output_file> \
  -s <file_suffix> \
  [-l <lab>] \
  [-c <chromosomes>] \
  [-d <depth_cutoff>]
```

#### Arguments

| Argument | Type | Required | Default | Description |
|----------|------|----------|---------|-------------|
| `-i, --input-dir` | str | Yes | - | Directory containing formatted methylation files (parquet/csv/tsv) |
| `-ro, --raw-output` | str | Yes | - | Output file for raw conversion rates (CSV) |
| `-so, --stat-output` | str | Yes | - | Output file for conversion rate statistics (CSV) |
| `-s, --suffix` | str | Yes | - | File suffix to match (e.g., `.parquet.lz4`) |
| `-l, --lab` | str | No | - | Filter by lab prefix (e.g., `BS1`). If not provided, all files are used |
| `-c, --chromosomes` | str | No | `chrM,lambda,pUC19` | Comma-separated control chromosomes to analyze |
| `-d, --depth-cutoff` | int | No | 1 | Minimum depth threshold for inclusion |

#### Output

- **Raw output**: CSV with conversion rates per sample
- **Stat output**: CSV with min/max/mean/median statistics grouped by lab

#### Example

```bash
python cr.py \
  -i /path/to/methylation/files \
  -ro conversion_rates.csv \
  -so conversion_stats.csv \
  -s .parquet.lz4 \
  -l BS1 \
  -d 5
```

---

### 2. `dmc_prepare.py` - DMC Data Preparation

Prepares methylation data in the format required by DMC tools (methylKit, methylSig, cpgtools).

#### Usage

```bash
python dmc_prepare.py \
  -i <input_dir> \
  -r <ref_dir> \
  -o <output_dir> \
  -l <labs> \
  -s <labels> \
  -m <method> \
  [-d <depth>] \
  [-c <chrom>] \
  [-a <action>] \
  [--ref-only] \
  [--beta-col <col_name>]
```

#### Arguments

| Argument | Type | Required | Default | Description |
|----------|------|----------|---------|-------------|
| `-i, --input-dir` | str | Yes | - | Directory with formatted parquet files |
| `-r, --ref-dir` | str | Yes | - | Directory with reference parquet files |
| `-o, --output-dir` | str | Yes | - | Output directory for prepared files |
| `-l, --labs` | str | Yes | - | Comma-separated lab names (e.g., `BS1,BS2,BS3`) |
| `-s, --labels` | str | Yes | - | Two comma-separated sample labels (e.g., `D6,D5`) |
| `-m, --method` | str | Yes | - | DMC method: `methylkit`, `methylsig`, or `cpgtools` |
| `-d, --depth` | int | No | 5 | Minimum depth threshold |
| `-c, --chrom` | str | No | `chr17` | Chromosome to process |
| `-a, --action` | str | Yes | - | Action: `prepare` or `remove` |
| `--ref-only` | flag | No | - | Only generate input within reference locations |
| `--beta-col` | str | No | `beta` | Beta column name in input files |

#### Output

- **For methylKit/methylSig**: BED files with coverage and methylation frequency
- **For cpgtools**: BED file with methylation counts and group file

#### Example

```bash
# Prepare data for methylKit
python dmc_prepare.py \
  -i /path/to/cytosine \
  -r /path/to/reference \
  -o /path/to/output \
  -l BS1,BS2,BS3 \
  -s D6,D5 \
  -m methylkit \
  -d 5 \
  -a prepare

# Remove prepared files
python dmc_prepare.py \
  -i /path/to/cytosine \
  -r /path/to/reference \
  -o /path/to/output \
  -l BS1,BS2,BS3 \
  -s D6,D5 \
  -m methylkit \
  -a remove
```

---

### 3. `methylkit.R` - methylKit DMC Analysis

Performs differential methylation calling using the methylKit R package.

#### Usage

```bash
Rscript methylkit.R \
  -l <lab> \
  [-d <depth>] \
  [--lo_perc <percentile>] \
  [--hi_perc <percentile>] \
  [--hi_count <count>] \
  [-p <ncores>]
```

#### Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `-l, --lab` | str | - | Lab identifier (e.g., `BS1`) |
| `-d, --depth` | int | 5 | Minimum coverage depth for filtering |
| `--lo_perc` | numeric | 0.1 | Lower percentile for coverage filtering |
| `--hi_perc` | numeric | 99.9 | Upper percentile for coverage filtering |
| `--hi_count` | int | NULL | Maximum count for filtering |
| `-p, --ncores` | int | 4 | Number of CPU cores to use |

#### Input

Expects BED files in `/data/input/` directory with format:
```
chrBase  chrom  start  strand  coverage  freqC  freqT
```

#### Output

TSV files in `/data/output/` with columns:
```
chrom  start  pvalue  qvalue
```

#### Docker Usage

```bash
docker run --rm \
  -v /path/to/input:/data/input \
  -v /path/to/output:/data/output \
  methylkit:1.33.3 \
  Rscript methylkit.R -l BS1 -d 5 -p 8
```

---

### 4. `methylsig.R` - methylSig DMC Analysis

Performs differential methylation calling using the methylSig R package with multiple statistical models.

#### Usage

```bash
Rscript methylsig.R \
  -l <lab> \
  [-d <min_depth>] \
  [-m <max_depth>] \
  [-p <ncores>]
```

#### Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `-l, --lab` | str | - | Lab identifier (e.g., `BS1`) |
| `-d, --depth` | int | 5 | Minimum coverage depth |
| `-m, --max` | int | 1000 | Maximum coverage depth |
| `-p, --ncores` | int | 4 | Number of CPU cores to use |

#### Input

Expects BED files in `/data/input/` directory with format:
```
chrom  start  strand  M  Cov
```

#### Output

Three TSV files per sample pair in `/data/output/`:
- `{lab}_{pair}.beta_binomial.methylsig` - Beta-binomial model results
- `{lab}_{pair}.binomial.methylsig` - Binomial model results
- `{lab}_{pair}.dss.methylsig` - DSS model results

Each file contains columns:
```
chrom  start  pvalue  fdr
```

#### Docker Usage

```bash
docker run --rm \
  -v /path/to/input:/data/input \
  -v /path/to/output:/data/output \
  methylsig:1.19.0 \
  Rscript methylsig.R -l BS1 -d 5 -m 1000 -p 8
```

---

### 5. `dmc_evaluate.py` - DMC Results Evaluation

Evaluates DMC results against reference truth sets using various genomic features.

#### Usage

```bash
python dmc_evaluate.py \
  -i <input_dir> \
  -r <ref_dir> \
  -c <cytosine_dir> \
  -o <output_dir> \
  -m <method> \
  [--p-type <p_or_q>] \
  [--cutoff <pvalue>] \
  [-g, --global-only] \
  [-v <log_level>]
```

#### Arguments

| Argument | Type | Required | Default | Description |
|----------|------|----------|---------|-------------|
| `-i, --input-dir` | str | Yes | - | Directory with DMC results |
| `-r, --ref-dir` | str | Yes | - | Directory with reference DMC files |
| `-c, --cytosine-dir` | str | Yes | - | Directory with formatted cytosine files |
| `-o, --output-dir` | str | Yes | - | Output directory for evaluation results |
| `-m, --method` | str | Yes | - | DMC method: `methylkit`, `methylsig-beta-binomial`, `methylsig-binomial`, `methylsig-DSS`, `cpgtools` |
| `--p-type` | str | No | `p` | Use p-values (`p`) or q-values (`q`) |
| `--cutoff` | float | No | 0.05 | P-value cutoff for DMC detection |
| `-g, --global-only` | flag | No | - | Only evaluate global feature (skip feature-specific evaluation) |
| `-v, --log-level` | str | No | `info` | Logging level: `debug`, `info`, `warning`, `error`, `critical` |

#### Output

CSV file with evaluation metrics:
```
sample_pair  lab  p_type  feature  fgroup  count  tp  tn  fn  fp  recall  precision  f1_score  mcc  ppv  npv  roc_auc
```

#### Example

```bash
python dmc_evaluate.py \
  -i /path/to/dmc/output \
  -r /path/to/reference/dmc \
  -c /path/to/cytosine \
  -o /path/to/evaluation \
  -m methylkit \
  --p-type p \
  --cutoff 0.05 \
  -v info
```

---

### 6. `evaluate.py` - Comprehensive Quality Metrics

Calculates RMSE, Spearman correlation, and SNR for methylation data evaluation.

#### Usage

```bash
python evaluate.py \
  -i <input_dir> \
  -o <output_dir> \
  [-l <level>] \
  [--rmse] \
  [--spearmanr] \
  [--snr] \
  [--gather] \
  [--pyro] \
  [--feature] \
  [--curated] \
  [--compact] \
  [-d <depth>] \
  [--beta-col <col>] \
  [-v <log_level>]
```

#### Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `-i, --input-dir` | Path | - | Input directory with parquet files |
| `-o, --output-dir` | Path | - | Output directory for results |
| `-l, --level` | str | `c` | Reference level: `c` or `cpg` |
| `--rmse` | flag | - | Calculate RMSE metrics |
| `--spearmanr` | flag | - | Calculate Spearman correlation |
| `--snr` | flag | - | Calculate Signal-to-Noise Ratio (requires ≥1TB RAM) |
| `--gather` | flag | - | Gather all results into single files |
| `--pyro` | flag | - | Use pyro beta values |
| `--feature` | flag | - | Include feature-specific metrics |
| `--curated` | flag | - | Use curated feature set |
| `--compact` | flag | - | Use compact feature set |
| `-d, --depth` | int | 1 | Minimum depth threshold |
| `--beta-col` | str | `beta` | Beta column name |
| `-v, --log-level` | str | `info` | Logging level |

#### Output

- **RMSE**: `{level}_rmse/` directory with per-sample RMSE values
- **Spearman**: `{level}_spearmanr/` directory with correlation values
- **SNR**: `{level}_snr/` directory with signal-to-noise ratios

#### Example

```bash
python evaluate.py \
  -i /path/to/bedgraph \
  -o /path/to/output \
  -l c \
  --rmse \
  --spearmanr \
  --feature \
  --gather \
  -d 5 \
  -v info
```

---

### 7. `snr.py` - Signal-to-Noise Ratio Calculator

Calculates SNR using PCA on methylation data across multiple samples.

#### Usage

```bash
python snr.py \
  -i <input_dir> \
  -t <tmp_dir> \
  -o <output_dir> \
  -j <job> \
  [-T <type>] \
  [-l <lab>] \
  [-ld <min_depth>] \
  [-hd <max_depth>] \
  [-c <beta_col>] \
  [-rm]
```

#### Arguments

| Argument | Type | Required | Default | Description |
|----------|------|----------|---------|-------------|
| `-i, --input-dir` | str | Yes | - | Directory with parquet files |
| `-t, --tmp-dir` | str | Yes | - | Temporary directory for intermediate files |
| `-o, --output-dir` | str | Yes | - | Output directory for SNR results |
| `-j, --job` | str | Yes | - | Job type: `prepare`, `calculate`, or `merge` |
| `-T, --type` | str | No | `c` | Data type: `c` or `cpg` |
| `-l, --lab` | str | No | - | Lab identifier (required for `prepare` and `calculate`) |
| `-ld, --min-depth` | int | No | -1 | Minimum depth threshold |
| `-hd, --max-depth` | int | No | 0 | Maximum depth threshold |
| `-c, --column` | str | No | `beta` | Beta column name |
| `-rm, --remove` | flag | No | - | Remove temporary files after processing |

#### Workflow

1. **Prepare**: Merge sample data for a lab
   ```bash
   python snr.py -i /path/to/input -t /path/to/tmp -o /path/to/output -j prepare -l BS1
   ```

2. **Calculate**: Compute SNR from prepared data
   ```bash
   python snr.py -i /path/to/input -t /path/to/tmp -o /path/to/output -j calculate -l BS1
   ```

3. **Merge**: Combine all SNR results
   ```bash
   python snr.py -i /path/to/input -t /path/to/tmp -o /path/to/output -j merge -rm
   ```

#### Output

- `{lab}.snr.csv` - SNR value and PCA statistics
- `{lab}.pcs.csv` - PC1 and PC2 values for each sample

---

## Complete Workflow Example

### Step 1: Prepare Data for methylKit

```bash
python dmc_prepare.py \
  -i /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE \
  -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/calibrated \
  -o /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylkit/input \
  -l BS1,BS2,BS3,BS4 \
  -s D6,D5 \
  -m methylkit \
  -d 5 \
  -a prepare
```

### Step 2: Run methylKit DMC Analysis

```bash
docker run --rm \
  -v /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylkit/input:/data/input \
  -v /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylkit/output:/data/output \
  -v /path/to/benchmark:/data/bin \
  methylkit:1.33.3 \
  Rscript /data/bin/methylkit.R -l BS1 -d 5 -p 8
```

### Step 3: Evaluate DMC Results

```bash
python dmc_evaluate.py \
  -i /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylkit/output \
  -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/dmc \
  -c /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE \
  -o /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylkit \
  -m methylkit \
  --p-type p \
  -v info
```

### Step 4: Calculate Quality Metrics

```bash
python evaluate.py \
  -i /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE \
  -o /mnt/eqa/zhangyuanfeng/methylation/evaluated/c \
  -l c \
  --rmse \
  --spearmanr \
  --feature \
  --gather \
  -d 5 \
  -v info
```

---

## Data Formats

### Input Parquet Format

Standard methylation data format with columns:
- `chrom`: Chromosome name
- `start`: Position (0-based)
- `beta`: Methylation level (0-100)
- `depth`: Coverage depth
- `m`: Methylated count
- `u`: Unmethylated count

### DMC Input Formats

**methylKit BED format:**
```
chrBase  chrom  start  strand  coverage  freqC  freqT
chr17.1  chr17  1      +       10        50     50
```

**methylSig BED format:**
```
chrom  start  strand  M    Cov
chr17  1      +       5    10
```

---

## Dependencies

### Python Packages
- polars
- pandas
- numpy
- scikit-learn
- scipy

### R Packages
- methylKit
- methylSig
- BiocParallel
- optparse

### External Tools
- Docker (for R scripts)
- samtools
- bedtools

---

## Notes

- **SNR Calculation**: Requires significant RAM (≥1TB recommended) for large datasets
- **Parallel Processing**: Most scripts support multi-core processing via `-p` or `-t` arguments
- **Reference Data**: Ensure reference files are in the correct format and location
- **Depth Filtering**: Adjust depth cutoffs based on your sequencing depth distribution
- **Feature Evaluation**: Feature-specific metrics require additional reference annotations

---

## Troubleshooting

### Common Issues

1. **File not found errors**: Verify input directory paths and file suffixes
2. **Memory errors**: Reduce dataset size or increase available RAM
3. **Docker errors**: Ensure Docker is installed and images are available
4. **Missing columns**: Check input file format matches expected schema

### Debug Mode

Enable debug logging for detailed error messages:
```bash
python script.py ... -v debug
```

---

## References

- methylKit: https://bioconductor.org/packages/methylKit/
- methylSig: https://bioconductor.org/packages/methylSig/