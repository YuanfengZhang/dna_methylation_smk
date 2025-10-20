suppressPackageStartupMessages({
  library(BiocParallel)
  library(bsseqData)
  library(bsseq)
  library(data.table)
  library(glue)
  library(methylSig)
  library(optparse)
})

option_list <- list(
  make_option(
    c("-l", "--lab"), type = "character"),
  make_option(
    c("-d", "--depth"), type = "integer", default = 5,
    help = "Minimum coverage depth for filtering"),
  make_option(
    c("-m", "--max"), type = "integer", default = 1000,
    help = "Maximum coverage depth for filtering"),
  make_option(
    c("-p", "--ncores"), type = "integer", default = 4,
    help = "Number of cores processes to use",
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
lab <- opt$lab
min_count <- opt$depth
max_count <- opt$max
ncores <- opt$ncores

cat(
  glue(
    "Runtime Config:
     lab\t\t{lab}
     min_count\t{min_count}
     max_count\t{max_count}
     ncores\t{ncores}",
    .trim = FALSE
  )
)

setDTthreads(threads = ncores)

bioc_param <- MulticoreParam(workers = opt$ncores)
register(bioc_param)
bplog(bioc_param) <- TRUE
bpthreshold(bioc_param) <- "INFO"


pairs <- list(
  c("D6", "D5"),
  c("D6", "F7"),
  c("D6", "M8"),
  c("BC", "BL")
)


log_msg <- function(group_pair, stage, msg, total = 9) {
  message(
    glue(
      "\n\n\033[1;32m{format(Sys.time(), '%Y-%m-%d-%H:%M:%S')}: {group_pair[1]} vs {group_pair[2]} ",
      "[{stage}/{total}] {msg}\033[0m\n"))
}


read_for_methylsig <- function(run_name, label) {
  dat <- fread(
    glue("data/input/{run_name}_{label}.bed"),
    colClasses = c("character", "integer", "character", "numeric", "numeric"))

  dat[strand == "-", start := start - 1]

  bsseq_object <- BSseq(
    pos         = dat$start,
    chr         = dat$chrom,
    M           = as.matrix(dat$M, ncol = 1),
    Cov         = as.matrix(dat$Cov, ncol = 1),
    sampleNames = label)

  strand(rowRanges(bsseq_object)) <- dat$strand

  return(bsseq_object)
}

format_methylsig_result <- function(test_result) {
  df <- as.data.frame(test_result)
  df$start <- ifelse(df$strand == "-", df$start + 1, df$start)
  colnames(df)[colnames(df) == "seqnames"] <- "chrom"
  df <- df[c("chrom", "start", "pvalue", "fdr")]
  df <- df[!is.na(df$pvalue) & !is.na(df$fdr), ]

  return(df)
}


dm_analysis <- function(lab, group_pair, min_count, max_count, ncores) {
  log_msg(group_pair = group_pair, stage = 1, msg = "constructing file paths")
  sample_ids <- c(
    glue("{group_pair[1]}_1"),
    glue("{group_pair[1]}_2"),
    glue("{group_pair[2]}_1"),
    glue("{group_pair[2]}_2"))

  run_name <- glue("{lab}_{group_pair[1]}_vs_{group_pair[2]}")

  output_beta_binomial <- file.path(glue("/data/output/{run_name}.beta_binomial.methylsig"))
  output_binomial <- file.path(glue("/data/output/{run_name}.binomial.methylsig"))
  output_dss <- file.path(glue("/data/output/{run_name}.dss.methylsig"))

  log_msg(group_pair = group_pair, stage = 2, msg = "reading files as BSseq object")
  bsseq_objects <- lapply(sample_ids, function(label) {
    read_for_methylsig(run_name = run_name, label = label)
  })
  bsseq_merged <- do.call(combine, bsseq_objects)
  pData(bsseq_merged)$group <- c("treated", "treated", "control", "control")

  log_msg(group_pair = group_pair, stage = 3, msg = "filtering & smoothing")
  bs <- filter_loci_by_coverage(bsseq_merged, min_count = min_count, max_count = max_count)
  bs <- BSmooth(bs)

  log_msg(group_pair = group_pair, stage = 4, msg = "Find DMCs with beta-binomial model")
  beta_binomial_result <- diff_methylsig(
    bs                = bs,
    group_column      = "group",
    comparison_groups = c("case" = "treated", "control" = "control"),
    disp_groups       = c("case" = TRUE, "control" = TRUE),
    local_window_size = 0,
    t_approx          = TRUE,
    n_cores           = ncores)

  log_msg(group_pair = group_pair, stage = 5, msg = "writing output for beta binomial model")
  write.table(
    format_methylsig_result(beta_binomial_result),
    file      = output_beta_binomial,
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE)

  system(glue("chmod 777 {output_beta_binomial}"))

  log_msg(group_pair = group_pair, stage = 6, msg = "Find DMCs with binomial model")
  binomial_result <- diff_binomial(
    bs                = bs,
    group_column      = "group",
    comparison_groups = c("case" = "treated", "control" = "control")
  )

  log_msg(group_pair = group_pair, stage = 7, msg = "writing output for binomial model")
  write.table(
    format_methylsig_result(binomial_result),
    file      = output_binomial,
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE)

  system(glue("chmod 777 {output_binomial}"))

  log_msg(group_pair = group_pair, stage = 8, msg = "Find DMCs with DML model")
  dss_result <- diff_dss_test(
    bs                       = bs,
    methylation_group_column = "group",
    methylation_groups       = c("case" = "treated", "control" = "control"),
    contrast = matrix(c(0, 1), ncol = 1),
    diff_fit                 = diff_dss_fit(bs      = bs,
                                            design  = bsseq::pData(bs),
                                            formula = as.formula("~ group")))

  log_msg(group_pair = group_pair, stage = 9, msg = "writing output for DML model")
  write.table(
    format_methylsig_result(dss_result),
    file      = output_dss,
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE)

  system(glue("chmod 777 {output_dss}"))
}

bplapply(pairs, function(group_pair) {
  tryCatch({
    log_msg(group_pair = group_pair, stage = 0, msg = "beginning")
    dm_analysis(
      lab        = lab,
      group_pair = group_pair,
      min_count  = min_count,
      max_count  = max_count,
      ncores     = ncores
    )
  },
  error = function(e) {
    log_msg(group_pair = group_pair, stage = -1, msg = "FATAL ERROR")
    message(glue("error message:\n{e$message}"))
    NULL})}, BPPARAM = bioc_param)
