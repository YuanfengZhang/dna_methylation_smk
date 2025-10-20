suppressPackageStartupMessages({
  library(BiocParallel)
  library(glue)
  library(methylKit)
  library(optparse)
})


option_list <- list(
  make_option(
    c("-l", "--lab"), type = "character"),
  make_option(
    c("-d", "--depth"), type = "integer", default = 5,
    help = "Minimum coverage depth for filtering"),
  make_option(
    c("--lo_perc"), type = "numeric", default = 0.1,
    help = "Lower percentile for filtering"),
  make_option(
    c("--hi_perc"), type = "numeric", default = 99.9,
    help = "Upper percentile for filtering"),
  make_option(
    c("--hi_count"), type = "integer", default = NULL,
    help = "Maximum count for filtering"),
  make_option(
    c("-p", "--ncores"), type = "integer", default = 4,
    help = "Number of cores processes to use",
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
lab <- opt$lab
depth <- opt$depth
lo_perc <- opt$lo_perc
hi_perc <- opt$hi_perc
hi_count <- opt$hi_count

bioc_param <- MulticoreParam(workers = opt$ncores)
register(bioc_param)
bplog(bioc_param) <- TRUE
bpthreshold(bioc_param) <- "TRACE"


cat(
  glue(
    "Runtime Config:
     lab\t\t{opt$lab}
     depth\t{opt$depth}
     lo_perc\t{opt$lo_perc}
     hi_perc\t{opt$hi_perc}
     hi_count\t{opt$hi_count %||% 'Null'}
     ncores\t{opt$ncores}",
    .trim = FALSE
  )
)

pairs <- list(
  c("D6", "D5"),
  c("D6", "F7"),
  c("D6", "M8"),
  c("BC", "BL")
)


log_msg <- function(group_pair, stage, msg, total = 5) {
  message(
    glue(
      "\n\n\033[1;32m{format(Sys.time(), '%Y-%m-%d-%H:%M:%S')}: {group_pair[1]} vs {group_pair[2]} ",
      "[{stage}/{total}] {msg}\033[0m\n"))
}


dm_analysis <- function(lab, group_pair, mincov, lo_perc, hi_perc, hi_count){
  # construct file paths
  log_msg(group_pair = group_pair, stage = 1, msg = "constructing file paths")
  sample_ids <- list(
    glue("{lab}_{group_pair[1]}_vs_{group_pair[2]}_{group_pair[1]}_1"),
    glue("{lab}_{group_pair[1]}_vs_{group_pair[2]}_{group_pair[1]}_2"),
    glue("{lab}_{group_pair[1]}_vs_{group_pair[2]}_{group_pair[2]}_1"),
    glue("{lab}_{group_pair[1]}_vs_{group_pair[2]}_{group_pair[2]}_2")
  )
  run_name <- glue("{lab}_{group_pair[1]}_vs_{group_pair[2]}")
  output_file <- file.path(glue("/data/output/{run_name}.methylkit"))

  # read files
  log_msg(group_pair = group_pair, stage = 2, msg = "reading files as methylRawList object")
  myobj <- methRead(
    as.list(glue("/data/input/{sample_ids}.bed")),
    sample.id = sample_ids,
    assembly  = "GRCh38.p14",
    treatment = c(1, 1, 0, 0),
    context   = "CpG",
    mincov    = mincov)

  # filter and normalize
  log_msg(group_pair = group_pair, stage = 3, msg = "filtering & normalizing")
  if (!(lab %in% c("BS0", "EM0"))) {
    filtered_myobj <- filterByCoverage(
      myobj,
      hi.count = hi_count,
      lo.perc  = lo_perc,
      hi.perc  = hi_perc)} else {
    filtered_myobj <- myobj}
  normalized_myobj <- normalizeCoverage(filtered_myobj)
  # DMC
  log_msg(group_pair = group_pair, stage = 4, msg = "finding DMCs")
  meth <- unite(normalized_myobj, destrand=FALSE)
  mydiff <- calculateDiffMeth(meth)
  mydata <- getData(mydiff)[, c("chr", "start", "pvalue", "qvalue")]
  colnames(mydata)[colnames(mydata) == "chr"] <- "chrom"

  log_msg(group_pair = group_pair, stage = 5, msg = "writing output")
  write.table(
    mydata,
    file      = output_file,
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE)
  system(glue("chmod 777 {output_file}"))
}

bplapply(pairs, function(group_pair) {
  tryCatch({
    log_msg(group_pair = group_pair, stage = 0, msg = "beginning")
    dm_analysis(
      lab        = lab,
      group_pair = group_pair,
      mincov     = depth,
      lo_perc    = lo_perc,
      hi_perc    = hi_perc,
      hi_count   = hi_count
    )
  },
  error = function(e) {
    log_msg(group_pair = group_pair, stage = -1, msg = "FATAL ERROR")
    message(glue("error message:\n{e$message}"))
    NULL
  })}, BPPARAM = bioc_param)
