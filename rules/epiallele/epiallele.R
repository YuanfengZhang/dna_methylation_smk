# load libs
suppressPackageStartupMessages({
  library(optparse)
  library(epialleleR)
  library(glue)
})

option_list <- list(
  make_option(
    c("--basename"), type = "character",
    dest = "basename", help = "BaseName of bam file"),
  make_option(
    c("--ref-name"), type = "character", dest = "ref_name",
    help = "file name of reference fasta file"),
  make_option(
    c("--threads"), type = "integer", dest = "threads",
    default = 1, help = "threads to use")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)
basename <- args$basename
ref_name <- args$ref_name
threads <- args$threads

required_options <- c("basename", "ref_name")
for (opt in required_options) {
  if (is.null(args[[opt]])) {
    cat("required option", opt, "not found\n")
    quit('no', status = 1, runLast = FALSE)
  }
}

copied_bam <- glue("/data/process/{basename}.bam")
preprocessed_bam <- glue("/data/process/{basename}.preprocessed.bam")
sortn_bam <- glue("/data/process/{basename}.sortn.bam")
report_file <- glue("/data/process/{basename}.raw_report")
ref_path <- glue("/data/ref/{ref_name}")

log_message <- function(message) {
  cat(glue("[{Sys.time()}] {message}\n\n"))
}

# Begin with showing the config
log_message(glue("Running config:"))
log_message(glue("\tbasename={basename}, ref_name={ref_name}, threads={threads}"))

if (!file.exists(sortn_bam)) {
  log_message("sortn bam missing, start processing from the original bam file...")

  log_message(glue("preprocessGenome"))
  genome <- preprocessGenome(ref_path, nthreads = threads)

  log_message(glue("callMethylation: {copied_bam} -> {sortn_bam}"))
  callMethylation(copied_bam, preprocessed_bam, genome, nthreads = threads)
  log_message("bam preprocessed")

  system(glue("samtools sort -n -@ {threads} -o {sortn_bam} {preprocessed_bam}"))
  log_message(glue("bam sorted"))

} else {
  log_message(glue("using existing sortn.bam: {sortn_bam}"))
}

log_message(glue("Generating cytosine report: {report_file}"))
generateCytosineReport(
  sortn_bam,
  report.file = report_file,
  threshold.reads = TRUE,
  threshold.context = "CX",
  report.context = "CX",
  paired = TRUE,
  min.mapq = 5,
  min.baseq = 10,
  skip.duplicates = TRUE,
  gzip = FALSE,
  nthreads = args$threads,
  verbose = TRUE
)

system(glue("chmod 777 -R /data/process"))
log_message("Done")
