install.packages(
  c("optparse", "glue"),
  repos = "https://mirrors.sjtug.sjtu.edu.cn/cran/",
  quiet = TRUE
)

library(optparse)
option_list <- list(
  make_option(
    "--package",
    type = "character",
    help = "Package to install"),
  make_option(
    "--cores",
    type = "integer",
    default = 4,
    help = "Number of cores to use for installation")
)

parser <- OptionParser(
    option_list = option_list
)
args <- parse_args(parser)

options(
    repos = c(
        CRAN = "https://mirrors.sjtug.sjtu.edu.cn/cran/"),
    Ncpus = args$cores
)
install.packages(
  c(
    "devtools",
    "testthat",
    "BiocManager")
  )
BiocManager::install(version='devel', ask = FALSE)

if (args$package == "") {

  cat("No package specified. Skipping installation.\n")

} else {
  package_lower <- tolower(args$package)

  if (package_lower == "methylsig") {

    BiocManager::install(c("methylSig", "bsseqData", "BiocParallel"))

  } else if (package_lower == "methylkit") {

    BiocManager::install(c(
        "methylKit", "genomation",
        "GenomicFeatures", "BiocParallel"
        ))

  } else if (package_lower == "dmrcate") {

    BiocManager::install(c(
        "DMRcate",
        "BiocParallel",
        "tissueTreg",  # example data for WGBS analysis
        "DMRcatedata",  # example data for EPICv2 analysis
        "FlowSorted.Blood.EPIC"  # example data for EPICv1 analysis
        ))
    system("mkdir -p /root/.cache/R/ExperimentHub")
    system("mkdir -p /root/.cache/R/AnnotationHub")

  } else if (package_lower == "annotatr") {

    BiocManager::install(
      c(
        "annotatr", "BiocStyle", "BiocParallel",
        "org.Hs.eg.db",
        "TxDb.Hsapiens.UCSC.hg38.knownGene"
      )
    )
    system("mkdir -p /root/.cache/R/AnnotationHub")
    library(annotatr)
    repeat {
      tryCatch({
        hg38_annotations <- build_annotations(
          genome      = "hg38",
          annotations = builtin_annotations()[grepl("hg38", builtin_annotations())]
        )
        break
      }, error = function(e) {
        message("Error: ", e)
        message("Retrying in 5 seconds...")
        Sys.sleep(5)
      })
    }
    system("mkdir -p /opt")
    saveRDS(hg38_annotations, file = "/opt/hg38_annotations.rds")

  } else if (package_lower == "epialleler") {

    install.packages(c("RUnit", "knitr", "ggplot2", "rmarkdown"))
    BiocManager::install(c(
      "epialleleR",
      "BiocParallel",
      "GenomeInfoDb",
      "SummarizedExperiment",
      "VariantAnnotation"))

  } else {
    BiocManager::install(args$package)
  }
}

cat("Installation process completed.\n")
