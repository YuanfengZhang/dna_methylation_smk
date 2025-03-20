install.packages(
    "optparse",
    repos = "https://mirrors.sjtug.sjtu.edu.cn/cran/",
    quiet = TRUE)

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
        CRAN  = "https://mirrors.sjtug.sjtu.edu.cn/cran/"),
        Ncpus = args$cores
        )
install.packages(
    c(
        "devtools",
        "testthat",
        "BiocManager")
    )
BiocManager::install(version='devel')

if (args$package == "") {

    cat("No package specified. Skipping installation.\n")

} else {
    package_lower <- tolower(args$package)

    if (package_lower == "methylsig") {

        BiocManager::install(c("methylSig", "bsseqData"))

    } else if (package_lower == "methylkit") {

        BiocManager::install(c(
            "methylKit", "genomation", "GenomicFeatures"
            ))

    } else if (package_lower == "dmrcate") {

        BiocManager::install(c(
            "DMRcate",
            "tissueTreg",  # example data for WGBS analysis
            "DMRcatedata",  # example data for EPICv2 analysis
            "FlowSorted.Blood.EPIC"  # example data for EPICv1 analysis
            ))
        system("mkdir -p /root/.cache/R/ExperimentHub")
        system("mkdir -p /root/.cache/R/AnnotationHub")

    } else {
        BiocManager::install(args$package)
    }
}

cat("Installation process completed.\n")
