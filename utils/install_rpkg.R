install.packages(
    "optparse",
    repos = "https://mirrors.sjtug.sjtu.edu.cn/cran/",
    quiet = TRUE)

library(optparse)
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
BiocManager::install(args$package)

