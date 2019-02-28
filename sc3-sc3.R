#!/usr/bin/env Rscript 
# Script to run all steps of SC3 in one go.

# Load optparse we need to check inputs
suppressPackageStartupMessages(require(optparse))

# Load common functions
suppressPackageStartupMessages(require(workflowscriptscommon))

# parse options
option_list <- list(
  make_option(
    c("-i", "--input-object-file"),
    action = "store",
    default = NA,
    type = "character",
    help = "File name in which a serialized R SingleCellExperiment object where object matrix found"
  ),
    make_option(
    c("-c", "--ks"),
    action = "store",
    default = NULL,
    type = "character",
    help = "A comma-separated string or single value representing the number of clusters k to be used for SC3 clustering."
  ),
  # SingleCellExperiment object without a counts slot will fail to run sc3()
  # SC3 publication states "in general the gene filter did not affect the accuracy of clustering"
  # setting default to FALSE allows more SingleCellExperiment object inputs to successfully run sc3()
  make_option(
    c("-f", "--gene-filter"),
    action = "store",
    default = FALSE,
    type = "logical",
    help = "A boolean variable which defines whether to perform gene filtering before SC3 clustering."
  ),
  make_option(
    c("-p", "--pct-dropout-min"),
    action = "store",
    default = 10,
    type = "integer",
    help = "If gene_filter = TRUE, then genes with percent of dropouts smaller than pct_dropout_min are filtered out before clustering."
  ),
  make_option(
    c("-q", "--pct-dropout-max"),
    action = "store",
    default = 90,
    type = "integer",
    help = "If gene_filter = TRUE, then genes with percent of dropouts larger than pct_dropout_max are filtered out before clustering."
  ),
    make_option(
    c("-d", "--d-region-min"),
    action = "store",
    default = 0.04,
    type = "double",
    help = "Defines the minimum number of eigenvectors used for kmeans clustering as a fraction of the total number of cells. Default is 0.04. See SC3 paper for more details."
  ),
  make_option(
    c("-e", "--d-region-max"),
    action = "store",
    default = 0.07,
    type = "double",
    help = "Defines the maximum number of eigenvectors used for kmeans clustering as a fraction of the total number of cells. Default is 0.07. See SC3 paper for more details."
  ),
  make_option(
    c("-n", "--svm-num-cells"),
    action = "store",
    default = NULL,
    type = "integer",
    help = "Number of randomly selected training cells to be used for SVM prediction. Default is NULL."
  ),
  make_option(
    c("-r", "--svm-train-inds"),
    action = "store",
    default = NULL,
    type = "numeric",
    help = "Text file with one integer per line. Will be used to create a numeric vector defining indices of training cells that should be used for SVM training. The default is NULL."
  ),
  make_option(
    c("-m", "--svm-max"),
    action = "store",
    default = 5000,
    type = "integer",
    help = "Define the maximum number of cells below which SVM are not run."
  ),
  make_option(
    c("-t", "--n-cores"),
    action = "store",
    default = NULL,
    type = "integer",
    help = "Number of threads/cores to be used in the user's machine."
  ),
  make_option(
    c("-k", "--kmeans-nstart"),
    action = "store",
    default = NULL,
    type = "integer",
    help = "nstart parameter passed to kmeans function. Default is 1000 for up to 2000 cells and 50 for more than 2000 cells."
  ),
    make_option(
    c("-a", "--kmeans-iter-max"),
    action = "store",
    default = 1e+09,
    type = "double",
    help = "iter.max parameter passed to kmeans function. Default is 1e+09."
  ),
    make_option(
    c("-g", "--k_estimator"),
    action = "store",
    default = FALSE,
    type = "logical",
    help = "A boolean variable which defines whether to estimate an optimal number of clusters. If user has already defined the ks parameter the estimation does not affect the user's paramater."
  ),
    make_option(
    c("-b", "--biology"),
    action = "store",
    default = FALSE,
    type = "logical",
    help = "A boolean variable which defines whether to compute differentially expressed genes, marker genes and cell outliers."
  ),
    make_option(
    c("-s", "--rand-seed"),
    action = "store",
    default = 1,
    type = "integer",
    help = "sets the seed of the random number generator. SC3 is a stochastic method, so setting the rand_seed to a fixed values can be used for reproducibility purposes."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = "character",
    help = "File name in which to store serialized R object of type 'SingleCellExperiment'."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c("input_object_file", "output_object_file"))

# Check parameter values defined
if ( ! file.exists(opt$input_object_file)){
  stop(paste("File object or matrix", opt$input_object_file, "does not exist"))
}

# See if the user has specified a file for svm_train_inds, in which case we need to read it and create a vector

if (! is.null(opt$svm_train_inds)){
  if (! file.exists(opt$svm_train_inds)){
    stop(paste("Supplied svm_train_inds file", opt$svm_train_inds, "does not exist"))
  }else{
    svm_train_inds <- readLines(opt$svm_train_inds)
  }
}else{
  svm_train_inds <- NULL
}


# Once arguments are satisfactory, load SC3 package
suppressPackageStartupMessages(require(SC3))
suppressPackageStartupMessages(require(SingleCellExperiment))

# Input from serialized R object
SingleCellExperiment <- readRDS(opt$input_object_file)

# SC3 needs this
rowData(SingleCellExperiment)$feature_symbol <- rownames(SingleCellExperiment)

# Temporary fix to deal with https://github.com/hemberg-lab/SC3/issues/53
assignInNamespace(x = "rowSums", value = Matrix::rowSums, ns = "base")

# Determine value of ks to pass to sc3_kmeans()
if (! is.null(opt$ks)){
  ## If the user has specified ks, parse the ks to a vector and use as ks value
  ks <- wsc_parse_numeric(opt, "ks")
} else {
  ## run sc3_estimate_k() and extract resulting estimate to pass to sc3_kmeans() as 'ks' parameter
  SingleCellExperiment <- sc3_estimate_k(SingleCellExperiment)
  ks <- metadata(SingleCellExperiment)$sc3$k_estimation
}
# Create SCE object from data and run SC3
SingleCellExperiment <- sc3(SingleCellExperiment, ks = ks, gene_filter = opt$gene_filter,
                            pct_dropout_min = opt$pct_dropout_min, pct_dropout_max = opt$pct_dropout_max,
                            d_region_min = opt$d_region_min, d_region_max = opt$d_region_max,
                            svm_num_cells = opt$svm_num_cells, svm_train_inds = opt$svm_train_inds,
                            svm_max = opt$svm_max, n_cores = opt$n_cores, kmeans_nstart = opt$kmeans_nstart,
                            kmeans_iter_max = opt$kmeans_iter_max, k_estimator = opt$k_estimator,
                            biology = opt$biology, rand_seed = opt$rand_seed)

if (is.null(SingleCellExperiment)){
  stop("sc3() failed")
}

# Output to a serialized R object
saveRDS(SingleCellExperiment, file = opt$output)
