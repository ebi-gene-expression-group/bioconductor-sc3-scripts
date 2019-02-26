#!/usr/bin/env Rscript 

# Load optparse we need to check inputs
suppressPackageStartupMessages(require(optparse))

# Load common functions
suppressPackageStartupMessages(require(workflowscriptscommon))

# parse options
option_list = list(
  make_option(
    c("-i", "--input-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which a SC3 'SingleCellExperiment' object has been stored after kmeans clustering."
  ),
  make_option(
    c("-t", "--output-text-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "A text file file to write marker matrices to. A 'k' column will defined from which value of 'k' the markers are derived."
  ),
  make_option(
    c("-k", "--ks"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "A comma-separated string or single value representing the number of clusters k to be used for SC3 clustering."
  ),
  make_option(
    c("-r", "--regime"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "defines what biological analysis to perform. \"marker\" for marker genes, \"de\"\ for differentiall expressed genes and \"outl\" for outlier cells."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name for R object of type 'SingleCellExperiment' from SC3 in which to store the consensus matrix."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file'))

# Check parameter values defined
if ( ! file.exists(opt$input_object_file)){
  stop((paste('SC3 file object', opt$input_object_file, 'does not exist.')))
}

# Parse the ks to a vector
if ( is.null(opt$ks)){
  stop((paste('Please provide a k.')))
}else{
  ks <- wsc_parse_numeric(opt, 'ks')
}

# Once arguments are satisfcatory, load Scater package
suppressPackageStartupMessages(require(SC3))
suppressPackageStartupMessages(require(SingleCellExperiment))

# Input from serialized R object
SingleCellExperiment <- readRDS(opt$input_object_file)

# Calculate consensus matrix
SingleCellExperiment  <- sc3_calc_biology(object = SingleCellExperiment, ks = ks, regime = opt$regime)

if ( ! is.na(opt$output_text_file)){

    # Get results by k
    results_cols <- grep(opt$regime, colnames(rowData(SingleCellExperiment)), value = TRUE)

    # Make the markers matrix

    markers <- do.call(rbind, lapply(ks, function(k){
      k_cols <- grep(paste0('_', k, '_'), results_cols, value = TRUE)
      mark <- cbind(k = k, data.frame(gene = rownames(SingleCellExperiment)), rowData(SingleCellExperiment)[, k_cols, drop = FALSE])
      colnames(mark)[c(3,4,5)] <- c('cluster', 'padj', 'auroc')
      mark
    }))

    # Some rows have NAs in cluster etc. Not sure why, let's just take them out

    markers <- markers[! is.na(markers$cluster), ]

    # Output text file

    write.table(markers, file = opt$output_text_file, sep = "\t", row.names = FALSE, quote = FALSE, na='')
}

# Print introspective information
cat(capture.output(SingleCellExperiment), sep='\n')

# Output to a serialized R object
saveRDS(SingleCellExperiment, file = opt$output_object_file)
