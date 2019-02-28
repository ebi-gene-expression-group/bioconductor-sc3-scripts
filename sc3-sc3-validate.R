#!/usr/bin/env Rscript 
# Script validate SC3 docker container.

#####
# Data sources
#####

## input data for demo and validation
#validation_input_source <- paste0("https://github.com/hemberg-lab/",
#  "scRNA.seq.course/raw/master/deng/deng-reads.rds")
validation_input_source <- paste0("https://scrnaseq-public-datasets",
  ".s3.amazonaws.com/scater-objects/deng-reads.rds")
validation_input <- "deng-reads.rds"

## reference output for validation
## TEMPORARY LOCATION for script development
validation_reference_source <- paste0("https://raw.githubusercontent.com/",
  "HumanCellAtlas/analysis-tools-registry-supplemental/master/",
  "sc3/reference_output.csv")

# Download reference test data for validation 
cat("***Loading test data\n")
download.file(validation_input_source, validation_input)

# Make sure the reference test data exists 
if (!file.exists(validation_input)){
    stop(paste0("Reference test data, ", validation_input,
                ", failed to download. Please try again."))
}

#####
# Load test data
#####
sce <- readRDS(validation_input)
#check reference test data is expected dimensions
if (! all( dim(sce) == c(22431, 268))) {
    stop(paste0("Unexpected content in reference test data file ",
                 validation_input, " - validation aborted."))
}

ref <- read.table(validation_reference_source, header = FALSE,
                  row.names = 1, sep = ",")
#check reference test data is expected dimensions
if (! all( dim(ref) == c(268, 19))) {
    stop(paste0("Unexpected content in reference output file at ",
                validation_reference_source, " - validation aborted."))
}

#####
# Load libraries
#####
suppressPackageStartupMessages(library("SC3"));
suppressPackageStartupMessages(library("SingleCellExperiment"));

#####
# Calculate results for validation
#####
#Validation occurs on output from last deterministic step of SC3 (sc3_calc_transfs)
#set seed and limit cores to 1 to generate reproducible output
sce <- sc3_prepare(sce, rand_seed = 1, n_cores = 1)
sce <- sc3_calc_dists(sce)
sce <- sc3_calc_transfs(sce)

# Clean up downloaded input file
file.remove(validation_input)

write_output_files <- function(result, reference) {
  transfs_result <- cbind(colnames(result),
                    metadata(result)$sc3$transformations$pearson_laplacian)
  write.table(transfs_result, file = "validation_output.csv",
            row.names = FALSE, col.names = FALSE,
            quote = FALSE, sep = ",")
  write.table(reference, file = "reference.csv",
            row.names = TRUE, col.names = FALSE,
            quote = FALSE, sep = ",")
  if ( file.exists("validation_output.csv") & file.exists("reference.csv") ) {
    return(0)
  } else {
    stop("Failed to write output files")
  }
}

#####
# Perform validation
#####

cat("***Validate cell names in output\n")
###confirm row names match between reference and validation output
#if cell name mis-match, make validation & reference file available and exit
if (! isTRUE(all.equal(colnames(sce), rownames(ref)))) {
  write_output_files(sce, ref)
  stop("Validation failed - output cells do not match reference")
}

cat("***Validate output values\n")
###check if values in reference and validation output pass similarity check
#coerce the reference data in preparation for calculating relative error
test_matrix <- as.matrix(ref)
#check whether largest relative error is larger than acceptable threshold
if (max ( abs(as.matrix(metadata(sce)$sc3$transformations$pearson_laplacian)
      - test_matrix) / abs(test_matrix)) < 1.0e-8) {
  cat("***Successful validation - output passes similarity check\n")
#if max relative error exceeds threshold,
#make validation and reference files available and exit
} else {
  write_output_files(sce, ref)
  stop("Validation failed - max relative error exceeds threshold")
}