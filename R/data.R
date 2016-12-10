#' @import S4Vectors Biobase GenomicRanges methods BiocGenerics GenomeInfoDb
#' SummarizedExperiment
#' @importFrom readr read_delim
NULL

#' TCGA Cancer Disease Codes Table
#'
#' A dataset for obtaining the cancer codes in TCGA for about 13 different
#' types of cancers.
#'
#' @format A data frame with 37 rows and 2 variables:
#' \describe{
#'      \item{Study.Abbreviation}{Disease Code used in TCGA}
#'      \item{Study.Name}{The full length study name (i.e., type of cancer)}
#' }
#' @return The TCGA `diseaseCodes` table
#' @source \url{https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations}
"diseaseCodes"

#' Barcode Sample Type Table
#'
#' A dataset that contains the mappings for sample codes in the TCGA
#' barcodes.
#' @format A data frame with 19 rows and 3 variables:
#' \describe{
#'      \item{Code}{Two digit code number found in the barcode}
#'      \item{Definition}{Long name for the sample type}
#'      \item{Short.Letter.Code}{Letter code for the sample type}
#' }
#' @return The TCGA `sampleTypes` table
#' @source \url{https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes}
"sampleTypes"

