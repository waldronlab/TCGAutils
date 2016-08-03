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
#' @source \url{https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm}
"diseaseCodes"
