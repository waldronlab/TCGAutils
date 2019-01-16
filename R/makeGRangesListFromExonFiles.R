#' Read Exon level files and create a GRangesList
#'
#' This function serves to read exon-level expression data. It works for exon
#' quantification (raw counts and RPKM) and junction quantification
#' (raw counts) files paths and represent such data as a
#' \linkS4class{GRangesList}. The data can be downloaded
#' via the TCGA Legacy Archive. File name and structure requirements are as
#' follows: The third position delimited by dots (".") in the file name should
#' be the universally unique identifier (UUID). The column containing the
#' ranged information is labeled "exon."
#'
#' @param filepaths A \code{character} vector of valid exon data file paths
#' @param sampleNames A \code{character} vector of TCGA barcodes to be applied
#' if not present in the data (default NULL)
#' @param fileNames A \code{character} vector of file names as downloaded from
#' the Genomic Data Commons Legacy archive (default NULL)
#' @param rangesColumn (default "exon") A single string indicating the name of
#' the column in the data containing the ranges information
#'
#' @return A \linkS4class{GRangesList} object
#'
#' @author M. Ramos
#'
#' @examples
#'
#' ## Load example file found in package
#' pkgDir <- system.file("extdata", package = "TCGAutils", mustWork = TRUE)
#' exonFile <- list.files(pkgDir, pattern = "cation\\.txt$", full.names = TRUE)
#'
#' filePrefix <- "unc.edu.32741f9a-9fec-441f-96b4-e504e62c5362.1755371."
#'
#' ## Add actual file name manually (due to Windows OS restriction)
#' makeGRangesListFromExonFiles(exonFile,
#'     fileNames = paste0(filePrefix, basename(exonFile)),
#'     sampleNames = "TCGA-AA-3678-01A-01R-0905-07")
#'
#' @export makeGRangesListFromExonFiles
makeGRangesListFromExonFiles <- function(filepaths, sampleNames = NULL,
    fileNames = NULL, rangesColumn = "exon")
{
    if (!is.null(sampleNames)) {
        if (length(filepaths) != length(sampleNames))
            stop("Inconsistent sample names obtained from file names")
    } else {
        queryNames <-
            if (!is.null(fileNames)) fileNames else basename(filepaths)
        sampleNames <-
            filenameToBarcode(queryNames, TRUE)[["aliquots.submitter_id"]]
    }
    btData <- lapply(filepaths, function(file) {
        if (requireNamespace("readr", quietly = TRUE))
            readr::read_delim(file, delim = "\t")
        else
            read.delim(file, sep = "\t")
    })

    if (!length(sampleNames))
        sampleNames <- NULL

    names(btData) <- sampleNames

    allrowdata <- if (requireNamespace("dplyr", quietly = TRUE))
        dplyr::bind_rows(btData)
    else
        do.call(rbind, btData)

    newGRanges <- GRanges(allrowdata[[rangesColumn]])
    mcols(newGRanges) <- allrowdata[, names(allrowdata) != rangesColumn]
    GRangesList(relist(newGRanges, btData))
}
