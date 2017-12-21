.parseFileName <- function(filepath) {
    fileName <- basename(filepath)
    splitFileName <- strsplit(fileName, "\\.")
    fileuuid <- sapply(strsplit(fileName, "\\."), "[", 3L)
    if (length(strsplit(fileuuid, "-")[[1]]) != 5L)
        stop("Inconsistent UUID in file name")
    bcode <- TCGAtranslateID(fileuuid, type = "entity_id")
    # extraInfo <- TCGAbiospec(bcodes$barcode)
    bcode
}

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
#' @param filepaths A vector of valid exon data file paths
#' @param sampleNames A vector of TCGA barcodes to be applied if not present in
#' the data
#' @return A \linkS4class{GRangesList} object
#'
#' @importFrom GenomicRanges GRanges GRangesList
#'
#' @author Marcel Ramos
#'
#' @examples
#'
#' pkgDir <- system.file("extdata", package = "TCGAutils", mustWork = TRUE)
#' exonFile <- list.files(pkgDir, pattern = "cation.txt$", full.names = TRUE)
#' makeGRangesListFromExonFiles(exonFile)
#'
#' @export makeGRangesListFromExonFiles
makeGRangesListFromExonFiles <-
    function(filepaths, sampleNames = NULL, rangeCol = "exon") {
    btData <- lapply(filepaths, function(file) {
        read_delim(file, delim = "\t")
    })
    if (!is.null(sampleNames)) {
        if (length(filepaths) != length(sampleNames))
            stop("Inconsistent sample names obtained from file names")
    } else {
        sampleNames <- lapply(filepaths, .parseFileName)
    }
    names(btData) <- sampleNames
    GRangesList(
        lapply(btData, function(range) {
            newGRanges <- GRanges(as.character(range[[rangeCol]]))
            mcols(newGRanges) <- range[, names(range) != rangeCol]
            newGRanges
        })
    )
}
