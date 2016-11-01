.parseFileNames <- function(filepaths) {
    fileuuids <- basename(dirname(filepaths))
    # fileNames <- basename(filepaths)
    # fileuuids <- sapply(strsplit(fileNames, "\\."), "[", 3)
    bcodes <- TCGAtranslateID(fileuuids)
    # extraInfo <- TCGAbiospec(bcodes$barcode)
    bcodes
}

#' Read Exon level files and create a GRangesList
#'
#' This function serves to read exon level data from a vector of file paths
#' and to create a \linkS4class{GRangesList} object for representation. File
#' name and structure requirements are as follows: The third atomic position
#' of the parsed file name should be the universally unique identifier (UUID).
#' File should use default separators (i.e., periods). The column containing
#' ranged information should be named "exon."
#'
#' @param filepaths A vector of valid exon data file paths
#' @param filenames logical (default FALSE) whether to attempt to parse and
#' translate the file names from exon files
#' @return A \linkS4class{GRangesList} object
#'
#' @author Marcel Ramos \email{mramos09@gmail.com}
#'
#' @export TCGAexonToGRangesList
TCGAexonToGRangesList <- function(filepaths, filenames=FALSE) {
    btData <- lapply(filepaths, function(file) {
        readr::read_delim(file, delim = "\t")
    })
    if (filenames) {
    sampNames <- try(.parseFileNames(filepaths))
    if (!is(sampNames, "try-error"))
        names(btData) <- sampNames
    }
    newGRL <- GenomicRanges::GRangesList(lapply(btData, function(range) {
        newGRanges <- GenomicRanges::GRanges(as.character(range[["exon"]]))
        mcols(newGRanges) <- range[, -(which(names(range) == "exon"))]
        newGRanges
    }))
    return(newGRL)
}
