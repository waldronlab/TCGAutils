.FileNamesToBarcodes <- function(filepaths) {
    fnames <- basename(filepaths)
    filesres <- files(legacy = TRUE)
    info <- results_all(
        select(filter(filesres, ~ file_name %in% fnames),
            "cases.samples.portions.analytes.aliquots.submitter_id")
    )
    id_list <- lapply(info[["cases"]], function(a) {
        a[[1L]][[1L]][[1L]]
    })
    # so we can later expand to a data.frame of the right size
    barcodes_per_file <- sapply(id_list,length)
    # And build the data.frame
    data.frame(file_name = rep(fnames, barcodes_per_file),
        file_id = rep(ids(info), barcodes_per_file),
        aliquots.submitter_id = unlist(id_list), row.names = NULL,
        stringsAsFactors = FALSE)
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
#' @author M. Ramos
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
    readr_avail <- requireNamespace("readr", quietly = TRUE)
    btData <- lapply(filepaths, function(file) {
        if (readr_avail)
            readr::read_delim(file, delim = "\t")
        else
            read.delim(file, sep = "\t")
    })
    if (!is.null(sampleNames)) {
        if (length(filepaths) != length(sampleNames))
            stop("Inconsistent sample names obtained from file names")
    } else {
        sampleNames <- .FileNamesToBarcodes(filepaths)[["aliquots.submitter_id"]]
    }
    if (!length(sampleNames))
        sampleNames <- NULL
    names(btData) <- sampleNames
    GRangesList(
        lapply(btData, function(range) {
            newGRanges <- GRanges(as.character(range[[rangeCol]]))
            mcols(newGRanges) <- range[, names(range) != rangeCol]
            newGRanges
        })
    )
}
