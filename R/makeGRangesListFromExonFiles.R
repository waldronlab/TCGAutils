#' Read exon-level expression files and create a `GRangesList`
#'
#' This function serves to read exon-level expression data. It works for exon
#' quantification (raw counts and RPKM) and junction quantification
#' (raw counts) file paths and represents such data as a
#' \linkS4class{GRangesList}. The data files can be downloaded
#' via the Genomic Data Commons (GDC) Legacy Archive.
#'
#' @details The `rangesColumn` name in the GDC data files is usually "exon"
#'   but can be changed with the `rangesColumn` argument, if different.
#'   To avoid programmatically obtaining TCGA barcodes from the GDC
#'   API, set the `getBarcodes` to `FALSE`. When `getBarcodes` is set to
#'   `FALSE`, the file names are used to name the elements of the `GRangesList`
#'   output.
#'
#' @param filepaths character() vector of file paths containing TCGA exon
#'     data usually obtained from the GDC
#'
#' @param sampleNames character() vector of TCGA barcodes to be used as
#'     names for the `GRangesList` output (default NULL)
#'
#' @param fileNames character() vector of file names as downloaded from
#'     the Genomic Data Commons Legacy archive (default `basename(filepaths)`)
#'
#' @param getBarcodes logical(1). Whether to query the GDC API with the
#'     `filenameToBarcode` and obtain the TCGA barcodes from the file names
#'     (default TRUE); see details.
#'
#' @param rangesColumn character(1). The name of the column in the data
#'     containing the ranges information (default "exon"); see details.
#'
#' @param nrows numeric(1). The number of rows to return from each of the files
#'     read in (all rows by default; default Inf)
#'
#' @md
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
    fileNames = basename(filepaths), getBarcodes = TRUE, rangesColumn = "exon",
    nrows = Inf)
{
    if (is.null(sampleNames) && getBarcodes) {
        sampleNames <-
            filenameToBarcode(filenames = fileNames)[[
                "cases.samples.portions.analytes.aliquots.submitter_id"
            ]]
    } else if (is.null(sampleNames)) {
        sampleNames <- fileNames
    }

    if (!identical(length(filepaths), length(sampleNames)))
        stop("'sampleNames' length is inconsistent with 'fileNames'")

    btData <- lapply(filepaths, function(file) {
        if (requireNamespace("readr", quietly = TRUE)) {
            readr::local_edition(1)
            readr::read_delim(file, delim = "\t", n_max = nrows)
        } else
            read.delim(file, sep = "\t",
                nrows = if (is.infinite(nrows)) -1 else nrows)
    })

    names(btData) <- sampleNames

    allrowdata <-
        if (requireNamespace("dplyr", quietly = TRUE))
            dplyr::bind_rows(btData)
        else
            do.call(rbind, btData)

    newGRanges <- GenomicRanges::GRanges(allrowdata[[rangesColumn]])
    mcols(newGRanges) <- allrowdata[, names(allrowdata) != rangesColumn]

    splitIndx <- rep(names(btData), vapply(btData, nrow, integer(1L)))
    S4Vectors::splitAsList(newGRanges, splitIndx)
}
