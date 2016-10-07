#' Take a MultiAssayExperiment and include curated variables
#'
#' This function works on the \code{pData} of a
#' \code{\link{MultiAssayExperiment}} object to add curated variable columns.
#' Endomorphic operation!
#'
#' @param MultiAssayExperiment
#'
#' @return A \code{\link{MultiAssayExperiment}} object
enhanceClinical <- function(MultiAssayExperiment = MultiAssayExperiment(),
                            cancerCode, pkgLocation = ".") {
    stopifnot(inherits(MultiAssayExperiment, "MultiAssayExperiment"))
    stopifnot(S4Vectors::isSingleString(cancerCode))
    stopifnot(file.exists(file.path(pkgLocation, "MultiAssayExperiment-TCGA")))

    clinicalDF <- pData(MultiAssayExperiment)
    curatedSubtypes <-
        readr::read_csv(
            file.path(pkgLocation, "MultiAssayExperiment-TCGA",
                      "inst/extdata/curatedSubtypes/",
                      paste0(toupper(cancerCode), "_clinical.csv")))
}
