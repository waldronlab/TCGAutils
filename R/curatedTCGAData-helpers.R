## Helpers for downloaded objects

#' @name curatedTCGAData-helpers
#'
#' @title Helper functions for managing MultiAssayExperiment from
#' curatedTCGAData
#'
#' @description
#' Additional helper functions for cleaning and uncovering metadata
#' within a downloaded \code{MultiAssayExperiment} from \code{curatedTCGAData}
#'
#' @aliases getSubtypeMap
#'
#' @param multiassayexperiment A \linkS4class{MultiAssayExperiment} object
#' @export
getSubtypeMap <- function(multiassayexperiment) {
    frameMap <- metadata(colData(multiassayexperiment))[["subtypes"]]
    if (is.null(frameMap)) {
        message("No data available")
    } else frameMap
}

#' @rdname curatedTCGAData-helpers
#'
#' @param diseaseCode A TCGA cancer code (e.g., "BRCA")
#' @export
getClinicalNames <- function(diseaseCode) {
    stopifnot(S4Vectors::isSingleString(diseaseCode))
    env <- new.env(parent = emptyenv())
    data("clinicalNames", envir = env)
    clinNames <- env[["clinicalNames"]]
    clinNames[[diseaseCode]]
}
