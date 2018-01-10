## Helpers for downloaded objects

#' @name curatedTCGAData-helpers
#'
#' @title Helper functions for managing MultiAssayExperiment from
#' curatedTCGAData
#'
#' @aliases getSubtypeMap
#'
#' @description
#' Additional helper functions for cleaning and uncovering metadata
#' within a downloaded \code{MultiAssayExperiment} from \code{curatedTCGAData}.
#' The \code{getSubtypeMap} function provides a 2 column \code{data.frame}
#' with in-data variable names and an interpreted names. The
#' \code{getClinicalNames} function provides a vector of variable names that
#' exist in the \code{colData} slot of a downloaded \code{MultiAssayExperiment}
#' object. These variables are obtained from
#' \link[RTCGAToolbox]{getFirehoseData} by default and tend to be present
#' across most cancer codes.
#'
#' @param multiassayexperiment A \linkS4class{MultiAssayExperiment} object
#'
#' @examples
#' \dontrun{
#' library(curatedTCGAData)
#'
#' coad <- curatedTCGAData(diseaseCode = "COAD",
#'     assays = "CNA*", dry.run = FALSE)
#' getSubtypeMap(coad)
#' }
#'
#' @return \itemize{
#'     \item{getSubtypeMap}: A \code{data.frame} with columns representing
#'     actual data variables and explanatory names
#'     \item{getClinicalNames}: A \code{vector} of names that correspond to
#'     a particular disease code.
#' }
#'
#' @export
getSubtypeMap <- function(multiassayexperiment) {
    if (!is(multiassayexperiment, "MultiAssayExperiment"))
        stop("Provide a 'MultiAssayExperiment' object")
    frameMap <- metadata(colData(multiassayexperiment))[["subtypes"]]
    if (is.null(frameMap)) {
        message("No subtype data available")
    } else frameMap
}

#' @rdname curatedTCGAData-helpers
#'
#' @param diseaseCode A TCGA cancer code (e.g., "BRCA")
#' @examples
#' getClinicalNames("COAD")
#'
#' @export
getClinicalNames <- function(diseaseCode) {
    stopifnot(S4Vectors::isSingleString(diseaseCode))
    env <- new.env(parent = emptyenv())
    data("clinicalNames", envir = env)
    clinNames <- env[["clinicalNames"]]
    clinNames[[diseaseCode]]
}
