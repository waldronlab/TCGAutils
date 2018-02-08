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

#' @rdname curatedTCGAData-helpers
#'
#' @param sampleCodes A string of sample type codes
#' (refer to \code{data(sampleTypes)}; default "01", "11")
#'
#' @section splitAssays:
#'     Separates samples by indicated sample codes into different assays
#'     in a \code{MultiAssayExperiment}. Refer to the \code{sampleTypes}
#'     data object for a list of available codes. This operation generates
#'     \strong{n} times the number of assays based on the number of sample codes
#'     entered. By default, primary solid tumors ("01") and solid tissue
#'     normals ("11") are seperated out.
#' @export
splitAssays <- function(multiassayexperiment, sampleCodes = c("01", "11")) {
    if (!is(multiassayexperiment, "MultiAssayExperiment"))
        stop("Provide a 'MultiAssayExperiment' object")

    env <- new.env(parent = emptyenv())
    data("sampleTypes", envir = env)
    sampleTypes <- env[["sampleTypes"]]
    if (!sampleCodes %in% sampleTypes[["Code"]] || !is.character(sampleCodes))
        stop("Provide valid sample types string")

    cnames <- colnames(multiassayexperiment)
    exps <- experiments(multiassayexperiment)

    egroups <- lapply(sampleCodes, function(scode) {
        logitype <- BiocGenerics::relist(TCGAsampleSelect(
            unlist(cnames, use.names = FALSE), scode), cnames)
        explist <- subsetByColumn(exps, logitype)
        names(explist) <- paste0(scode, "_", names(explist))
        explist
    })
    egroups <- do.call(c, egroups)
    sampmap <- generateMap(egroups, colData(multiassayexperiment),
        idConverter = TCGAbarcode)
    BiocGenerics:::replaceSlots(multiassayexperiment,
        ExperimentList = egroups,
        sampleMap = sampmap)
}

#' @rdname curatedTCGAData-helpers
#' @param vial (logical default FALSE) whether to display vials in the
#' table output
#'
#' @section sampleTables:
#'     Display all the available samples in each of the assays
#' @export
sampleTables <- function(multiassayexperiment, vial = FALSE) {
    lapply(colnames(multiassayexperiment), function(x) {
        scodes <- TCGAbarcode(x, participant = FALSE, sample = TRUE)
        if (!vial)
            scodes <- substr(scodes, 1L, 2L)
        table(unname(scodes))
   })
}
