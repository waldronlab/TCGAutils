#' Take a MultiAssayExperiment and include curated variables
#'
#' This function works on the \code{colData} of a
#' \code{\linkS4class{MultiAssayExperiment}} object to merge curated variable
#' columns or other clinical variables that would like to be added. It is
#' recommended that the user run the scripts in the
#' \code{MultiAssayExperiment-TCGA} repository that build the "enhanced" type
#' of data but not necessary if using different clinical data.
#' Please see the repository's README for more information.
#'
#' @param MultiAssayExperiment A \linkS4class{MultiAssayExperiment} object
#' @param colData A \code{DataFrame} or {data.frame} to merge with
#' clinical data in the MultiAssayExperiment object
#'
#' @return A \code{\link{MultiAssayExperiment}} object
#'
#' @examples
#'
#' library(MultiAssayExperiment)
#'
#' mergeColData(MultiAssayExperiment(), S4Vectors::DataFrame())
#'
#' @export mergeColData
mergeColData <- function(MultiAssayExperiment, colData) {
    if (!is(MultiAssayExperiment, "MultiAssayExperiment"))
        stop("Provide a valid MultiAssayExperiment object")
    if (!is(colData, "DataFrame") && !is.data.frame(colData))
        stop("'colData' must be 'DataFrame' or 'data.frame'")
    if (is.null(rownames(colData)) && length(colData))
        stop("'colData' data must have rownames")

    maeClinical <- colData(MultiAssayExperiment)
    mergedClin <- merge(maeClinical, colData,
        by = c("row.names", intersect(names(maeClinical), names(colData))),
        all = TRUE, sort = FALSE, stringsAsFactors = FALSE)

    rownames(mergedClin) <- mergedClin[["Row.names"]]
    mergedClin <- mergedClin[, names(mergedClin) != "Row.names", drop = FALSE]
    colData(MultiAssayExperiment) <- as(mergedClin, "DataFrame")
    MultiAssayExperiment
}

#' Minimize the number of variables in colData
#'
#' This function removes variables that have a high number of missing data
#' and contain keywords.
#'
#' @param multiassayexperiment A \linkS4class{MultiAssayExperiment} object with colData
#' @param maxNAfrac (numeric default 0.2) A decimal between 0 and 1 to indicate
#' the amount of NA values allowed per column
#' @param keystring (character) A vector of keywords to match and remove
#' variables
#'
#' @return A \code{\linkS4class{MultiAssayExperiment}} object
#'
#' @examples
#'
#' library(curatedTCGAData)
#'
#' gbm <- curatedTCGAData("GBM", "CNASNP", FALSE)
#' (gbm_trimmed <- trimColData(gbm))
#' head(colData(gbm_trimmed))[1:5]
#'
#' @export trimColData
trimColData <- function(multiassayexperiment, maxNAfrac = 0.2,
    keystring = c("portion", "analyte")) {
    if (!is(multiassayexperiment, "MultiAssayExperiment"))
        stop("Provide a 'MultiAssayExperiment' input")
    DF <- colData(multiassayexperiment)
    keystring <- na.omit(keystring)

    NAabove <- vapply(DF, function(x) mean(is.na(x)) >= maxNAfrac, logical(1L))

    keymat <- vapply(keystring, function(string)
        grepl(string, names(DF)), logical(length(DF)))
    keymatch <- apply(keymat, 1L, any)

    todrop <- NAabove | keymatch
    colData(multiassayexperiment) <- DF[, !todrop]

    multiassayexperiment
}
