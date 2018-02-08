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
        stop("Clinical data must be 'DataFrame' or 'data.frame'")

    maeClinical <- colData(MultiAssayExperiment)
    mergedClin <- merge(maeClinical, colData,
        by = "row.names", all = TRUE, sort = FALSE, stringsAsFactors = FALSE)
    rownames(mergedClin) <- mergedClin[["Row.names"]]
    mergedClin <- mergedClin[, names(mergedClin) != "Row.names", drop = FALSE]
    colData(MultiAssayExperiment) <- as(mergedClin, "DataFrame")
    MultiAssayExperiment
}
