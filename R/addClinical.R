#' Take a MultiAssayExperiment and include curated variables
#'
#' This function works on the \code{colData} of a
#' \code{\linkS4class{MultiAssayExperiment}} object to add curated variable
#' columns. It is recommended that the user run the scripts in the
#' \code{MultiAssayExperiment-TCGA} repository that build the "enhanced" type
#' of data. Please see the repository's README for more information.
#'
#' @param MultiAssayExperiment A \linkS4class{MultiAssayExperiment} object
#' @param clinicalData A \code{DataFrame} or {data.frame} to merge with
#' clinical data in the MultiAssayExperiment object
#'
#' @return A \code{\link{MultiAssayExperiment}} object
#'
#' @export addClinical
addClinical <- function(MultiAssayExperiment, clinicalData) {
    if (!is(MultiAssayExperiment, "MultiAssayExperiment"))
        stop("Provide a valid MultiAssayExperiment object")
    if (!is(clinicalData, "DataFrame") && !is.data.frame(clinicalData))
        stop("Clinical data must be 'DataFrame' or 'data.frame'")

    maeClinical <- colData(MultiAssayExperiment)
    mergedClin <- merge(maeClinical, clinicalData,
        by = "row.names", all = TRUE, sort = FALSE, stringsAsFactors = FALSE)
    rownames(mergedClin) <- mergedClin[["Row.names"]]
    mergedClin <- mergedClin[, names(mergedClin) != "Row.names", drop = FALSE]
    colData(MultiAssayExperiment) <- as(mergedClin, "DataFrame")
    MultiAssayExperiment
}
