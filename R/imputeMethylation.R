#' @name imputeMethylation
#'
#' @title This function imputes Methylation values inside a \code{MultiAssayExperiment}
#'
#' @description These function allow the user to enter a \code{multiAssayExperiment} and impute all
#' the NA values inside Methylation.
#'
#' @param multiassayexperiment A \code{MultiAssayExperiment} with genes in the rows, samples in the columns
#' @param i A \code{vector} of indices indicating the position inside \code{MultiAssayExperiment} of the
#' Methylation experiments.
#' @inheritDotParams impute::impute.knn
#'
#' @return MultiAssayExperiment with imputed Methylation values
#'
#' @examples
#' library(curatedTCGAData)
#'
#' mae <- curatedTCGAData("GBM", "Methylation", FALSE)
#'
#' for (i in seq_along(experiments(mae))) {
#'   newmat <- assay(mae[[i]])
#'   mode(newmat) <- "numeric"
#'   mae[[i]] <- newmat
#' }
#'
#' results <- imputeMethylation(mae, c(1,2))
#'
#' @export
imputeMethylation <- function(multiassayexperiment, i, ...) {
    if (!requireNamespace("impute"))
        stop("Install the 'impute' package to run 'imputeMethylation'")

    if (!is(multiassayexperiment, "MultiAssayExperiment"))
        stop("Methylation values has to be a matrix")
    if (!any(is.character(i), is.numeric(i), is.logical(i)))
        stop("'i' has to be character or integer or logical")

    sub.multiassayexperiment <- multiassayexperiment[,,i]
    assays <- assays(sub.multiassayexperiment)
    data.imputed <- lapply(assays, function(mat) {impute.knn(mat, ...)$data})

    for (x in i) {
        multiassayexperiment[[x]] <- data.imputed[[x]]
    }

    return(multiassayexperiment)
}
