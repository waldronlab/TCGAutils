#' @name imputeAssay
#'
#' @title This function imputes assays values inside a \code{MultiAssayExperiment}
#'
#' @description These function allow the user to enter a \code{MultiAssayExperiment} and impute all
#' the NA values inside assays.
#'
#' @param multiassayexperiment A \code{MultiAssayExperiment} with genes in the rows, samples in the columns
#' @param i A \code{vector} of indices indicating the position inside \code{MultiAssayExperiment} of the
#' assays experiments to impute, default = 1.
#' @inheritDotParams impute::impute.knn
#'
#' @return MultiAssayExperiment with imputed assays values
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
#' results <- imputeAssay(mae, c(1,2))
#'
#' @export
imputeAssay <- function(multiassayexperiment, i = 1, ...) {
    if (!requireNamespace("impute"))
        stop("Install the 'impute' package to run 'imputeAssay'")

    if (!is(multiassayexperiment, "MultiAssayExperiment"))
        stop("Input has to be a MultiAssayExperiment")
    if (!any(is.character(i), is.numeric(i), is.logical(i)))
        stop("'i' has to be character or integer or logical")

    sub.multiassayexperiment <- multiassayexperiment[,,i]
    assays <- assays(sub.multiassayexperiment)
    data.imputed <- lapply(assays, function(mat) {impute::impute.knn(mat, ...)$data})

    for (x in i) {
        multiassayexperiment[[x]] <- data.imputed[[x]]
    }

    return(multiassayexperiment)
}
