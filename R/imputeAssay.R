#' @name imputeAssay
#'
#' @title This function imputes assays values inside a
#' \code{MultiAssayExperiment}
#'
#' @description These function allow the user to enter a
#' \code{MultiAssayExperiment} and impute all the NA values inside assays.
#'
#' @param multiassayexperiment A \code{MultiAssayExperiment} with genes in the
#' rows, samples in the columns
#' @param i A numeric, logical, or character \code{vector} indicating the
#' assays to perform imputation on (default 1L)
#' @inheritDotParams impute::impute.knn
#'
#' @return MultiAssayExperiment with imputed assays values
#'
#' @examples
#'
#' example(getSubtypeMap)
#'
#' ## convert data to matrix and add as experiment
#' gbm <-
#'   c(gbm, RPPA_matrix = data.matrix(assay(gbm[["GBM_RPPAArray-20160128"]])))
#'
#' imputeAssay(gbm, i = "RPPA_matrix")
#'
#' @export
imputeAssay <- function(multiassayexperiment, i = 1, ...) {
    if (!requireNamespace("impute"))
        stop("Install the 'impute' package to run 'imputeAssay'")

    if (!is(multiassayexperiment, "MultiAssayExperiment"))
        stop("Input has to be a MultiAssayExperiment")
    if (!any(is.character(i), is.numeric(i), is.logical(i)))
        stop("'i' has to be character or numeric or logical")

    sub.multiassayexperiment <- multiassayexperiment[,,i]
    assays <- assays(sub.multiassayexperiment)
    assayclasses <- vapply(assays, is.matrix, logical(1L))
    if (!all(assayclasses))
        stop("Only matrix assay(s) can be imputed")
    data.imputed <- lapply(assays, function(mat) {
        impute::impute.knn(mat, ...)$data
    })

    for (x in i) {
        multiassayexperiment[[x]] <- data.imputed[[x]]
    }

    return(multiassayexperiment)
}
