#' Create a SummarizedExperiment from FireHose GISTIC
#'
#' @description Use the output of \code{getFirehoseData} to create a
#' \linkS4class{SummarizedExperiment}. This can be done for three types of
#' data, G-scores thresholded by gene, copy number by gene, and copy number by
#' peak regions.
#'
#' @param gistic A \link[RTCGAToolbox]{FirehoseGISTIC-class} object
#' @param dataType Either one of "ThresholdedByGene", "AllByGene", "Peaks"
#' @param ... Additional arguments passed to 'RTCGAToolbox::getGISTICPeaks'.
#'
#' @author L. Geistlinger, M. Ramos
#'
#' @examples
#'
#' library(RTCGAToolbox)
#' co <- getFirehoseData("COAD", clinical = FALSE, GISTIC = TRUE,
#'     destdir = tempdir())
#' makeSummarizedExperimentFromGISTIC(co, "AllByGene")
#'
#' @return A \code{SummarizedExperiment} object
#' @export
makeSummarizedExperimentFromGISTIC <- function(gistic, dataType, ...) {
    .Deprecated("RTCGAToolbox::makeSummarizedExperimentFromGISTIC")
    if (!requireNamespace("RTCGAToolbox"))
        stop("Please install 'RTCGAToolbox' to use this function")
    RTCGAToolbox::makeSummarizedExperimentFromGISTIC(gistic, dataType, ...)
}
