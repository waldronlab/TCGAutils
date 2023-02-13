#' @rdname TCGAutils-defunct
#'
#' @title Defunct functions in TCGAutils
#' 
#' @description These functions have been moved to either `MultiAssayExperiment`
#'   or `RTCGAToolbox`
#' 
#' @aliases splitAssays makeSummarizedExperimentFromGISTIC
#'
#' @return Defunct functions return errors
#'
#' @export
splitAssays <- function(...) {
    .Defunct("TCGAsplitAssays")
}

#' @name TCGAutils-defunct
#' 
#' @details `splitAssays` was moved to `MultiAssayExperiment` and the TCGA
#'   specific implementation was named `TCGAsplitAssays` which resided in this
#'   package. `makeSummarizedExperimentFromGISTIC` now resides in
#'   `RTCGAToolbox`. See there for details.
#'
#' @param gistic An `FirehoseGISTIC` object from `RTCGAToolbox`
#' @param dataType Either one of "ThresholdedByGene", "AllByGene", "Peaks"
#' @param ... Additional arguments passed to 'RTCGAToolbox::getGISTICPeaks'.
#'
#' @author L. Geistlinger, M. Ramos
#'
#' @export
makeSummarizedExperimentFromGISTIC <- function(gistic, dataType, ...) {
    .Defunct("RTCGAToolbox::makeSummarizedExperimentFromGISTIC")
}