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
    if (!requireNamespace("RTCGAToolbox"))
        stop("Please install 'RTCGAToolbox' to use this function")
    gist <- RTCGAToolbox::getData(gistic, "GISTIC", dataType)
    rel.cols <- grepl("^TCGA", colnames(gist))
    gistData <- as.matrix(gist[, rel.cols])
    if (dataType == "Peaks") {
        gist <- RTCGAToolbox::getGISTICPeaks(gistic, ...)
        rel.cols <- grepl("^TCGA", colnames(gist))
        gistData <- as.matrix(gist[, rel.cols])
        rowranges <- gist[["rowRanges"]]
        rowranges <- as(rowranges, "GRanges")
        # get the peak type (amplification / deletion)
        peak.type <- vapply(strsplit(gist[["Unique.Name"]], " "),
            function(x) x[[1L]], character(1L))
        # create the SE
        rowdata <- cbind.data.frame(gist[, !rel.cols], type = peak.type,
            stringsAsFactors = FALSE)
        gisticSE <- SummarizedExperiment(gistData,
            rowRanges = rowranges)
        rowData(gisticSE) <- rowdata
    } else if (dataType %in% c("ThresholdedByGene", "AllByGene")) {
        colnames(gistData) <- .standardBarcodes(colnames(gistData))
        gisticSE <- SummarizedExperiment(gistData)
    }
    gisticSE
}
