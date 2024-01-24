#' Select primary tumors from TCGA datasets
#'
#' Tumor selection is decided using the `sampleTypes` data. For 'LAML' datasets,
#' the primary tumor code used is "03" otherwise, "01" is used.
#'
#' @param multiassayexperiment A
#'   [`MultiAssayExperiment`][MultiAssayExperiment::MultiAssayExperiment-class]
#'   with TCGA data as obtained from [curatedTCGAData::curatedTCGAData()]
#'
#' @return A `MultiAssayExperiment` containing only primary tumor samples
#'
#' @examples
#'
#' example(getSubtypeMap)
#'
#' TCGAprimaryTumors(gbm)
#'
#' @export TCGAprimaryTumors
TCGAprimaryTumors <- function(multiassayexperiment) {
    if (!is(multiassayexperiment, "MultiAssayExperiment"))
        stop("Provide a 'MultiAssayExperiment' object as input")

    exptnames <- names(experiments(multiassayexperiment))
    dcodes <- vapply(strsplit(exptnames, "_"), `[[`, character(1L), 1L)

    primaries <- ifelse(dcodes == "LAML", "03", "01")
    primaries <- setNames(primaries, dcodes)

    logisub <- Map(function(barcodes, tumorcode) {
        TCGAsampleSelect(barcodes, tumorcode)
    }, colnames(multiassayexperiment), primaries)

    multiassayexperiment[, logisub, ]
}
