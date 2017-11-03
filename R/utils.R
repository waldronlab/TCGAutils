## Helper for finding barcode column
## **Takes the first result!**
.findBarcodeCol <- function(DF) {
    apply(DF, 2, function(column) {
        logicBCode <- grepl("^TCGA", column)
        logicBCode
    }) %>% apply(., 2, all) %>% Filter(isTRUE, .) %>% names %>% `[[`(1L)
}

## Standardize barcode format
standardBarcodes <- function(sampleBarcode) {
    if (!length(sampleBarcode)) {
        stop("<internal> Barcode must be of positive length")
    }
    sampleBC <- base::sample(sampleBarcode, 10L, replace = TRUE)
    bcodeTest <- grepl("\\.", sampleBC)
    if (all(bcodeTest))
        sampleBarcode <- gsub("\\.", "-", sampleBarcode)
    toupper(sampleBarcode)
}

## Find columns that are all NA
.findNAColumns <- function(dataset) {
    apply(dataset, 2, function(column) {
        all(is.na(column))
    })
}

## Helpers for downloaded objects

#' @name TCGA-helpers
#'
#' @title Helper functions for managing MultiAssayExperiment from
#' curatedTCGAData
#'
#' @description
#' Additional helper functions for cleaning and uncovering metadata
#' within a downloaded \code{MultiAssayExperiment} from \code{curatedTCGAData}
#'
#' @aliases getSubtypeMap
#'
#' @param multiassayexperiment A \linkS4class{MultiAssayExperiment} object
#' @export
getSubtypeMap <- function(multiassayexperiment) {
    frameMap <- metadata(colData(multiassayexperiment))[["subtypes"]]
    if (is.null(frameMap)) {
        message("No data available")
    } else frameMap
}

#' @rdname TCGA-helpers
#'
#' @param diseaseCode A TCGA cancer code (e.g., "BRCA")
#' @export
getClinicalNames <- function(diseaseCode) {
    stopifnot(S4Vectors::isSingleString(diseaseCode))
    env <- new.env(parent = emptyenv())
    data("clinicalNames", envir = env)
    clinNames <- env[["clinicalNames"]]
    clinNames[[diseaseCode]]
}
