## Helper for finding barcode column
## **Takes the first result!**
.findBarcodeCol <- function(DF) {
    cnames <- names(DF)
    containsBC <- vapply(head(DF), function(column) {
        all(grepl("^TCGA", column))
    }, logical(1L))
    names(containsBC) <- cnames
    bcIdx <- which(containsBC)
    stopifnot(S4Vectors::isSingleInteger(which(containsBC)))
    names(containsBC)[bcIdx]
}

.translateBuild <- function(fromBuild, toBuild = "UCSC") {
    buildDF <- S4Vectors::DataFrame(
        Date = c("July 2004", "May 2004", "March 2006", "February 2009",
            "December 2013"),
        NCBI = c("34", "35", "36", "37", "38"),
        UCSC = c("hg16", "hg17", "hg18", "hg19", "hg38")
    )
    matchBuild <- switch (toBuild, UCSC = "NCBI", NCBI = "UCSC" )
    buildIndex <- match(fromBuild, buildDF[[matchBuild]])
    if (is.na(buildIndex)) {
        warning("build could not be matched")
        return(NA_character_)
    }
    buildDF[[toBuild]][buildIndex]
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
