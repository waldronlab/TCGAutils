#' Select samples from barcodes from lookup table
#'
#' The TCGA barcode contains several pieces of information which can
#' be parsed by the \link{TCGAbarcode} function. To select a specific type of
#' sample, enter the appropriate sampleCode argument from the lookup table.
#' See lookup table in \code{data("sampleTypes")}.
#'
#' @param barcodes Standard TCGA barcodes containing patient identifiers,
#' sample, portion, plate, center codes.
#' @param sampleCode Either a character or numeric vector of length one. See
#' the \code{sampleType} dataset.
#'
#' @return A logical vector of the same length as 'barcodes' indicating matches
#'
#' @examples
#'
#' example("TCGAbarcode")
#' TCGAsampleSelect(barcodes, 11)
#'
#' @export TCGAsampleSelect
TCGAsampleSelect <- function(barcodes, sampleCode) {
    stopifnot(
        is.character(sampleCode) || is.numeric(sampleCode),
        identical(length(sampleCode), 1L),
        !is.na(sampleCode), !is.logical(sampleCode)
    )
    bcodeCharLen <- unique(nchar(barcodes))
    if (!S4Vectors::isSingleNumber(bcodeCharLen))
        warning("Inconsistent barcode lengths: ",
            paste(bcodeCharLen, collapse = ", "))
     if (any(bcodeCharLen < 15L))
        stop("'barcodes' should be at least 15 characters ",
                "with sample information")
    local_data_store <- new.env(parent = emptyenv())
    data("sampleTypes", envir = local_data_store)
    sampleTypes <- local_data_store[["sampleTypes"]]
    if (is.numeric(sampleCode)) {
        if (identical(nchar(as.character(sampleCode)), 1L))
            sampleCode <- paste0("0", sampleCode)
    }
    if (!sampleCode %in% sampleTypes[["Code"]])
        stop("'sampleCode' not in look up table")
    LongName <- as.character(
        sampleTypes[match(sampleCode, sampleTypes[["Code"]]), "Definition"]
    )
    message("Selecting '", LongName, "' samples")
    sampleSnippet <- TCGAbarcode(barcodes, sample = TRUE, participant = FALSE)
    barcodeSamples <- substr(sampleSnippet, 1L, 2L)
    return(barcodeSamples == sampleCode)
}
