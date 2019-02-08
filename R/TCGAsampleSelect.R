#' Select samples from barcodes from lookup table
#'
#' The TCGA barcode contains several pieces of information which can
#' be parsed by the \link{TCGAbarcode} function. To select a specific type of
#' sample, enter the appropriate sampleCode argument from the lookup table.
#' See lookup table in \code{data("sampleTypes")}.
#'
#' @param barcodes Standard TCGA barcodes containing patient identifiers,
#' sample, portion, plate, center codes.
#' @param sampleCodes Either a character or numeric vector of TCGA sample codes.
#' See the \code{sampleType} dataset.
#'
#' @return A logical vector of the same length as 'barcodes' indicating matches
#'
#' @examples
#'
#' example("TCGAbarcode")
#' TCGAsampleSelect(barcodes, c(11, 01))
#'
#' @export TCGAsampleSelect
TCGAsampleSelect <- function(barcodes, sampleCodes) {
    stopifnot(
        is.character(sampleCodes) || is.numeric(sampleCodes),
        !is.na(sampleCodes), !is.logical(sampleCodes)
    )
    bcodeCharLen <- unique(nchar(barcodes))
    if (!S4Vectors::isSingleNumber(bcodeCharLen))
        warning("Inconsistent barcode lengths: ",
            paste(bcodeCharLen, collapse = ", "))
     if (any(bcodeCharLen < 15L))
        stop("'barcodes' should be at least 15 characters ",
                "with sample information")
    local_data_store <- new.env(parent = emptyenv())
    data("sampleTypes", envir = local_data_store, package = "TCGAutils")
    sampleTypes <- local_data_store[["sampleTypes"]]

    ncodes <- nchar(as.character(sampleCodes))
    singles <- ncodes == 1L
    sampleCodes[singles] <- paste0("0", sampleCodes[singles])

    if (any(!sampleCodes %in% sampleTypes[["Code"]]))
        stop("At least one sample code not found in look-up table")

    sampleSnippet <- TCGAbarcode(barcodes, sample = TRUE, participant = FALSE)
    barcodeSamples <- substr(sampleSnippet, 1L, 2L)
    return(setNames(barcodeSamples %in% sampleCodes, barcodeSamples))
}
