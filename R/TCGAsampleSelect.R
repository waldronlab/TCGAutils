#' Select samples from barcodes from lookup table
#'
#' The TCGA barcode contains several pieces of information which can
#' be parsed by the \link{TCGAbarcode} function. To select a specific type of
#' sample, enter the appropriate sampleCode argument from the lookup table.
#' See lookup table in `data("sampleTypes")`. Barcode inputs can be a
#' character vector or a \linkS4class{CharacterList} object.
#'
#' @param barcodes Either a TCGA barcode vector or \linkS4class{CharacterList}
#'   containing patient identifiers, sample, portion, plate, and center codes.
#' @param sampleCodes Either a character or numeric vector of TCGA sample codes.
#'   See the `sampleType` dataset.
#'
#' @return A logical vector or \linkS4class{LogicalList} of the same length as
#'   'barcodes' indicating sample type matches
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
    if (clist <- is(barcodes, "CharacterList")) {
        bcodes <- barcodes
        barcodes <- unlist(barcodes, use.names = FALSE)
    }

    .checkBarcodes(barcodes, check.sample = TRUE)

    sampleCodes <- .addLeadingZero(sampleCodes)
    .checkSampleCodes(sampleCodes, strict = TRUE)

    sampleSnippet <- TCGAbarcode(barcodes, sample = TRUE, participant = FALSE)
    barcodeSamples <- substr(sampleSnippet, 1L, 2L)
    barc <- setNames(barcodeSamples %in% sampleCodes, barcodeSamples)
    if (exists("clist") && isTRUE(clist))
        barc <- BiocGenerics::relist(barc, bcodes)
    return(barc)
}
