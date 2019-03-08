#' Select samples from barcodes from lookup table
#'
#' The TCGA barcode contains several pieces of information which can
#' be parsed by the \link{TCGAbarcode} function. To select a specific type of
#' sample, enter the appropriate sampleCode argument from the lookup table.
#' See lookup table in \code{data("sampleTypes")}. Barcode inputs can be a
#' character vector or a \linkS4class{CharacterList} object.
#'
#' @param barcodes Either a TCGA barcode vector or \linkS4class{CharacterList}
#'   containing patient identifiers, sample, portion, plate, and center codes.
#' @param sampleCodes Either a character or numeric vector of TCGA sample codes.
#'   See the \code{sampleType} dataset.
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
    barc <- setNames(barcodeSamples %in% sampleCodes, barcodeSamples)
    if (exists("clist") && isTRUE(clist))
        barc <- utils::relist(barc, bcodes)
    return(barc)
}
