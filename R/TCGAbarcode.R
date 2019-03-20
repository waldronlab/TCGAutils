.uniqueDelim <- function(ids) {
    unique(unlist(
        strsplit(gsub("[a-zA-Z0-9]", "", ids), "")
    ))
}

.checkBarcodes <- function(barcodes, check.sample = FALSE, minIndex = NULL) {
    if (!all(startsWith(toupper(barcodes), "TCGA")))
        stop("Barcodes must start with 'TCGA'")
    filler <- .uniqueDelim(barcodes)
    if (length(filler) != 1L)
        stop("Barcode delimiters not consistent")
    bcodelens <- unique(nchar(barcodes))
    if (length(bcodelens) > 1L)
        warning("Inconsistent barcode lengths: ",
            paste(bcodelens, collapse = ", "))
    if (check.sample) {
        if (any(bcodelens < 15L))
        stop("'barcodes' should be at least 15 characters ",
                "with sample information")
    }
    if (!is.null(minIndex)) {
        splitCodes <- strsplit(barcodes, filler)
        obsIdx <- unique(lengths(splitCodes))
        if (obsIdx < minIndex)
            stop("Minimum barcode fields required: ", minIndex,
                "; first three are 'project-TSS-participant'")
    }
}

#' Parse data from TCGA barcode
#'
#' This function returns the specified snippet of information obtained from
#' the TCGA barcode.
#'
#' @param barcodes A character vector of TCGA barcodes
#' @param participant Logical (default TRUE) participant identifier chunk
#' @param sample Logical (default FALSE) includes the numeric sample code of
#' the barcode and the vial letter
#' @param portion Logical (default FALSE) includes the portion and analyte
#' codes of the barcode
#' @param plate Logical (default FALSE) returns the plate value
#' @param center Logical (default FALSE) returns a matrix with the plate and
#' center codes
#' @param index A numerical vector of TCGA barcode positions desired when
#' split by the delimiter (i.e., hyphen '-')
#'
#' @return A character vector or data matrix of TCGA barcode information
#'
#' @author M. Ramos
#'
#' @examples
#' barcodes <- c("TCGA-B0-5117-11A-01D-1421-08",
#' "TCGA-B0-5094-11A-01D-1421-08",
#' "TCGA-E9-A295-10A-01D-A16D-09")
#'
#' ## Patient identifiers
#' TCGAbarcode(barcodes)
#'
#' ## Sample identifiers
#' TCGAbarcode(barcodes, sample = TRUE)
#'
#' @export TCGAbarcode
TCGAbarcode <- function(barcodes, participant = TRUE, sample = FALSE,
    portion = FALSE, plate = FALSE, center = FALSE, index = NULL)
{
    .checkBarcodes(barcodes)
    filler <- .uniqueDelim(barcodes)
    barcodeMat <- do.call(rbind, strsplit(barcodes, filler))
    if (is.null(index)) {
        if (participant) index <- c(index, 1L:3L)
        if (sample) index <- c(index, 4L)
        if (portion) index <- c(index, 5L)
        if (plate) index <- c(index, 6L)
        if (center) index <- c(index, 7L)
    }
    apply(barcodeMat[, index, drop = FALSE], 1L, paste, collapse = filler)
}

