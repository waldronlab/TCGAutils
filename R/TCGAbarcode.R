#' @import S4Vectors Biobase GenomicRanges
#' @importFrom plyr mapvalues
#' @importFrom readr read_delim
NULL

#' Parse data from TCGA barcode
#'
#' This function returns the specified snippet of information obtained from
#' the TCGA barcode.
#'
#' @param barcodes A character vector of TCGA barcodes
#' @param participant Logical (default TRUE) participant identifier chunk
#' @param sample Logical (default FALSE) includes the numeric sample code of
#' the barcode
#' @param vial Logical (default FALSE) includes the sample vial label
#' @param portion Logical (default FALSE) includes the portion and analyte
#' codes of the barcode
#' @param plate Logical (default FALSE) returns the plate value
#' @param center Logical (default FALSE) returns a matrix with the plate and
#' center codes
#' @param index A numerical vector of TCGA barcode positions desired
#'
#' @return A character vector or data matrix of TCGA barcode information
#'
#' @author Marcel Ramos \email{mramos09@gmail.com}
#'
#' @examples
#' library(TCGAbiolinks)
#' query <- TCGAquery(tumor = "OV", level = 3)
#' barcodes <- unlist(strsplit(query$barcode[1], ","))
#' TCGAsamplebc <- TCGAbarcode(barcodes)
#'
#' @export TCGAbarcode
TCGAbarcode <- function(barcodes, participant = TRUE,
                        sample = TRUE, vial = FALSE,
                        portion = FALSE, plate = FALSE,
                        center = FALSE, index = NULL)
{
    if (!all(nchar(barcodes) == 28L)) {
        stop("inconsistent barcode lengths")
    }
    stopifnot(all(startsWith(toupper(barcodes), "TCGA")))

    filler <- unique(substr(barcodes, 5, 5))
    if (length(filler) != 1L)  stop("barcode delimiters not consistent")

    barcodeMat <- do.call(rbind, strsplit(barcodes, "-"))
    if (!vial) {
        barcodeMat <- cbind(barcodeMat, substr(barcodeMat[, 4], 3, 3))
        barcodeMat[, 4] <- substr(barcodeMat[, 4], 1, 2)
    }
    if (is.null(index)) {
        if (participant) index <- c(index, 1:3)
        if (sample) index <- c(index, 4)
        if (portion) index <- c(index, 5)
        if (plate) index <- c(index, 6)
        if (center) index <- c(index, 7)
    }
    apply(barcodeMat[, index], 1, paste, collapse = filler)
}
