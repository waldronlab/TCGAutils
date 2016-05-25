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
#' @param position A numerical vector of TCGA barcode positions desired
#' @param part Logical (default TRUE) participant identifier chunk
#' @param sample Logical (default FALSE) includes the numeric sample code of
#' the barcode
#' @param vial Logical (default FALSE) includes the sample vial label
#' @param portion Logical (default FALSE) includes the portion and analyte
#' codes of the barcode
#' @param plate Logical (default FALSE) returns the plate value
#' @param center Logical (default FALSE) returns a matrix with the plate and
#' center codes
#' @param collapse Logical (default FALSE) concatenates requested barcode
#' chunks
#' @return A character vector or data matrix of TCGA barcode information
#'
#' @author Marcel Ramos \email{mramos09@gmail.com}
#'
#' @examples
#' library(TCGAbiolinks)
#' query <- TCGAquery(tumor = "OV", level = 3)
#' barcodes <- unlist(strsplit(query$barcode[1], ","))
#' TCGAsamplebc <- TCGAbarcode(barcodes, sample = TRUE, collapse = TRUE)
#'
#' @export TCGAbarcode
TCGAbarcode <- function(barcodes, position = NULL, part = TRUE, sample=FALSE,
                        vial = FALSE, portion = FALSE, plate = FALSE,
                        center=FALSE, collapse=FALSE) {
    if(!is.character(barcodes)){
        stop("Barcodes must be a character vector!")
    } else {
        randIndices <- sample(length(barcodes), 10, replace = TRUE)
        if (!all(grepl("^TCGA", barcodes[randIndices], ignore.case = TRUE))) {
            stop("All barcodes must start with TCGA!")
        }
    }
    filler <- substr(barcodes[1], 5,5)
    if (filler != "-") {
        barcodes <- gsub(paste0("\\", filler), "-", barcodes)
    }
    if(is.null(position)) {
        index <- NULL
        if(part) index <- c(index, 1:3)
        if(sample) index <- c(index, 4)
        if(portion) index <- c(index, 5)
        if(plate) index <- c(index, 6)
        if(center) index <- c(index, 7)
        if(is.null(index)){
            stop("Enter a barcode portion to output!")
        }
    } else {
        if(any(position > 7)) {
            stop("Please enter a valid number or range!")
        }
        index <- position
    }
    if(identical(index, 1:3)){
        bcc <- lapply(strsplit(barcodes, "\\-"), "[", index)
        bcc <- tolower(sapply(bcc, paste, collapse="-"))
        return(bcc)
    }
    if(collapse) {
        bcc <- lapply(strsplit(barcodes, "\\-"), "[", index)
        bcc <- tolower(sapply(bcc, paste, collapse = "-"))
        if(!vial & max(index)==4){
            if(part){
                lettNo <- 15
            } else {
                lettNo <- 2
            }
            bcc <- substr(bcc, 1, lettNo)
        }
    } else {
        bcc <- sapply(strsplit(barcodes, "\\-"), "[", index)
        if(length(index) > 1){
            bcc <- tolower(t(bcc))
        }
        if(!vial & max(index)==4 & !part){
            bcc <- substr(bcc, 1, 2)
        }
    }
    return(bcc)
}
