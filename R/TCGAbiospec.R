#' Extract biospecimen data from the TCGA barcode
#'
#' This function uses the full TCGA barcode to return a data frame of the
#' data pertinent to laboratory variables such as vials, portions, analytes,
#' plates and the center.
#'
#' @param barcodes A character vector of TCGA barcodes
#' @return A \code{dataframe} with sample type, sample code, portion, plate,
#' and center columns.
#'
#' @author M. Ramos
#'
#' @examples
#' example("TCGAbarcode")
#' TCGAbiospec(barcodes)
#'
#' @export TCGAbiospec
TCGAbiospec <- function(barcodes) {
    .checkBarcodes(barcodes)
    filler <- .uniqueDelim(barcodes)
    local_data_store <- new.env(parent = emptyenv())
    data("sampleTypes", envir = local_data_store)
    sampleTypes <- local_data_store[["sampleTypes"]]
    sampCode <- substr(TCGAbarcode(barcodes, FALSE, TRUE), 1L, 2L)
    sample_type <- sampleTypes[["Definition"]][
        match(sampCode, sampleTypes[["Code"]])]
    newBiospec <-
        data.frame(
            patient_id = TCGAbarcode(barcodes),
            sample_type,
            sample_code = TCGAbarcode(barcodes, FALSE, TRUE),
            stringsAsFactors = FALSE
        )
    portPlateCent <- do.call(rbind.data.frame,
        args = c(strsplit(TCGAbarcode(barcodes, index = 5L:7L), filler),
            list(stringsAsFactors = FALSE)))
    names(portPlateCent) <- c("portion", "plate", "center")
    cbind.data.frame(newBiospec, portPlateCent, stringsAsFactors = FALSE)
}
