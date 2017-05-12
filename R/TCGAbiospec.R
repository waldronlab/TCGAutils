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
#' @author Marcel Ramos \email{marcel.ramos@roswellpark.org}
#'
#' @examples
#' example("TCGAbarcode")
#' TCGAbiospec(barcodes)
#'
#' @export TCGAbiospec
TCGAbiospec <- function(barcodes) {
    if (!all(nchar(barcodes) == 28L)) {
        stop("inconsistent barcode lengths")
    }
    stopifnot(all(startsWith(toupper(barcodes), "TCGA")))
    filler <- unique(substr(barcodes, 5, 5))
    if (length(filler) != 1L) {
        stop("barcode delimiters not consistent")
    }
    local_data_store <- new.env(parent = emptyenv())
    data("sampleTypes", envir = local_data_store)
    sampleTypes <- local_data_store[["sampleTypes"]]
    sample_type <- sampleTypes[["Definition"]][
        match(substr(TCGAbarcode(barcodes, participant = FALSE, sample = TRUE), 1, 2),
              sampleTypes[["Code"]])]
    newBiospec <- data.frame(patient_id = TCGAbarcode(barcodes),
                             sample_type,
                             sample_code = TCGAbarcode(barcodes,
                                                       participant = FALSE,
                                                       sample = TRUE),
                             stringsAsFactors = FALSE)
    portPlateCent <- data.frame(
        do.call(rbind,
                strsplit(TCGAbarcode(barcodes, index = 5:7), filler)),
        stringsAsFactors = FALSE)
    names(portPlateCent) <- c("portion", "plate", "center")
    newBiospec <- cbind(newBiospec, portPlateCent)
    return(newBiospec)
}
