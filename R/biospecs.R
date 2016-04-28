#' Extract biospecimen data from the TCGA barcode
#'
#' This function uses the full TCGA barcode to return a data frame of the
#' data pertinent to laboratory variables such as vials, portions, analytes,
#' plates and the center.
#'
#' @param barcodes A character vector of TCGA barcodes
#' @return A \code{dataframe} with sample type, sample code, vial, portion, analyte,
#' plate, and center columns.
#'
#' @author Marcel Ramos \email{mramos09@gmail.com}
#'
#' @examples
#'
#' \dontrun{
#' barcodes <- colnames(getElement(a2, "RNASeqGene"))
#' right_bctable <- biospecs(barcodes)
#' }
#'
#' @export biospecs
 biospecs <- function(barcodes) {
  sample_type <- as.character(
    sampleTypes[,"Definition"][
      match(as.numeric(barcode(barcodes, sample = TRUE, part = FALSE)),
            sampleTypes[,"Code"])])
  tb <- data.frame(patientids = barcode(barcodes),
                   sample_type,
                   sample_code = as.character(
                     barcode(barcodes, part=FALSE, sample = TRUE)),
                   vial = as.character(
                     substr(barcode(barcodes, part=FALSE, sample = TRUE,
                                  vial = TRUE), 3,3)),
                   portion = as.character(
                     substr(barcode(barcodes, part = FALSE, portion = TRUE),
                            1,2)),
                   analyte = as.character(
                     substr(barcode(barcodes, part = FALSE, portion = TRUE),
                            3,3)),
                   plate = as.character(
                     barcode(barcodes, part = FALSE, plate = TRUE)),
                   center = as.character(
                     barcode(barcodes, part = FALSE, center=TRUE)),
                   stringsAsFactors = FALSE)
  return(tb)
}
