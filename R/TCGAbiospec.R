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
#' example("TCGAbarcode")
#' right_bctable <- TCGAbiospec(barcodes)
#'
#' @export TCGAbiospec
TCGAbiospec <- function(barcodes) {
  sample_type <- sampleTypes[["Definition"]][
      match(TCGAbarcode(barcodes, sample = TRUE, part = FALSE),
            sampleTypes[["Code"]])]
  vialSlice <- TCGAbarcode(barcodes, part=FALSE, sample = TRUE, vial = TRUE)
  tb <- data.frame(patientids = TCGAbarcode(barcodes),
                   sample_type,
                   sample_code = as.character(
                     TCGAbarcode(barcodes, part=FALSE, sample = TRUE)),
                   vial = ifelse(nchar(vialSlice) < 3, NA_character_,
                                 substr(vialSlice, 3,3)),
                   portion = as.character(
                     substr(TCGAbarcode(barcodes, part = FALSE, portion = TRUE),
                            1,2)),
                   analyte = as.character(
                     substr(TCGAbarcode(barcodes, part = FALSE, portion = TRUE),
                            3,3)),
                   plate = as.character(
                     TCGAbarcode(barcodes, part = FALSE, plate = TRUE)),
                   center = as.character(
                     TCGAbarcode(barcodes, part = FALSE, center=TRUE)),
                   stringsAsFactors = FALSE, row.names = NULL)
  tb <- tb[, apply(!is.na(tb), 2, all)]
  return(tb)
}
