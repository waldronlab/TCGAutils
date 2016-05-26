#' Set TCGA data parameters
#'
#' This function is for use inside the \link{makeGRangesList} function for
#' converting raw data. The parameters in this function indicate properties
#' of the TCGA data.
#'
#' @param primary The column name for primary samples/specimen, defaults to
#' "Tumor_Sample_Barcode" (case ignored).
#' @param standard logical (default TRUE) whether to use standard range column
#' indicators found in the tcga data, these include "chromosome",
#' "start_position", "end_position", and "strand"
#' @param rangeID An optional \code{character} vector of length four indicating
#' genomic positions when standard argument is FALSE
#' @param idFUN A function that helps in the parsing of barcode identifiers
#' (TCGAbarcode used by default)
#'
#' @return A \code{list} of arguments for use by \link{makeGRangesList}
#'
#' @author Marcel Ramos \email{mramos09@gmail.com}
#' @export tcga
tcga <- function (primary = NULL, standard = TRUE,
                  rangeID = character(), idFUN = NULL) {
    if (is.null(primary)) {
        primary <- "Tumor_Sample_Barcode"
    }
    if (standard) {
        rangeID <- c("chromosome", "start_position", "end_position", "strand")
    } else {
        if (!is.character(rangeID) || !(length(rangeID) %in% 3L:4L)) {
            stop("range indicator vector length should be between 3 and 4")
        }
    }
    if (is.null(idFUN)) {
        idFUN <- function(x, ...) {
            TCGAbarcode(x, ...)
        }
    }
    newArgs <- list(primary = primary, rangeID = rangeID, idFUN = idFUN)
    return(newArgs)
}
