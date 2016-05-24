#' Set TCGA data parameters
#'
#' This function is for use inside the \link{makeGRangesList} function for
#' converting raw data. The parameters in this function indicate properties
#' of the TCGA data.
#'
#' @param standard logical (default TRUE) whether to use standard range column
#' indicators found in the tcga data, these include "chromosome",
#' "start_position", "end_position", and "strand".
#' @param single logical (default FALSE) whether to use a single column for
#' the conversion of ranged values. Such column may look like: "chr1:100-200:+"
#' @param targetCol A \code{character} string that indicates the name of the
#' column to use when a single column is used (i.e., TRUE)
#' @param idFUN A function that helps in the parsing of barcode identifiers
#' (barcode used by default)
#'
#' @return A \code{list} of arguments for use by \link{makeGRangesList}
#'
#' @author Marcel Ramos \email{mramos09@gmail.com}
#' @export tcga
tcga <- function (standard = TRUE, single = FALSE, targetCol = character(), idFUN = NULL) {
    if (standard) {
        rangeID <- c("chromosome", "start_position", "end_position", "strand")
    }
    if (single) {
        rangeID <- NULL
        if (length(targetCol) == 0L) {
            stop("No target column specified")
        }
    }
    if (is.null(idFUN)) {
        idFUN <- function(x, ...) {
            barcode(x, ...)
        }
    }
    newArgs <- list(rangeID = rangeID,
                    targetCol = targetCol, idFUN = idFUN)
    return(newArgs)
}
