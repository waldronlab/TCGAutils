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
