ccle <- function (standard = TRUE, single = FALSE, targetCol = character(), idFUN = NULL) {
    if (standard) {
        rangeID <- c("chrom", "start", "end", "strand")
    }
    if (is.null(idFUN)) {
    cutCLN <- function(CLN) {
        gsub("^[^_]+_", "", CLN, perl = TRUE)
    }
    idFUN <- cutCLN
    }
    newArgs <- list(rangeID = rangeID, targetCol = targetCol, idFUN = idFUN)
    return(newArgs)
}
