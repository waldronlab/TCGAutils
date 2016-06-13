makeGRangesListFromTCGA <- function(x, primary, idFUN = NULL, ...)
{
    if (is.list(x) && !inherits(x, "data.frame"))
        x <- do.call(rbind, x)

    partitioning <- x[[primary]]
    x <- x[, -match(primary, names(x)), drop = FALSE]

    twoMeta <- all(c("num_probes", "segment_mean") %in% tolower(names(x)))
    hugo <- "hugo_symbol" %in% tolower(names(x))

    grl <- makeGRangesListFromDataFrame(df, idFUN(partitioning), ...)

    if (hugo) {
        hugoName <- names(x)[which("hugo_symbol" %in% tolower(names(x)))]
        names(unlist(grl)) <- x[[hugoName]]
    }
    if (twoMeta) {
        numProb <- names(x)[which("num_probes" %in% tolower(names(x)))]
        segMean <- names(x)[which("segment_mean" %in% tolower(names(x)))]
        browser()
        mcols(grl) <- cbind(mcols(grl), DataFrame(num_probes = numProb,
                                                  segment_mean = segMean))
    }
    grl
}
