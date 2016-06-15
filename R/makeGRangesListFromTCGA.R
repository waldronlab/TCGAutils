.setHGBuild <- function(hgbuild) {
    buildDF <- DataFrame(Date = c("July 2004", "May 2004", "March 2006",
                                  "February 2009"),
                         NCBI = c("34", "35", "36", "37"),
                         UCSC = c("hg16", "hg17", "hg18", "hg19"))
    buildIndex <- match(hgbuild, buildDF[["NCBI"]])
    if (is.na(buildIndex)) {
        warning("build could not be matched")
        return(NA_character_)
    } else {
        ucscBuild <- buildDF$UCSC[match(hgbuild, buildDF[["NCBI"]])]
        return(ucscBuild)
    }
}

makeGRangesListFromTCGA <-
    function(x, primary,
             feature.field = "hugo_symbol",
             idFUN = NULL, ...)
{
    if (is.list(x) && !inherits(x, "data.frame"))
        x <- do.call(rbind, x)

    twoMeta <- all(c("num_probes", "segment_mean") %in% tolower(names(x)))
    hugo <- feature.field %in% tolower(names(x))
    ncbi <- "ncbi_build" %in% tolower(names(x))

    if (is.null(idFUN)) {
        idFUN <- function(bcode) {
            TCGAbarcode(bcode, sample = TRUE, collapse = TRUE)
        }}

    grl <- makeGRangesListFromDataFrame(x, primary, ...)
    names(grl) <- idFUN(names(grl))

    if (hugo) {
        hugoName <- names(x)[which(feature.field %in% tolower(names(x)))]
        # names with unlist assignment not working TODO
        # names(BiocGenerics::unlist(grl, use.names = FALSE)) <- x[[hugoName]]
    }

    if (twoMeta) {
        numProb <- names(x)[which("num_probes" %in% tolower(names(x)))]
        segMean <- names(x)[which("segment_mean" %in% tolower(names(x)))]
        mcols(grl) <- cbind(mcols(grl), DataFrame(num_probes = numProb,
                                                  segment_mean = segMean))
    }
    if (ncbi) {
        ncbi_build <- names(x)[which("ncbi_build" %in% tolower(names(x)))]
        build_name <- unique(x[[ncbi_build]])
        if (length(build_name) != 1L) {
            warning("inconsistent ncbi_build values in data")
        } else {
            ucscBuild <- .setHGBuild(build_name)
            GenomeInfoDb::genome(grl) <- ucscBuild
        }
    }
    grl
}
