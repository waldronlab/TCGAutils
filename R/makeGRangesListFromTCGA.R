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
        ucscBuild <- buildDF$UCSC[buildIndex]
        return(ucscBuild)
    }
}

#' Make a GRangesList from TCGA data
#'
#' \code{makeGRangesListFromTCGA} allows the user to convert objects of class
#' data.frame or DataFrame to a \link{GRangesList}. It includes additional
#' features specific to TCGA data such as, hugo symbols, probe numbers,
#' segment means, and ucsc build (if available).
#'
#' @param x A \code{data.frame} or \code{DataFrame} class object. \code{list}
#' class objects are coerced to data.frame or DataFrame.
#' @param primary A \code{character} vector of length one indicating the
#' column to be used as sample identifiers
#' @param feature.field A \code{character} vector of length one indicating the
#' column to be used as names for each of the ranges in the data
#' @param ... Additional arguments to pass on to
#' \link{makeGRangesListFromDataFrame}
#'
#' @return A \link{GRangesList} class object
#'
#' @export makeGRangesListFromTCGA
makeGRangesListFromTCGA <-
    function(x, primary,
             feature.field = "hugo_symbol", ...)
    {
        if (is.list(x) && !inherits(x, "data.frame"))
            x <- do.call(rbind, x)

        if (!S4Vectors::isSingleString(primary))
            stop("'primary' must be a single sting")

        twoMeta <- all(c("num_probes", "segment_mean") %in% tolower(names(x)))
        hugo <- feature.field %in% tolower(names(x))
        ncbi <- "ncbi_build" %in% tolower(names(x))

        grl <- makeGRangesListFromDataFrame(x, primary, ...)

        if (hugo) {
            hugoName <- names(x)[which(feature.field %in% tolower(names(x)))]
            tempGRL <- BiocGenerics::unlist(grl)
            names(tempGRL) <- x[[hugoName]]
            grl <- BiocGenerics::relist(tempGRL, grl)
        }

        if (twoMeta) {
            numProb <- names(x)[match("num_probes", tolower(names(x)))]
            segMean <- names(x)[match("segment_mean", tolower(names(x)))]
            mcols(grl) <- cbind(mcols(grl), DataFrame(num_probes = numProb,
                                                      segment_mean = segMean))
        }
        if (ncbi) {
            ncbi_build <- names(x)[match("ncbi_build", tolower(names(x)))]
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
