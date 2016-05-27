#' Obtain minimum necessary names for the creation of a GRangesList object
#'
#' This function attempts to match chromosome, start position, end position and
#' strand names in the given character vector.
#'
#' @param namesVector A \code{character} vector or names in a dataset
#' @param regEx A \code{character} vector of regular expressions to
#' pass to \code{grepl} for detecting chromosome, start range, end range, and
#' strand column names, in that order. If empty, it will default to common use
#' cases of the names i.e.,\code{c("^chr", "start", "end", "strand")}.
#'
#' @return The 3 to 4 length \code{character} vector of minimum necessary names to
#' create a \link{GRangesList} object
#'
#' @examples
#' myDataColNames <- c("txStart_position", "txEnd_position", "strand",
#'                  "chromo", "num_probes", "segment_mean")
#' getRangeNames(myDataColNames)
#'
#' @export getRangeNames
getRangeNames <- function(namesVector, regEx = character()) {
    if (!is.character(namesVector)) {
        stop("provide a character vector of names")
    }
    if (length(regEx) == 0L) {
        rangeNames <- c("^chr", "start", "end", "strand")
    } else if (length(regEx) != 4L) {
        stop("regEx length is not 4")
    } else {
        rangeNames <- regEx
    }
    AvailNames <- lapply(rangeNames, function(nam) {
        grepl(nam, namesVector, ignore.case = TRUE)
    })
    names(AvailNames) <- c("chr", "start", "end", "strand")
    timesMatched <- vapply(AvailNames, sum, numeric(1L))
    multipleMatches <- which(timesMatched > 1)
    multiMatched <- names(AvailNames)[multipleMatches]
    grangeNames <- sapply(rangeNames, function(name) {
        grep(name, namesVector, value = TRUE, ignore.case = TRUE)
    })
    nonMatch <- vapply(grangeNames, function(name) {
        length(name) == 0L
    }, logical(1L))
    if (any(multipleMatches)) {
        stop("'", multiMatched, "' had more than one range indicator match")
    } else if (any(nonMatch)) {
        warning("not all ranged indicators found")
        if (all(timesMatched[1:3] == 1L)) {
            grangeNames <- BiocGenerics::Filter(function(string) {
                length(string) != 0L}, grangeNames)
            return(unlist(grangeNames))
        } else {
            stop("minimum necessary range indicators not found")
        }
    } else {
        return(grangeNames)
    }
}
