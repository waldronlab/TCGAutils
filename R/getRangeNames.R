#' Obtain minimum necessary names for the creation of a GRangesList object
#'
#' This function attempts to match chromosome, start position, end position and
#' strand names in the given character vector.
#'
#' @param namesVector A \code{character} vector or names in a dataset
#'
#' @return The 3 to 4 length \code{character} vector of minimum necessary names to
#' create a \link{GRangesList} object
#'
#' @examples
#' myDataColNames <- c("txStart_position", "txEnd_position", "strand", "chromo",
#' "num_probes", "segment_mean")
#' getRangeNames(myDataColNames)
#'
#' @export getRangeNames
getRangeNames <- function(namesVector) {
    if (!is.character(namesVector)) {
        stop("provide a character vector")
    }
    rangeNames <- c("^chr", "start", "end", "strand")
    AvailNames <- lapply(rangeNames, function(nam) {
        grepl(nam, namesVector, ignore.case = TRUE)
    })
    names(AvailNames) <- c("chr", "start", "end", "strand")
    timesMatched <- sapply(AvailNames, function(x) {sum(x)})
    multipleMatches <- which(timesMatched > 1)
    multiMatched <- names(AvailNames)[multipleMatches]
    grangeNames <- vapply(rangeNames, function(name) {
        grep(name, namesVector, value = TRUE, ignore.case = TRUE)
    }, FUN.VALUE = character(1L))
    if (any(multipleMatches)) {
        stop("'", multiMatched, "' had more than one range indicator match")
    } else if (length(grangeNames) != 4L) {
        warning("not all ranged indicators found")
        if (all(timesMatched[1:3] == 1L)) {
            return(grangeNames)
        } else {
            stop("minimum necessary range indicators not found")
        }
    } else {
        return(grangeNames)
    }
}
