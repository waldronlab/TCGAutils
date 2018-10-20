#' @rdname builds
#'
#' @title Utilities for working with build numbers
#'
#' @description A few functions are available to search for build versions,
#' either from NCBI or UCSC.
#'
#' \itemize{
#'   \item \code{translateBuild}: translates between UCSC and NCBI build
#'   versions
#'   \item \code{extractBuild}: use grep patterns to find the first build
#'   within the string input
#'   \item \code{uniformBuilds}: replace build occurrences below a threshold
#'   level of occurence with the alternative build
#' }
#'
#' @param from A build version name
#' @param to The name of the desired version
#'
#' @examples
#'
#' translateBuild("GRCh35", "UCSC")
#'
#' @export
translateBuild <- function(from, to = "UCSC") {
    if (!S4Vectors::isSingleString(to) && !S4Vectors::isSingleString(from))
        stop("Enter a single valid genomic build")
    if (!to %in% c("UCSC", "NCBI"))
        stop ("Only UCSC and NCBI supported")

    buildDF <- S4Vectors::DataFrame(
        Date = c("July 2004", "May 2004", "March 2006", "February 2009",
            "December 2013"),
        NCBI = c("34", "35", "36", "37", "38"),
        UCSC = c("hg16", "hg17", "hg18", "hg19", "hg38")
    )
    if (to == "UCSC")
        from <- gsub("[GgRrCcHh]", "", from)
    matchBuild <- switch(to, UCSC = "NCBI", NCBI = "UCSC")
    buildIndex <- match(from, buildDF[[matchBuild]])
    if (is.na(buildIndex)) {
        warning("build could not be matched")
        return(NA_character_)
    }
    if (to == "NCBI")
        paste0("GRCh", buildDF[[to]][buildIndex])
    else
        buildDF[[to]][buildIndex]
}

#' @rdname builds
#'
#' @param string A single character string
#' @param build A vector of build version names (default UCSC, NCBI)
#'
#' @examples
#'
#' extractBuild(
#' "SCENA_p_TCGAb29and30_SNP_N_GenomeWideSNP_6_G05_569110.nocnv_grch38.seg.txt"
#' )
#'
#' @export
extractBuild <- function(string, build = c("UCSC", "NCBI")) {
    if (!S4Vectors::isSingleString(string))
        stop("Provide a single string for build search")
    builds <- vector(mode = "character", length(build))
    names(builds) <- build
    for (i in build) {
        pattrn <- switch(i, UCSC = "[Hh][Gg][0-9]{2}",
            NCBI = "[Gg][Rr][Cc][Hh][0-9]{2}")
        builds[[i]] <- stringr::str_extract(string, pattrn)
    }
    builds <- Filter(function(x) !is.na(x), builds)
    if (!length(builds))
        NA_character_
    else if (length(builds))
        builds[1L]
}

#' @rdname builds
#'
#' @param builds A character vector of builds
#' @param cutoff A threshold value for translating builds below the threshold
#'
#' @examples
#'
#' buildvec <- rep(c("GRCh37", "hg19"), times = c(5, 1))
#' uniformBuilds(buildvec)
#'
#' @export uniformBuilds
uniformBuilds <- function(builds, cutoff = 0.2) {
    if (length(unique(builds)) == 1L)
        return(builds)
    wbuilds <- tolower(builds)
    tots <- length(wbuilds)
    ubuilds <- unique(tolower(wbuilds))
    if (length(ubuilds) > 2)
        stop("Only two build types at a time can be used")
    names(ubuilds) <- unique(builds)
    props <- vapply(ubuilds, function(biobuild) {
        sum(biobuild == wbuilds)/tots
    }, numeric(1L))
    props <- props[props < cutoff]
    offbuild <- names(props)
    results <- Filter(function(x) !is.na(x),
        lapply(c("NCBI", "UCSC"), function(buildfmt) {
            suppressWarnings(translateBuild(offbuild, buildfmt))
        })
    )
    builds[wbuilds == offbuild] <- unlist(results)
    builds
}
