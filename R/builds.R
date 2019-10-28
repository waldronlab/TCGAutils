#' @name builds
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
#' @param from character() A vector of build versions typically from `genome()`
#'     (e.g., "37"). The build vector must be homogenous (i.e.,
#'     `length(unique(x)) == 1L`).
#'
#' @param to character(1) The name of the desired build version (either "UCSC"
#'     or "NCBI")
#'
#' @examples
#'
#' translateBuild("GRCh35", "UCSC")
#'
#' @return
#'     translateBuild: A character vector of translated genome builds
#'
#'     extractBuild: A character string of the build information available
#'
#'     uniformBuilds: A character vector of builds where all builds are
#'         identical `identical(length(unique(build)), 1L)`
#' @export
translateBuild <- function(from, to = "UCSC") {
    lfro <- length(from)
    from <- unique(from)

    if (!identical(length(from), 1L))
        stop("Enter a consistent vector of genomic builds")
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

    build <-
        if (to == "NCBI")
            paste0("GRCh", buildDF[[to]][buildIndex])
        else
            buildDF[[to]][buildIndex]
    rep(build, lfro)
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
#' navec <- c(rep(c("GRCh37", "hg19"), times = c(5, 1)), "NA")
#' uniformBuilds(navec)
#'
#' @export uniformBuilds
uniformBuilds <- function(builds, cutoff = 0.2) {
    if (length(unique(builds)) == 1L)
        return(builds)
    wbuilds <- tolower(builds)
    nabuilds <- wbuilds == "na"
    wbuilds[nabuilds] <- NA_character_
    tots <- length(wbuilds)
    nabuilds <- is.na(wbuilds)
    propna <- sum(nabuilds) / tots
    if (propna > cutoff)
        stop("Frequency of NA values higher than the cutoff tolerance")
    wbuilds <- na.omit(wbuilds)
    ubuilds <- unique(tolower(wbuilds))

    if (identical(length(ubuilds), 1L)) {
        builds[nabuilds] <- unique(builds[!nabuilds])
        return(builds)
    } else if (length(ubuilds) > 2)
        stop("Only two build types at a time can be used")

    names(ubuilds) <- unique(wbuilds)
    props <- vapply(ubuilds, function(biobuild) {
        sum(biobuild == wbuilds) / tots
    }, numeric(1L))

    offbuild <- names(props[props < cutoff])
    mainbuild <- names(props[props > cutoff])
    mainbuild <- builds[match(mainbuild, tolower(builds))]
    if (any(nabuilds))
        builds[nabuilds] <- mainbuild

    results <- Filter(function(x) !is.na(x),
        lapply(c("NCBI", "UCSC"), function(buildfmt) {
            suppressWarnings(translateBuild(offbuild, buildfmt))
        })
    )
    builds[wbuilds == offbuild] <- unlist(results)
    builds
}

