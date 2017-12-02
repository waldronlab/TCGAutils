.find_with_xfix <- function(df_colnames, xfix1, xfix2,
        start.field, end.field, xfixType = "pre") {
    fixint <- intersect(xfix1, xfix2)
    fixint <- fixint[fixint != ""]
    if (!S4Vectors::isSingleString(fixint))
        stop("'start.field' and 'end.field' ", xfixType, "fixes do not match")
    names(fixint) <- xfixType

    fixFUN <- switch(xfixType, pre = I, suf = rev)
    start.field <- paste(fixFUN(c(fixint, start.field)), collapse = "")
    validEnd <- vapply(end.field, function(efield)
        paste(fixFUN(c(fixint, efield)), collapse = "") %in% df_colnames,
        logical(1L))
    stopifnot(sum(validEnd) == 1L)
    end.field <- paste(fixFUN(c(fixint, end.field[validEnd])), collapse = "")
    if (!length(start.field) && !length(end.field))
        list(c(start.field = "", end.field = ""), "")
    else
    list(c(start.field = start.field, end.field = end.field), fixint)
}

## Helper functions
.find_start_end_cols <- function (df_colnames, start.field, end.field) {
    idx1 <- which(df_colnames %in% start.field)
    idx2 <- which(df_colnames %in% end.field)
    prefixes1 <- .collect_prefixes(df_colnames, start.field)
    prefixes2 <- .collect_prefixes(df_colnames, end.field)
    suffixes1 <- .collect_suffixes(df_colnames, start.field)
    suffixes2 <- .collect_suffixes(df_colnames, end.field)
    if (length(idx1) == 1L && length(idx2) == 1L) {
        return(list(c(start = idx1, end = idx2), ""))
    }
    if (length(idx1) != 1L && length(prefixes1) ||
        length(idx2) != 1L && length(prefixes2)) {
    startend.fields <- .find_with_xfix(df_colnames, prefixes1, prefixes2,
        start.field, end.field, "pre")
    idx1 <- which(df_colnames %in% startend.fields[[1L]][["start.field"]])
    idx2 <- which(df_colnames %in% startend.fields[[1L]][["end.field"]])
    }
    if (!length(idx1) && !length(idx2)) {
    startend.fields <- .find_with_xfix(df_colnames, suffixes1, suffixes2,
        start.field, end.field, "suf")
    idx1 <- which(df_colnames %in% startend.fields[[1L]][["start.field"]])
    idx2 <- which(df_colnames %in% startend.fields[[1L]][["end.field"]])
    }
    if (length(idx1) == 1L && length(idx2) == 1L) {
        list(c(start = idx1, end = idx2), startend.fields[2L])
    } else {
        list(c(start = NA_integer_, end = NA_integer_), "")
    }
}

.collect_prefixes <- function (df_colnames, field) {
    df_colnames_nc <- nchar(df_colnames)
    prefixes <- lapply(field, function(suf) {
        pref_nc <- df_colnames_nc - nchar(suf)
        idx <- which(substr(df_colnames, pref_nc + 1L, df_colnames_nc) ==
                         suf)
        substr(df_colnames[idx], 1L, pref_nc[idx])
    })
    unique(unlist(prefixes))
}

.collect_suffixes <- function(df_colnames, field) {
    suffixes <- lapply(field, function(pre) {
        idx <- which(startsWith(df_colnames, pre))
        substr(df_colnames[idx], nchar(field) + 1L,
               nchar(df_colnames[idx]))
    })
    unique(unlist(suffixes))
}

.find_strands_col <- function(df_colnames, strand.field, xfix) {
    fixFUN <- switch(names(xfix[[1]]), pre = I, suf = rev)
    idx <- which(df_colnames %in%
        paste(fixFUN(c(xfix, strand.field)), collapse = ""))
    if (length(idx) == 0L)
        idx <- which(df_colnames %in% strand.field)
    if (length(idx) == 0L)
        return(NA_integer_)
    if (length(idx) >= 2L) {
        warning("Multiple strand measurements detected, taking first one")
        idx <- idx[[1L]]
    }
    idx
}

.find_seqnames_col <- function (df_colnames, seqnames.field, xfix) {
    fixFUN <- switch(names(xfix[[1]]), pre = I, suf = rev)
    idx <- which(df_colnames %in%
        paste(fixFUN(c(xfix, seqnames.field)), collapse = ""))
    if (length(idx) == 0L)
        idx <- which(df_colnames %in% seqnames.field)
    if (length(idx) == 0L)
        return(NA_integer_)
    if (length(idx) >= 2L)
        warning("cannnot determine seqnames column unambiguously")
        return(idx[[1L]])
    idx
}

.find_width_col <- function (df_colnames, width.field, xfix) {
    fixFUN <- switch(names(xfix[[1]]), pre = I, suf = rev)
    idx <- which(df_colnames %in%
        paste(fixFUN(c(xfix, width.field)), collapse = ""))
    if (length(idx) == 0L)
        idx <- which(df_colnames %in% width.field)
    if (length(idx) == 0L)
        return(NA_integer_)
    if (length(idx) >= 2L) {
        warning("cannnot determine width column unambiguously")
        return(idx[[1L]])
    }
    idx
}

#' Obtain minimum necessary names for the creation of a GRangesList object
#'
#' This function attempts to match chromosome, start position, end position and
#' strand names in the given character vector. Modified helper from the
#' \code{GenomicRanges} package.
#'
#' @param df_colnames A \code{character} vector of names in a dataset
#' @param seqnames.field A \code{character} vector of the chromosome name
#' @param start.field A \code{character} vector that indicates the column name
#' of the start positions of ranged data
#' @param end.field A \code{character} vector that indicates the end position
#' of ranged data
#' @param strand.field A \code{character} vector of the column name that
#' indicates the strand type
#' @param ignore.strand logical (default FALSE) whether to ignore the strand
#' field in the data
#' @return Index positions vector indicating columns with appropriate names
#'
#' @examples
#' myDataColNames <- c("Start_position", "End_position", "strand",
#'                  "chromosome", "num_probes", "segment_mean")
#' findGRangesCols(myDataColNames)
#'
#' @export findGRangesCols
findGRangesCols <- function (df_colnames,
    seqnames.field = c("seqnames", "seqname", "chromosome",
                    "chrom", "chr", "chromosome_name", "seqid"),
    start.field = "start",
    end.field = c("end", "stop"),
    strand.field = "strand",
    ignore.strand = FALSE) {

    df_colnames0 <- tolower(df_colnames)
    seqnames.field0 <- GenomicRanges:::.normarg_field(seqnames.field, "seqnames")
    start.field0 <- GenomicRanges:::.normarg_field(start.field, "start")
    end.field0 <- GenomicRanges:::.normarg_field(end.field, "end")
    start_end_cols <- .find_start_end_cols(df_colnames0, start.field0,
                                           end.field0)
    xfix <- start_end_cols[[2L]]
    width_col <- .find_width_col(df_colnames0, "width", xfix)
    seqnames_col <- .find_seqnames_col(df_colnames0, seqnames.field0, xfix)
    if (ignore.strand) {
        strand_col <- NA_integer_
    } else {
        strand.field0 <- GenomicRanges:::.normarg_field(strand.field, "strand")
        strand_col <- .find_strands_col(df_colnames0, strand.field0, xfix)
    }
    c(seqnames = seqnames_col, start_end_cols[[1L]], width = width_col,
      strand = strand_col)
}

.find_col <- function(df_colnames, field, xfix = "pre") {
    FUN <- switch(xfix, pre = I, suf = rev)
    name.field <- paste0(FUN(c(field, suf)), collapse = "")
    idx <- which(df_colnames %in% name.field)
    if (!length(idx))
        idx <-  which(df_colnames %in% field)
    if (!length(idx))
        stop("cannnot find ", field, " column")
    if (length(idx) >= 2L)
        stop("cannnot determine seqnames column unambiguously")
    idx
}
