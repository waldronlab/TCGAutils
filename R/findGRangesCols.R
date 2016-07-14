#' Obtain minimum necessary names for the creation of a GRangesList object
#'
#' This function attempts to match chromosome, start position, end position and
#' strand names in the given character vector.
#'
#' @param namesVector A \code{character} vector or names in a dataset
#'
#' @return Index positions vector indicating columns with appropriate names
#'
#' @examples
#' myDataColNames <- c("Start_position", "End_position", "strand",
#'                  "chromosome", "num_probes", "segment_mean")
#' findGRangesCols(myDataColNames)
#'
#' @export findGRangesCols
findGRangesCols <- function (df_colnames,
                             seqnames.field = c("seqnames", "seqname",
                                                "chromosome", "chrom", "chr",
                                                "chromosome_name", "seqid"),
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
    prefix <- start_end_cols[[2L]]
    width_col <- GenomicRanges:::.find_width_col(df_colnames0, "width", prefix)
    seqnames_col <- .find_seqnames_col(df_colnames0, seqnames.field0,
                                       prefix)
    if (ignore.strand) {
        strand_col <- NA_integer_
    }
    else {
        strand.field0 <- GenomicRanges:::.normarg_field(strand.field, "strand")
        strand_col <- GenomicRanges:::.find_strand_col(df_colnames0, strand.field0,
                                                       prefix)
    }
    c(seqnames = seqnames_col, start_end_cols[[1L]], width = width_col,
      strand = strand_col)
}

.find_start_end_cols <-function (df_colnames, start.field, end.field) {
    idx1 <- which(df_colnames %in% start.field)
    idx2 <- which(df_colnames %in% end.field)
    if (length(idx1) == 1L && length(idx2) == 1L)
        return(list(c(start = idx1, end = idx2), ""))
    if (length(idx1) == 0L && length(idx2) == 0L) {
        prefixes1 <- .collect_prefixes(df_colnames, start.field)
        prefixes2 <- .collect_prefixes(df_colnames, end.field)
        suffixes1 <- .collect_suffixes(df_colnames, start.field)
        suffixes2 <- .collect_suffixes(df_colnames, end.field)
        if (length(prefixes1) != 0L && length(prefixes2) != 0L) {
            if (length(prefixes1) >= 2L && length(prefixes2) >= 2L) {
                warning("multiple prefixes found, using first match")
                if (prefixes1[[1L]] == prefixes2[[1L]])
                    prefix <- prefixes1[[1L]]
            } else if (length(prefixes1) == 1L && length(prefixes2) == 1L &&
                       prefixes1 == prefixes2) {
                prefix <- prefixes1
            }
            idx1 <- which(df_colnames %in% paste0(prefix, start.field))
            idx2 <- which(df_colnames %in% paste0(prefix, end.field))
            if (length(idx1) == 1L && length(idx2) == 1L)
                return(list(c(start = idx1, end = idx2), prefix))
        } else if (length(suffixes1) != 0L && length(suffixes2) != 0L) {
            if (length(suffixes1) >= 2L && length(suffixes2) >= 2L) {
                warning("multiple suffixes found, using first match")
                if (suffixes1[[1L]] == suffixes2[[1L]])
                    suffix <- suffixes1[[1L]]
            } else if (length(suffixes1) == 1L && length(suffixes2) == 1L &&
                       suffixes1 == suffixes2) {
                suffix <- suffixes1
            }
            idx1 <- which(df_colnames %in% paste0(start.field, suffix))
            idx2 <- which(df_colnames %in% paste0(end.field, suffix))
            if (length(idx1) == 1L && length(idx2) == 1L)
                return(list(c(start = idx1, end = idx2), ""))
        } else {
            stop("cannnot determine start/end columns")
        }
    }
}

.find_seqnames_col <- function (df_colnames, seqnames.field, prefix) {
    idx <- which(df_colnames %in% paste0(prefix, seqnames.field))
    if (length(idx) == 0L)
        idx <- which(df_colnames %in% seqnames.field)
    if (length(idx) == 0L)
        stop("cannnot find seqnames column")
    if (length(idx) >= 2L)
        warning("cannnot determine seqnames column unambiguously")
        return(idx[[1L]])
    idx
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
