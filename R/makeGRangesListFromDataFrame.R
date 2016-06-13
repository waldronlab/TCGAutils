#' Convert DataFrame to GRangesList
#'
#' \code{makeGRangesListFromDataFrame} extends the
#' \link[GenomicRanges]{makeGRangesFromDataFrame} functionality from
#' \code{GenomicRanges}. It returns a \link{GRangesList} object.
#'
#' @param df A \code{DataFrame} class object
#' @param partitioning.field Typically a \code{factor} vector that defines the
#' grouping order of the \code{DataFrame} rows
#' @inheritParams GenomicRanges::makeGRangesFromDataFrame
#'
#' @return A \linkS4class{GRangesList} class object
#'
#' @seealso GenomicRanges::makeGRangesFromDataFrame
#'
#' @export makeGRangesListFromDataFrame
makeGRangesListFromDataFrame <-
    function(df, partitioning.field,
             keep.extra.columns = FALSE,
             ignore.strand = FALSE,
             seqinfo = NULL,
             seqnames.field = c(
                 "seqnames", "seqname", "chromosome",
                 "chrom", "chr", "chromosome_name", "seqid"),
             start.field = "start",
             end.field = c("end", "stop"),
             strand.field = "strand",
             starts.in.df.are.0based = FALSE)
    {
        if (!S4Vectors::isSingleString(partitioning.field))
            stop("'partitioning.field' must be a single string")

        partitioningIdx <- match(partitioning.field, names(df))
        if (is.na(partitioningIdx))
            stop("'partitioning.field' is not in 'names(df)'")

        gr <- GenomicRanges::makeGRangesFromDataFrame(
            df[, -match(partitioning.field, names(df))],
            keep.extra.columns=keep.extra.columns,
            ignore.strand=ignore.strand,
            seqinfo = seqinfo,
            seqnames.field = seqnames.field,
            start.field = start.field,
            end.field = end.field,
            strand.field = strand.field,
            starts.in.df.are.0based = starts.in.df.are.0based)

        S4Vectors::split(gr, df[[partitioning.field]])
    }
