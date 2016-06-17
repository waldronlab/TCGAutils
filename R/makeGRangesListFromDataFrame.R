#' Make a GRangesList object from a data.frame or DataFrame
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
#' @section Value:
#' A \linkS4class{GRangesList} class object.
#'
#' Refer to the \link{makeGRangesFromDataFrame} documentation for more details.
#'
#' @examples
#' ##----
#' ## BASIC EXAMPLES
#' ##----
#'
#' df <- data.frame(chr="chr1", start=11:15, end=12:16,
#'                  strand=c("+","-","+","*","."), score=1:5,
#'                  specimen = c("a", "a", "b", "b", "c"))
#' df
#' makeGRangesListFromDataFrame(df, partitioning.field = "specimen")
#'
#' ##----
#' ## Adding feature names
#' ##----
#' rownames(df) <- paste0("GENE", 1:5)
#' makeGRangesListFromDataFrame(df, partitioning.field = "specimen")
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
