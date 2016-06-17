#' Make a RangedSummarizedExperiment from a data.frame or DataFrame
#'
#' \code{makeRangedSummarizedExperimentFromDataFrame} uses \code{data.frame} or
#' \code{DataFrame} column names to create a \link{GRanges} object for the
#' \code{rowRanges} of the resulting \link{SummarizedExperiment} object.
#' It requires that non-range data columns be coercible into a \code{numeric}
#' \code{matrix} for the \link{SummarizedExperiment} constructor.
#'
#' @inheritParams GenomicRanges::makeGRangesFromDataFrame
#'
#' @return A \link{SummarizedExperiment} object with rowRanges
#'
#' @export makeRangedSummarizedExperimentFromDataFrame
makeRangedSummarizedExperimentFromDataFrame <-
    function(df,
             seqnames.field = c(
                 "seqnames", "seqname", "chromosome",
                 "chrom", "chr", "chromosome_name", "seqid"),
             start.field = "start",
             end.field = c("end", "stop"),
             strand.field = "strand",
             ignore.strand = FALSE,
             seqinfo = NULL,
             starts.in.df.are.0based = FALSE)
    {
        ## feature.field not needed if DF has rownames
        rowRanges <- makeGRangesFromDataFrame(
            df,
            keep.extra.columns = FALSE,
            ignore.strand = ignore.strand,
            seqinfo = seqinfo,
            seqnames.field = seqnames.field,
            start.field = start.field,
            end.field = end.field,
            strand.field = strand.field,
            starts.in.df.are.0based = starts.in.df.are.0based)

        granges_cols <-
            GenomicRanges:::.find_GRanges_cols(names(df),
                                               seqnames.field = seqnames.field,
                                               start.field = start.field,
                                               end.field = end.field,
                                               strand.field = strand.field,
                                               ignore.strand = ignore.strand)

        droppedColumns <- names(df)[na.omit(granges_cols)]
        idx <- match(droppedColumns, names(df))
        counts <- as.matrix(df[, -idx, drop = FALSE])

        if (!is(as.vector(counts), "numeric"))
            stop("failed to coerce non-range columns to 'numeric'")

        SummarizedExperiment::SummarizedExperiment(
            assays=SimpleList(counts=counts), rowRanges=rowRanges)
    }
