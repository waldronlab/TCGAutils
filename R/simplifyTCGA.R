#' @importFrom GenomicFeatures genes microRNAs
#' @importFrom GenomeInfoDb keepStandardChromosomes seqlevelsStyle
#' seqlevelsStyle<-
NULL

.checkHas <-
    function(x, pattern = c("^hsa", "^cg", "symbols"), threshold = 0.9) {
    if (identical(pattern, "symbols"))
        pattern <- "^[A-Z0-9]{1,6}|^C[0-9]orf[0-9]{1,4}"
    mean(c(FALSE, grepl(pattern, rownames(x))), na.rm = TRUE) > 0.9
}

.isSummarizedExperiment <- function(x) {
    is(x, "SummarizedExperiment") & !is(x, "RangedSummarizedExperiment")
}

.convertTo <- function(x, which, FUN, keep, unmap) {
    for (i in which(which)) {
        lookup <- FUN(rownames(x[[i]]))
        ranges <- lookup[["mapped"]]
        rse <- x[[i]][names(ranges), ]
        # rowData not merged with mcols of RHS in `rowRanges<-` method
        mcols(ranges) <-
            S4Vectors::DataFrame(rowData(rse), S4Vectors::mcols(ranges))
        SummarizedExperiment::rowRanges(rse) <- ranges
        x <- c(x, setNames(S4Vectors::List(rse),
            paste0(names(x)[i], "_ranged")))
        if (length(lookup[["unmapped"]]) && unmap) {
            se <- x[[i]][lookup[["unmapped"]], ]
            x <- c(x, setNames(S4Vectors::List(se),
                paste0(names(x)[i], "_unranged")))
        }
    }
    if (!keep & any(which))
        x <- x[, , -match(names(which(which)), names(x))]
    x
}

#' @name hidden-helpers
#' @title A small document for helper functions
#' @param x A character vector
#' @param gn A GRanges object with some of its names found in x
#' @return A list of length 2: unmapped (character vector) and mapped (GRanges)
#' @keywords internal
.makeListRanges <- function(x, gn) {
    res <- list(unmapped = x[!x %in% names(gn)])
    x <- x[x %in% names(gn)]
    gn <- gn[match(x, names(gn))]
    res[["mapped"]] <- gn
    return(res)
}

.getGN <- function(gen, FUN) {
    stopifnot(is.character(gen), length(gen) == 1L)

    fun <- switch(FUN,
        genes = GenomicFeatures::genes,
        microRNAs = GenomicFeatures::microRNAs
    )

    txdb <- if (identical(gen, "hg18"))
        TxDb.Hsapiens.UCSC.hg18.knownGene::TxDb.Hsapiens.UCSC.hg18.knownGene
    else if (identical(gen, "hg19"))
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene

    gn <- keepStandardChromosomes(fun(txdb), pruning.mode = "coarse")
    seqlevelsStyle(gn) <- "NCBI"

    if (identical(FUN, "genes"))
        names(gn) <- AnnotationDbi::mapIds(
            org.Hs.eg.db::org.Hs.eg.db,
            names(gn),
            keytype = "ENTREZID",
            column = "SYMBOL"
        )
    else if (identical(FUN, "microRNAs"))
        names(gn) <- mcols(gn)[["mirna_id"]]

    gn
}

#' @rdname hidden-helpers
#' @return list of length 2: "unmapped" is a character vector providing
#' unmapped symbols, "mapped" is a GRanges object with ranges of mapped symbols
#' @keywords internal
.getRangesOfSYMBOLS <- function(x) {
    gn <- .getGN("hg19", "genes")
    .makeListRanges(x, gn)
}

#' @rdname hidden-helpers
#' @param x A SummarizedExperiment containing hsa miR IDs as rownames
#' @keywords internal
.getRangesOfMir <- function(x) {
    mr <- .getGN("hg19", "microRNAs")
    .makeListRanges(x, mr)
}

.checkPkgsAvail <- function(pkgnames) {
    vapply(pkgnames, function(pkgname) {
    if (!requireNamespace(pkgname, quietly = TRUE)) {
        func <- as.character(sys.call(1L)[[1L]])
        func <- func[!(func %in% c("::", "TCGAutils"))]
        stop("Install the '", pkgname, "' package to use '", func, "'",
            call. = FALSE)
    } else
        TRUE
    }, logical(1L))
}

.getRangesOfCpG <- function(x) {
    local_data_store <- new.env(parent = emptyenv())
    data("Locations", envir = local_data_store,
        package = "IlluminaHumanMethylation450kanno.ilmn12.hg19")
    Locations <- local_data_store[["Locations"]]

    clist <- list(seqnames = "chr", pos = "pos", strand = "strand")
    gps <- do.call(GenomicRanges::GPos,
        lapply(clist, function(x) Locations[, x]))
    names(gps) <- rownames(Locations)
    seqlevelsStyle(gps) <- "NCBI"

    .makeListRanges(x, gps)
}

#' @rdname simplifyTCGA
#'
#' @title Functions to convert rows annotations to ranges and RaggedExperiment
#' to RangedSummarizedExperiment
#'
#' @description This group of functions will convert row annotations as
#' either gene symbols or miRNA symbols to row ranges based on database
#' resources 'TxDB' and 'org.Hs' packages. It will also simplify the
#' representation of \linkS4class{RaggedExperiment} objects to
#' \linkS4class{RangedSummarizedExperiment}.
#'
#' @details The original SummarizedExperiment containing either gene symbol
#'   or miR annotations is replaced or supplemented by a
#'   \linkS4class{RangedSummarizedExperiment} for those that could be mapped to
#'   \linkS4class{GRanges}, and optionally another
#'   \linkS4class{SummarizedExperiment} for annotations that
#'   could not be mapped to \linkS4class{GRanges}.
#'
#' RaggedExperiment mutation objects become a genes by patients
#' RangedSummarizedExperiment object containing '1' if there is a non-silent
#' mutation somewhere in the gene, and '0' otherwise as obtained from the
#' `Variant_Classification` column in the data.
#'
#' "CNA" and "CNV" segmented copy number are reduced using a weighted mean in
#' the rare cases of overlapping (non-disjoint) copy number regions.
#'
#' These functions rely on 'TxDb.Hsapiens.UCSC.hg19.knownGene' and
#' 'org.Hs.eg.db' to map to the 'hg19' NCBI build. Users should use the
#' `liftOver` procedure for datasets that are provided against a different
#' reference genome (usually 'hg18'). An example of this procedure is provided
#' in the vignette.
#'
#' @param obj A MultiAssayExperiment object obtained from curatedTCGAData
#'
#' @param keep.assay logical (default FALSE) Whether to keep the
#'   `SummarizedExperiment` assays that have been converted to
#'   `RangedSummarizedExperiment`
#'
#' @param unmapped logical (default TRUE) Include an assay of data that was
#'   not able to be mapped in reference database
#'
#' @param suffix character (default "_simplified") A character string to append
#' to the newly modified assay for `qreduceTCGA`.
#'
#' @return A \linkS4class{MultiAssayExperiment} with any gene expression, miRNA,
#'   copy number, and mutations converted to RangedSummarizedExperiment objects
#'
#' @author L. Waldron
#'
#' @md
#'
#' @examples
#'
#' library(curatedTCGAData)
#' library(GenomeInfoDb)
#'
#' accmae <-
#'     curatedTCGAData(diseaseCode = "ACC",
#'     assays = c("CNASNP", "Mutation", "miRNASeqGene", "GISTICT"),
#'     version = "1.1.38",
#'     dry.run = FALSE)
#'
#' ## update genome annotation
#' rex <- accmae[["ACC_Mutation-20160128"]]
#'
#' ## Translate build to "hg19"
#' tgenome <- vapply(genome(rex), translateBuild, character(1L))
#' genome(rex) <- tgenome
#'
#' accmae[["ACC_Mutation-20160128"]] <- rex
#'
#' simplifyTCGA(accmae)
#'
#' @export
simplifyTCGA <- function(obj, keep.assay = FALSE, unmapped = TRUE) {
    obj <- qreduceTCGA(obj, keep.assay)
    obj <- mirToRanges(obj, keep.assay, unmapped)
    symbolsToRanges(obj, keep.assay, unmapped)
}

#' @name simplifyTCGA
#' @aliases symbolsToRanges
#' @export
symbolsToRanges <- function(obj, keep.assay = FALSE, unmapped = TRUE) {
    can.fix <- vapply(experiments(obj), function(y) {
        .checkHas(y, "symbols") & .isSummarizedExperiment(y)
    }, logical(1L))

    .checkPkgsAvail(c("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"))
    .convertTo(
        x = obj,
        which = can.fix,
        FUN = .getRangesOfSYMBOLS,
        keep = keep.assay,
        unmap = unmapped
    )
}

#' @name simplifyTCGA
#' @aliases mirToRanges
#' @export
mirToRanges <- function(obj, keep.assay = FALSE, unmapped = TRUE) {
    can.fix <- vapply(experiments(obj), function(y) {
        .checkHas(y, "^hsa") & .isSummarizedExperiment(y)
    }, logical(1L))

    .checkPkgsAvail(c("TxDb.Hsapiens.UCSC.hg19.knownGene", "mirbase.db"))
    .convertTo(
        x = obj,
        which = can.fix,
        FUN = .getRangesOfMir,
        keep = keep.assay,
        unmap = unmapped
    )
}

#' @name simplifyTCGA
#' @aliases CpGtoRanges
#' @export
CpGtoRanges <- function(obj, keep.assay = FALSE, unmapped = TRUE) {
    can.fix <- vapply(experiments(obj), function(y) {
        .checkHas(y, "^cg") & .isSummarizedExperiment(y)
    }, logical(1L))

    .checkPkgsAvail(c("IlluminaHumanMethylation450kanno.ilmn12.hg19"))

    .convertTo(
        x = obj,
        which = can.fix,
        FUN = .getRangesOfCpG,
        keep = keep.assay,
        unmap = unmapped
    )
}

#' @name simplifyTCGA
#' @aliases qreduceTCGA
#' @export
qreduceTCGA <- function(obj, keep.assay = FALSE, suffix = "_simplified") {
    .checkPkgsAvail(c("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"))
    gn <- genes(
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
    gn <- keepStandardChromosomes(GenomicRanges::granges(gn),
        pruning.mode = "coarse")
    seqlevelsStyle(gn) <- "NCBI"
    names(gn) <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, names(gn),
        keytype = "ENTREZID", column = "SYMBOL")

    weightedmean <- function(scores, ranges, qranges) {
        isects <- GenomicRanges::pintersect(ranges, qranges)
        sum(scores * BiocGenerics::width(isects)) /
            sum(BiocGenerics::width(isects))
    }

    nonsilent <- function(scores, ranges, qranges)
        any(scores != "Silent")

    isRE <- function(x) vapply(experiments(x), function(y)
          is(y, "RaggedExperiment"), logical(1L))

    isMut <- function(x) grepl("Mutation", names(x))

    for (i in which(isMut(obj))) {
        sqls <- seqlevelsStyle(obj[[i]])
        seqlevelsStyle(gn) <- sqls
        mutations <- RaggedExperiment::qreduceAssay(obj[[i]], gn, nonsilent,
            "Variant_Classification")
        rownames(mutations) <- names(gn)
        mutations[is.na(mutations)] <- 0
        remove.rows <- is.na(rownames(mutations))
        mutations <- SummarizedExperiment(mutations[!remove.rows,],
            rowRanges = gn[!remove.rows])
        el <- ExperimentList(x = mutations)
        names(el) <- paste0(names(obj)[i], suffix)
        obj <- c(obj, el)
    }
    for (i in which(isRE(obj) & !isMut(obj))) {
        sqls <- seqlevelsStyle(obj[[i]])
        seqlevelsStyle(gn) <- sqls
        suppressWarnings(
            cn <- RaggedExperiment::qreduceAssay(obj[[i]], gn,
                weightedmean, "Segment_Mean")
        )
        rownames(cn) <- names(gn)
        remove.rows <- is.na(rownames(cn))
        cn <- SummarizedExperiment(cn[!remove.rows, ], rowRanges = gn[!remove.rows])
        el <- ExperimentList(x = cn)
        names(el) <- paste0(names(obj)[i], suffix)
        obj <- c(obj, el)
    }
    if (!keep.assay) {
        obj <- obj[, , !isRE(obj)]
    }
    return(obj)
}
