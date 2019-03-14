#' @importFrom GenomicFeatures genes microRNAs
#' @importFrom GenomeInfoDb keepStandardChromosomes seqlevelsStyle
#' seqlevelsStyle<-
NULL

.cMAE <- function(mae, x, name) {
    if (missing(name))
        stop("<internal> Provide a name for the 'ExperimentList' element")
    newlist <- setNames(list(x), name)
    c(mae, ExperimentList(newlist))
}

.checkHas <- function(x, what = c("^hsa", "^cg", "symbols"), threshold = 0.9) {
    if (identical(what, "symbols"))
        what <- "^[A-Z0-9]{1,6}|^C[0-9]orf[0-9]{1,4}"
    mean(c(FALSE, grepl(what, rownames(x))), na.rm = TRUE) > 0.9
}

.hasMir <- function(x) {
    mean(c(FALSE, grepl("^hsa", rownames(x))), na.rm = TRUE) > 0.9
}

.hasCpG <- function(x) {
    mean(c(FALSE, grepl("^cg", rownames(x))), na.rm = TRUE) > 0.9
}

.hasSymbols <- function(x) {
    mean(c(
        FALSE,
        grepl("^[A-Z0-9]{1,6}|^C[0-9]orf[0-9]{1,4}", rownames(x))
    ), na.rm = TRUE) > 0.9
}

.isSummarizedExperiment <- function(x) {
    is(x, "SummarizedExperiment") & !is(x, "RangedSummarizedExperiment")
}

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

#' @return list of length 2: "unmapped" is a character vector providing
#' unmapped symbols, "mapped" is a GRanges object with ranges of mapped symbols
#' @keywords internal
.getRangesOfSYMBOLS <- function(x) {
    entrez <-
        AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, x, keytype = "SYMBOL",
        column = "ENTREZID")
    gn <- GenomicFeatures::genes(
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
    names(gn) <-
        AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
            names(gn),
            keytype = "ENTREZID",
            column = "SYMBOL")
    gn <- keepStandardChromosomes(GenomicRanges::granges(gn),
        pruning.mode = "coarse")
    seqlevelsStyle(gn) <- "NCBI"

    return(.makeListRanges(x, gn))
}

#' @param x A SummarizedExperiment containing hsa miR IDs as rownames
#' @keywords internal
.getRangesOfMir <- function(x) {
    mr <- microRNAs(
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
    names(mr) <- mr$mirna_id
    mr <- keepStandardChromosomes(granges(mr),
        pruning.mode = "coarse")
    seqlevelsStyle(mr) <- "NCBI"
    return(.makeListRanges(x, mr))
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
    annote450k <- minfi::getAnnotation(
        IlluminaHumanMethylation450kanno.ilmn12.hg19::
            IlluminaHumanMethylation450kanno.ilmn12.hg19)
    clist <- list(seqnames = "chr", ranges = "pos", strand = "strand")
    gps <- do.call(GRanges, lapply(clist, function(x) annote450k[, x]))
    names(gps) <- rownames(annote450k)

    res <- split(x, x %in% names(gps))
    names(res) <- c("unmapped", "mapped")[seq_along(res)]
    gps <- gps[match(res[["mapped"]], names(gps)), ]
    res[["mapped"]] <- gps
    return(res)
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
#' mutation somewhere in the gene, and '0' otherwise. "CNA" and "CNV" segmented
#' copy number are reduced using a weighted mean in the rare cases of
#' overlapping (non-disjoint) copy number regions.
#'
#' These functions rely on 'TxDb.Hsapiens.UCSC.hg19.knownGene' and
#' 'org.Hs.eg.db' to map to the 'hg19' NCBI build.
#'
#' @param obj A MultiAssayExperiment object obtained from curatedTCGAData
#' @param keep.assay logical (default FALSE) Whether to keep the
#'   SummarizedExperiment assays that have been converted to
#'   RangedSummarizedExperiment
#' @param unmapped logical (default TRUE) Include an assay of data that was
#'   not able to be mapped in reference database
#' @param suffix character (default "_simplified") A character string to append
#' to the newly modified assay for `qreduceTCGA`.
#'
#' @return A \linkS4class{MultiAssayExperiment} with any gene expression, miRNA,
#'   copy number, and mutations converted to RangedSummarizedExperiment objects
#'
#' @author L. Waldron
#'
#' @examples
#' library(curatedTCGAData)
#' library(GenomeInfoDb)
#'
#' accmae <-
#'     curatedTCGAData(diseaseCode = "ACC",
#'     assays = c("CNASNP", "Mutation", "miRNASeqGene", "GISTICT"),
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
        .hasSymbols(y) & .isSummarizedExperiment(y)
    }, logical(1L))

    .checkPkgsAvail(c("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"))
    for (i in which(can.fix)) {
        lookup <- .getRangesOfSYMBOLS(rownames(obj[[i]]))
        rse <- obj[[i]][names(lookup$mapped),]
        SummarizedExperiment::rowRanges(rse) <- lookup$mapped
        obj <- .cMAE(obj, rse, name = paste0(names(obj)[i], "_ranged"))
        if (length(lookup$unmapped) && unmapped) {
            se <- obj[[i]][lookup$unmapped, ]
            obj <- .cMAE(obj, se, name = paste0(names(obj)[i], "_unranged"))
        }
    }
    if (!keep.assay && any(can.fix))
        obj <- obj[, ,-which(can.fix)]
    return(obj)
}

#' @name simplifyTCGA
#' @aliases mirToRanges
#' @export
mirToRanges <- function(obj, keep.assay = FALSE, unmapped = TRUE) {
    can.fix <- vapply(experiments(obj), function(y) {
        .hasMir(y) & .isSummarizedExperiment(y)
    }, logical(1L))

    .checkPkgsAvail(c("TxDb.Hsapiens.UCSC.hg19.knownGene", "mirbase.db"))
    for (i in which(can.fix)) {
        lookup <- .getRangesOfMir(rownames(obj[[i]]))
        rse <- obj[[i]][names(lookup$mapped), ]
        SummarizedExperiment::rowRanges(rse) <- lookup$mapped
        obj <- .cMAE(obj, rse, paste0(names(obj)[i], "_ranged"))
        if (length(lookup$unmapped) && unmapped) {
            se <- obj[[i]][lookup$unmapped, ]
            obj <- .cMAE(obj, se, paste0(names(obj)[i], "_unranged"))
        }
        if (!keep.assay & any(can.fix))
            obj <- obj[, , -which(can.fix)]
    }
    return(obj)
}

CpGtoRanges <- function(obj, keep.assay = FALSE, unmapped = TRUE) {
    can.fix <- vapply(experiments(obj), function(y) {
        .hasCpG(y) & .isSummarizedExperiment(y)
    }, logical(1L))

    .checkPkgsAvail(c("IlluminaHumanMethylation450kanno.ilmn12.hg19", "minfi"))

    browser()
    for (i in which(can.fix)) {
        lookup <- .getRangesOfCpG(rownames(obj[[i]]))
        rse <- obj[[i]][names(lookup[["mapped"]]), ]
        SummarizedExperiment::rowRanges(rse) <- lookup[["mapped"]]
        obj <- .cMAE(obj, rse, paste0(names(obj)[i], "_ranged"))
        if (length(lookup[["unmapped"]]) && unmapped) {
            se <- obj[[i]][lookup[["unmapped"]], ]
            obj <- .cMAE(obj, se, paste0(names(obj)[i], "_unranged"))
        }
        if (!keep.assay & any(can.fix))
            obj <- obj[, , -which(can.fix)]
    }
    return(obj)
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

    weightedmean <- function(scores, ranges, qranges)
    ## weighted average score per query range
    sum(scores * BiocGenerics::width(ranges)) / sum(BiocGenerics::width(ranges))

    nonsilent <- function(scores, ranges, qranges)
        any(scores != "Silent")

    isRE <- function(x) vapply(experiments(x), function(y)
          is(y, "RaggedExperiment"), logical(1L))

    isMut <- function(x) grepl("Mutation", names(x))

    for (i in which(isMut(obj))) {
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
