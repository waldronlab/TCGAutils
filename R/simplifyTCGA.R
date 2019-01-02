.cMAE <- function(mae, x, name) {
    if (missing(name))
        stop("<internal> Provide a name for the 'ExperimentList' element")
    newlist <- setNames(list(x), name)
    c(mae, ExperimentList(newlist))
}

.hasMir <- function(x) {
    mean(c(FALSE, grepl("^hsa", rownames(x))), na.rm = TRUE) > 0.9
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

.makeListRanges <- function(x, gn) {
    ## x is a character vector
    ## gn is a GRanges object with some of its names found in x
    res <- list(unmapped = x[!x %in% names(gn)])
    x <- x[x %in% names(gn)]
    gn <- gn[match(x, names(gn))]
    res$mapped <- gn
    ## Returns a list of length 2:
    ## unmapped = character vector of x not found in gn
    ## mapped = GRanges of gn that are found in x
    return(res)
}

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
    ## returns a list of length 2: "unmapped" is a character vector providing
    ## unmapped symbols, "mapped" is a GRanges object with ranges of mapped symbols.
    return(.makeListRanges(x, gn))
}

.getRangesOfMir <- function(x) {
## x is a SummarizedExperiment containing hsa miR IDs as rownames
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

#' Simplify curatedTCGAData objects by replacing RaggedExperiment objects
#'
#' @param obj A MultiAssayExperiment object obtained from curatedTCGAData
#' @param keep logical (default FALSE) whether to remove the original
#'   RaggedExperiment objects from the returned MultiAssayExperiment
#' @param suffix character (default "_simplified") A character string to append
#' to the newly modified assay.
#'
#' @return
#' A MultiAssayExperiment object with RaggedExperiments converted to
#' RangedSummarizedExperiment with rows corresponding to gene symbol.
#'
#' "Mutations" objects become a genes x patients RangedSummarizedExperiment
#' containing 1 if there is a non-silent mutation somewhere in the gene, and 0
#' otherwise. "CNA" and "CNV" segmented copy number are reduced using a weighted
#' mean in the rare cases of overlapping (non-disjoint) copy number regions.
#' @details
#' Relies on TxDb.Hsapiens.UCSC.hg19.knownGene and org.Hs.eg.db to map to hg19
#' NCBI build.
#' @examples
#' library(curatedTCGAData)
#' library(GenomeInfoDb)
#'
#' accmae <-
#'     curatedTCGAData("ACC", c("CNASNP", "Mutation"), dry.run = FALSE)
#'
#' qreduceTCGA(accmae)
#'
#' @importFrom GenomicFeatures genes microRNAs
#' @importFrom GenomeInfoDb keepStandardChromosomes seqlevelsStyle
#' seqlevelsStyle<-
#'
#' @author L. Waldron
#'
#' @export
qreduceTCGA <- function(obj, keep = FALSE, suffix = "_simplified") {
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
    if (!keep) {
        obj <- obj[, , !isRE(obj)]
    }
    return(obj)
}

#' Convert SummarizedExperiment elements with gene symbols to
#' RangedSummarizedExperiment
#'
#' @param obj A MultiAssayExperiment object obtained from curatedTCGAData
#' @param keep If FALSE (default), remove the SummarizedExperiment objects
#' that have been converted to RangedSummarizedExperiment
#'
#' @return a MultiAssayExperiment where any of the original SummarizedExperiment
#'   containing gene symbols as rownames have been replaced or supplemented by a
#'   RangedSummarizedExperiment for miR that could be mapped to GRanges, and
#'   another SummarizedExperiment for miR that could not be mapped to GRanges.
#'
#' @seealso mirToRanges
#'
#' @details Any SummarizedExperiment elements with gene symbols as rownames will
#'   have ranges added. Symbols where ranges can't be found are put in a new
#'   SummarizedExperiment.
#'
#' @author L. Waldron
#'
#' @examples
#' library(MultiAssayExperiment)
#'
#' data(miniACC)
#'
#' symbolsToRanges(miniACC)
#'
#' @export
symbolsToRanges <- function(obj, keep = FALSE) {
    can.fix <- vapply(experiments(obj), function(y) {
        .hasSymbols(y) & .isSummarizedExperiment(y)
    }, TRUE)

    .checkPkgsAvail(c("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"))
    for (i in which(can.fix)) {
        lookup <- .getRangesOfSYMBOLS(rownames(obj[[i]]))
        rse <- obj[[i]][names(lookup$mapped),]
        SummarizedExperiment::rowRanges(rse) <- lookup$mapped
        obj <- .cMAE(obj, rse, name = paste0(names(obj)[i], "_ranged"))
        if (length(lookup$unmapped)) {
            se <- obj[[i]][lookup$unmapped,]
            obj <- .cMAE(obj, se, name = paste0(names(obj)[i], "_unranged"))
        }
    }
    if (!keep && any(can.fix))
        obj <- obj[, ,-which(can.fix)]
    return(obj)
}

#' Convert SummarizedExperiment elements with microRNA to
#' RangedSummarizedExperiment
#'
#' @param obj A MultiAssayExperiment object obtained from curatedTCGAData
#' @param keep If FALSE (default), remove the SummarizedExperiment
#'   objects that have been converted to RangedSummarizedExperiment
#' @return a MultiAssayExperiment where any of the original SummarizedExperiment
#'   containing gene symbols as rownames have been replaced or supplemented by a
#'   RangedSummarizedExperiment for miR that could be mapped to GRanges, and
#'   another SummarizedExperiment for miR that could not be mapped to GRanges.
#' @seealso symbolsToRanges
#'
#' @author L. Waldron
#'
#' @examples
#' library(curatedTCGAData)
#'
#' accmae <- curatedTCGAData("ACC", "miRNASeqGene", dry.run = FALSE)
#'
#' mirToRanges(accmae)
#' @export
mirToRanges <- function(obj, keep = FALSE) {
    can.fix <- vapply(experiments(obj), function(y) {
        .hasMir(y) & .isSummarizedExperiment(y)
    }, TRUE)

    .checkPkgsAvail(c("TxDb.Hsapiens.UCSC.hg19.knownGene", "mirbase.db"))
    for (i in which(can.fix)) {
        lookup <- .getRangesOfMir(rownames(obj[[i]]))
        rse <- obj[[i]][names(lookup$mapped), ]
        SummarizedExperiment::rowRanges(rse) <- lookup$mapped
        obj <- .cMAE(obj, rse, paste0(names(obj)[i], "_ranged"))
        if (length(lookup$unmapped > 0)) {
            se <- obj[[i]][lookup$unmapped, ]
            obj <- .cMAE(obj, se, paste0(names(obj)[i], "_unranged"))
        }
        if (!keep & any(can.fix))
            obj <- obj[, , -which(can.fix)]
    }
    return(obj)
}

#' All-in-one simplification of curatedTCGAData objects
#'
#' @param obj A MultiAssayExperiment from curatedTCGAData
#' @param keep If FALSE (default), remove the original
#'   MultiAssayExperiment elements that have simplified versions in the output.
#'
#' @return a MultiAssayExperiment with any gene expression, miRNA, copy number,
#'   and mutations converted to RangedSummarizedExperiment objects
#'
#' @author L. Waldron
#'
#' @seealso mirToRanges, symbolsToRanges, qreduceTCGA
#' @examples
#' library(curatedTCGAData)
#' library(GenomeInfoDb)
#'
#' accmae <- curatedTCGAData("ACC",
#'     c("CNASNP", "Mutation", "miRNASeqGene", "GISTICT"),
#'     dry.run = FALSE)
#'
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
simplifyTCGA <- function(obj, keep = FALSE) {
    obj <- qreduceTCGA(obj, keep)
    obj <- mirToRanges(obj, keep)
    symbolsToRanges(obj, keep)
}
