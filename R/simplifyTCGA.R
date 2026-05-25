#' @importFrom GenomicFeatures genes
#' @importFrom GenomeInfoDb keepStandardChromosomes seqlevelsStyle
#' @importFrom GenomeInfoDb seqlevelsStyle<-
#' @importFrom BiocBaseUtils isScalarCharacter isCharacter
NULL

.checkHas <-
    function(x, pattern, threshold = 0.9)
{
    if (identical(pattern, "symbols"))
        pattern <- "^[A-Z0-9]{1,6}|^C[0-9]orf[0-9]{1,4}"
    mean(c(FALSE, grepl(pattern, rownames(x))), na.rm = TRUE) > threshold
}

.isFixable <- function(mae, pattern = c("^hsa", "^cg", "symbols")) {
    if (missing(pattern) || !isScalarCharacter(pattern))
        stop("<internal> Provide a single 'pattern' value to search for")
    vapply(
        experiments(mae),
        function(y) {
            .checkHas(x = y, pattern = pattern) &&
            (
                is(y, "SummarizedExperiment") &&
                    !is(y, "RangedSummarizedExperiment")
            )
        },
        logical(1L)
    )
}

.convertTo <- function(x, which, FUN, keep, unmap) {
    for (i in which(which)) {
        assay <- x[[i]]
        lookup <- FUN(rownames(assay))
        ranges <- lookup[["mapped"]]
        rowidx <- mcols(ranges)[["rowIdx"]]
        rowidx <- Filter(Negate(is.na), rowidx)
        if (!is.null(rowidx) && length(rowidx))
            rse <- `rownames<-`(assay[rowidx, ], names(ranges))
        else
            rse <- assay[names(ranges), ]
        # rowData not merged with mcols of RHS in `rowRanges<-` method
        mcols(ranges) <-
            S4Vectors::DataFrame(rowData(rse), S4Vectors::mcols(ranges))
        SummarizedExperiment::rowRanges(rse) <- ranges
        x <- c(x, setNames(S4Vectors::List(rse),
            paste0(names(x)[i], "_ranged")))
        if (length(lookup[["unmapped"]]) && unmap) {
            se <- assay[lookup[["unmapped"]], ]
            x <- c(x, setNames(S4Vectors::List(se),
                paste0(names(x)[i], "_unranged")))
        }
    }
    if (!keep & any(which))
        x <- x[, , -match(names(which(which)), names(x))]
    x
}

#' @name hidden-helpers
#' @title Helper functions to get genomic ranges from different identifiers
#' @param x `character()` vector of either micro RNA identifiers
#'   (`.getRangesOfMir`) or gene symbols (`.getRangesOfSYMBOLS`) or CpG probe
#'   identifiers (`.getRangesOfCpG`)
#' @param gn `GRanges()` with some of its names found in x or translated x
#' @return list of length 2: "unmapped" `character()` symbols, "mapped" is a
#'   `GRanges()` object with ranges of mapped symbols
#' @keywords internal
.makeListRanges <- function(x, gn) {
    res <- list(unmapped = x[!x %in% names(gn)])
    x <- x[x %in% names(gn)]
    gn <- gn[match(x, names(gn))]
    res[["mapped"]] <- gn
    res
}

#' @name hidden-helpers
#' @keywords internal
.makeMiRNAListRanges <- function(x, gn) {
    checkInstalled("miRNAmeConverter")
    nc <- miRNAmeConverter::MiRNANameConverter()
    mirna_version <-
        miRNAmeConverter::assessVersion(nc, names(gn))[1L, "version"]
    trout <- miRNAmeConverter::translateMiRNAName(
        nc, x, versions = mirna_version
    )
    new_x <- trout[[paste0("v", mirna_version, ".0")]]
    res <- list(unmapped = setdiff(x, trout[["input"]]))
    rowIdx <- match(tolower(trout[["input"]]), tolower(x))
    gn <- gn[match(new_x, names(gn))]
    mcols(gn)[["rowIdx"]] <- rowIdx
    res[["mapped"]] <- gn
    res
}

.getGN <- function(gen) {
    stopifnot(isScalarCharacter(gen))

    txdb <- if (identical(gen, "hg18"))
        TxDb.Hsapiens.UCSC.hg18.knownGene::TxDb.Hsapiens.UCSC.hg18.knownGene
    else if (identical(gen, "hg19"))
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene

    gn <- keepStandardChromosomes(
        GenomicFeatures::genes(txdb), pruning.mode = "coarse"
    )
    seqlevelsStyle(gn) <- "NCBI"

    names(gn) <- AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        names(gn),
        keytype = "ENTREZID",
        column = "SYMBOL"
    )
    gn
}

#' @rdname hidden-helpers
#' @keywords internal
.getRangesOfMir <- function(x) {
    stopifnot(isCharacter(x))

    mirnas_gr <- .get_hsa_gff3("hg19")

    miR <- mirnas_gr[
        mcols(mirnas_gr)[["type"]] %in% c("miRNA", "microRNA", "tRNA")
    ]
    miR <- keepStandardChromosomes(miR, pruning.mode = "coarse")
    seqlevelsStyle(miR) <- "NCBI"

    names(miR) <- mcols(miR)[["Name"]]
    .makeMiRNAListRanges(x, miR)
}

#' @rdname hidden-helpers
#' @keywords internal
.getRangesOfSYMBOLS <- function(x) {
    gn <- .getGN("hg19")
    .makeListRanges(x, gn)
}

#' @rdname hidden-helpers
#' @keywords internal
.getRangesOfCpG <- function(x) {
    local_data_store <- new.env(parent = emptyenv())
    data(
        "Locations",
        envir = local_data_store,
        package = "IlluminaHumanMethylation450kanno.ilmn12.hg19"
    )
    Locations <- local_data_store[["Locations"]]

    clist <- list(seqnames = "chr", pos = "pos", strand = "strand")
    gps <- do.call(
        GenomicRanges::GPos,
        lapply(clist, function(x) Locations[, x])
    )
    names(gps) <- rownames(Locations)
    seqlevelsStyle(gps) <- "NCBI"

    .makeListRanges(x, gps)
}

#' @title Functions to convert rows annotations to ranges and RaggedExperiment
#' to RangedSummarizedExperiment
#'
#' @description This group of functions will convert row annotations as
#' either gene symbols or miRNA symbols to row ranges based on database
#' resources 'TxDB' and 'org.Hs' packages. It will also simplify the
#' representation of
#' [RaggedExperiment][RaggedExperiment::RaggedExperiment-class] objects to
#' [RangedSummarizedExperiment][SummarizedExperiment::RangedSummarizedExperiment-class].
#'
#' @details The original `SummarizedExperiment` containing either gene symbol
#'   or miR annotations is replaced or supplemented by a
#'   [RangedSummarizedExperiment][SummarizedExperiment::RangedSummarizedExperiment-class]
#'   for those that could be mapped to
#'   [GRanges][GenomicRanges::GRanges-class], and optionally another
#'   [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment-class]
#'   for annotations that could not be mapped to
#'   [GRanges][GenomicRanges::GRanges-class].
#'
#' @section qreduceTCGA:
#'
#' Using `TxDb.Hsapiens.UCSC.hg19.knownGene` as the reference, `qreduceTCGA`
#' reduces the data by applying either the `weightedmean` or `nonsilent`
#' function (see below) to non-mutation or mutation data, respectively.
#' Internally, it uses [RaggedExperiment::qreduceAssay()] to reduce the ranges
#' to the gene-level.
#'
#' `qreduceTCGA` will update `genome(x)` based on the NCBI reference annotation
#' which includes the patch number, e.g., GRCh37.p14, as provided by the
#' `seqlevelsStyle` setter, `seqlevelsStyle(gn) <- "NCBI"`. `qreduceTCGA`
#' uses the NCBI genome annotation as the default reference.
#'
#'     nonsilent <- function(scores, ranges, qranges)
#'         any(scores != "Silent")
#'
#' `RaggedExperiment` mutation objects become a genes by patients
#' `RangedSummarizedExperiment` object containing '1' if there is a non-silent
#' mutation somewhere in the gene, and '0' otherwise as obtained from the
#' `Variant_Classification` column in the data.
#'
#'     weightedmean <- function(scores, ranges, qranges) {
#'         isects <- GenomicRanges::pintersect(ranges, qranges)
#'         sum(scores * BiocGenerics::width(isects)) /
#'             sum(BiocGenerics::width(isects))
#'     }
#'
#' "CNA" and "CNV" segmented copy number are reduced using a weighted mean in
#' the rare cases of overlapping (non-disjoint) copy number regions.
#'
#' These functions rely on `TxDb.Hsapiens.UCSC.hg19.knownGene` and
#' `org.Hs.eg.db` to map to the 'hg19' NCBI build. Use the `liftOver` procedure
#' for datasets that are provided against a different reference genome (usually
#' 'hg18'). See an example in the vignette.
#'
#' @param obj A `MultiAssayExperiment` object obtained from `curatedTCGAData`
#'
#' @param keep.assay logical (default FALSE) Whether to keep the
#'   `SummarizedExperiment` assays that have been converted to
#'   `RangedSummarizedExperiment`
#'
#' @param unmapped logical (default TRUE) Include an assay of data that was
#'   not able to be mapped in reference database
#'
#' @param suffix character (default "_simplified") A character string to append
#'   to the newly modified assay for `qreduceTCGA`.
#'
#' @return A
#'   [`MultiAssayExperiment`][MultiAssayExperiment::MultiAssayExperiment-class]
#'   with any gene expression, miRNA, copy number, and mutations converted to
#'   [`RangedSummarizedExperiment`][SummarizedExperiment::RangedSummarizedExperiment-class]
#'   objects
#'
#' @author L. Waldron, M. Ramos
#'
#' @examples
#' library(curatedTCGAData)
#' library(GenomeInfoDb)
#'
#' accmae <- curatedTCGAData(
#'     diseaseCode = "ACC",
#'     assays = c("CNASNP", "Mutation", "miRNASeqGene", "GISTICT"),
#'     version = "1.1.38",
#'     dry.run = FALSE
#' )
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
#' @export
simplifyTCGA <- function(obj, keep.assay = FALSE, unmapped = TRUE) {
    obj <- qreduceTCGA(obj, keep.assay)
    obj <- mirToRanges(obj, keep.assay, unmapped)
    symbolsToRanges(obj, keep.assay, unmapped)
}

#' @rdname simplifyTCGA
#' @importFrom BiocBaseUtils checkInstalled
#' @export
symbolsToRanges <- function(obj, keep.assay = FALSE, unmapped = TRUE) {
    checkInstalled(c("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"))

    can.fix <- .isFixable(mae = obj, pattern = "symbols")

    .convertTo(
        x = obj,
        which = can.fix,
        FUN = .getRangesOfSYMBOLS,
        keep = keep.assay,
        unmap = unmapped
    )
}

#' @rdname simplifyTCGA
#' @export
mirToRanges <- function(obj, keep.assay = FALSE, unmapped = TRUE) {
    checkInstalled("Bioc.gff")

    can.fix <- .isFixable(mae = obj, pattern = "^hsa")

    .convertTo(
        x = obj,
        which = can.fix,
        FUN = .getRangesOfMir,
        keep = keep.assay,
        unmap = unmapped
    )
}

#' @rdname simplifyTCGA
#' @export
CpGtoRanges <- function(obj, keep.assay = FALSE, unmapped = TRUE) {
    checkInstalled("IlluminaHumanMethylation450kanno.ilmn12.hg19")

    can.fix <- .isFixable(mae = obj, pattern = "^cg")

    .convertTo(
        x = obj,
        which = can.fix,
        FUN = .getRangesOfCpG,
        keep = keep.assay,
        unmap = unmapped
    )
}

#' @rdname simplifyTCGA
#' @export
qreduceTCGA <- function(obj, keep.assay = FALSE, suffix = "_simplified") {
    checkInstalled(c("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"))
    gn <- genes(
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    )
    gn <- keepStandardChromosomes(
        GenomicRanges::granges(gn),
        pruning.mode = "coarse"
    )
    seqlevelsStyle(gn) <- "NCBI"
    names(gn) <- AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        names(gn),
        keytype = "ENTREZID",
        column = "SYMBOL"
    )

    weightedmean <- function(scores, ranges, qranges) {
        isects <- GenomicRanges::pintersect(ranges, qranges)
        sum(scores * BiocGenerics::width(isects)) /
            sum(BiocGenerics::width(isects))
    }

    nonsilent <- function(scores, ranges, qranges)
        any(scores != "Silent")

    isRE <-
        function(x) vapply(experiments(x), is, logical(1L), "RaggedExperiment")

    isMut <- function(x) grepl("Mutation", names(x))

    for (i in which(isMut(obj))) {
        sqls <- seqlevelsStyle(obj[[i]])
        seqlevelsStyle(gn) <- sqls
        ## remove patch release info
        gname <- genome(gn)
        genome(gn) <- gsub("\\.p[0-9]{1,2}$", "", genome(gn))
        mutations <- RaggedExperiment::qreduceAssay(
            obj[[i]],
            gn,
            nonsilent,
            "Variant_Classification"
        )
        rownames(mutations) <- names(gn)
        mutations[is.na(mutations)] <- 0
        remove.rows <- is.na(rownames(mutations))
        mut_ranges <- gn[!remove.rows]
        ## replace patch release info
        genome(mut_ranges) <- gname
        mutations <- SummarizedExperiment(
            mutations[!remove.rows, ], rowRanges = mut_ranges
        )
        el <- ExperimentList(x = mutations)
        names(el) <- paste0(names(obj)[i], suffix)
        obj <- c(obj, el)
    }
    for (i in which(isRE(obj) & !isMut(obj))) {
        sqls <- seqlevelsStyle(obj[[i]])
        seqlevelsStyle(gn) <- sqls
        suppressWarnings(
            cn <- RaggedExperiment::qreduceAssay(
                obj[[i]],
                gn,
                weightedmean,
                "Segment_Mean"
            )
        )
        rownames(cn) <- names(gn)
        remove.rows <- is.na(rownames(cn))
        cn <- SummarizedExperiment(
            cn[!remove.rows, ], rowRanges = gn[!remove.rows]
        )
        el <- ExperimentList(x = cn)
        names(el) <- paste0(names(obj)[i], suffix)
        obj <- c(obj, el)
    }
    if (!keep.assay) {
        obj <- obj[, , !isRE(obj)]
    }
    return(obj)
}
