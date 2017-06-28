#' @importFrom psygenet2r extract
NULL

.getDataMatrix <- function(object) {
    getElement(object, "DataMatrix")
}

.getFilenames <- function(object) {
    getElement(object, "Filename")
}

.getGISTIC <- function(x, type) {
    x <- getElement(x, type)
    annoteCols <- !grepl("TCGA", names(x))
    annoteRowDF <- x[, annoteCols]
    rownames(annoteRowDF) <-
        annoteRowDF[, grepl("gene", names(annoteRowDF), ignore.case = TRUE)]
    x <- x[, !annoteCols]
    x <- vapply(x, type.convert, numeric(nrow(x)))
    colnames(x) <- .stdIDs(colnames(x))
    SummarizedExperiment(SimpleList(x), rowData = annoteRowDF)
}

.removeShells <- function(x, type) {
    dataTypes <- c("Clinical", "RNAseq_Gene", "miRNASeq_Gene",
    "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP", "CNA_Seq", "CNA_CGH",
    "Methylation", "Mutation", "mRNA_Array", "miRNA_Array", "RPPA_Array",
    "GISTIC_A", "GISTIC_T")
    dataTypes <- gsub("_", "", dataTypes)
    type <- match.arg(type, dataTypes)
    type <- gsub("Seq$", "seq", type)
    if (is(x, "list")) {
        if (length(x) == 1L) {
            x <- x[[1L]]
        if (is(x, "FirehoseCGHArray") || is(x, "FirehosemRNAArray")) {
            fname <- .getFilenames(x)
            x <- SummarizedExperiment(.getDataMatrix(x))
            metadata(x) <- list(filename = fname)
            if (is(x, "data.frame")) {
                x <- DataFrame(x)
                metadata(x) <- .getFilenames(x)
            }
        }
        } else {
            x <- List(lapply(x, .getDataMatrix))
            metadata(x) <- lapply(x, .getFilenames)
        }
    } else {
        x <- getElement(x, type)
    }
    return(x)
}

.fileSelect <- function() {
    g <- readline(
        paste0("The selected data type has more than one",
               "file available.\nPlease select the desired file.",
       "\n(Enter 0 for the first file with the most number of samples)\n_"))
    g <- suppressWarnings(as.integer(g))
    if(is.na(g)){
        stop("Your selection must be an integer!")
    } else {
        return(g)
    }
}

.getBuildCol <- function(x) {
    ncbi <- tolower(names(x)) %in% "ncbi_build"
    if (sum(ncbi) > 1L)
        stop("Multiple ncbi_build columns detected")
    names(x)[ncbi]
}

.hasBuildInfo <- function(x) {
    buildInfo <- .getBuildCol(x)
    as.logical(length(buildInfo))
}

.getBuild <- function(x) {
    binf <- .hasBuildInfo(x)
    if (binf) {
        BCOL <- .getBuildCol(x)
        build <- unique(x[[BCOL]])
        if (length(build) > 1L)
            stop("Inconsistent genome build column")
        build <- as.character(build)
        return(.getHGBuild(build))
    } else {
        stop("Build not available")
    }
}

.ansRangeNames <- function(x) {
    if (is(x, "list")) { return(list()) }
    granges_cols <- findGRangesCols(names(x), seqnames.field = "Chromosome",
        start.field = c("Start", "Start_position"),
        end.field = c("End", "End_position"))
    fielders <- list(seqnames.field = "seqnames", start.field = "start",
        end.field = "end", strand.field = "strand")
    Fargs <- lapply(fielders, function(name) { names(x)[granges_cols[[name]]] })
    Fargs[["ignore.strand"]] <- is.na(Fargs[["strand.field"]])
    Filter(function(g) {!is.na(g)}, Fargs)
}

.findSampleCol <- function(x) {
    tsb <- match("tumor_sample_barcode", tolower(names(x)))
    if (length(tsb) == 1L && !is.na(tsb)) {
        primary <- names(x)[tsb]
    } else if (is.na(tsb)) {
        primary <- names(x)[tolower(names(x)) == "sample"]
    } else {
        stop("'split.field' could not be found")
    }
    return(primary)
}

## TODO: All filename metadata

## TODO: Filter and create by equal lengths per sample (RSE vs RE)

.hasConsistentRanges <- function(object) {
    primary <- .findSampleCol(object)
    S4Vectors::isSingleInteger(unique(
        vapply(base::split(object, object[[primary]]), nrow, integer(1L))
        ))
}

.hasRangeNames <- function(x) {
    if (!any(is.data.frame(x), is(x, "DataFrame"), is.matrix(x)))
        stop("(internal) 'x' must be rectangular")
    !all(is.na(findGRangesCols(names(x), seqnames.field = "Chromosome",
        start.field = c("Start", "Start_position"),
        end.field = c("End", "End_position"))))
}

## Safe to assume equal number of ranges == equal ranges (?)

.makeRangedSummarizedExperimentFromDataFrame <- function(df, ...,
        seqinfo = NULL, starts.in.df.are.0based = FALSE) {
    args <- list(...)
    if (!is.null(args[["build"]]))
        build <- args[["build"]]
    if (!.hasConsistentRanges(df))
        stop("All ranges must be equal in number by 'split.field'")
    split.field <- .findSampleCol(df)
    ansRanges <- .ansRangeNames(df)
    strictRanges <- Filter(function(x) !is.logical(x), ansRanges)
    RangeInfo <- c(strictRanges, list(split.field = split.field))
    numInfo <- df[, !(names(df) %in% RangeInfo)]
    numAssays <- ncol(numInfo)
    nameAssays <- names(numInfo)
    numInfo <- base::split(numInfo, df[, split.field])
    countList <- vector(mode = "list", length = numAssays)
    for (i in seq_len(numAssays)) {
        countList[[i]] <- do.call(cbind, lapply(numInfo,
            function(smalldf) { smalldf[[i]] }))
    }
    names(countList) <- nameAssays
    rowRanges <- makeGRangesListFromDataFrame(df[, unlist(RangeInfo)],
                                              split.field = split.field)
    if (exists("build"))
        GenomeInfoDb::genome(rowRanges) <- build
    return(SummarizedExperiment(assays = SimpleList(countList),
        rowRanges = rowRanges))
}

.makeRaggedExperimentFromDataFrame <- function(df, ...) {
    args <- list(...)
    if (!is.null(args[["build"]]))
        build <- args[["build"]]
    split.field <- .findSampleCol(df)
    ansRanges <- .ansRangeNames(df)
    rangeInfo <- c(ansRanges, list(split.field = split.field))
    dropIdx <- .omitAdditionalIdx(object, ansRanges)
    if (length(dropIdx))
        df <- df[, -dropIdx]
    newGRL <- do.call(makeGRangesListFromDataFrame,
            args = c(list(df = df, keep.extra.columns = TRUE), rangeInfo))
    if (exists("build"))
        GenomeInfoDb::genome(newGRL) <- build
    RaggedExperiment::RaggedExperiment(newGRL)
}

.omitAdditionalIdx <- function(object, rangeNames) {
    rangeNames <- Filter(function(x) !is.logical(x), rangeNames)
    omitAdditional <- c("seqnames", "seqname", "chromosome", "chrom",
        "chromosome_name", "ranges", "seqlevels", "seqlengths", "seq_id",
        "iscircular", "start", "end", "width", "element", "chr")
    diffNames <- setdiff(omitAdditional, tolower(rangeNames))
    which(tolower(names(object)) %in% diffNames)
}

## TODO: Handle strand conditionally

.extractRanged <- function(object, rangeNames) {
    primary <- .findSampleCol(object)
    dropIdx <- .omitAdditionalIdx(object, rangeNames)
    if (length(dropIdx))
        object <- object[, -dropIdx]
    mygrl <- makeGRangesListFromTCGA(df = object,
        split.field = primary,
        seqnames.field = rangeNames[["seqnames.field"]],
        start.field = rangeNames[["start.field"]],
        end.field = rangeNames[["end.field"]],
        strand.field =
            if (is.null(rangeNames[["strand"]]))
                { "strand" } else { rangeNames[["strand"]] },
        keep.extra.columns = TRUE,
        ignore.strand = rangeNames[["ignore.strand"]])
    return(mygrl)
}

setGeneric("extract", getGeneric("extract", package = "psygenet2r"))

#' @export
setMethod("extract", "List", function(object, ...) {
    args <- list(...)
    type <- args[["type"]]
    for (i in seq_along(object))
    object[[i]] <- TCGAextract(object[[i]], type)
})

#' Extract data from \code{FirehoseData} object into \code{ExpressionSet} or
#' \code{GRangesList} object
#'
#' This function processes data from a \code{\linkS4class{FirehoseData}}
#' object from the \code{RTCGAToolbox} package. Raw data is converted to
#' conventional Bioconductor objects. The function returns an
#' \linkS4class{ExpressionSet} or \linkS4class{GRangesList} class object. Note:
#' this function works best with the modifications found in the github fork:
#' \code{LiNk-NY/RTCGAToolbox}. In cases where range data are found
#' (i.e., "mutations") the default extraction method is used
#' (see makeGRangesListFromTCGA).
#'
#' @section type:
#' Choices include: "RNAseq_Gene",
#' "Clinic", "miRNASeq_Gene", "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP",
#' "CNA_Seq", "CNA_CGH", "Methylation", "Mutation", "mRNA_Array",
#' "miRNA_Array", "RPPA_Array", "GISTIC_A", "GISTIC_T". The "GISTIC_A" type of
#' dataset represents GISTIC data by all genes. "GISTIC_T" represents data
#' thresholded by genes. Lowercase entries and entries without the "underscore"
#' character are also valid inputs.
#'
#' @param object A \code{FirehoseData} object from which to extract data.
#' @param type The type of data to extract from the "FirehoseData" object,
#' see type section.
#' @return Either an \linkS4class{ExpressionSet} object or a
#' \linkS4class{GRangesList} object.
#'
#' @author Marcel Ramos \email{marcel.ramos@roswellpark.org}
#'
#' @examples \dontrun{
#' library(RTCGAToolbox)
#' dataFolder <- normalizePath("~/Documents/data")
#' coadmut <- getFirehoseData("COAD", runDate = "20151101", Mutation = TRUE,
#'                          destdir = dataFolder)
#' cm <- TCGAextract(coadmut, "mutations")
#' }
#' @importClassesFrom RTCGAToolbox FirehoseData FirehosemRNAArray
#' FirehoseCGHArray FirehoseMethylationArray
#' @export TCGAextract
TCGAextract <- function(object, type = c("Clinical", "RNAseq_Gene",
    "miRNASeq_Gene", "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP", "CNA_Seq",
    "CNA_CGH", "Methylation", "Mutation", "mRNA_Array", "miRNA_Array",
    "RPPA_Array", "GISTIC_A", "GISTIC_T"), ...) {
    sNames <- c("Clinical", "RNASeqGene", "RNASeq2GeneNorm", "miRNASeqGene",
                "CNASNP", "CNVSNP", "CNAseq", "CNACGH", "Methylation",
                "mRNAArray", "miRNAArray", "RPPAArray", "Mutations", "GISTIC")
    object <- .removeShells(object, type)
    if (type == "Clinical") { return(object) }
    if (is(object, "matrix")) {
        return(SummarizedExperiment(assays = SimpleList(object)))
    }
    if (is(object, "SummarizedExperiment")) { return(object) }
    hasRanged <- .hasRangeNames(object)
    if (hasRanged) {
        if (.hasBuildInfo(object)) {
            build <- .getBuild(object)
        }
        if (.hasConsistentRanges(object)) {
            object <- .makeRangedSummarizedExperimentFromDataFrame(object,
                build = if (exists("build")) { build } else { NULL })
        } else {
            object <- .makeRaggedExperimentFromDataFrame(object)
        }
        return(object)
    }
    if (is(object, "List")) {
        return(extract(object, type = type, ...))
    }
    rangeslots <- c("CNVSNP", "CNASNP", "CNAseq", "CNACGH", "Mutations")
    slotreq <- grep(paste0("^", type) , sNames, ignore.case=TRUE, value=TRUE)
    gisticType <- grepl("^GISTIC", type, ignore.case = TRUE)
    if (gisticType) {
        slotreq <- switch(type, GISTICA = "AllByGene",
                          GISTICT = "ThresholdedByGene")
        type <- gsub("A$|T$", "", type)
        object <- getElement(object, type)
        result <- .getGISTIC(object, slotreq)
        return(result)
    }
    if (slotreq %in% "Methylation") {
        object <- getElement(object, "DataMatrix")
        headers <- names(object)
        annote <- object[, !grepl("TCGA", headers)]
        isNumRow <- all(grepl("^[0-9]*$",
            sample(rownames(object), size = 100L, replace = TRUE)))
        if (isNumRow) {
            geneSymbols <- annote[, grep("symbol", names(annote),
                                         ignore.case = TRUE, value = TRUE)]
            rNames <- geneSymbols
        } else { rNames <- rownames(object) }
        dm <- data.matrix(object[, grepl("TCGA", names(object))])
        rownames(dm) <- rNames
        colnames(dm) <- .stdIDs(colnames(dm))
        newSE <- SummarizedExperiment::SummarizedExperiment(
            assays = SimpleList(dm), rowData = annote)
        return(newSE)
        } else if (slotreq %in% rangeslots) {
            object <- .extractRanged(object, rangeNames)
        }
    return(object)
}
