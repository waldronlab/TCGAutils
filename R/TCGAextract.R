#' @importFrom psygenet2r extract
NULL

.getDataMatrix <- function(object) {
    getElement(object, "DataMatrix")
}

.getFilenames <- function(object) {
    getElement(object, "Filename")
}

.standardizeBC <- function(x) {
    colnames(x) <- .stdIDs(colnames(x))
    return(x)
}

.getGISTIC <- function(x, type) {
    x <- getElement(x, type)
    annoteCols <- !grepl("TCGA", names(x))
    annoteRowDF <- x[, annoteCols]
    rownames(annoteRowDF) <-
        annoteRowDF[, grepl("gene", names(annoteRowDF), ignore.case = TRUE)]
    x <- x[, !annoteCols]
    x <- vapply(x, type.convert, numeric(nrow(x)))
    colnames(x) <- .standardizeBC(colnames(x))
    SummarizedExperiment(SimpleList(x), rowData = annoteRowDF)
}

.removeShell <- function(x, type) {
    dataTypes <- c("Clinical", "RNASeqGene", "miRNASeqGene",
    "RNASeq2GeneNorm", "CNASNP", "CNVSNP", "CNASeq", "CNACGH",
    "Methylation", "Mutation", "mRNAArray", "miRNAArray", "RPPAArray",
    "GISTICA", "GISTICT")
    type <- match.arg(type, dataTypes)
    type <- gsub("A$|T$", "", type)
    x <- getElement(x, type)
    return(x)
}

.searchBuild <- function(x) {
    gsub("(^.+)_(hg[0-9]{2})_(.+$)", "\\2",
         x = x, ignore.case = TRUE)
}


## TODO: Fix this

.searchPlatform <- function(x) {
    gsub("(^.+)_(cgh_[0-9]{3}[ak])_(.+$)", "\\2", ffname )
}

.unNestList <- function(x) {
    suppclasses <- all(vapply(x, function(y) {
        any(is(y, "FirehosemRNAArray"), is(y, "FirehoseCGHArray")) },
        logical(1L)))
    if (suppclasses) {
        x <- lapply(x, function(y) {
            fname <- .getFilenames(y)
            if (!.hasBuildInfo(y))
                 build <- .searchBuild(fname)
            y <- .getDataMatrix(y)
            y <- DataFrame(y)
            metadata(y) <- list(filename = fname, build = build)
            return(y)
        })
    }
    if (length(x) == 1L) {
        x <- x[[1L]]
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

.hasConsistentRanges <- function(object) {
    primary <- .findSampleCol(object)
    if (is(object, "DataFrame"))
        asListData <- IRanges::splitAsList(object, object[[primary]])
    else
        asListData <- base::split(object, object[[primary]])
    S4Vectors::isSingleInteger(unique( vapply(asListData, nrow, integer(1L)) ))
}

.hasRangeNames <- function(x) {
    if (is(x, "list")) { return(FALSE) }
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
    metadat <- metadata(df)
    if (!.hasConsistentRanges(df))
        stop("All ranges must be equal in number by 'split.field'")
    split.field <- .findSampleCol(df)
    ansRanges <- .ansRangeNames(df)
    strictRanges <- Filter(function(x) !is.logical(x), ansRanges)
    RangeInfo <- c(strictRanges, list(split.field = split.field))
    numInfo <- df[, !(names(df) %in% RangeInfo)]
    numAssays <- ncol(numInfo)
    nameAssays <- names(numInfo)
    if (is(df, "DataFrame"))
        numInfo <- IRanges::splitAsList(numInfo, df[[split.field]])
    else
        numInfo <- base::split(numInfo, df[[split.field]])
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
    newSE <- SummarizedExperiment(assays = SimpleList(countList),
        rowRanges = rowRanges)
    metadata(newSE) <- metadat
    return(newSE)
}

.makeRaggedExperimentFromDataFrame <- function(df, ...) {
    args <- list(...)
    if (!is.null(args[["build"]]))
        build <- args[["build"]]
    metadat <- metadata(df)
    split.field <- .findSampleCol(df)
    ansRanges <- .ansRangeNames(df)
    rangeInfo <- c(ansRanges, list(split.field = split.field))
    dropIdx <- .omitAdditionalIdx(df, ansRanges)
    if (length(dropIdx))
        df <- df[, -dropIdx]
    newGRL <- do.call(makeGRangesListFromDataFrame,
            args = c(list(df = df, keep.extra.columns = TRUE), rangeInfo))
    if (exists("build"))
        GenomeInfoDb::genome(newGRL) <- build
    newRE <- RaggedExperiment::RaggedExperiment(newGRL)
    metadata(newRE) <- metadat
    return(newRE)
}

.omitAdditionalIdx <- function(object, rangeNames) {
    rangeNames <- Filter(function(x) !is.logical(x), rangeNames)
    omitAdditional <- c("seqnames", "seqname", "chromosome", "chrom",
        "chromosome_name", "ranges", "seqlevels", "seqlengths", "seq_id",
        "iscircular", "start", "end", "width", "element", "chr")
    diffNames <- setdiff(omitAdditional, tolower(rangeNames))
    which(tolower(names(object)) %in% diffNames)
}

## Genome build from FILENAME
## RSE helper function from genome symbols to build (RNASeq, ExpSets)

.extractList <- function(object, ...) {
    args <- list(...)
    type <- args[["type"]]
    for (i in seq_along(object))
    object[[i]] <- TCGAextract(object[[i]], type)
    return(object)
}

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
#' Choices include: "RNAseqGene",
#' "Clinic", "miRNASeqGene", "RNASeq2GeneNorm", "CNASNP", "CNVSNP",
#' "CNASeq", "CNACGH", "Methylation", "Mutation", "mRNAArray",
#' "miRNAArray", "RPPAArray", "GISTICA", "GISTICT". The "GISTICA" type of
#' dataset represents GISTIC data by all genes. "GISTICT" represents data
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
#' cm <- TCGAextract(coadmut, "Mutation")
#' }
#' @importClassesFrom RTCGAToolbox FirehoseData FirehosemRNAArray
#' FirehoseCGHArray FirehoseMethylationArray
#' @export TCGAextract
TCGAextract <- function(object, type = c("Clinical", "RNASeqGene",
    "miRNASeqGene", "RNASeq2GeneNorm", "CNASNP", "CNVSNP", "CNASeq",
    "CNACGH", "Methylation", "Mutation", "mRNAArray", "miRNAArray",
    "RPPAArray", "GISTICA", "GISTICT"), ...) {
    if (length(type) != 1L)
        stop("Please specify a single data type")
    rangeslots <- c("CNVSNP", "CNASNP", "CNAseq", "CNACGH", "Mutation")
    if (!is(object, "DataFrame") && !is.data.frame(object))
        object <- .removeShell(object, type)
    if (!length(object)) { return(object) }
    if (is.list(object) && !is.data.frame(object)) {
        object <- .unNestList(object)
    }
    if (type == "Clinical") { return(object) }
    if (is(object, "matrix")) {
        return(SummarizedExperiment(assays = SimpleList(object)))
    }
    if (is(object, "list") && !is(object, "DataFrame")) {
        return(.extractList(object, type = type, ...))
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
            object <- .makeRaggedExperimentFromDataFrame(object,
                build = if (exists("build")) { build } else { NULL })
        }
        return(object)
    } else {
        object <- .standardizeBC(object)
        metadat <- metadata(object)
        object <- SummarizedExperiment(assays = SimpleList(object))
        metadata(object) <- metadat
    }
    gisticType <- grepl("^GISTIC", type, ignore.case = TRUE)
    if (gisticType) {
        slotreq <- switch(type, GISTICA = "AllByGene",
                          GISTICT = "ThresholdedByGene")
        result <- .getGISTIC(object, slotreq)
        return(result)
    }
    if (type == "Methylation") {
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
        dm <- .standardizeBC(dm)
        object <- SummarizedExperiment::SummarizedExperiment(
            assays = SimpleList(dm), rowData = annote)
        }
    return(object)
}
