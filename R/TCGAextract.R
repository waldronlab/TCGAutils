#' @importFrom psygenet2r extract
NULL

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

.removeShell <- function(x, type) {
    dataTypes <- c("RNAseq_Gene", "miRNASeq_Gene",
    "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP", "CNA_Seq", "CNA_CGH",
    "Methylation", "Mutation", "mRNA_Array", "miRNA_Array", "RPPA_Array",
    "GISTIC_A", "GISTIC_T")
    dataTypes <- gsub("_", "", dataTypes)
    type <- match.arg(type, dataTypes)
    getElement(x, type)
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

.ansRangeNames <- function(x) {
    granges_cols <- findGRangesCols(names(x), seqnames.field = "Chromosome",
        start.field = c("Start", "Start_position"),
        end.field = c("End", "End_position"))
    fielders <- list(seqnames.field = "seqnames", start.field = "start",
        end.field = "end", strand.field = "strand")
    Fargs <- lapply(fielders, function(name) { names(x)[granges_cols[[name]]] })
    Fargs[["ignore.strand"]] <- is.na(Fargs[["strand.field"]])
    Filter(function(g) {!is.na(g)}, Fargs)
}

setClassUnion("RTCGAArray", c("FirehosemRNAArray", "FirehoseCGHArray"))

setGeneric("extract", getGeneric("extract", package = "psygenet2r"))

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

#' @export
setMethod("extract", "RTCGAArray", function(object, ...) {
    dataMat <- getElement(object, "DataMatrix")
    headers <- names(dataMat)
    rangeNames <- .ansRangeNames(dataMat)
    if (length(rangeNames)) {
        sampID <- .findSampleCol(dataMat)
        rowRanges <- do.call(makeGRangesListFromDataFrame,
            args = c(list(df = dataMat, split.field = sampID,
                          keep.extra.columns = TRUE), rangeNames))
        rangeNames <- rangeNames[-match("ignore.strand", names(rangeNames))]
        rangeNames <- c(rangeNames, split.field = sampID)
        rse <- SummarizedExperiment(assays = SimpleList(dataMat[,
        -which(!names(dataMat) %in% rangeNames)]), rowRanges = rowRanges)
        return(rse)
    }
})

setMethod("extract", "ANY", function(object, ...) {
    object
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
    sNames <- slotNames(object)
    object <- .removeShell(object, type)
    if (is(object, "list")  && length(object) == 1L)
        object <- object[[1L]]
    if (is(object, "matrix"))
        return(SummarizedExperiment::SummarizedExperiment(
                assays = SimpleList(object)))
    object <- extract(object, ...)

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
            primary <- .findSampleCol(dm)
            granges_cols <-
                findGRangesCols(names(dm),
                                seqnames.field = "Chromosome",
                                start.field = c("Start", "Start_position"),
                                end.field = c("End", "End_position"))
            ans_seqnames <- names(dm)[granges_cols[["seqnames"]]]
            ans_start <- names(dm)[granges_cols[["start"]]]
            ans_end <- names(dm)[granges_cols[["end"]]]
            ans_strand <- names(dm)[granges_cols[["strand"]]]
            omitAdditional <- c("seqnames", "ranges", "seqlevels",
                                "seqlengths", "iscircular", "start", "end",
                                "width", "element", "chr")
            diffNames <- setdiff(omitAdditional,
                                 tolower(names(dm)[na.omit(granges_cols)]))
            dropIdx <- which(tolower(names(dm)) %in% diffNames)
            if (length(dropIdx)) {
                dm <- dm[, -dropIdx]
            }
            ignore.strand <- is.na(ans_strand)
            mygrl <- makeGRangesListFromTCGA(df = dm,
                                             split.field = primary,
                                             seqnames.field = ans_seqnames,
                                             start.field = ans_start,
                                             end.field = ans_end,
                                             strand.field = ans_strand,
                                             keep.extra.columns = TRUE,
                                             ignore.strand = ignore.strand)
            if(exists("sourceName")) {
                mygrl@metadata <- c(mygrl@metadata,
                                    list("fileName" = sourceName[fileNo]))
            }
            return(mygrl)
        }
        eset <- ExpressionSet(dm)
        if (exists("annote")) {
            featureData(eset) <- AnnotatedDataFrame(annote)
        }
        return(eset)
    }

