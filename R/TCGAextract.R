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
#' @author Marcel Ramos \email{mramos09@@gmail.com}
#'
#' @examples \dontrun{
#' library(RTCGAToolbox)
#' dataFolder <- normalizePath("~/Documents/data")
#' coadmut <- getFirehoseData("COAD", runDate = "20151101", Mutation = TRUE,
#'                          destdir = dataFolder)
#' cm <- TCGAextract(coadmut, "mutations")
#' }
#'
#' @seealso makeGRangesListFromTCGA()
#' @export TCGAextract
TCGAextract <- function(object, type = NULL) {
    if (!is.null(type)) {
        if (is.character(type)) {
            type <- tolower(gsub("_", "", type))
            type <- gsub("s$", "", type)
        } else {
            stop("Data type must be a character string")
        }
    } else {
        stop("Specify type")
    }
    choices <- tolower(
        gsub("_", "",
             c("RNAseq_Gene", "miRNASeq_Gene", "RNAseq2_Gene_Norm",
               "CNA_SNP", "CNV_SNP", "CNA_Seq", "CNA_CGH", "Methylation",
               "Mutation", "mRNA_Array", "miRNA_Array", "RPPA_Array")))
    rangeslots <- c("CNVSNP", "CNASNP", "CNAseq", "CNACGH", "Mutations")
    if (type %in% choices) {
        slotreq <- grep(paste0("^", type) , methods::slotNames(object),
                        ignore.case=TRUE, perl=TRUE, value=TRUE)
        if (methods::is(getElement(object, slotreq), "list")) {
            elemlength <- length(getElement(object, slotreq))
            if (elemlength > 1L) {
                if (interactive()) {
                    sourceName <- sapply(getElement(object, slotreq),
                                         function(FHarray) {
                                             getElement(FHarray, "Filename")
                                         })
                    dimensions <- sapply(lapply(getElement(object, slotreq),
                                                function(tmp) {
                                                    getElement(tmp,
                                                               "DataMatrix")
                                                }), dim)
                    cat(paste0("[", seq(length(sourceName)), "] ",
                               sourceName, paste0("\n\tNumber of rows: ",
                                                  dimensions[1,],
                                                  "\tNumber of columns: ",
                                                  dimensions[2,]) ),
                        fill = TRUE, sep = "\n")
                    fileNo <- .fileSelect()
                    if (fileNo == 0) {
                        fileNo <- which.max(sapply(
                            lapply(getElement(object, slotreq),
                                   function(tmp) {
                                       getElement(tmp, "DataMatrix")
                                   }), ncol)
                        )
                    }
                    message("Selecting file: [", fileNo, "] ",
                            sourceName[fileNo])
                    dm <- getElement(object, slotreq)[[fileNo]]@DataMatrix
                } else {
                    dm <- lapply(getElement(object, slotreq),
                                 function(tmp) {
                                     getElement(tmp, "DataMatrix")
                                 })
                    keeplist <- which.max(sapply(dm, ncol))
                    dm <- dm[[keeplist]]
                    warning(paste("Taking the array platform with",
                                  "the greatest number of samples:", keeplist))
                }
            } else if(elemlength == 1L) {
                dm <- getElement(object, slotreq)[[1]]@DataMatrix
            } else if(elemlength == 0L) {
                dm <- matrix(NA, nrow=0, ncol=0)
            }
        } else {
            dm <- getElement(object, slotreq)
        }
    } else  if (type %in% c("gistica", "gistict")) {
        if(type=="gistica"){
            slotreq <- "AllByGene"
        } else {
            slotreq <- "ThresholdedByGene"
        }
        dm <- getElement(object@GISTIC, slotreq)
    } else {
        stop(paste("Data type not yet supported or could not be matched."))
    }
    if (dim(dm)[1] == 0 | dim(dm)[2] == 0) {
        stop(type, " does not contain any data!")
    } else {
        if (slotreq %in% c("Methylation", "AllByGene", "ThresholdedByGene")) {
            annote <- dm[, !grepl("TCGA", names(dm))]
            isNumRow <- all(grepl("^[0-9]*$", rownames(dm)))
            if (isNumRow) {
                geneSymbols <- annote[, grep("symbol", names(annote),
                                             ignore.case = TRUE, value = TRUE)]
                rNames <- geneSymbols
            } else {
                rNames <- rownames(dm)
            }
            dm <- apply(dm[grepl("TCGA", names(dm))], 2, as.numeric, as.matrix)
            rownames(dm) <- rNames
            filler <- substr(colnames(dm)[1], 5, 5)
            if (filler != "-") {
                colnames(dm) <- gsub(paste0("\\", filler), "-", colnames(dm))
            }
            newSE <- SummarizedExperiment::SummarizedExperiment(
                assays = SimpleList(dm), rowData = annote)
            return(newSE)
        } else if (slotreq %in% rangeslots) {
            tsb <- match("tumor_sample_barcode", tolower(names(dm)))
            if (length(tsb) == 1L && !is.na(tsb)) {
                primary <- names(dm)[tsb]
            } else if (is.na(tsb)) {
                primary <- names(dm)[tolower(names(dm)) == "sample"]
            } else {
                stop("'partitioning.field' could not be found")
            }
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
            ignore.strand <- ifelse(is.na(ans_strand), TRUE, FALSE)
            mygrl <- makeGRangesListFromTCGA(df = dm,
                                             partitioning.field = primary,
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
}
