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
#' This function obtains and processes data from a
#' \code{\linkS4class{FirehoseData}} object. Processing involves
#' converting raw data to structured S4 Bioconductor objects. The function
#' returns an \linkS4class{ExpressionSet} or \linkS4class{GRangesList} class
#' object.
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
#' @examples
#'
#' \dontrun{
#' b2 <- extract(a2, "Methylation")
#' }
#'
#' @export extract
extract <- function(object, type = NULL) {
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
        slotreq <- grep(paste0("^", type) , slotNames(object),
                        ignore.case=TRUE, perl=TRUE, value=TRUE)
        if (is(getElement(object, slotreq), "list")) {
            elemlength <- length(getElement(object, slotreq))
            if (elemlength > 1L) {
                if (interactive()) {
                    sourceName <- sapply(getElement(object, slotreq),
                                         function(FHarray) {
                                             getElement(FHarray, "Filename")
                                         })
                    dimensions <- sapply(lapply(getElement(object, slotreq),
                                                function(tmp) {
                                                    getElement(tmp, "DataMatrix")
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
                    message("Selecting file: [", fileNo, "] ", sourceName[fileNo])
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
        stop("There is no data for that data type!")
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
        } else if (slotreq %in% rangeslots) {
            colnames(dm) <- tolower(colnames(dm))
            mygrl <- makeGRangesList(dm)
            if(exists("sourceName")) {
                mygrl@metadata <- list("fileName" = sourceName[fileNo])
            }
            return(mygrl)
        }
        colnames(dm) <- barcode(colnames(dm), sample=TRUE, collapse=TRUE)
        eset <- ExpressionSet(dm)
        if (exists("annote")) {
            featureData(eset) <- AnnotatedDataFrame(annote)
        }
        return(eset)
    }
}
