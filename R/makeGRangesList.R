#' Convert raw TCGA Mutation data to GRangesList
#'
#' This function takes the data.frame of raw data from the output of a TCGA
#' data pipeline and converts it to a \linkS4class{GRangesList} class.
#' Input data can be entered as either a \code{data.frame} with a sample
#' indicator titled as "tumor_sample_barcode" or "sample." If input data is a
#' entered as a list, all element names are expected to be TCGA
#' sample identifiers.
#'
#' @param inputData A \code{data.frame} or \code{list} class of TCGA
#' mutation data
#' @return A \linkS4class{GRangesList} class object
#'
#' @author Marcel Ramos \email{mramos09@gmail.com}
#'
#' @importFrom GenomeInfoDb genome
#' @export makeGRangesList
makeGRangesList <- function(inputData) {
    if (is(inputData, "data.frame")) {
    names(inputData) <- tolower(names(inputData))
    sampleIndicator <- ifelse(is.null(inputData$sample),
                              "tumor_sample_barcode", "sample")
    ## Convert data to list for GRangesList
    inputData <- split(inputData, barcode(as.character(
        inputData[, sampleIndicator]),
        sample = TRUE, collapse = TRUE))
    }
    inputData <- lapply(inputData, function(elements) {
        names(elements) <- tolower(names(elements))
        elements
    })
    longNames <- c("chromosome", "start_position", "end_position")
    shortNames <- c("chrom", "start", "end")
    twoMeta <- ifelse(all(c("num_probes", "segment_mean") %in%
                              names(inputData[[1]])), TRUE, FALSE)
    hugo <- ifelse("hugo_symbol" %in% names(inputData[[1]]), TRUE, FALSE)
    if ("ncbi_build" %in% names(inputData[[1]])) {
        ncbi_build <- Reduce(intersect, lapply(inputData,
                                               function(x)
                                                   { x[, "ncbi_build"] }))
        if (length(ncbi_build) == 1L) {
            ncbi <- ncbi_build
        } else {
            message("NCBI build was not consistent")
        }
    }
    inputData <- lapply(inputData, function(elements) {
        names(elements) <- plyr::mapvalues(names(elements),
                                           longNames, shortNames,
                                           warn_missing = FALSE)
        elements
    })
    inputData <- lapply(inputData, function(elements) {
        elements[, c("start", "end")] <-
            sapply(elements[, c("start", "end")], as.numeric)
        elements
    })
    if (!all(grepl("chr", inputData[[1]]$chrom[1:5], ignore.case = TRUE))) {
        inputData <- lapply(inputData, function(elements) {
            elements[, "chrom"] <- paste0("chr", elements[, "chrom"])
            elements
        })
    }
    metadats <- lapply(inputData, FUN = function(mydata) {
        mydata <- mydata[, !names(mydata)
                         %in%
                             c("seqnames", "ranges", "strand", "seqlevels",
                               "seqlengths", "isCircular", "start", "end",
                               "width", "element", "dataset", "chrom",
                               "center", "gene", "type", "chr",
                               "ref_allele", "tum_allele1", "tum_allele2")]
        return(mydata)
    })
    mygrl <- GRangesList(lapply(inputData, FUN = function(gr){
        NewGR <- GRanges(gr[, shortNames])
        if (twoMeta) {
            mcols(NewGR)$num_probes <- gr[, "num_probes"]
            mcols(NewGR)$segment_mean <- gr[, "segment_mean"]
        }
        if (hugo) { names(NewGR) <- gr[, "hugo_symbol"] }
        return(NewGR)
    }
    ))
    if (exists("ncbi")) {
        GenomeInfoDb::genome(mygrl) <- ncbi
    }
    metadata(mygrl) <- metadats
    return(mygrl)
}
