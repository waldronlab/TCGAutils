#' Convert raw Mutation data to GRangesList
#'
#' This function takes the data.frame or list of raw data and converts it to
#' a \linkS4class{GRangesList} class. Input data can be entered along with a
#' function that indicates the parameters of the data such as, range column
#' names, sample/specimen identifier column name and an ID parsing function
#' (such as for TCGA barcodes). If input data is entered as a list, all
#' element names are expected to be sample/specimen identifiers. See data
#' specific functions for more details.
#'
#' @param inputData A \code{data.frame} or \code{list} class of mutation data
#' @param dataparam A function for specifying a list of parameters to pass to
#' the function depending on the type of data
#' @param ... Additional arguments passed to the identifier parser function
#'
#' @return A \linkS4class{GRangesList} class object
#'
#' @examples \dontrun{
#' makeGRangesList(meso, tcga(primary = "Tumor_Sample_Barcode",
#'                              standard = TRUE)
#'          sample = TRUE, collapse = TRUE)
#' }
#'
#' @author Marcel Ramos \email{mramos09@gmail.com}
#'
#' @seealso tcga(), ccle()
#'
#' @importFrom GenomeInfoDb genome
#' @export makeGRangesList
makeGRangesList <- function(inputData, dataparam = NULL, ...) {
    if (is(inputData, "data.frame")) {
        names(inputData) <- tolower(names(inputData))
        if (!is.null(dataparam)) {
            sampleIndicator <- tolower(dataparam$primary)
            if (length(sampleIndicator) == 0L) {
                stop("please indicate a target sample column the data parameters")
            }
        } else {
            warning("trying to obtain target column...")
            sampleIndicator <- ifelse(is.null(inputData$sample),
                                      "tumor_sample_barcode", "sample")
        }
        ## Convert data to list for GRangesList
        f <- dataparam$idFUN
        inputData <- split(inputData, f(as.character(
            inputData[, sampleIndicator]), ...))
    }
    inputData <- lapply(inputData, function(elements) {
        names(elements) <- tolower(names(elements))
        elements
    })
    longNames <- dataparam$rangeID
    shortNames <- c("chrom", "start", "end", "strand")
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
