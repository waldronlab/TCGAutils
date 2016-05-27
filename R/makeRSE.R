#' Create a RangedSummarizedExperiment from a data frame
#'
#' This function uses input data and some regular expression to guess or
#' indicate the method for matching ranged column names. Such column names are
#' used to create \code{rowRanges} for the \link{SummarizedExperiment}.
#'
#' @param inputData A \code{data.frame} of ranged and expression data
#' @param geneCol The name of the gene symbol column as character
#' @param regEx A vector of regular expressions for "grepping" each of the 4
#' columns needed to create the necessary \link{GRanges} object
#'
#' @return A \link{SummarizedExperiment} object with rowRanges
#'
#' @author Marcel Ramos \email{mramos09@gmail.com}
#'
#' @export makeRSE
makeRSE <- function(inputData, geneCol = character(), regEx) {
    if (length(geneCol) == 0L) {
        geneCol <- "gene"
    } else if (!is.character(geneCol)) {
        stop("provide a valid gene column name")
    }
    genName <- grep(geneCol, names(inputData), value = TRUE, ignore.case = TRUE)
    grangeNames <- getRangeNames(names(inputData), regEx = regEx)
    DFranges <- inputData[, grangeNames]
    chrName <- grangeNames[1]
    if (all(grepl("^[0-9]+$", sample(DFranges[, chrName], size = 5)))) {
        DFranges[, chrName] <- paste0("chr", DFranges[, chrName])
    }
    if (!chrName %in% c("seqnames", "chr", "chrom")) {
        names(DFranges)[which(names(DFranges) == chrName)] <- "chrom"
    }
    RowRanges <- as(DFranges, "GRanges")
    if(length(unique(inputData[, genName])) == dim(inputData)[1]) {
        names(RowRanges) <- inputData[, genName]
        rNames <- inputData[, genName]
    }
    dm <- inputData[, !(names(inputData) %in% c(grangeNames, genName))]
    if(exists("rNames")) {
        rownames(dm) <- rNames
    }
    newSE <- SummarizedExperiment::SummarizedExperiment(
        assays = SimpleList(counts = dm), rowRanges = RowRanges)
    return(newSE)
}
