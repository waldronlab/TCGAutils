#' Create a sampleMap from an experiment list and phenoData dataframe
#'
#' This function helps create a sampleMap in preparation of a
#' \code{MultiAssayExperiment} object
#'
#' @param experiments A named \code{list} of experiments compatible with the
#' MultiAssayExperiment API
#' @param pData A \code{data.frame} of clinical data with patient identifiers
#' as rownames
#' @param idConverter A function to be used against the sample or specimen
#' identifiers to match those in the rownames of the \code{pData} (default NULL)
#' @param ... Additonal arguments to pass to the 'idConverter' function.
#'
#' @return A \code{DataFrame} class object of mapped samples and patient
#' identifiers including assays
#'
#' @author Marcel Ramos \email{mramos09@gmail.com}, Lucas Schiffer
#'
#' @examples \dontrun{
#' ## For TCGA data
#' newMap <- generateMap(myExpList, myPheno, TCGAbarcode)
#' }
#'
#' @export generateMap
generateMap <- function(experiments, pData, idConverter = NULL, ...) {
    if (requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
        if (!inherits(experiments, "ExperimentsList")) {
            experiments <- MultiAssayExperiment::ExperimentList(experiments)
        }
        samps <- colnames(experiments)
    } else {
        samps <- lapply(experiments, colnames)
        warning("attempting to use colnames on each experiment")
    }
    assay <- factor(rep(names(samps), lengths(samps)), levels=names(samps))
    colname <- unlist(samps, use.names=FALSE)
    if (is.null(idConverter)) {
        matches <- match(colname, rownames(pData))
    } else {
        matches <- match(idConverter(colname, ...), rownames(pData))
    }
    if (length(matches) && all(is.na(matches)))
        stop("no way to map pData to ExperimentList")
    primary <- rownames(pData)[matches]
    autoMap <- S4Vectors::DataFrame(
        assay=assay, primary=primary, colname=colname)

    if (nrow(autoMap) && any(is.na(autoMap$primary))) {
        notFound <- autoMap[is.na(autoMap$primary), ]
        warning("Data from rows:",
                sprintf("\n %s - %s", notFound[, 2], notFound[, 3]),
                "\ndropped due to missing phenotype data")
        autoMap <- autoMap[!is.na(autoMap$primary), ]
    }
    autoMap
}
