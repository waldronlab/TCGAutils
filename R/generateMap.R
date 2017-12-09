#' Create a sampleMap from an experiment list and phenoData dataframe
#'
#' This function helps create a sampleMap in preparation of a
#' \code{MultiAssayExperiment} object. This especially useful when the
#' sample identifiers are not very different, as in the case of TCGA barcodes.
#' An \code{idConverter} function can be provided to truncate such sample
#' identifiers and obtain patient identifiers.
#'
#' @param experiments A named \code{list} of experiments compatible with the
#' MultiAssayExperiment API
#' @param colData A \code{data.frame} of clinical data with patient identifiers
#' as rownames
#' @param idConverter A function to be used against the sample or specimen
#' identifiers to match those in the rownames of the \code{colData}
#' (default NULL)
#' @param ... Additonal arguments to pass to the 'idConverter' function.
#'
#' @return A \code{DataFrame} class object of mapped samples and patient
#' identifiers including assays
#'
#' @importFrom MultiAssayExperiment ExperimentList
#' @author Marcel Ramos, Martin Morgan, Lucas Schiffer
#'
#' @examples \dontrun{
#' ## For TCGA data
#' newMap <- generateMap(myExpList, myPheno, TCGAbarcode)
#' }
#'
#' @export generateMap
generateMap <- function(experiments, colData, idConverter = NULL, ...) {
    if (!is(experiments, "ExperimentList"))
    experiments <- ExperimentList(experiments)
    samps <- colnames(experiments)
    assay <- factor(rep(names(samps), lengths(samps)), levels=names(samps))
    colname <- unlist(samps, use.names=FALSE)
    if (is.null(idConverter)) {
        matches <- match(colname, rownames(colData))
    } else {
        matches <- match(idConverter(colname, ...), rownames(colData))
    }
    if (length(matches) && all(is.na(matches)))
        stop("no way to map colData to ExperimentList")
    primary <- rownames(colData)[matches]
    autoMap <- S4Vectors::DataFrame(assay=assay,
                                    primary=primary,
                                    colname=colname)
    missingPrimary <- is.na(autoMap[["primary"]])
    if (nrow(autoMap) && any(missingPrimary)) {
        notFound <- autoMap[missingPrimary, ]
        warning("[", sum(missingPrimary), "] rows dropped, no phenotype data:",
                sprintf("\n %s - %s", notFound[, 2], notFound[, 3]))
        autoMap <- autoMap[!is.na(autoMap$primary), ]
    }
    autoMap
}
