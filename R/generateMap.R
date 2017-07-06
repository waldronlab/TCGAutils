#' Create a sampleMap from an experiment list and phenoData dataframe
#'
#' This function helps create a sampleMap in preparation of a
#' \code{MultiAssayExperiment} object
#'
#' @param experiments A named \code{list} of experiments compatible with the
#' MultiAssayExperiment API
#' @param colData A \code{data.frame} of clinical data with patient identifiers
#' as rownames
#' @param idConverter A function to be used against the sample or specimen
#' identifiers to match those in the rownames of the \code{colData}
#' (default NULL)
#' @param force Force use of \code{colnames} to all the data list elements for
#' map creation
#' @param ... Additonal arguments to pass to the 'idConverter' function.
#'
#' @return A \code{DataFrame} class object of mapped samples and patient
#' identifiers including assays
#'
#' @author Marcel Ramos \email{marcel.ramos@roswellpark.org}, Martin Morgan, Lucas Schiffer
#'
#' @examples \dontrun{
#' ## For TCGA data
#' newMap <- generateMap(myExpList, myPheno, TCGAbarcode)
#' }
#'
#' @export generateMap
generateMap <- function(experiments, colData, idConverter = NULL, force = FALSE, ...) {
    if (is.null(names(experiments)))
        stop("experiments list/List must be named")
    if (!is(experiments, "ExperimentList") && !force) {
        experiments <- MultiAssayExperiment::ExperimentList(experiments)
        samps <- colnames(experiments)
    } else { samps <- lapply(experiments, colnames) }
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
    if (nrow(autoMap) && any(is.na(autoMap$primary))) {
        notFound <- autoMap[is.na(autoMap$primary), ]
        warning("Data from rows:",
                sprintf("\n %s - %s", notFound[, 2], notFound[, 3]),
                "\ndropped due to missing phenotype data")
        autoMap <- autoMap[!is.na(autoMap$primary), ]
    }
    autoMap
}
