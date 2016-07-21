#' Create a sampleMap from an experiment list and phenoData dataframe
#'
#' This function helps create a sampleMap in preparation of a
#' \code{MultiAssayExperiment} object
#'
#' @param exlist A named \code{list} of experiments compatible with the
#' MultiAssayExperiment API
#' @param mPheno A \code{data.frame} of clinical data with patient identifiers
#' as rownames
#' @param idConverter A function to be used against the sample or specimen
#' identifiers to match those in the rownames of the \code{pData} (default NULL)
#' @param ... Additonal arguments to pass to the 'idConverter' function.
#'
#' @return A \code{DataFrame} class object of mapped samples and patient
#' identifiers including assays
#'
#' @author Marcel Ramos \email{mramos09@gmail.com}
#'
#' @examples \dontrun{
#' ## For TCGA data
#' newMap <- generateMap(myExpList, myPheno, TCGAbarcode)
#' }
#'
#' @export generateMap
generateMap <- function(exlist, mPheno, idConverter = NULL, ...) {
    if (requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    exlist <- MultiAssayExperiment::ExperimentList(exlist)
    samps <- as.list(colnames(exlist))
    } else {
    samps <- lapply(exlist, colnames)
    warning("attempting to use colnames on each experiment")
    }
    listM <- lapply(seq_along(samps), function(i, x) {
        S4Vectors::DataFrame(colname = x[[i]], assay = names(x)[i])
    }, x = samps)
    full_map <- do.call(S4Vectors::rbind, listM)
    if (is.null(idConverter)) {
        matches <- match(full_map$colname, rownames(mPheno))
    } else {
        matches <- match(idConverter(full_map$colname, ...), rownames(mPheno))
    }
    if (all(is.na(matches))) {
        stop("no way to map pData to ExperimentList")
    }
    assay <- full_map$assay
    primary <- rownames(mPheno)[matches]
    colname <- full_map$colname
    autoMap <- S4Vectors::DataFrame(assay, primary, colname)
    if (any(is.na(autoMap$primary))) {
        notFound <- autoMap[is.na(autoMap$primary), ]
        warning("Data dropped from rows:",
                sprintf("\n %s - %s", notFound[, 2], notFound[, 3]),
                "\ndue to missing phenotype data")
    }
    autoMap <- autoMap[!is.na(autoMap$primary), ]
    return(autoMap)
}
