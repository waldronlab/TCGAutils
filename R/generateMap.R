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
#' @author M. Ramos, M. Morgan, L. Schiffer
#'
#' @examples
#' ## Minimal example
#' expList <- list(assay1 = matrix(1:6, ncol = 2L,
#'         dimnames = list(paste0("feature", 1:3), c("A-J", "B-J"))),
#'     assay2 = matrix(1:4, ncol = 2,
#'         dimnames = list(paste0("gene", 1:2), c("A-L", "B-L"))))
#'
#' ## A ficticious colData
#' myPheno <- data.frame(var1 = c("Yes", "No"), var2 = c("High", "Low"),
#'     row.names = c("a", "b"))
#'
#' ## Provide an idConverter function to take corresponding letter codes
#' ## for each biological unit
#' generateMap(expList, myPheno,
#'     idConverter = function(x) substr(tolower(x), 1L, 1L))
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
        warning("Data from rows:",
                sprintf("\n %s - %s", notFound[, 2], notFound[, 3]),
                "\ndropped due to missing phenotype data")
        autoMap <- autoMap[!is.na(autoMap$primary), ]
    }
    autoMap
}
