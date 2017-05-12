#' Remove unmatched samples from a list of experiments
#'
#' This function is intended to drop any unmatched samples from all of the
#' listed experiments that are not present in the rownames of the pheno data.
#'
#' @param experiments A named \link[MultiAssayExperiment]{ExperimentList}
#' of experiments compatible with the \code{MultiAssayExperiment} API
#' @param colData A \code{data.frame} of clinical data with patient identifiers
#' as rownames
#' @return A named \code{list} or \link[MultiAssayExperiment]{ExperimentList}
#' if the \code{MultiAssayExperiment} package is available
#'
#' @author Marcel Ramos \email{marcel.ramos@roswellpark.org}
#'
#' @export TCGAcleanExpList
TCGAcleanExpList <- function(experiments, colData) {
    if (requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    experiments <- MultiAssayExperiment::ExperimentList(experiments)
    sampNames <- as.list(colnames(experiments))
    } else {
    sampNames <- lapply(experiments, colnames)
    warning("attempting to use colnames on each experiment")
    }
    patientIDS <- rownames(colData)
    filler <- substr(patientIDS[1], 5, 5)
    if (filler != "-") {
        patientIDS <- gsub(paste0("\\", filler), "-", patientIDS)
    }
    logicSub <- lapply(sampNames, function(assay) {
        TCGAbarcode(assay) %in% patientIDS
    })
    newElist <- mapply(function(x, y) {
        x[, y, drop = FALSE]
    },
    x = experiments, y = logicSub, SIMPLIFY = FALSE)
    if (requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    newElist <- MultiAssayExperiment::ExperimentList(newElist)
    }
    return(newElist)
}

