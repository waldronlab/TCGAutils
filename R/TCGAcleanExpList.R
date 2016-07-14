#' Remove unmatched samples from a list of experiments
#'
#' This function is intended to drop any unmatched samples from all of the
#' listed experiments that are not present in the rownames of the pheno data.
#'
#' @param exlist A named \code{\link[MultiAssayExperiment]{ExperimentList}} of
#' experiments compatible with the \code{MultiAssayExperiment} API
#' @param mPheno A \code{data.frame} of clinical data with patient identifiers
#' as rownames
#' @return A named \code{list} of experiments
#'
#' @author Marcel Ramos \email{mramos09@@gmail.com}
#'
#' @export TCGAcleanExpList
TCGAcleanExpList <- function(exlist, mPheno) {
    if (requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    exlist <- MultiAssayExperiment::ExperimentList(exlist)
    sampNames <- as.list(colnames(exlist))
    } else {
    sampNames <- lapply(exlist, colnames)
    warning("attempting to use colnames on each experiment")
    }
    patientIDS <- rownames(mPheno)
    filler <- substr(patientIDS[1], 5, 5)
    if (filler != "-") {
        patientIDS <- gsub(paste0("\\", filler), "-", patientIDS)
    }
    logicSub <- lapply(sampNames, function(assay) {
        TCGAbarcode(assay) %in% patientIDS
    })
    newElist <- mapply(function(x, y) {x[, y]}, x = exlist, y = logicSub,
                       SIMPLIFY = FALSE)
    if (requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    newElist <- MultiAssayExperiment::ExperimentList(newElist)
    }
    return(newElist)
}
