#' Remove unmatched samples from a list of experiments
#'
#' This function is intended to drop any unmatched samples from all of the
#' listed experiments that are not present in the rownames of the pheno data.
#'
#' @param exlist A named \linkS4class{Elist} of experiments compatible with the
#' MultiAssayExperiment API
#' @param mPheno A \code{data.frame} of clinical data with patient identifiers
#' as rownames
#' @return A named \code{list} of experiments
#'
#' @author Marcel Ramos \email{mramos09@@gmail.com}
#'
#' @export cleanExpList
cleanExpList <- function(exlist, mPheno) {
    exlist <- MultiAssayExperiment::Elist(exlist)
    sampNames <- as.list(colnames(exlist))
    patientIDS <- tolower(rownames(mPheno))
    filler <- substr(patientIDS[1], 5, 5)
    if (filler != "-") {
        patientIDS <- gsub(paste0("\\", filler), "-", patientIDS)
    }
    logicSub <- lapply(sampNames, function(assay) {
        barcode(assay) %in% patientIDS
    })
    newElist <- mapply(function(x, y) {x[, y]}, x = exlist, y = logicSub,
                       SIMPLIFY = FALSE)
    newElist <- MultiAssayExperiment::Elist(newElist)
    return(newElist)
}
