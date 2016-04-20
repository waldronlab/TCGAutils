#' Remove unmatched samples from a list of experiments
#'
#' This function is intended to drop any unmatched samples from all of the
#' listed experiments that are not present in the rownames of the pheno data.
#'
#' @param exlist A named \code{list} of experiments compatible with the
#' MultiAssayExperiment API
#' @param mPheno A \code{data.frame} of clinical data with patient identifiers
#' as rownames
#' @return A named \code{list} of experiments
#'
#' @author Marcel Ramos \email{mramos09@@gmail.com}
#'
#' @export cleanExpList
cleanExpList <- function(exlist, mPheno) {
    sampNames <- lapply(exlist, colnames)
    PatientsInSamp <- lapply(sampNames,
                             function(bcode) {
                                 unique(barcode(bcode))
                             })
    patientID <- rownames(mPheno)
    validIDs <- Reduce(intersect, PatientsInSamp, patientID)
    logicSub <- lapply(exlist, function(elem) {
        barcode(colnames(elem)) %in% validIDs
    })
    return(Map(function(x, y) {x[, y]}, x = exlist, y = logicSub))
}
