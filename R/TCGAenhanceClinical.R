#' Take a MultiAssayExperiment and include curated variables
#'
#' This function works on the \code{colData} of a
#' \code{\linkS4class{MultiAssayExperiment}} object to add curated variable
#' columns. It is recommended that the user run the scripts in the
#' \code{MultiAssayExperiment-TCGA} repository that build the "enhanced" type
#' of data. Please see the repository's README for more information.
#'
#' @param MultiAssayExperiment A \linkS4class{MultiAssayExperiment} object
#' @param cancerCode A single string indicating the TCGA cancer code
#' (e.g., "PRAD")
#' @param repoLocation The directory location where the
#' "MultiAssayExperiment-TCGA" repository can be found
#'
#' @return A \code{\link{MultiAssayExperiment}} object
#' @export TCGAenhanceClinical
TCGAenhanceClinical <- function(MultiAssayExperiment = MultiAssayExperiment(),
                            cancerCode, repoLocation = ".") {
    stopifnot(inherits(MultiAssayExperiment, "MultiAssayExperiment"))
    stopifnot(S4Vectors::isSingleString(cancerCode))
    stopifnot(file.exists(file.path(repoLocation, "MultiAssayExperiment-TCGA")))

    clinicalDF <- colData(MultiAssayExperiment)
    enhancedDataset <-
        readr::read_csv(
            file.path(repoLocation, "MultiAssayExperiment-TCGA",
                      "inst/extdata/Clinical/enhanced/",
                      paste0(toupper(cancerCode), ".csv")))
    merge(clinicalDF, enhancedDataset, by.x = rownames(clinicalDF),
          by.y = "patientID")
}
