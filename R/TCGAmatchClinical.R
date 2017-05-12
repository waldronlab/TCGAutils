#' Match Clinical with Experiment TCGA Data
#'
#' This function parses TCGA patient barcode identifiers to find matches in
#' both the experiment data and the clinical datasets. It creates a
#' corresponding clinical dataset for all of the unique patients in the
#' experiment dataset.
#'
#' @param expData Either a \code{list} or \code{data.frame} of data where
#' identifiers are element names or column names, respectively
#' @param phenoDat A \code{data.frame} of clinical variables where rownames
#' indicate patient identifiers to match on
#' @return A corresponding \code{data.frame} of clinical variables of the same
#' size as the number of unique patient identifiers with patient identifiers as
#' a "patientids" column
#'
#' @author Marcel Ramos \email{marcel.ramos@roswellpark.org}
#'
#' @export TCGAmatchClinical
TCGAmatchClinical <- function (expData, phenoDat) {
    if (methods::is(expData, "list")) {
        filler <- substr(names(expData)[1], 5, 5)
        if (filler != "-") {
            names(expData) <- gsub(paste0("\\", filler),
                                   "\\-", names(expData)) }
        getNames <- function(x){ names(x) }
    } else {
        filler <- substr(colnames(expData)[1], 5, 5)
        if (filler != "-") {
            names(expData) <- gsub(paste0("\\", filler),
                                   "\\-", colnames(expData)) }
        getNames <- function(x) { colnames(x) }
    }
    filler <- substr(rownames(phenoDat)[1], 5, 5)
    if (filler != "-") {
        rownames(phenoDat) <- gsub(paste0("\\", filler),
                                   "\\-", rownames(phenoDat)) }
    commonNames <- intersect(TCGAbarcode(getNames(expData)), TCGAbarcode(rownames(phenoDat)))
    namesRight <- getNames(expData)[match(commonNames, TCGAbarcode(getNames(expData)))]
    righttab <- TCGAbiospec(namesRight)
    clindup <- matrix(NA, nrow = length(commonNames))
    ## Rownames here use samples
    rownames(clindup) <- namesRight
    ## Match on patient identifiers
    clindup <- phenoDat[match(commonNames, rownames(phenoDat)),]
    clindup <- cbind(clindup,
                     righttab[, -which(names(righttab) == "patientids")])
    return(clindup)
}
