#' Create a sampleMap from an experiment list and phenoData dataframe
#'
#' This function helps create a sampleMap in preparation of a
#' \code{MultiAssayExperiment} object
#'
#' @param exlist A named \code{list} of experiments compatible with the
#' MultiAssayExperiment API
#' @param mPheno A \code{data.frame} of clinical data with patient identifiers
#' as rownames
#' @return A \code{DataFrame} class object of mapped samples and patient
#' identifiers including assays
#'
#' @author Marcel Ramos \email{mramos09@gmail.com}
#'
#' @export TCGAgenerateMap
TCGAgenerateMap <- function(exlist, mPheno) {
    if (requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    exlist <- MultiAssayExperiment::Elist(exlist)
    samps <- as.list(colnames(exlist))
    } else {
    samps <- lapply(exlist, colnames)
    warning("attempting to use colnames on each experiment")
    }
    listM <- lapply(seq_along(samps), function(i, x) {
        S4Vectors::DataFrame(assay = x[[i]], assayname = names(x)[i])
    }, x = samps)
    full_map <- do.call(S4Vectors::rbind, listM)
    # matches <- match(full_map$assay, rownames(mPheno))
    matches <- match(TCGAbarcode(full_map$assay), rownames(mPheno))
    if (all(is.na(matches))) {
        stop("no way to map pData to Elist")
    }
    primary <- rownames(mPheno)[matches]
    autoMap <- S4Vectors::cbind(DataFrame(primary), full_map)
    if (any(is.na(autoMap$primary))) {
        notFound <- autoMap[is.na(autoMap$primary), ]
        warning("Data from rows:",
                sprintf("\n %s - %s", notFound[, 2], notFound[, 3]),
                "\ndropped due to missing phenotype data")
    }
    autoMap <- autoMap[!is.na(autoMap$primary), ]
    return(autoMap)
}
