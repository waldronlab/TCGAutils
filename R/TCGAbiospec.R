.strsep <- function(text, pos) {
    stopifnot(length(unique(nchar(text))) == 1L)
    lengthText <- unique(nchar(text))
    allIndx <- seq_len(lengthText)
    stopifnot(pos %in% allIndx)
    fgroup <- seq_len(pos)
    sgroup <- allIndx[!allIndx %in% fgroup]
    list(
        substr(text, min(fgroup), max(fgroup)),
        substr(text, min(sgroup), max(sgroup))
    )
}

#' Extract biospecimen data from the TCGA barcode
#'
#' This function uses the full TCGA barcode to return a data frame of the
#' data pertinent to laboratory variables such as vials, portions, analytes,
#' plates and the center.
#'
#' @param barcodes A character vector of TCGA barcodes
#' @return A \code{dataframe} with sample type, sample code, portion, plate,
#' and center columns.
#'
#' @author M. Ramos
#'
#' @examples
#' example("TCGAbarcode")
#' TCGAbiospec(barcodes)
#'
#' @export TCGAbiospec
TCGAbiospec <- function(barcodes) {
    .checkBarcodes(barcodes)
    filler <- .uniqueDelim(barcodes)
    maxIndx <- unique(lengths(strsplit(barcodes, filler)))
    if (maxIndx < 4L)
        stop("Provide a longer barcode")

    local_data_store <- new.env(parent = emptyenv())
    data("sampleTypes", envir = local_data_store, package = "TCGAutils")
    sampleTypes <- local_data_store[["sampleTypes"]]
    sampCode <- TCGAbarcode(barcodes, FALSE, TRUE)
    sampVial <- .strsep(sampCode, 2L)
    names(sampVial) <- c("sample", "vial")
    sample_definition <- sampleTypes[["Definition"]][
        match(sampVial[["sample"]], sampleTypes[["Code"]])]
    biospec <-
        data.frame(
            submitter_id = TCGAbarcode(barcodes),
            sample_definition,
            as.data.frame(sampVial, stringsAsFactors = FALSE),
            stringsAsFactors = FALSE
        )
    if (maxIndx == 4L)
        return(biospec)
    else
        splitDex <- seq(5L, maxIndx)

    tailBarcode <- strsplit(TCGAbarcode(barcodes, index = splitDex), filler)
    splitCol <- splitDex == 5L
    tailBarcode <- lapply(tailBarcode, function(x)
        c(unlist(.strsep(x[[1L]], 2L)), x[!splitCol]))

    portPlateCent <- do.call(rbind.data.frame,
        args = c(tailBarcode, list(stringsAsFactors = FALSE)))
    names(portPlateCent) <-
        c("portion", "analyte", "plate", "center")[seq_along(c(splitDex, 1L))]

    cbind.data.frame(biospec, portPlateCent, stringsAsFactors = FALSE)
}
