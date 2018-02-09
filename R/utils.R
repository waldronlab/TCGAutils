## Helper for finding barcode column
## **Takes the first result!**
.findBarcodeCol <- function(DF) {
    cnames <- names(DF)
    containsBC <- vapply(head(DF), function(column) {
        all(startsWith(column, "TCGA"))
    }, logical(1L))
    names(containsBC) <- cnames
    bcIdx <- which(containsBC)
    stopifnot(S4Vectors::isSingleInteger(which(containsBC)))
    names(containsBC)[bcIdx]
}

## Standardize barcode format
.standardBarcodes <- function(sampleBarcode) {
    if (!length(sampleBarcode)) {
        stop("<internal> Barcode must be of positive length")
    }
    sampleBC <- base::sample(sampleBarcode, 10L, replace = TRUE)
    bcodeTest <- grepl("\\.", sampleBC)
    if (all(bcodeTest))
        sampleBarcode <- gsub("\\.", "-", sampleBarcode)
    toupper(sampleBarcode)
}

## Find columns that are all NA
.findNAColumns <- function(dataset) {
    apply(dataset, 2L, function(column) {
        all(is.na(column))
    })
}
