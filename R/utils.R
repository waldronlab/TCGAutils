## Helper for finding barcode column
## **Takes the first result!**
.findBarcodeCol <- function(DF) {
    apply(DF, 2, function(column) {
        logicBCode <- grepl("^TCGA", column)
        logicBCode
    }) %>% apply(., 2, all) %>% Filter(isTRUE, .) %>% names %>% `[[`(1L)
}

## Standardize barcode format
.stdIDs <- function(sampleBarcode) {
    bcodeTest <- grepl("\\.", sample(sampleBarcode, 10L, replace = TRUE))
    if (all(bcodeTest))
        sampleBarcode <- gsub("\\.", "-", sampleBarcode)
    toupper(sampleBarcode)
}

## Find columns that are all NA
.findNAColumns <- function(dataset) {
    apply(dataset, 2, function(column) {
              all(is.na(column))
                              })
}
