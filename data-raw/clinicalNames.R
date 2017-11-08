# Locate Clinical datasets for each cancer
# Script used with https://github.com/waldronlab/MultiAssayExperiment-TCGA

library(TCGAutils)
library(RTCGAToolbox)
data(diseaseCodes)

TCGAcodes <- RTCGAToolbox::getFirehoseDatasets()
excludedCodes <- c("COADREAD", "GBMLGG", "KIPAN", "STES", "FPPP", "CNTL",
                   "LCML", "MISC")
TCGAcodes <- TCGAcodes[-which(TCGAcodes %in% excludedCodes)]

myDataDir <- "data/Clinical"

if (!dir.exists(myDataDir))
    dir.create(myDataDir, recursive = TRUE)

lapply(TCGAcodes, function(cancer) {
    if (!file.exists(file.path(myDataDir, cancer, "clinical.csv"))) {
        clinDat <- RTCGAToolbox::getFirehoseData(dataset = cancer,
                                                 destdir = tempfile())
        clinFrame <- RTCGAToolbox::getData(clinDat, "clinical")
        rownames(clinFrame) <- .standardBarcodes(rownames(clinFrame))

        dir.create(file.path(myDataDir, cancer))

        write.csv(clinFrame, file.path(myDataDir, cancer, "clinical.csv"))
        message(cancer, " clinical data saved.")
    } else {
        message(cancer, " clinical data already exists!")
    }
})

names(TCGAcodes) <- TCGAcodes

clinicalNames <- IRanges::CharacterList(lapply(TCGAcodes, function(cancer) {
    clinDat <- read.csv(file.path(myDataDir, cancer, "clinical.csv"),
        row.names = 1L)
    allNA <- vapply(clinDat, function(col) all(is.na(col)), logical(1L))
    clinDat <- clinDat[, !allNA]
    names(clinDat)[names(clinDat) != "Composite.Element.REF"]
}))

devtools::use_data(clinicalNames)
