# Locate Clinical datasets for each cancer
# Script used with https://github.com/waldronlab/MultiAssayExperiment-TCGA

if (!requireNamespace("RTCGAToolbox"))
    stop("Install `RTCGAToolbox` to generate 'clinicalNames' data")

TCGAcodes <- RTCGAToolbox::getFirehoseDatasets()

excludedCodes <- c("COADREAD", "GBMLGG", "KIPAN", "STES", "FPPP", "CNTL",
    "LCML", "MISC")
TCGAcodes <- TCGAcodes[-which(TCGAcodes %in% excludedCodes)]

myDataDir <- tempdir()

lapply(TCGAcodes, function(cancer) {
    if (!file.exists(file.path(myDataDir, cancer, "clinical.csv"))) {
        clinDat <- RTCGAToolbox::getFirehoseData(dataset = cancer,
            destdir = myDataDir)
        clinFrame <- RTCGAToolbox::getData(clinDat, "clinical")
        rownames(clinFrame) <-
            TCGAutils:::.standardBarcodes(rownames(clinFrame))

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

devtools::use_data(clinicalNames, overwrite = TRUE)
