## Download example dataset from legacy archive
system("gdc-client download --manifest ./inst/resources/gdc_manifest_20171221_000206.txt -d ./data-raw/")
exonFile <- list.files("data-raw", pattern = "cation.txt$", recursive = TRUE, full.names = TRUE)
exonEx <- read.delim(exonFile, nrows = 100)
write.table(exonEx, paste0("inst/extdata/", basename(exonFile)), sep = "\t",
    row.names = FALSE)
