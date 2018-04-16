## Download example dataset from legacy archive
if (!requireNamespace("GenomicDataCommons"))
    stop("Please install 'GenomicDataCommons' to update file")

manifile <- files(legacy = TRUE) %>%
    filter(~ file_id == "d56a5dec-cb55-457f-8d93-dd1f3911ae9f") %>%
        manifest()

write.table(manifile, file = "inst/resources/gdc_manifest_20171221_000206.txt")

gdcdata(manifile[["id"]], destination_dir = "inst/extdata",
    overwrite = TRUE, progress = FALSE)

flist <- list.files("inst/extdata", pattern = "cation.txt$",
    recursive = TRUE, full.names = TRUE)
flist <- flist[grepl("^unc", basename(flist))]

exonFile <- "inst/extdata/bt.exon_quantification.txt"
file.rename(flist, exonFile)

exonEx <- read.delim(exonFile, nrows = 100)

write.table(exonEx, paste0("inst/extdata/", basename(exonFile)), sep = "\t",
    row.names = FALSE)
