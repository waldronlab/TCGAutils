## Download example dataset from legacy archive
if (!requireNamespace("GenomicDataCommons"))
    stop("Please download 'GenomicDataCommons' to update file")

library(GenomicDataCommons)

manifile <- files() |>
    filter(~ file_id == "d56a5dec-cb55-457f-8d93-dd1f3911ae9f") |>
        manifest()

gdcdata(manifile[["id"]], use_cached = TRUE)

flist <- list.files(gdc_cache(), pattern = "cation.txt$", recursive = TRUE,
    full.names = TRUE)
flist <- flist[grepl("^unc", basename(flist))]

exonFile <- "bt.exon_quantification.txt"
file.rename(flist, exonFile)

exonEx <- read.delim(exonFile, nrows = 100)

write.table(exonEx, file.path("inst", "extdata", basename(exonFile)),
    sep = "\t", row.names = FALSE)
