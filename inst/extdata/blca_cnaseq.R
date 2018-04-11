## Generate blca_cnaseq data
if (!requireNamespace("RTCGAToolbox"))
    stop("Download package 'RTCGAToolbox' to regenerate data")

blca <- RTCGAToolbox::getFirehoseData("BLCA", clinical = FALSE,
    CNASeq = TRUE, destdir = tempdir())
bl <- RTCGAToolbox::getData(blca, "CNASeq")
blsplit <- lapply(split(bl, bl[["Sample"]]), function(x)
    x[sample(seq_len(nrow(x)), 2L), ])

blframe <- dplyr::bind_rows(blsplit)
blca_cnaseq <- blframe[c(TRUE, TRUE, FALSE, FALSE), ]

write.table(blca_cnaseq, file = "blca_cnaseq.txt")
