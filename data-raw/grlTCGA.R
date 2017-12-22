library(RTCGAToolbox)

blca <- getFirehoseData("BLCA", clinical = FALSE, CNASeq = TRUE)
bl <- getData(blca, "CNASeq")
blsplit <- lapply(split(bl, bl[["Sample"]]), function(x)
    x[sample(seq_len(nrow(x)), 2L), ])

blframe <- dplyr::bind_rows(blsplit)
blframe <- blframe[c(TRUE, TRUE, FALSE, FALSE), ]

write.table(blframe, file = "inst/extdata/grlTCGA.txt")
