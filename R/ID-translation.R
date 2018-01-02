## function to figure out exact endpoint based on TCGA barcode
.barcodeEndpoint <- function(sectionLimit = "participant") {
    startPoint <- "cases"
    p.analyte <- paste0(startPoint, ".samples.portions.analytes")
    switch(sectionLimit,
        participant = paste0(startPoint, ".submitter_id"),
        sample = paste0(startPoint, ".samples.submitter_id"),
        portion = p.analyte,
        analyte = p.analyte,
        plate = paste0(p.analyte, ".aliquots.submitter_id"),
        center = paste0(p.analyte, ".aliquots.submitter_id")
    )
}

.findBarcodeLimit <- function(barcode) {
    filler <- .uniqueDelim(barcode)
    maxIndx <- unique(lengths(strsplit(barcode, filler)))
    c(rep("participant", 3L), "sample", "portion", "plate", "center")[maxIndx]
}

#' @name ID-translation
#'
#' @title Translate study identifiers from barcode to UUID
#'
#' @description This function allows the user to enter a character vector of identifiers
#' and use the GDC API to translate from TCGA barcodes to Universally Unique
#' Identifiers (UUID) and vice versa.
#'
#' @param file_ids A \code{character} vector of UUIDs corresponding to
#' files
#' @param legacy (logical default FALSE) whether to search the legacy archives
#'
#' @return A \code{character} vector of TCGA barcode identifiers
#'
#' @examples
#' ## Translate UUIDs <--> TCGA Barcode
#' uuids <- c("0001801b-54b0-4551-8d7a-d66fb59429bf",
#' "002c67f2-ff52-4246-9d65-a3f69df6789e",
#' "003143c8-bbbf-46b9-a96f-f58530f4bb82")
#'
#' UUID.barcode(uuids)
#'
#' @author Sean Davis
#'
#' @export UUID.barcode
UUID.barcode <-  function(file_ids, legacy = FALSE) {
    filesres <- files(legacy = legacy)
    info <- results_all(
        select(filter(filesres, ~ file_id %in% file_ids),
        "cases.samples.portions.analytes.aliquots.submitter_id")
    )
    # The mess of code below is to extract TCGA barcodes
    # id_list will contain a list (one item for each file_id)
    # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
    id_list <- lapply(info[["cases"]], function(x) {
        x[[1]][[1]][[1]]
    })
    # so we can later expand to a data.frame of the right size
    barcodes_per_file <- sapply(id_list, length)
    # And build the data.frame
    data.frame(
        file_id = rep(ids(info), barcodes_per_file),
        submitter_id = unlist(id_list),
        row.names = NULL,
        stringsAsFactors = FALSE
    )
}

## Still in alpha stage!

# @rdname ID-translation
#
# @examples
# ## Translate TCGA Barcode <--> UUIDs
# barcodes <- c("TCGA-B0-5117-11A-01D-1421-08",
# "TCGA-B0-5094-11A-01D-1421-08",
# "TCGA-E9-A295-10A-01D-A16D-09")
# pt_identifiers <- TCGAbarcode(barcodes)
# barcode.UUID(pt_identifiers)
#
# ## Complex Example
# exampleCodes <- c("TCGA-CK-4948", "TCGA-D1-A17N",
# "TCGA-4V-A9QX", "TCGA-4V-A9QM")
#
# barcode.UUID(exampleCodes)
#
# @export
barcode.UUID <-  function(barcodes, legacy = FALSE) {
    .checkBarcodes(barcodes)
    filesres <- files(legacy = legacy)
    lastVal <- .findBarcodeLimit(barcodes)
    selector <- .barcodeEndpoint(lastVal)
    info <- results_all(
        select(filter(filesres, as.formula(
            paste("~ ", selector, "%in% barcodes")
        )),
        selector)
    )

    id_list <- lapply(info[["cases"]], function(x) {
        x[[1]][[1]][[1]]
    })

    barcodes_per_file <- sapply(id_list, length)

    data.frame(
        file_id = rep(ids(info), barcodes_per_file),
        barcode = unlist(id_list),
        row.names = NULL,
        stringsAsFactors = FALSE
    )
}
