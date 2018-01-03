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
#' @title Translate study identifiers from barcode to UUID and vice versa
#'
#' @description These functions allow the user to enter a character vector of
#' identifiers and use the GDC API to translate from TCGA barcodes to
#' Universally Unique Identifiers (UUID) and vice versa.
#'
#' @details
#' The \code{end_point} options reflect endpoints in the Genomic Data Commons
#' API. These are summarized as follows:
#' \itemize{
#'   \item{participant}: This default snippet of information includes project,
#'   tissue source site (TSS), and participant number
#'   (barcode format: TCGA-XX-XXXX)
#'   \item{sample}: This adds the sample information to the participant barcode
#'   (TCGA-XX-XXXX-11X)
#'   \item{portion, analyte}: Either of these options adds the portion and
#'   analyte information to the sample barcode (TCGA-XX-XXXX-11X-01X)
#'   \item{plate, center}: Additional plate and center information is returned,
#'   i.e., the full barcode (TCGA-XX-XXXX-11X-01X-XXXX-XX)
#' }
#' Only these keywords need to be used to target the specific barcode endpoint.
#'
#' @param file_ids A \code{character} vector of UUIDs corresponding to
#' files
#' @param end_point The cutoff point of the barcode that should be returned.
#' See details for options.
#' @param legacy (logical default FALSE) whether to search the legacy archives
#'
#' @return A \code{data.frame} of TCGA barcode identifiers and UUIDs
#'
#' @examples
#' ## Translate UUIDs >> TCGA Barcode
#'
#' uuids <- c("0001801b-54b0-4551-8d7a-d66fb59429bf",
#' "002c67f2-ff52-4246-9d65-a3f69df6789e",
#' "003143c8-bbbf-46b9-a96f-f58530f4bb82")
#'
#' UUIDtoBarcode(uuids, "sample")
#'
#' @author Sean Davis, M. Ramos
#'
#' @export UUIDtoBarcode
UUIDtoBarcode <-  function(file_ids, end_point = "participant", legacy = FALSE) {
    APIendpoint <- .barcodeEndpoint(end_point)
    filesres <- files(legacy = legacy)
    info <- results_all(
        select(filter(filesres, ~ file_id %in% file_ids),
        APIendpoint)
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
    resultFrame <- data.frame(
        file_id = rep(ids(info), barcodes_per_file),
        barcode = if (!length(ids(info))) character(0L) else unlist(id_list),
        row.names = NULL,
        stringsAsFactors = FALSE
    )
    names(resultFrame)[2L] <- APIendpoint
    resultFrame
}

#' @rdname ID-translation
#'
#' @param barcodes A \code{character} vector of TCGA barcodes
#'
#' @examples
#' ## Translate TCGA Barcode >> UUIDs
#'
#' fullBarcodes <- c("TCGA-B0-5117-11A-01D-1421-08",
#' "TCGA-B0-5094-11A-01D-1421-08",
#' "TCGA-E9-A295-10A-01D-A16D-09")
#'
#' sample_ids <- TCGAbarcode(fullBarcodes, sample = TRUE)
#'
#' barcodeToUUID(sample_ids)
#'
#' participant_ids <- c("TCGA-CK-4948", "TCGA-D1-A17N",
#' "TCGA-4V-A9QX", "TCGA-4V-A9QM")
#'
#' barcodeToUUID(participant_ids)
#'
#' @export barcodeToUUID
barcodeToUUID <-  function(barcodes, legacy = FALSE) {
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

    resultFrame <- data.frame(
        barcode = if (!length(ids(info))) character(0L) else unlist(id_list),
        file_id = rep(ids(info), barcodes_per_file),
        row.names = NULL,
        stringsAsFactors = FALSE
    )
    resultFrame <- resultFrame[resultFrame[["barcode"]] %in% barcodes, ]
    rownames(resultFrame) <- NULL
    resultFrame
}
