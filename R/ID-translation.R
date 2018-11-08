## function to figure out exact endpoint based on TCGA barcode
.barcode_files <- function(startPoint = "cases", submitter_id = TRUE) {
    keywords <- c("cases", "samples", "portions", "analytes", "aliquots")
    last <- match.arg(startPoint, keywords)
    indx <- seq_len(which(keywords == last))
    sub_id <- if (submitter_id) "submitter_id" else NULL
    paste(c(keywords[indx], sub_id), collapse = ".")
}

.subword_id <- function(keyword) {
    ret <- paste0(keyword, "_ids")
    setNames(paste0("submitter_", ret), ret)
}

.barcode_cases <- function(bcodeType = "case") {
    if (identical(bcodeType, "case"))
        setNames("submitter_id", "case_id")
    else
        .subword_id(bcodeType)
}

.findBarcodeLimit <- function(barcode) {
    filler <- .uniqueDelim(barcode)
    splitCodes <- strsplit(barcode, filler)
    maxIndx <- unique(lengths(splitCodes))
    if (!S4Vectors::isSingleInteger(maxIndx))
        stop("Inconsistent barcode lengths found")
    if (maxIndx < 3L)
        stop("Minimum barcode fields required: 'project-TSS-participant'")
    key <- c(rep("case", 3L), "sample", "analyte", "aliquot", "aliquot")[maxIndx]
    if (identical(key, "analyte")) {
        analyte_chars <- unique(
            vapply(splitCodes, function(x) nchar(x[[maxIndx]]), integer(1L))
        )
        if (!S4Vectors::isSingleInteger(analyte_chars))
            stop("Inconsistent '", key, "' barcodes")
        if (analyte_chars < 3)
            key <- "portion"
    } else if (identical(key, "aliquot")) {
        if (identical(maxIndx, 6L)) {
            ali_chars <- vapply(splitCodes, function(x)
                nchar(x[c(maxIndx-1L, maxIndx)]), integer(2L))
            if (identical(ali_chars, c(2L, 3L)))
                key <- "slide"
        }
    }
    key
}

.buildIDframe <- function(info, id_list) {
    barcodes_per_file <- lengths(id_list)
    # And build the data.frame
    data.frame(
        id = rep(ids(info), barcodes_per_file),
        barcode = if (!length(ids(info))) character(0L) else unlist(id_list),
        row.names = NULL,
        stringsAsFactors = FALSE
    )
}

#' @name ID-translation
#'
#' @title Translate study identifiers from barcode to UUID and vice versa
#'
#' @description These functions allow the user to enter a character vector of
#' identifiers and use the GDC API to translate from TCGA barcodes to
#' Universally Unique Identifiers (UUID) and vice versa. These relationships
#' are not one-to-one. Therefore, a \code{data.frame} is returned for all
#' inputs. The UUID to TCGA barcode translation only applies to file and case
#' UUIDs. Two-way UUID translation is available from 'file_id' to 'case_id'
#' and vice versa. Please double check any results before using these
#' features for analysis. Case / submitter identifiers are translated by
#' default, see the \code{id_type} argument for details. All identifiers are
#' converted to lower case.
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
#' These endpoints only apply to "file_id" type translations to TCGA barcodes
#' (see \code{id_type} argument).
#'
#' @param id_vector A \code{character} vector of UUIDs corresponding to
#' either files or cases (default assumes case_ids)
#' @param id_type Either \code{case_id} or \code{file_id} indicating the type of
#' \code{id_vector} entered (default "case_id")
#' @param end_point The cutoff point of the barcode that should be returned,
#' only applies to \code{file_id} type queries. See details for options.
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
#' UUIDtoBarcode(uuids, id_type = "file_id", end_point = "sample")
#'
#' UUIDtoBarcode("ae55b2d3-62a1-419e-9f9a-5ddfac356db4", id_type = "case_id")
#'
#' @author Sean Davis, M. Ramos
#'
#' @export UUIDtoBarcode
UUIDtoBarcode <-  function(id_vector, id_type = c("case_id", "file_id"),
    end_point = "cases", legacy = FALSE) {
    id_type <- match.arg(id_type)
    APIendpoint <- .barcode_files(end_point)
    if (id_type == "case_id") {
        targetElement <- APIendpoint <- "submitter_id"
    } else if (id_type == "file_id") {
        targetElement <- "cases"
    }
    funcRes <- switch(id_type,
        file_id = files(legacy = legacy),
        case_id = cases(legacy = legacy))
    info <- results_all(
        select(filter(funcRes, as.formula(
            paste("~ ", id_type, "%in% id_vector")
            )),
        APIendpoint)
    )
    # The mess of code below is to extract TCGA barcodes
    # id_list will contain a list (one item for each file_id)
    # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
    id_list <- lapply(info[[targetElement]], function(x) {
        x[[1]][[1]][[1]]
    })

    rframe <- .buildIDframe(info, id_list)
    names(rframe) <- c(id_type, APIendpoint)
    rframe
}

#' @rdname ID-translation
#'
#' @param to_type The desired UUID type to obtain, can either be "case_id" or
#' "file_id"
#'
#' @examples
#' ## Translate file UUIDs >> case UUIDs
#'
#' uuids <- c("0001801b-54b0-4551-8d7a-d66fb59429bf",
#' "002c67f2-ff52-4246-9d65-a3f69df6789e",
#' "003143c8-bbbf-46b9-a96f-f58530f4bb82")
#'
#' UUIDtoUUID(uuids)
#'
#' @export UUIDtoUUID
UUIDtoUUID <- function(id_vector, to_type = c("case_id", "file_id"),
    legacy = FALSE) {
    id_vector <- tolower(id_vector)
    type_ops <- c("case_id", "file_id")
    to_type <- match.arg(to_type)
    from_type <- type_ops[!type_ops %in% to_type]
    if (!length(from_type))
        stop("Provide a valid UUID type")

    endpoint <- switch(to_type,
        case_id = "cases.case_id",
        file_id = "files.file_id")
    apifun <- switch(to_type,
        file_id = cases(legacy = legacy),
        case_id = files(legacy = legacy))
    info <- results_all(
        select(filter(apifun, as.formula(
            paste("~ ", from_type, "%in% id_vector")
            )),
        endpoint)
    )
    targetElement <- gsub("(\\w+).*", "\\1", endpoint)
    id_list <- lapply(info[[targetElement]], function(x) {x[[1]]})

    rframe <- .buildIDframe(info, id_list)
    names(rframe) <- c(from_type, endpoint)
    rframe
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
barcodeToUUID <-
    function(barcodes, legacy = FALSE)
{
    .checkBarcodes(barcodes)
    bend <- .findBarcodeLimit(barcodes)
    endtargets <- .barcode_cases(bend)
    expander <- gsub("cases\\.", "", .barcode_files(paste0(bend, "s"), FALSE))
    expander <- if (identical(expander, "cases")) "samples" else expander

    info <- results_all(
        expand(filter(cases(), as.formula(
            paste("~ ", endtargets, "%in% barcodes")
        )),
        expander)
    )
    rframe <- cbind.data.frame(
        unlist(info[[endtargets]]), unlist(info[[names(endtargets)]]),
        stringsAsFactors = FALSE, row.names = NULL)
    names(rframe) <- c(endtargets, names(endtargets))
    rframe[na.omit(match(barcodes, rframe[[endtargets]])), , drop = FALSE]
}

#' @rdname ID-translation
#'
#' @param filenames A \code{character} vector of filenames obtained from
#' the GenomicDataCommons
#'
#' @examples
#' library(GenomicDataCommons)
#'
#' fquery <- files() %>%
#'     filter(~ cases.project.project_id == "TCGA-COAD" &
#'         data_category == "Copy Number Variation" &
#'         data_type == "Copy Number Segment")
#'
#' fnames <- results(fquery)$file_name[1:6]
#'
#' filenameToBarcode(fnames)
#'
#' @export filenameToBarcode
filenameToBarcode <- function(filenames, legacy = FALSE) {
    filesres <- files(legacy = legacy)
    info <- results_all(
        select(filter(filesres, ~ file_name %in% filenames),
            "cases.samples.portions.analytes.aliquots.submitter_id")
    )
    id_list <- lapply(info[["cases"]], function(a) {
        a[[1L]][[1L]][[1L]]
    })
    # so we can later expand to a data.frame of the right size
    barcodes_per_file <- lengths(id_list)
    # And build the data.frame
    data.frame(file_name = rep(filenames, barcodes_per_file),
        file_id = rep(ids(info), barcodes_per_file),
        aliquots.submitter_id = unlist(id_list), row.names = NULL,
        stringsAsFactors = FALSE)
}
