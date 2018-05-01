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
#' Universally Unique Identifiers (UUID) and vice versa. These relationships
#' are not one-to-one. Therefore, a \code{data.frame} is returned for all
#' inputs. The UUID to TCGA barcode translation only applies to file and case
#' UUIDs. API queries for this translation service with other types of UUIDS
#' are not fully supported. Please double check any results before using these
#' features for analysis. Case / submitter identifiers are translated by
#' default, see the \code{id_type} argument for details.
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
#' either files or cases
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
    end_point = "participant", legacy = FALSE) {
    id_type <- match.arg(id_type)
    APIendpoint <- .barcodeEndpoint(end_point)
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
    # so we can later expand to a data.frame of the right size
    barcodes_per_file <- lengths(id_list)
    # And build the data.frame
    resultFrame <- data.frame(
        id = rep(ids(info), barcodes_per_file),
        barcode = if (!length(ids(info))) character(0L) else unlist(id_list),
        row.names = NULL,
        stringsAsFactors = FALSE
    )
    names(resultFrame) <- c(id_type, APIendpoint)
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
barcodeToUUID <-  function(barcodes, id_type = c("case_id", "file_id"),
    legacy = FALSE) {
    .checkBarcodes(barcodes)
    id_type <- match.arg(id_type)
    if (id_type == "case_id") {
        targetElement <- APIendpoint <- "submitter_id"
        barcodes <- unique(TCGAbarcode(barcodes))
    } else if (id_type == "file_id") {
        targetElement <- "cases"
        APIendpoint <- .barcodeEndpoint(.findBarcodeLimit(barcodes))
    }
    funcRes <- switch(id_type,
        file_id = files(legacy = legacy),
        case_id = cases(legacy = legacy))
    info <- results_all(
        select(filter(funcRes, as.formula(
            paste("~ ", APIendpoint, "%in% barcodes")
        )),
        APIendpoint)
    )

    id_list <- lapply(info[[targetElement]], function(x) {
        x[[1]][[1]][[1]]
    })

    barcodes_per_file <- lengths(id_list)

    resultFrame <- data.frame(
        barcode = if (!length(ids(info))) character(0L) else unlist(id_list),
        ids = rep(ids(info), barcodes_per_file),
        row.names = NULL,
        stringsAsFactors = FALSE
    )
    names(resultFrame) <- c(APIendpoint, id_type)
    resultFrame <- resultFrame[
        na.omit(match(barcodes, resultFrame[[APIendpoint]])), ]
    rownames(resultFrame) <- NULL
    resultFrame
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

#' @rdname ID-translation
#'
#' @section builds:
#'
#' A couple of functions are available to search for build versions, either from
#' NCBI or UCSC. \code{translateBuild} will translate between UCSC and NCBI
#' build versions. \code{extractBuild} will use grep patterns to find the first
#' build within a string.
#'
#' @param from A build version name
#' @param to The name of the desired version
#'
#' @examples
#'
#' translateBuild("GRCh35", "UCSC")
#'
#' @export
translateBuild <- function(from, to = "UCSC") {
    if (!S4Vectors::isSingleString(to) && !S4Vectors::isSingleString(from))
        stop("Enter a single valid genomic build")
    if (!to %in% c("UCSC", "NCBI"))
        stop ("Only UCSC and NCBI supported")

    buildDF <- S4Vectors::DataFrame(
        Date = c("July 2004", "May 2004", "March 2006", "February 2009",
            "December 2013"),
        NCBI = c("34", "35", "36", "37", "38"),
        UCSC = c("hg16", "hg17", "hg18", "hg19", "hg38")
    )
    if (to == "UCSC")
        from <- gsub("[GgRrCcHh]", "", from)
    matchBuild <- switch (to, UCSC = "NCBI", NCBI = "UCSC" )
    buildIndex <- match(from, buildDF[[matchBuild]])
    if (is.na(buildIndex)) {
        warning("build could not be matched")
        return(NA_character_)
    }
    if (to == "NCBI")
        paste0("GRCh", buildDF[[to]][buildIndex])
    else
        buildDF[[to]][buildIndex]
}

#' @rdname ID-translation
#'
#' @param string A single character string
#' @param build A vector of build version names (default UCSC, NCBI)
#'
#' @examples
#'
#' extractBuild(
#' "SCENA_p_TCGAb29and30_SNP_N_GenomeWideSNP_6_G05_569110.nocnv_grch38.seg.txt"
#' )
#'
#' @export
extractBuild <- function(string, build = c("UCSC", "NCBI")) {
    if (!S4Vectors::isSingleString(string))
        stop("Provide a single string for build search")
    builds <- vector(mode = "character", length(build))
    names(builds) <- build
    for (i in build) {
        pattrn <- switch(i, UCSC = "[Hh][Gg][0-9]{2}",
            NCBI = "[Gg][Rr][Cc][Hh][0-9]{2}")
        builds[[i]] <- stringr::str_extract(string, pattrn)
    }
    builds <- Filter(function(x) !is.na(x), builds)
    if (!length(builds))
        NA_character_
    else if (length(builds))
        builds[1L]
}
