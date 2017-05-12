.build_json_request <- function(identifier, type="file_id",
                                dataType = NULL,
                                size = 1000000) {
    options(scipen = 999)
    if (grepl("^TCGA", sample(identifier, 1L), ignore.case = TRUE))
        field <- "cases.submitter_id"
    else
        field <- switch(type, file_id = "files.file_id",
                        file_name = "files.file_name",
                        entity_id = "entity_id")
    if (!is.null(dataType)) {
        internalList <- structure(list(
            op = "and",
            content = structure(list(
                op = c("in", "="),
                content = structure(list(
                    field = c(field, "files.data_type"),
                    value = list(identifier, dataType)),
                    .Names = c("field", "value"),
                    class = "data.frame", row.names = 1:2)),
                .Names = c("op", "content"),
                class = "data.frame", row.names = 1:2)),
            .Names = c("op", "content"))
    } else {
        internalList <- list(op = "in", content = list(field = field, value = identifier))
    }
    requestList <- list(filters = internalList, format = "TSV",
                        fields = paste("file_id", "file_name",
                                       "cases.samples.portions.analytes.aliquots.submitter_id",
                                       sep = ","),
                        size = as.character(size))
    request <- jsonlite::toJSON(
        requestList,
        pretty = TRUE, auto_unbox = TRUE)
    options(scipen = 0)
    list(request = request, field = field)
}

#' Translate study identifiers from barcode to UUID and vice versa
#'
#' This function allows the user to enter a character vector of identifiers
#' and use the GDC API to translate from TCGA barcodes to Universally Unique
#' Identifiers (UUID) and vice versa.
#'
#' @param identifier A \code{character} vector of identifiers
#' @param type A single \code{character} string indicating identifier type,
#' can use "file_id" (default) or "file_name" ignored if translating TCGA
#' barcodes
#' @param dataType A single \code{character} string indicating data type
#' used for searching (e.g., "Gene Expression Quantification")
#' @return A \code{character} vector of alternate identifiers
#'
#' @examples
#' ## Translate UUIDs <--> TCGA Barcode
#' uuids <- c("0001801b-54b0-4551-8d7a-d66fb59429bf",
#' "002c67f2-ff52-4246-9d65-a3f69df6789e",
#' "003143c8-bbbf-46b9-a96f-f58530f4bb82")
#'
#' TCGAtranslateID(uuids)
#'
#' ## Translate TCGA Barcode <--> UUIDs
#' barcodes <- c("TCGA-B0-5117-11A-01D-1421-08",
#' "TCGA-B0-5094-11A-01D-1421-08",
#' "TCGA-E9-A295-10A-01D-A16D-09")
#' pt_identifiers <- TCGAbarcode(barcodes)
#' TCGAtranslateID(pt_identifiers, type="file_id")
#'
#' ## Complex Example
#' exampleCodes <- c("TCGA-CK-4948",
#' "TCGA-D1-A17N",
#' "TCGA-4V-A9QX",
#' "TCGA-4V-A9QM")
#' exampleType <- "Gene Expression Quantification"
#' TCGAtranslateID(exampleCodes, dataType = exampleType)
#'
#' @author Marcel Ramos \email{marcel.ramos@roswellpark.org}
#'
#' @export TCGAtranslateID
#' @importFrom httr POST content_type_json http_status
#' @importFrom jsonlite toJSON
TCGAtranslateID <- function(identifier, type="file_id",
                            dataType = NULL) {
    stopifnot(is(identifier, "character"))
    resultingList <- .build_json_request(identifier, type=type,
                                         dataType = dataType)
    urlEndpoint <- strsplit(resultingList[["field"]], "\\.")[[1L]][1]
    if (resultingList[["field"]] == "entity_id")
        urlEndpoint <- "annotations"
    response <- httr::POST(paste0("https://gdc-api.nci.nih.gov/",
                                  urlEndpoint),
                           body = resultingList[["request"]],
                           encode = "json",
                           httr::content_type_json())
    if (httr::http_status(response)$category != "Success")
        stop("Unsuccessful request")
    result <- suppressMessages(httr::content(response,
                            type = "text/tab-separated-values",
                            encoding = "UTF-8"))
    result <- as.data.frame(result, stringsAsFactors = FALSE)
    if (identical(dim(result), c(0L, 0L)))
        stop("No identifier data returned")
    dataCol <- ifelse(urlEndpoint == "cases",
                      "case_id",
                      grep("submitter_id$",
                           names(result),
                           ignore.case = TRUE,
                           value = TRUE))
    result[[dataCol]]
}
