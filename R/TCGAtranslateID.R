.build_json_request <- function(identifier, size = 1000000) {
    options(scipen = 999)
    field <- ifelse(grepl("^TCGA", identifier[1], ignore.case = TRUE),
               "cases.submitter_id",
               "files.file_id")
    request <- jsonlite::toJSON(structure(
        list(filters = structure(
            list(op = "in",
                 content = structure(
                     list(
                         field = field,
                         value = identifier),
                     .Names = c("field", "value"))),
            .Names = c("op", "content")), format = "TSV",
            fields = paste("file_id", "file_name",
                "cases.samples.portions.analytes.aliquots.submitter_id",
                sep = ","),
            size = as.character(size)),
        .Names = c("filters", "format", "fields",
                   "size")), pretty = TRUE, auto_unbox = TRUE)
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
#'
#' TCGAtranslateID(barcodes)
#'
#' @author Marcel Ramos \email{mramos09@gmail.com}
#'
#' @export TCGAtranslateID
#' @importFrom httr POST content_type_json
#' @importFrom jsonlite toJSON
TCGAtranslateID <- function(identifier) {
    stopifnot(is(identifier, "character"))
    resultingList <- .build_json_request(identifier)
    urlEndpoint <- strsplit(resultingList[["field"]], "\\.")[[1L]][1]
    response <- httr::POST(paste0("https://gdc-api.nci.nih.gov/",
                                  urlEndpoint),
                           body = resultingList[["request"]],
                           encode = "json",
                           httr::content_type_json())
    if (http_status(response)$category != "Success")
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
