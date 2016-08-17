.build_json_request <- function(UUID, size = 100000) {
    options(scipen = 999)
    request <- jsonlite::toJSON(structure(
        list(filters = structure(
            list(op = "in",
                 content = structure(
                     list(
                         field = "files.file_id",
                         value = UUID),
                     .Names = c("field", "value"))),
            .Names = c("op", "content")), format = "TSV",
            fields =
                "cases.samples.portions.analytes.aliquots.submitter_id",
            size = as.character(size)),
        .Names = c("filters", "format", "fields",
                   "size")), pretty = TRUE, auto_unbox = TRUE)
    options(scipen = 0)
    return(request)
}

#' Take UUIDs and return TCGA Barcodes
#'
#' This function will use the GDC API to get metadata from Universally Unique
#' Identifiers (UUID)
#'
#' @param ids A character vector of UUIDs
#' @return A \code{character} vector of TCGA Barcode Identifiers
#'
#' @examples
#' uuids <- c("0004d251-3f70-4395-b175-c94c2f5b1b81",
#' "000d566c-96c7-4f1c-b36e-fa2222467983",
#' "0011a67b-1ba9-4a32-a6b8-7850759a38cf")
#' TCGAtranslateUUID(uuids)
#'
#' @author Marcel Ramos \email{mramos09@gmail.com}
#'
#' @export TCGAtranslateUUID
TCGAtranslateUUID <- function(ids) {
    newRequest <- .build_json_request(ids)
    response <- httr::POST("https://gdc-api.nci.nih.gov/files",
                     body = newRequest,
                     encode = "json", httr::content_type("application/json"))
    result <- httr::content(response,
                            type = "text/tab-separated-values",
                            encoding = "UTF-8")
    result <- as.data.frame(result, stringsAsFactors = FALSE)
    result[[1L]]
}
