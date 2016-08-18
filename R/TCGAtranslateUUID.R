.build_json_request <- function(UUID, size = 100) {
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
            fields = paste("file_id", "file_name",
                "cases.samples.portions.analytes.aliquots.submitter_id",
                sep = ","),
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
#' uuids <- c("0001801b-54b0-4551-8d7a-d66fb59429bf",
#' "002c67f2-ff52-4246-9d65-a3f69df6789e",
#' "003143c8-bbbf-46b9-a96f-f58530f4bb82")
#' TCGAtranslateUUID(uuids)
#'
#' @author Marcel Ramos \email{mramos09@gmail.com}
#'
#' @export TCGAtranslateUUID
TCGAtranslateUUID <- function(ids) {
    stopifnot(is(ids, "character"))
    newRequest <- .build_json_request(ids)
    response <- httr::POST("https://gdc-api.nci.nih.gov/files",
                           body = newRequest,
                           encode = "json",
                           httr::content_type_json())
    result <- suppressMessages(httr::content(response,
                            type = "text/tab-separated-values",
                            encoding = "UTF-8"))
    result <- as.data.frame(result, stringsAsFactors = FALSE)
    result[[2L]]
}
