#' Translate TCGA Barcode to Universally Unique Identifiers and vice versa
#' 
#' This function allows the user to enter a character vector of identifiers
#' to translate from barcodes to UUIDs and vice versa. 
#' 
#' @param identifier A \code{character} vector of either TCGA
#' barcodes or UUID identifiers
#' 
#' @author Marcel Ramos \email{mramos09@@gmail.com}
#' 
#' @return A \code{data.frame} of original and translated identifiers
#' @export translateIDS
#' @importFrom httr POST content content_type
translateIDS <- function(identifier) {
    identifier <- unique(identifier)
    if (length(identifier) > 500) {
        stop("enter 500 at a time for now")
    }
    keyword <- ifelse(grepl("TCGA", identifier[1], ignore.case = TRUE),
               "barcode",
               "uuid")
    queryURL <-
        paste0("https://tcga-data.nci.nih.gov/uuid/uuidws/mapping/json/",
                  keyword, "/batch")
    idQuery <- paste(identifier, collapse = ",")
    id_table <- httr::POST(queryURL, body = idQuery, encode = "json",
                         httr::content_type("text/plain"))
    id_table <- do.call(rbind,
                            lapply(
                                httr::content(id_table,
                                              "parsed")$uuidMapping, unlist))
    if (dim(id_table)[2] == 1L) {id_table <- t(id_table)}
    id_table <- as.data.frame(id_table, stringsAsFactors = FALSE)
    return(id_table)
}
