## Helper for finding barcode column
## **Takes the first result!**
.findBarcodeCol <- function(DF) {
    cnames <- names(DF)
    containsBC <- vapply(head(DF), function(column) {
        all(startsWith(column, "TCGA"))
    }, logical(1L))
    names(containsBC) <- cnames
    bcIdx <- which(containsBC)
    stopifnot(S4Vectors::isSingleInteger(which(containsBC)))
    names(containsBC)[bcIdx]
}

## Standardize barcode format
.standardBarcodes <- function(sampleBarcode) {
    if (!length(sampleBarcode)) {
        stop("<internal> Barcode must be of positive length")
    }
    sampleBC <- base::sample(sampleBarcode, 10L, replace = TRUE)
    bcodeTest <- grepl("\\.", sampleBC)
    if (all(bcodeTest))
        sampleBarcode <- gsub("\\.", "-", sampleBarcode)
    toupper(sampleBarcode)
}

## Find columns that are all NA
.findNAColumns <- function(dataset) {
    apply(dataset, 2L, function(column) {
        all(is.na(column))
    })
}

#' @importFrom BiocFileCache BiocFileCache bfcquery bfcnew
#' @importFrom httr cache_info HEAD
#' @importFrom rappdirs user_cache_dir
#'
#' @keywords internal
.cacheNeedsUpdate <- function(url) {
    message("updating resource from ", url)
    needsUpdate <- TRUE
    cache <- rappdirs::user_cache_dir(appname = "TCGAutils")

    tryCatch({
        bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
        query <- bfcquery(bfc, url, "rname")

        if (!nrow(query)) {
            file <- bfcnew(bfc, url)
            needsUpdate <- TRUE
        } else {
            file <- query$rpath
            id <- query$rid
            mtime <- file.mtime(query$rpath)
            expires <- httr::cache_info(httr::HEAD(url))$expires
            needsUpdate <- expires < Sys.Date()
        }
    }, error = function(err) {
        stop(
            "could not connect or cache url ", url,
            "\n reason: ", conditionMessage(err)
        )
    })

    setNames(needsUpdate, file)
}
