#' TCGA Cancer Disease Codes Table
#'
#' A dataset for obtaining the cancer codes in TCGA for about 13 different
#' types of cancers.
#'
#' @format A data frame with 37 rows and 2 variables:
#' \describe{
#'      \item{Study.Abbreviation}{Disease Code used in TCGA}
#'      \item{Available}{Cancer datasets available via curatedTCGAData}
#'      \item{SubtypeData}{Subtype curation data available via curatedTCGAData}
#'      \item{Study.Name}{The full length study name (i.e., type of cancer)}
#' }
#' @return The TCGA `diseaseCodes` table
#' @source \url{https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations}
"diseaseCodes"

.parseDiseaseCodes <- function(from, to = "./data/diseaseCodes.rda") {
    htcc <- xml2::read_html(from)
    diseaseCodes <- rvest::html_table(htcc, fill = TRUE)[[2L]]
    names(diseaseCodes) <- make.names(colnames(diseaseCodes))

    excludedCodes <- c("COADREAD", "GBMLGG", "KIPAN", "STES", "FPPP", "CNTL",
                       "LCML", "MISC")
    available <- !diseaseCodes[["Study.Abbreviation"]] %in% excludedCodes
    diseaseCodes[["Available"]] <- factor(available,  levels = c("TRUE", "FALSE"),
        labels = c("Yes", "No"))

    subtypeCodes <- c("ACC", "BLCA", "BRCA", "COAD", "GBM", "HNSC", "KICH",
        "KIRC", "KIRP", "LAML", "LGG", "LUAD", "LUSC", "OV", "PRAD", "SKCM",
        "STAD", "THCA", "UCEC")
    diseaseCodes[["SubtypeData"]] <- factor(
        diseaseCodes[["Study.Abbreviation"]] %in% subtypeCodes,
        levels = c("TRUE", "FALSE"), labels = c("Yes", "No"))

    diseaseCodes <- diseaseCodes[order(diseaseCodes[["Study.Abbreviation"]]), ]
    ## Rearrange column order
    diseaseCodes <- diseaseCodes[,
        c("Study.Abbreviation", "Available", "SubtypeData", "Study.Name")]
    rownames(diseaseCodes) <- NULL

    ## Coerce to standard data.frame (no tibble required)
    diseaseCodes <- as(diseaseCodes, "data.frame")

    ## For easy subsetting use:
    ## diseaseCodes[["Study.Abbreviation"]][diseaseCodes$Available == "Yes"]

    ## Save dataset for exported use
    save(diseaseCodes, file = to, compress = "bzip2")
    TRUE
}

.get_cache <- function() {
    cache <- rappdirs::user_cache_dir("TCGAutils")
    BiocFileCache::BiocFileCache(cache)
}

update_data_file <- function(fileURL, verbose = FALSE , resource) {
    bfc <- .get_cache()
    rid <- BiocFileCache::bfcquery(bfc, fileURL, "rname")$rid
    if (!length(rid)) {
        if (verbose)
            message( "Downloading ", resource, " file" )
        rid <- names(BiocFileCache::bfcadd(bfc, fileURL, download = FALSE))
    }
    if (!isFALSE(BiocFileCache::bfcneedsupdate(bfc, rid)))
        BiocFileCache::bfcdownload(bfc, rid, ask = FALSE, FUN = .parseDiseaseCodes)
    if (verbose)
        message(resource, " update complete")
}

# url1 <- "https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations"
# update_data_file(url1, verbose = FALSE, resource = "diseaseCodes")

#' Barcode Sample Type Table
#'
#' A dataset that contains the mappings for sample codes in the TCGA
#' barcodes.
#' @format A data frame with 19 rows and 3 variables:
#' \describe{
#'      \item{Code}{Two digit code number found in the barcode}
#'      \item{Definition}{Long name for the sample type}
#'      \item{Short.Letter.Code}{Letter code for the sample type}
#' }
#' @return The TCGA `sampleTypes` table
#' @source \url{https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes}
"sampleTypes"

.parseSampleTypes <- function(from, to = "./data/sampleTypes.rda") {
    stcc <- read_html(from)

    sampleTypes <- html_table(stcc, fill = TRUE)[[2L]]
    names(sampleTypes) <- make.names(colnames(sampleTypes))

    ## Coerce to standard data.frame (no tibble required)
    sampleTypes <- as(sampleTypes, "data.frame")

    ## Save dataset for exported use
    save(sampleTypes, file = to, compress = "bzip2")
    TRUE
}

# url2 <- "https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes"
# update_data_file(url2, verbose = FALSE, resource = "sampleTypes")

#' Clinical dataset names in TCGA
#'
#' A dataset of names for each of the TCGA cancer codes available.
#' These names were obtained by the clinical datasets from
#' \link[RTCGAToolbox]{getFirehoseData}. They serve to subset the current
#' datasets provided by \code{curatedTCGAData}.
#'
#' @format A \linkS4class{CharacterList} of names for 33 cancer codes
#'
#' @return The clinical dataset column names in TCGA as provided by the
#' \code{RTCGAToolbox}
"clinicalNames"
