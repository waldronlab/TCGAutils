.getLinks <- function(keyWord1, keyWord2, datasetLink = NULL, doc) {
    # Function from RTCGAToolbox
    keyWord <- keyWord1
    keyWord <- paste0("//a[contains(@href, '",keyWord,"')]")
    plinks <- rvest::html_nodes(doc, xpath = keyWord)
    plinks <- rvest::html_attr(plinks, "href")
    if (is.null(datasetLink))
        plinks <- plinks[grepl(keyWord2,plinks)]
    else
        plinks <- plinks[grepl(paste0("*.",datasetLink,keyWord2),plinks)]
    return(plinks)
}

#' Find the file names used in RTCGAToolbox
#'
#' Part of this function is from the RTCGAToolbox. It aims to extract the file
#' name used inside of the \link{getFirehoseData} function. The arguments of the
#' function parallel those in the \link{getFirehoseData} function. It is only
#' available for select data types.
#'
#' @param disease The TCGA cancer disease code, e.g., "COAD"
#' @param runDate The single \code{string} used in the \code{getFirehoseData}
#' function (default "20160128")
#' @param CNASNP A \code{logical} (default = FALSE) vector indicating whether
#' to get the file name from this data type
#' @param CNVSNP A \code{logical} (default = FALSE) vector indicating whether
#' to get the file name from this data type
#' @param CNASeq A \code{logical} (default = FALSE) vector indicating whether
#' to get the file name from this data type
#' @param CNACGH A \code{logical} (default = FALSE) vector indicating whether
#' to get the file name from this data type
#'
#' @return A \code{character} vector of length one indicating the file name
#'
#' @importFrom xml2 read_html
#' @importFrom rvest html_nodes html_attr
#'
#' @examples
#'
#' getFileNames("COAD", CNVSNP = TRUE)
#'
#' @export getFileNames
getFileNames <- function(disease, runDate = "20160128", CNASNP = FALSE,
    CNVSNP = FALSE, CNASeq = FALSE, CNACGH = FALSE) {
    if (!any(c(CNASNP, CNVSNP, CNASeq, CNACGH)))
        stop("Set a data type to TRUE")

    fh_url <- "http://gdac.broadinstitute.org/runs/stddata__"
    fh_url <- paste0(fh_url, substr(runDate,1,4), "_",
                     substr(runDate,5,6), "_",
                     substr(runDate,7,8), "/data/")
    fh_url <- paste0(fh_url, disease, "/", runDate, "/")
    doc <- xml2::read_html(fh_url)

    plinks <- vector(mode = "character", length = 1L)
    names(plinks) <- disease

    if (CNASNP)
        plinks <- .getLinks(
            "Level_3__segmented_scna_hg19__seg.Level_3",
            paste0("[.]Merge_snp__.*.__Level_3__segmented",
                   "_scna_hg19__seg.Level_3.*.tar[.]gz$"),
            disease, doc)
    if (CNVSNP)
        plinks <- .getLinks(
            "Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3",
            paste0("[.]Merge_snp__.*.__Level_3__segmented_scna_",
                   "minus_germline_cnv_hg19__seg.Level_3.*.tar[.]gz$"),
            disease, doc)
    if (CNASeq)
        plinks <- .getLinks("__Level_3__segmentation__seg.Level_3",
                            paste0("[.]Merge_cna__.*.dnaseq.*.__Level_3__",
                                   "segmentation__seg.Level_3.*.tar[.]gz$"),
                            disease, doc)
    if (CNACGH)
        plinks <- .getLinks("__Level_3__segmentation__seg.Level_3",
                            paste0("[.]Merge_cna__.*.cgh.*.__Level_3__",
                                   "segmentation__seg.Level_3.*.tar[.]gz$"),
                            disease, doc)

    if(S4Vectors::isSingleString(plinks)) {
        return(unname(plinks))
    } else {
        return(NULL)
    }
}
