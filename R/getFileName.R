.getLinks <- function(keyWord1, keyWord2, datasetLink = NULL, doc) {
    # Function from RTCGAToolbox
    keyWord <- keyWord1
    keyWord <- paste0("//a[contains(@href, '",keyWord,"')]")
    plinks <- rvest::html_nodes(doc, xpath = keyWord)
    plinks <- rvest::html_attr(plinks, "href")

    if (is.null(datasetLink))
        plinks[grepl(keyWord2,plinks)]
    else
        plinks[grepl(paste0("*.",datasetLink,keyWord2),plinks)]
}

#' Find the file names used in RTCGAToolbox
#'
#' Part of this function is from the RTCGAToolbox. It aims to extract the file
#' name used inside of the \link[RTCGAToolbox]{getFirehoseData} function.
#' The arguments of the function parallel those in the
#' \link[RTCGAToolbox]{getFirehoseData} function. It is only available for
#' select data types.
#'
#' @param disease The TCGA cancer disease code, e.g., "COAD"
#' @param runDate The single \code{string} used in the \code{getFirehoseData}
#' function (default "20160128")
#' @param dataType A single character vector (default "CNASNP") indicating the
#' data type for which to get the source file name
#'
#' @return A single \code{character} file name
#'
#' @examples
#'
#' getFileName("COAD", dataType = "CNASNP")
#'
#' @export getFileName
getFileName <- function(disease, runDate = "20160128",
    dataType = c("CNASNP", "CNVSNP", "CNAseq", "CNACGH", "Mutation")) {

    dataType <- match.arg(dataType,
        c("CNASNP", "CNVSNP", "CNAseq", "CNACGH", "Mutation"))

    fh_url <- "https://gdac.broadinstitute.org/runs/stddata__"
    fh_url <- paste0(fh_url, substr(runDate,1,4), "_",
        substr(runDate,5,6), "_", substr(runDate,7,8), "/data/")
    fh_url <- paste0(fh_url, disease, "/", runDate, "/")
    doc <- xml2::read_html(fh_url)

    switch(dataType,
        CNASNP = .getLinks(
            "Level_3__segmented_scna_hg19__seg.Level_3",
            paste0("[.]Merge_snp__.*.__Level_3__segmented",
                "_scna_hg19__seg.Level_3.*.tar[.]gz$"),
            disease, doc),
        CNVSNP = .getLinks(
            "Level_3__segmented_scna_minus_germline_cnv_hg19__seg.Level_3",
            paste0("[.]Merge_snp__.*.__Level_3__segmented_scna_",
                "minus_germline_cnv_hg19__seg.Level_3.*.tar[.]gz$"),
            disease, doc),
        CNASeq = .getLinks("__Level_3__segmentation__seg.Level_3",
            paste0("[.]Merge_cna__.*.dnaseq.*.__Level_3__",
                "segmentation__seg.Level_3.*.tar[.]gz$"),
            disease, doc),
        CNACGH = .getLinks("__Level_3__segmentation__seg.Level_3",
            paste0("[.]Merge_cna__.*.cgh.*.__Level_3__",
                "segmentation__seg.Level_3.*.tar[.]gz$"),
            disease, doc),
        Mutation = .getLinks("Mutation_Packager_Calls",
            "[.]Mutation_Packager_Calls[.]Level_3[.].*.tar[.]gz$",
            disease, doc)
    )
}
