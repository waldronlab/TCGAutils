## Update internally stored datasets when URL is updated
## write to data folder

.parseSampleTypes <- function(url) {
    stcc <- read_html(url)

    sampleTypes <- html_table(stcc, fill = TRUE)[[2L]]
    names(sampleTypes) <- make.names(colnames(sampleTypes))

    ## Coerce to standard data.frame (no tibble required)
    sampleTypes <- as(sampleTypes, "data.frame")

    ## Save dataset for exported use
    devtools::use_data(sampleTypes, internal = FALSE, overwrite = TRUE)
}

.updateSampleTypes <-
    function(url = "https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes")
{
    needsUpdate <- .cacheNeedsUpdate(url)
    if (needsUpdate) {
        if (!requireNamespace("rvest") || !requireNamespace("devtools"))
            stop ("Please download 'rvest' and/or 'devtools' to",
                "\n  update web resources")
        .parseSampleTypes(url)
    }
}

.parseDiseaseCodes <- function(url) {
    htcc <- read_html(url)
    diseaseCodes <- html_table(htcc, fill = TRUE)[[2L]]
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
    devtools::use_data(diseaseCodes, internal = FALSE, overwrite = TRUE)
}

.updateDiseaseCodes <-
    function(url = "https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations")
{
    needsUpdate <- .cacheNeedsUpdate(url)
    if (needsUpdate) {
        if (!requireNamespace("rvest") || !requireNamespace("devtools"))
            stop ("Please download 'rvest' and/or 'devtools' to",
                "\n  update web resources")
        .parseDiseaseCodes(url)
    }
}
