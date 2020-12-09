## Extract sample types table from TCGA website
.parseSampleTypes <- function(from, to) {
    stcc <- xml2::read_html(from)

    sampleTypes <- rvest::html_table(stcc, fill = TRUE)[[2L]]

    ## convert code column to character
    codeCol <- sampleTypes[["Code"]]
    singleDigit <- codeCol < 10L
    sampleTypes[["Code"]][singleDigit] <-
        paste0("0", sampleTypes[["Code"]][singleDigit])

    names(sampleTypes) <- make.names(colnames(sampleTypes))

    ## Coerce to standard data.frame (no tibble required)
    sampleTypes <- as(sampleTypes, "data.frame")

    ## Save dataset for exported use
    save(sampleTypes, file = to, compress = "bzip2")
    TRUE
}

url2 <-
"https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes"
## update_data_file in data-raw/diseaseCodes.R
update_data_file(url2, verbose = FALSE, resource = "sampleTypes",
    FUN = .parseSampleTypes)
