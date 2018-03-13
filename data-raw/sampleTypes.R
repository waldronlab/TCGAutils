## Extract sample types table from TCGA website
updateSampleTypes <- function() {
if (!requireNamespace("rvest") || !requireNamespace("devtools"))
    stop ("Please download 'rvest' to update web resources")

stTableLink <- "https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes"

needsUpdate <- .needsCacheUpdate(stTableLink)

if (needsUpdate) {

stcc <- read_html(stTableLink)

sampleTypes <- html_table(stcc, fill = TRUE)[[2L]]
names(sampleTypes) <- make.names(colnames(sampleTypes))

## Coerce to standard data.frame (no tibble required)
sampleTypes <- as(sampleTypes, "data.frame")

## Save dataset for exported use
devtools::use_data(sampleTypes, internal = FALSE, overwrite = TRUE)
}

}
