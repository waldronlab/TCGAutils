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

.get_hsa_url <- function(gen) {
    ## handle cases where genome is specified as "GRCh37.p5"
    if (startsWith(gen, "GRCh37"))
        gen <- "GRCh37"
    else if (startsWith(gen, "GRCh38"))
        gen <- "GRCh38"

    release <- switch(
        gen,
        GRCh37 = ,
        hg19 = "20",
        GRCh38 = ,
        hg38 = "22",
        stop("Unsupported genome: ", gen)
    )
    glue::glue(
        "https://www.mirbase.org/download_version_genome_files/",
        "{release}/hsa.gff3"
    )
}

.get_hsa_genome <- function(file) {
    gff_lines <- readLines(file, n = 50)
    genome_line <- grepv("genome-build-id", gff_lines)
    gnm <- strsplit(genome_line, ":\\s+")[[1L]] |>
        utils::tail(n = 1L) |>
        trimws()
    if (!length(gnm))
        NA_character_
    else
        gnm
}

.get_hsa_gff3 <- function(gen, redownload) {
    url <- .get_hsa_url(gen)
    gff_local <- .cache_url_file(url, redownload = redownload)
    res <- Bioc.gff::import(gff_local)
    # res <- res[mcols(res)[["type"]] %in% c("miRNA", "microRNA", "tRNA"), ]
    gnm <- .get_hsa_genome(gff_local)
    genome(res) <- gnm
    if (identical(gen, "hg18")) {
        ## perform liftOver operation from GRCh38 to hg18
        checkInstalled(c("AnnotationHub", "rtracklayer"))
        chain <- AnnotationHub::AnnotationHub()[["AH14221"]]
        ranges18 <- rtracklayer::liftOver(res, chain)
        res <- res[as.logical(lengths(ranges18))]
    }
    res
}

.cache_url_file <- function(url, redownload) {
    checkInstalled("BiocFileCache")
    bfc <- BiocFileCache::BiocFileCache()
    bquery <- BiocFileCache::bfcquery(bfc, url, "rname", exact = TRUE)
    ## only re-download manually b/c bfcneedsupdate always returns TRUE
    if (identical(nrow(bquery), 1L) && redownload)
        tryCatch({
            BiocFileCache::bfcdownload(
                x = bfc, rid = bquery[["rid"]], ask = FALSE
            )
        }, error = function(e) {
            msg <- conditionMessage(e)
            if (grepl("download failed", msg, TRUE))
                warning(msg, call. = FALSE)
            else
                stop(msg, call. = FALSE)
            invisible()
        })

    BiocFileCache::bfcrpath(
        bfc, rnames = url, exact = TRUE, download = TRUE, rtype = "web"
    )
}
