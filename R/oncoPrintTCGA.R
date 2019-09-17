#' OncoPrint for TCGA Mutation Assays
#'
#' @param multiassayexperiment A MultiAssayExperiment preferably from
#'    `curatedTCGAData``
#'
#' @param mutation character(1) The name of the assay containing mutation data
#'
#' @param variantCol character(1) The name of the metadata column containing
#'     the mutation categories
#'
#' @param brewer_pal character(1) The name of the `RColorBrewer` palette for
#'     the different mutation types
#'
#' @examples
#'
#' acc <- curatedTCGAData("ACC", "Mutation", FALSE)
#'
#' oncoPrintTCGA(acc)
#'
#' @return
oncoPrintTCGA <-
    function(multiassayexperiment, mutation = "*_Mutation-*",
        variantCol = "Variant_Classification", brewer_pal = "Set3", ntop = 25)
{
    stopifnot(
        is.character(mutation) && length(mutation) == 1L,
        is(multiassayexperiment, "MultiAssayExperiment"),
        is.character(variantCol) && length(variantCol) == 1L,
        is.character(brewer_pal) && length(brewer_pal) == 1L
    )

    .checkPkgsAvail(c("TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"))

    mutname <-
        grep(glob2rx(mutation), names(multiassayexperiment), value = TRUE)

    if (length(mutname) > 1)
        stop("Only one mutation assay supported at this time")

    Variants <- mcols(multiassayexperiment[[mutname]])[[variantCol]]
    Variants <- gsub("_", " ", Variants)

    mcols(multiassayexperiment[[mutname]])[[variantCol]] <- Variants
    ## remove single Translation_Site_* mutation
    TranslationSite <-
        mcols(multiassayexperiment[[mutname]])[[variantCol]] ==
            "Translation Start Site"
    mcols(multiassayexperiment[[mutname]])[[variantCol]][TranslationSite] <-
        NA_character_

    variants <- mcols(multiassayexperiment[[mutname]])[[variantCol]]

    mutclass <- table(variants)
    vclass <- setNames(names(mutclass), names(mutclass))

    gn <- sort(.getGN())

    ragex <- multiassayexperiment[[mutname]]
    stopifnot(is(ragex, "RaggedExperiment"))

    gnomel <- length(genome(ragex))
    genome(ragex) <- rep(translateBuild(genome(ragex)[1L]), gnomel)

    simplify_fun <- function(scores, ranges, qranges)
        { sum(S4Vectors::`%in%`(scores, vclass)) }

    res <- qreduceAssay(
        ragex, gn, simplify_fun, "Variant_Classification", background = 0
    )

    rownames(resto) <- names(gn)
    resto <- resto[!is.na(rownames(resto)), ]

    topgenes <- head(sort(rowSums(resto), decreasing = TRUE), ntop)

    qualcolors <- RColorBrewer::brewer.pal(n = length(vclass), brewer_pal)
    colors <- setNames(qualcolors, names(mutclass))

    mutfuns <- lapply(colors, function(couleur) {
        args <- alist(x =, y =, w =, h =)
        args <- as.pairlist(args)
        body <- substitute({
            grid::grid.rect(x, y, w, h, gp = gpar(fill = z, col = NA))
        }, list(z = couleur))
        eval(call("function", args, body))
    })

    background <- function(x, y, w, h)
        grid::grid.rect(x, y, w, h,
            gp = grid::gpar(fill = "#FFFFFF", col = "#FFFFFF"))
    mutfuns2 <- c(background = background, mutfuns)

    resto <- resto[topgenes, ]

    return(
        ComplexHeatmap::oncoPrint(
            resto, alter_fun = mutfuns2, col = colors, show_pct = FALSE
        )
    )
}
