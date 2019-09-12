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
        variantCol = "Variant_Classification", brewer_pal = "Set3")
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

    simplify_funs <- lapply(vclass,
        function(variant) {
            args <- alist(scores =, ranges =, qranges =)
            args <- as.pairlist(args)
            body <- substitute({
                any(S4Vectors::`%in%`(scores, z))
            }, list(z = variant))
            eval(call("function", args, body))
        }
    )

    ragex <- multiassayexperiment[[mutname]]
    stopifnot(is(ragex, "RaggedExperiment"))

    gnomel <- length(genome(ragex))
    genome(ragex) <- rep(translateBuild(genome(ragex)[1L]), gnomel)

    list_mats <- lapply(simplify_funs, function(variant_fun) {
        res <-
            RaggedExperiment::qreduceAssay(ragex, gn, variant_fun, variantCol)
        rownames(res) <- names(gn)
        res[is.na(res)] <- 0
        res[!is.na(rownames(res)), ]
    })

    nomutgenes <- Reduce(intersect,
        lapply(list_mats, function(x) rownames(x[rowSums(x) == 0, ]))
    )

    updated_mats <- lapply(list_mats, function(genemat) {
        genemat[!rownames(genemat) %in% nomutgenes, ]
    })

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

    tops <- ceiling(nrow(updated_mats[[1]]) * 0.0005)
    mgenes <- Reduce(union, lapply(updated_mats, function(x)
        names(sort(rowSums(x), decreasing = TRUE)[seq_len(tops)]))
    )
    ug <- length(unique(mgenes))
    upmats2 <- lapply(updated_mats, function(mat) mat[mgenes, ])

    return(
        ComplexHeatmap::oncoPrint(
            upmats2, alter_fun = mutfuns2, col = colors, show_pct = FALSE
        )
    )
}
