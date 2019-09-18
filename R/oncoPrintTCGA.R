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
#' library(curatedTCGAData)
#' acc <- curatedTCGAData("ACC", "Mutation", FALSE)
#'
#' oncoPrintTCGA(acc)
#'
#' @return An oncoPrint plot of mutations
#'
#' @export
oncoPrintTCGA <-
    function(multiassayexperiment, mutation = "*_Mutation-*",
        variantCol = "Variant_Classification", brewer_pal = "Set3", ntop = 25,
        incl.thresh = 0.01)
{
    stopifnot(
        is.character(mutation) && length(mutation) == 1L,
        is(multiassayexperiment, "MultiAssayExperiment"),
        is.character(variantCol) && length(variantCol) == 1L,
        is.character(brewer_pal) && length(brewer_pal) == 1L
    )

    .checkPkgsAvail("org.Hs.eg.db")

    mutname <-
        grep(glob2rx(mutation), names(multiassayexperiment), value = TRUE)

    if (length(mutname) > 1)
        stop("Only one mutation assay supported at this time")

    ragex <- multiassayexperiment[[mutname]]
    stopifnot(is(ragex, "RaggedExperiment"))

    Variants <- mcols(ragex)[[variantCol]]
    Variants <- gsub("_", " ", Variants)
    mcols(ragex)[[variantCol]] <- Variants

    types <- table(Variants)
    tottypes <- sum(types)
    incl <- (types/tottypes) > incl.thresh & names(types) != "Silent"
    types <- types[incl]
    validvariants <- names(types)

    ragex <- ragex[mcols(ragex)[[variantCol]] %in% validvariants, ]
    rowRanges(ragex) <- unstrand(rowRanges(ragex))

    genomeannot <- unique(genome(ragex))

    if (length(genomeannot) > 1)
        stop("'genome' annotation is not consistent")

    if (!startsWith(genomeannot, "hg")) {
        genomeannot <- translateBuild(genomeannot)
        genome(ragex) <- rep(genomeannot, length(genome(ragex)))
    }

    .checkPkgsAvail(paste0("TxDb.Hsapiens.UCSC.", genomeannot, ".knownGene"))

    gn <- sort(.getGN(genomeannot))
    gn <- unstrand(gn)

    simplify_fun <- function(scores, ranges, qranges)
        { any(scores != "Silent") }

    res <- qreduceAssay(
        ragex, gn, simplify_fun, "Variant_Classification", background = FALSE
    )
    rownames(res) <- names(gn)
    res <- res[!is.na(rownames(res)), ]

    topgenes <- head(sort(rowSums(res), decreasing = TRUE), ntop)
    gn2 <- gn[names(gn) %in% names(topgenes)]

    qualcolors <-
        RColorBrewer::brewer.pal(n = length(validvariants), brewer_pal)
    colors <- setNames(qualcolors, validvariants)

    mutfuns <- lapply(colors, function(couleur) {
        args <- alist(x =, y =, w =, h =)
        args <- as.pairlist(args)
        body <- substitute({
            grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = z, col = NA))
        }, list(z = couleur))
        eval(call("function", args, body))
    })

    background <- function(x, y, w, h)
        grid::grid.rect(x, y, w, h,
            gp = grid::gpar(fill = "#FFFFFF", col = "#FFFFFF"))
    mutfuns2 <- c(background = background, mutfuns)

    validvariants <- setNames(validvariants, validvariants)
    simplify_funs <- lapply(validvariants,
        function(variant) {
            args <- alist(scores =, ranges =, qranges =)
            args <- as.pairlist(args)
            body <- substitute({
                as.numeric(any(S4Vectors::`%in%`(scores, z)))
            }, list(z = variant))
            eval(call("function", args, body))
        }
    )

    list_mats <- lapply(simplify_funs, function(variant_fun) {
        res <- qreduceAssay(ragex, gn2, variant_fun,
            "Variant_Classification", background = 0)
        rownames(res) <- names(gn2)
        res <- res[!is.na(rownames(res)), ]
    })

    return(
        ComplexHeatmap::oncoPrint(
            list_mats, alter_fun = mutfuns2, col = colors, show_pct = FALSE
        )
    )
}

