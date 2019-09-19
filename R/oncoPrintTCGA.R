#' OncoPrint for TCGA Mutation Assays
#'
#' @param multiassayexperiment A MultiAssayExperiment preferably from
#'    `curatedTCGAData``
#'
#' @param matchassay character(1) The name of the assay containing mutation
#'     data, this can be a pattern (e.g., "*_Mutation-*", the default)
#'
#' @param variantCol character(1) The name of the metadata column containing
#'     the mutation categories, usually "Variant_Classification" in TCGA
#'
#' @param brewerPal character(1) The name of the `RColorBrewer::brewer.pal`
#'     palette, (default: "Set3")
#'
#' @param ntop integer(1) The number of the top N genes for displaying based
#'     on per-sample mutation frequency
#'
#' @param incl.thresh double(1) The inclusion threshold for empirical mutations,
#'     mutations less frequent than this value will not be included
#'
#' @param rowcol character(1) The name of the column in the metadata to annotate
#'     the rows with either "Hugo_Symbol" (default) or
#'
#' @examples
#'
#' library(curatedTCGAData)
#' acc <- curatedTCGAData("ACC", "Mutation", FALSE)
#'
#' oncoPrintTCGA(acc)
#'
#' @return An oncoPrint plot of mutations
#'
#' @export
oncoPrintTCGA <-
    function(multiassayexperiment, matchassay = "*_Mutation-*",
        variantCol = "Variant_Classification", brewerPal = "Set3", ntop = 25,
        incl.thresh = 0.01, rowcol = "Hugo_Symbol")
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

    rownames(ragex) <- mcols(ragex)[[rowcol]]
    somaticnonsilent <- mcols(ragex)[["Mutation_Status"]] == "Somatic" &
        mcols(ragex)[[variantCol]] != "Silent"
    ragex <- ragex[somaticnonsilent, ]

    Variants <- mcols(ragex)[[variantCol]]
    Variants <- gsub("_", " ", Variants)
    mcols(ragex)[[variantCol]] <- Variants

    types <- table(Variants)
    tottypes <- sum(types)
    incl <- (types/tottypes) > incl.thresh
    types <- types[incl]
    validvariants <- setNames(names(types), names(types))

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
    gn <- gn[!is.na(names(gn))]

    simplify_fun <- function(scores, ranges, qranges)
        { any(scores != "Silent") }

    res <- qreduceAssay(
        ragex, gn, simplify_fun, "Variant_Classification", background = FALSE
    )
    rownames(res) <- names(gn)

    topgenes <- head(sort(rowSums(res), decreasing = TRUE), ntop)
    gn2 <- gn[match(names(topgenes), names(gn))]

    qualcolors <-
        RColorBrewer::brewer.pal(n = length(validvariants), brewer_pal)
    colors <- setNames(qualcolors, validvariants)

    colfuns <- lapply(colors, function(couleur) {
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
    mutfuns <- c(background = background, colfuns)

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
    })

    return(
        ComplexHeatmap::oncoPrint(
            list_mats, alter_fun = mutfuns, col = colors, show_pct = FALSE
        )
    )
}

