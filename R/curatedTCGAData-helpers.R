#' @import methods
#' @importFrom xml2 read_html
#' @importFrom rvest html_nodes html_attr
#' @importFrom GenomicRanges GRanges GRangesList makeGRangesListFromDataFrame
#' granges
#' @importFrom GenomeInfoDb genome
#' @importFrom MultiAssayExperiment ExperimentList colData colData<- metadata
#' subsetByColumn experiments
#' @importFrom utils data head read.delim
#' @importFrom stats as.formula na.omit setNames
#' @importFrom stringr str_extract
#' @importFrom SummarizedExperiment SummarizedExperiment mcols mcols<- rowData
#'   rowData<-
#' @importFrom GenomicDataCommons files results_all select filter ids cases
#'   expand
#' @importFrom S4Vectors isSingleNumber isSingleInteger isSingleString
#' DataFrame
NULL

## Helpers for downloaded objects

#' @name curatedTCGAData-helpers
#'
#' @title Helper functions for managing MultiAssayExperiment from
#' curatedTCGAData
#'
#' @aliases getSubtypeMap
#'
#' @description
#' Additional helper functions for cleaning and uncovering metadata
#' within a downloaded \code{MultiAssayExperiment} from \code{curatedTCGAData}.
#' The \code{getSubtypeMap} function provides a 2 column \code{data.frame}
#' with in-data variable names and an interpreted names. The
#' \code{getClinicalNames} function provides a vector of variable names that
#' exist in the \code{colData} slot of a downloaded \code{MultiAssayExperiment}
#' object. These variables are obtained from
#' \link[RTCGAToolbox]{getFirehoseData} by default and tend to be present
#' across most cancer codes.
#'
#' @param multiassayexperiment A \linkS4class{MultiAssayExperiment} object
#'
#' @examples
#' \dontrun{
#' library(curatedTCGAData)
#'
#' coad <- curatedTCGAData(diseaseCode = "COAD",
#'     assays = "CNA*", dry.run = FALSE)
#' getSubtypeMap(coad)
#' }
#'
#' @return \itemize{
#'     \item{getSubtypeMap}: A \code{data.frame} with columns representing
#'     actual data variables and explanatory names
#'     \item{getClinicalNames}: A \code{vector} of names that correspond to
#'     a particular disease code.
#' }
#'
#' @export
getSubtypeMap <- function(multiassayexperiment) {

    if (!is(multiassayexperiment, "MultiAssayExperiment"))
        stop("Provide a 'MultiAssayExperiment' object")

    frameMap <- metadata(colData(multiassayexperiment))[["subtypes"]]
    frameMap[] <- lapply(frameMap, as.character)

    if (is.null(frameMap))
        return(message("No subtype data available"))

    subColIdx <- grep("subtype", names(frameMap))

    pats <-
        frameMap[[subColIdx]] %in% c("patient", "SAMPLE", "Complete TCGA ID")

    frameMap[pats, subColIdx] <- "patientID"
    frameMap
}

#' @rdname curatedTCGAData-helpers
#'
#' @param diseaseCode A TCGA cancer code (e.g., "BRCA")
#' @examples
#' getClinicalNames("COAD")
#'
#' @export
getClinicalNames <- function(diseaseCode) {
    stopifnot(S4Vectors::isSingleString(diseaseCode))
    env <- new.env(parent = emptyenv())
    data("clinicalNames", envir = env)
    clinNames <- env[["clinicalNames"]]
    clinNames[[diseaseCode]]
}

.samplesInData <- function(mae) {
    IRanges::CharacterList(lapply(sampleTables(mae), names))
}

.validSamples <- function(samplist, sampleCodes) {
    if (!is.null(sampleCodes)) {
        .checkSamplesInData(samplist, sampleCodes)
        samplist[S4Vectors::`%in%`(samplist, sampleCodes)]
    }
    else
        samplist
}

.checkSamplesInData <- function(samplist, sampleCodes) {
    invalidCodes <- IRanges::LogicalList(lapply(samplist,
        function(acode) !sampleCodes %in% acode))

    if (all(all(invalidCodes)))
        stop("'sampleCodes' not found in assay data, check 'sampleTables()'",
            "\n    and see the 'data(\"sampleTypes\")' table", call. = FALSE)

    if (any(any(invalidCodes))) {
        missingcodes <-
            IRanges::CharacterList(lapply(invalidCodes[any(invalidCodes)],
                function(inv) sampleCodes[inv]))
        warning("Some 'sampleCodes' not found in assays", call. = FALSE)
        missingcodes
    }
}

#' @rdname curatedTCGAData-helpers
#'
#' @param sampleCodes character (default NULL) A string of sample type codes
#' (refer to \code{data(sampleTypes)}; splitAssays section)
#'
#' @section splitAssays:
#'     Separates samples by indicated sample codes into different assays
#'     in a \code{MultiAssayExperiment}. Refer to the \code{sampleTypes}
#'     data object for a list of available codes. This operation generates
#'     \strong{n} times the number of assays based on the number of sample codes
#'     entered. By default, all assays will be split by samples present in
#'     the data.
#' @export
splitAssays <- function(multiassayexperiment, sampleCodes = NULL) {
    if (!is(multiassayexperiment, "MultiAssayExperiment"))
        stop("Provide a 'MultiAssayExperiment' object")

    if (!is.null(sampleCodes)) {
        env <- new.env(parent = emptyenv())
        data("sampleTypes", envir = env, package = "TCGAutils")
        sampleTypes <- env[["sampleTypes"]]
        if (!sampleCodes %in% sampleTypes[["Code"]] ||
            !is.character(sampleCodes))
            stop("Provide valid sample types string")
    }

    sampList <- .samplesInData(multiassayexperiment)
    sampList <- .validSamples(sampList, sampleCodes)

    egroups <- unlist(Map(function(exps, sampcodes, enames) {
        expnames <- setNames(sampcodes, paste0(sampcodes, "_", enames))
        lapply(expnames, function(code) {
            logitype <- TCGAsampleSelect(colnames(exps), code)
            exps[, logitype]
        })
    }, exps = experiments(multiassayexperiment), sampcodes = sampList,
    enames = names(multiassayexperiment), USE.NAMES = FALSE))

    sampmap <- generateMap(egroups, colData(multiassayexperiment),
        idConverter = TCGAbarcode)

    BiocGenerics:::replaceSlots(multiassayexperiment,
        ExperimentList = ExperimentList(egroups),
        sampleMap = sampmap)
}

#' @rdname curatedTCGAData-helpers
#' @param vial (logical default FALSE) whether to display vials in the
#' table output
#'
#' @section sampleTables:
#'     Display all the available samples in each of the assays
#' @export
sampleTables <- function(multiassayexperiment, vial = FALSE) {
    lapply(colnames(multiassayexperiment), function(x) {
        scodes <- TCGAbarcode(x, participant = FALSE, sample = TRUE)
        if (!vial)
            scodes <- substr(scodes, 1L, 2L)
        table(unname(scodes))
   })
}
