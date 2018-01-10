#' @import methods
#' @importFrom xml2 read_html
#' @importFrom rvest html_nodes html_attr
#' @importFrom GenomicRanges GRanges GRangesList makeGRangesListFromDataFrame
#' @importFrom MultiAssayExperiment ExperimentList colData colData<- metadata
#' @importFrom utils data head read.delim
#' @importFrom stats as.formula na.omit
#' @importFrom IRanges CharacterList
#' @importFrom SummarizedExperiment SummarizedExperiment mcols mcols<- rowData
#' rowData<-
#' @importFrom GenomicDataCommons files results_all select filter ids cases
#' @importFrom S4Vectors isSingleNumber isSingleInteger isSingleString
#' DataFrame
#' @importFrom RTCGAToolbox getData
NULL

#' TCGAutils: Helper fucntions for working with TCGA and MultiAssayExperiment
#' data
#'
#' TCGAutils is a toolbox to work with TCGA specific datasets. It allows the
#' user to manipulate and translate TCGA barcodes, conveniently convert a list
#' of data files to \linkS4class{GRangesList}. Take datasets from GISTIC and
#' return a \linkS4class{SummarizedExperiment} class object. The package also
#' provides functions for working with data from the \code{curatedTCGAData}
#' experiment data package. It provides convenience functions for extracting
#' subtype metadata data and adding clinical data to existing
#' \linkS4class{MultiAssayExperiment} objects.
"_PACKAGE"
