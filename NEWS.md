## Changes in version 1.20.0

### New features

* `makeSummarizedExperimentFromGISTIC` and `splitAssays` are now defunct.

## Changes in version 1.18.0

### Minor changes and bug fixes

* Use https instead of http in `getFileName` helper.
* Warn when column names in assays are not mappable and subsequently dropped in
`generateMap`.
* Updated `qreduceTCGA` documentation for clarity.

## Changes in version 1.16.0

### New features

* The `UUIDhistory` function allows users to map old UUIDs to new UUIDs
according to the latest data release for UUIDs that were affected and no longer
query-able.
* The `slides` argument has been added to the `filenameToBarcode` function for
translating slide file names into barcodes. Currently, the API returns all
barcodes of the associated case ID.
* Add sections in the vignette regarding GDC Data Updates and UUID history
lookup

### Minor changes and bug fixes

* Update examples in package to new GDC Data Release, see vignette.
* Use `AnnotationHub` to download chain file in main vignette.
* Slide file names now resolve to a single TCGA barcode in `filenameToBarcode`
(Thanks @hermidalc)
* Improved error messages and documentation for `makeGRangesListFromExonFiles`

## Changes in version 1.14.0

### Minor changes and bug fixes

* `UUIDtoBarcode` with the `from_type = "file_id"` argument now returns the IDs
in the proper order when more than one `UUID` is input.
* Update `makeGRangesListFromCopyNumber` examples with new names from API e.g.,
'associated_entities.entity_submitter_id'

## Changes in version 1.12.0

### New features

* `makeSummarizedExperimentFromGISTIC` has been moved to `RTCGAToolbox`.
* `splitAssays` now deprecated for `TCGAsplitAssays` to avoid conflict with
`MultiAssayExperiment::splitAssays`

### Minor changes and bug fixes

* Properly identifies genome annotation (`hg*`) in `oncoPrintTCGA`
* `qreduceTCGA` now works with updates to `seqlevelsStyle` where genome
annotation include patch versions when available

## Changes in version 1.10.0

### New features

* `correctBuild` attempts to provide the official name of a particular human
genome build to agree with changes in `GenomeInfoDb`
* `isCorrect` checks that the build name matches the official name

### Minor changes and bug fixes

* Documentation improvements to `simplifyTCGA`
* Improvements to `findGRangesCols` to locate ranged columns in a `DataFrame`
* Fixed a bug in `UUIDtoBarcode` where only the first record was returned
(#26, @DarioS)
* Fixed a bug in `filenameToBarcode` when multiple inputs were used (#22,
@DarioS)

## Changes in version 1.8.0

### New features

* `README.md` now includes a cheat sheet for reference
* `mergeColData` and `oncoPrintTCGA` sections updated/included in the vignette

### Minor changes and bug fixes

* `translateBuild` more robust to consistent inputs
* `translateBuild` returns vector output instead of single string as before
* `makeSummarizedExperimentFromGISTIC` now has a more open interface with
`...` input to `RTCGAToolbox::getGISTICPeaks`
* `oncoPrintTCGA` now uses `seqlevels` from input throughout

## Changes in version 1.6.0

### New features

* `oncoPrintTCGA`: Create an `oncoPrint` visualization for mutation data
* Support `aliquot_ids` as input to `UUIDtoBarcode` function
* Additional sections in the vignette: `CpGtoRanges`, `UUIDtoBarcode` for
`aliquot_ids`
* `TCGAprimaryTumors` allows users to select all primary tumors for a given
`curatedTCGAData` `MultiAssayExperiment` object (suggested by @vjcitn)

### Minor changes and bug fixes

* Now merging clinical data using both rows and columns in `mergeColData`
* Added informative error when query results are empty in `UUIDtoBarcode`
* Updates to `makeGRangesListFromExonFiles` to use `S4Vectors::splitAsList`
(@hpages)

## Changes in version 1.4.0

### New features

* `trimColData` added to remove any extra columns from the `colData` slot
(thanks to @vjcitn)
* `CpGtoRanges` translates CpG islands to genomic positions using an annotation
package and `minfi`
* Overhaul of the barcode translation services allows accurate translation
of identifiers
* `splitAssays` now separates all assays by sample codes contained therein
by default, previous behavior had default values
* Documentation for `simplifyTCGA` was modified to include similar operations,
such as, `symbolsToRanges`, `mirToRanges`, `CpGtoRanges`, etc.
* Vignette includes comprehensive examples of new functionality

### Minor changes and bug fixes

* `getFileNames` renamed to `getFileName`
* `TCGAsampleSelect` now allows multiple sample type inputs as the
`sampleCodes` argument
* `getSubtypeMap` updates column names to accurately represent patient
identifiers
* More robust checks were added to `splitAssays` to ensure valid sample codes
in the input and provided as arguments
* `makeGRangesListFromExonFiles` is optimized to use `dplyr` when available
and fast operations from `IRanges`
* Various enhancements to `*toRanges` functions, including re-using underlying
common helper function
* The internal `weightedmean` function in `qreduceTCGA` has been updated for
correctness
* The `keep` arugment in `qreduceTCGA` and related functions was changed
to `keep.assay`

## Changes in version 1.2.0

### New features

* `imputeAssay` added to impute data for MultiAssayExperiment assays
* `UUIDtoUUID` translation available to translate from file to case IDs
* A suite of functions is available to enhance existing MultiAssayExperiment
datasets: `qreduceTCGA`, `mirToRanges`, `symbolsToRanges`. Thanks to @lwaldron

### Minor changes and bug fixes

* Various changes to examples for compatibility with RaggedExperiment
* Bug fix to internal functions for finding GRanges columns

## Changes in version 1.1.5

* `uniformBuilds` cleans up a vector of inconsistently labelled builds by
changing the build with the lowest frequency

## Changes in version 1.1.4

### New features

* The `UUIDtoUUID` function can translate from case to file UUIDs and vice
versa
* `imputeAssay` allows missing data imputation using KNN for
`MultiAssayExperiment` assays

## Changes in version 1.1.1

### New features

* exported the internal helper, `filenameToBarcode`. See examples

## Changes in version 0.99.68

### Minor changes and bug fixes

* Minor changes in response to review, avoid switching from logical to numeric
index, updated vignette introduction
* Fix examples to updated `GenomicDataCommons` interface
* Move `RTCGAToolbox` to `Suggests` field in DESCRIPTION
* Removed `BiocFileCache` from `Imports` field

## Changes in version 0.99.51

### New features

* Added a group of ID translation helper functions (see ?ID-translation)
* Added a group of helper functions that work with `curatedTCGAData`
* `UUIDtoBarcode` function added thanks to @seandavi
* Renamed `makeGRangesListFromTCGA` to `makeGRangesListFromCopyNumber`
* `makeSummarizedExperimentFromGISTIC` is now available to convert
`RTCGAToolbox`'s `FirehoseGISTIC` data class to `SummarizedExperiment`
* Created a function to merge external `colData` to a `MultiAssayExperiment`
`colData` slot
* Revamped vignette documentation

### Minor changes and bug fixes

* Improvements to `TCGAbiospec` and `TCGAbarcode`
* Updated `sampleTypes` and `clinicalNames` datasets
* Updated DESCRIPTION file with appropriate imports and exports
* Various improvements to `findGRangesCols`
* `generateMap` is now updated to the recent `MultiAssayExperiment` API with
improved example
* Updated `getFileNames` to most recent `RTCGAToolbox` API
* Various updates to data generating scripts in `data-raw` folder
* Format updates to NEWS file
* Added tests

## Changes in version 0.2.0

### New features

* Package renamed to `TCGAutils` for working with TCGA data

## Changes in version 0.1.0

### New features

* `TCGAtranslateID` now works with GDC API

### Minor changes and bug fixes

* Code cleaned up
* Added proper import directives

## Changes in version 0.0.70

### New features

* `makeGRangesListFromDataFrame` now moved to `GenomicRanges`
* `makeSummarizedExperimentFromDataFrame` now moved to `SummarizedExperiment`
* `getFileNames` function will obtain filenames used in `RTCGAToolbox`
* Improved `getFileNames` with `xml2` and `rvest` dependencies and removes the
`XML` dependency

### Minor changes and bug fixes

* `TCGAextract` now uses the `findGRangesCols` to automatically detect ranged
data columns
* Arguments in functions `TCGA*` now renamed to match `MultiAssayExperiment`
conventions
* Informative errors in `TCGAextract`

## Changes in version 0.0.60

* `makeGRangesListFromTCGA` data builds on `makeGRangesListFromDataFrame`
* `makeGRangesListFromDataFrame` and
`makeRangedSummarizedExperimentFromDataFrame` will be moving to standard
Bioconductor packages soon.
* `tcga` and `ccle` functions soon to be deprecated.
* Upcoming: `TCGAbarcode` will be modified for efficiency

## Changes in version 0.0.50

* Add your own identifier parsing function for generating a `sampleMap` in
`generateMap`!
* Add proper genome build to ranged based objects.
* Return `SummarizedExperiment` class for certain data types.
* Fix genome build bugs

## Changes in version 0.0.46

* `makeRSE` function for creating a `RangedSummarizedExperiment` object from a
data frame.
* Bug fixes to `getRangeNames` including the option to enter a regular
expression vector for finding ranged column names.
* `matchClinical` renamed to `TCGAmatchClinical`

## Changes in version 0.0.44

* `getRangedNames` function will try to extract minimum necessary names for
creating ranges (works on a vector of names)
* minor bug fixes to `TCGAbiospec`, `TCGAextract`, `makeGRangesList`

## Changes in version 0.0.40

* Package renamed to `BiocInterfaces`!
* `TCGA` specific functions now start with the letters `TCGA`
* Included: more examples of use of the `TCGAbarcode` function
* Updated `makeGRangesList` function to work with `tcga` and `ccle` data
    parameter functions

## Changes in version 0.0.2

* Added a `NEWS.md` file to track changes to the package.
* TCGAmisc now a standalone package! (previously in `RTCGAToolbox`)
* Provides helper functions for converting raw data into S4 objects (e.g.,
`GRangesList`)
* Provides functions for creating a MultiAssayExperiment object such as:
    * `generateTCGAmap`
    * `cleanExpList`
