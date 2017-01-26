# TCGAutils 0.2.0

## New features

* Package renamed to `TCGAutils` for working with TCGA data

# BiocInterfaces 0.1.0

## New features

* `TCGAtranslateID` now works with GDC API

## Minor changes and bug fixes 

* Code cleaned up
* Added proper import directives

# BiocInterfaces 0.0.70

## New features

* `makeGRangesListFromDataFrame` now moved to `GenomicRanges`
* `makeSummarizedExperimentFromDataFrame` now moved to `SummarizedExperiment`
* `getFileNames` function will obtain filenames used in `RTCGAToolbox`
* Improved `getFileNames` with `xml2` and `rvest` dependencies and removes the
`XML` dependency

## Minor changes and bug fixes

* `TCGAextract` now uses the `findGRangesCols` to automatically detect ranged
data columns
* Arguments in functions `TCGA*` now renamed to match `MultiAssayExperiment`
conventions
* Informative errors in `TCGAextract`

# BiocInterfaces 0.0.60

* `makeGRangesListFromTCGA` data builds on `makeGRangesListFromDataFrame`
* `makeGRangesListFromDataFrame` and `makeRangedSummarizedExperimentFromDataFrame` will be
moving to standard Bioconductor packages soon. 
* `tcga` and `ccle` functions soon to be deprecated. 
* Upcoming: `TCGAbarcode` will be modified for efficiency

# BiocInterfaces 0.0.50

* Add your own identifier parsing function for generating a `sampleMap` in `generateMap`!
* Add proper genome build to ranged based objects.
* Return `SummarizedExperiment` class for certain data types.
* Fix genome build bugs

# BiocInterfaces 0.0.46

* `makeRSE` function for creating a `RangedSummarizedExperiment` object from a
data frame. 
* Bug fixes to `getRangeNames` including the option to enter a regular expression
vector for finding ranged column names. 
* `matchClinical` renamed to `TCGAmatchClinical`

# BiocInterfaces 0.0.44

* `getRangedNames` function will try to extract minimum necessary names for creating ranges 
(works on a vector of names)
* minor bug fixes to `TCGAbiospec`, `TCGAextract`, `makeGRangesList`

# BiocInterfaces 0.0.40

* Package renamed to `BiocInterfaces`!
* `TCGA` specific functions now start with the letters `TCGA`
* Included: more examples of use of the `TCGAbarcode` function
* Updated `makeGRangesList` function to work with `tcga` and `ccle` data
    parameter functions

# TCGAmisc 0.0.2

* Added a `NEWS.md` file to track changes to the package.
* TCGAmisc now a standalone package! (previously in `RTCGAToolbox`)
* Provides helper functions for converting raw data into S4 objects (e.g.,
`GRangesList`)
* Provides functions for creating a MultiAssayExperiment object such as:
    * `generateTCGAmap`
    * `cleanExpList`
