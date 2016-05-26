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
