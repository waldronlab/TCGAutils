# BiocInterfaces 

The `BiocInterfaces` package provides a suite of helper functions that aid in
the creation of several conventional Bioconductor classes such as,
`GRangesList`, `ExpressionSet`, and `MultiAssayExperiment`. Many of the
functions contained herein work on raw and derived data objects from The
Cancer Genome Atlas (TCGA). Functions like `TCGAextract` interact with
derived data classes (i.e., `FirehoseData`). Other functions like
`TCGAexonToGRangesList` work on text files to create a `GRangesList` object
from exon level data. Additional helper functions are available to clean up
datasets in preparation for integration to a `MultiAssayExperiment` class
object. For more on the `MultiAssayExperiment`, please see the [github repo][]
or download it from `bioc-devel`.

Please feel free to report bugs at our [github issue page][]

[github issue page]: https://github.com/waldronlab/BiocInterfaces
[github repo]: https://github.com/vjcitn/MultiAssayExperiment
