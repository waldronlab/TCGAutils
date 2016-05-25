# BiocInterfaces 

The `BiocInterfaces` package provides a suite of helper functions that aid in the
creation of several conventional Bioconductor classes such as, `GRangesList`, `ExpressionSet`,
and `MultiAssayExperiment`. Many of the functions contained herein work on raw and derived
data objects from The Cancer Genome Atlas (TCGA). Functions like `TCGAextract` interact with
derived data classes (i.e., `FirehoseData`). Other functions like `TCGAexonToGRangesList`
work on text files to create a `GRangesList` object from exon level data.

Please feel free to report bugs at our
[github page](https://github.com/waldronlab/TCGAmisc).
