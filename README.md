# TCGAutils

The `TCGAutils` package provides a suite of helper functions that aid in
the creation of several conventional Bioconductor classes such as,
`GRangesList`, `ExpressionSet`, `RaggedExperiment`, and `MultiAssayExperiment`.

Many of the functions contained herein work on raw and derived data objects
from The Cancer Genome Atlas (TCGA) and the `RTCGAToolbox` package.

Please make sure to download the latest and compatible version of the
`RTCGAToolbox` package by running: 

```
BiocInstaller::biocLite("LiNk-NY/RTCGAToolbox")
```

Functions like `TCGAextract` interact with derived data classes
(such as `FirehoseData`). Other functions like `TCGAexonToGRangesList` work on
text files to create a `GRangesList` object from exon level data.

Additional helper functions are available to clean up datasets in preparation
for integration to a `MultiAssayExperiment` object.

For more on the `MultiAssayExperiment`, please see the [github repo][]
or download it from `bioc-release`.

Please feel free to report bugs at our [github issue page][]

[github issue page]: https://github.com/waldronlab/TCGAutils
[github repo]: https://github.com/waldronlab/MultiAssayExperiment

