# TCGAutils

The `TCGAutils` package provides a suite of helper functions that aid in
the management and cleaning of data from The Cancer Genome Atlas (TCGA).

Many of the functions contained herein work on raw and derived data objects
from The Cancer Genome Atlas (TCGA), the `RTCGAToolbox` package and
`curatedTCGAData` experiment data package.

Please make sure to download the latest version of `RTCGAToolbox`
from Bioconductor.

```
BiocInstaller::biocLite("RTCGAToolbox")
```

Functions like `RTCGAToolbox::biocExtract` interact with derived data classes
(such as `FirehoseData`). Other functions like `TCGAexonToGRangesList` work on
text files to create a `GRangesList` object from exon level data.

Additional helper functions are available to clean up datasets in preparation
for integration to a `MultiAssayExperiment` object.

`MultiAssayExperiment` can be downloaded from both the release and development
versions of Bioconductor.

Please feel free to report bugs at our [github issue page][]

[github issue page]: https://github.com/waldronlab/TCGAutils/issues

