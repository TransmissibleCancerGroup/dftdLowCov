# dftdLowCov
R package of code relevant to DFTD project

`dftdLowCov` collects utility functions used in the analysis and figure-plotting code of the DFTD low coverage project.

## Installation

    require(devtools)
    devtools::install_github("TransmissibleCancerGroup/dftdLowCov")

## Data files

The package contains a number of data files. These are the also the files included as supplementary tables in the paper. The files can be listed with the command `list.files(system.file("extdata", "", package = "dftdLowCov"))`. To get the path of a specific file, e.g. "2020-05-07_TableS3_CNVs.xlsx", use the command `system.file("extdata", "2020-05-07_TableS3_CNVs.xlsx", package = "dftdLowCov")`.
