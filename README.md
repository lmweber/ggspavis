# ggspavis

[![R build status](https://github.com/lmweber/ggspavis/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/lmweber/ggspavis/actions)

The `ggspavis` package contains visualization functions for spatial transcriptomics data, including functions to generate several types of plots, including spot plots, feature (molecule) plots, reduced dimension plots, spot-level quality control (QC) plots, and feature-level QC plots, for datasets from the 10x Genomics Visium and other technological platforms. Datasets are assumed to be in either [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment) or [SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment) Bioconductor format.

The package is available from [Bioconductor](https://bioconductor.org/packages/ggspavis).


## Installation

The release version of the package can be installed from Bioconductor:

```
install.packages("BiocManager")
BiocManager::install("ggspavis")
```

The latest development version can be installed either from the development version of Bioconductor or from GitHub:

```
remotes::install_github("lmweber/ggspavis")
```


## Tutorial

A vignette containing a tutorial and examples is available from the [Bioconductor](https://bioconductor.org/packages/ggspavis) package page.


## Citation

- Righell D.\*, Weber L.M.\*, Crowell H.L.\*, Pardo B., Collado-Torres L., Ghazanfar S., Lun A.T.L., Hicks S.C.\*, and Risso D.\* (2022). *SpatialExperiment: infrastructure for spatially resolved transcriptomics data in R using Bioconductor.* [Bioinformatics](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac299/6575443), btac299.
