# ggspavis

[![R build status](https://github.com/lmweber/ggspavis/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/lmweber/ggspavis/actions)

The `ggspavis` package contains visualization functions for spatially resolved transcriptomics data stored in [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment) R/Bioconductor format.

The `ggspavis` package is available from [Bioconductor](https://bioconductor.org/packages/ggspavis).

A vignette containing examples and documentation is available from [Bioconductor](https://bioconductor.org/packages/ggspavis).


## Installation

The `ggspavis` package can be installed from Bioconductor. Note that Bioconductor follows a "release" and "development" schedule, where the release version is considered to be stable and updated every 6 months, and the development version contains latest updates (and then becomes the next release version every 6 months).


### Release version

To install the stable release version, install the latest release version of R from [CRAN](https://cran.r-project.org/), then install the Bioconductor package installer and the `ggspavis` package:

```
install.packages("BiocManager")
BiocManager::install("ggspavis")
```

### Development version

The development version can be installed from GitHub:

```
install.packages("remotes")
remotes::install_github("lmweber/ggspavis", build_vignettes = TRUE)
```


## Citation

Righell D.\*, Weber L.M.\*, Crowell H.L.\*, Pardo B., Collado-Torres L., Ghazanfar S., Lun A.T.L., Hicks S.C.<sup>+</sup>, and Risso D.<sup>+</sup> (2021), *SpatialExperiment: infrastructure for spatially resolved transcriptomics data in R using Bioconductor*, [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.01.27.428431v3).

