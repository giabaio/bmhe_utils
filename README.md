# BMHE 
## A utility package for post-processing of Bayesian models performed in `OpenBUGS` or `JAGS`

## Installation
The package is only available from this `GitHub` repository, for now. On Windows machines, you need to install a few dependencies, including [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first, e.g. by running
```r
pkgs <- c("MASS", "Rtools", "remotes")
repos <- c("https://cran.rstudio.com", "https://inla.r-inla-download.org/R/stable") 
install.packages(pkgs, repos=repos, dependencies = "Depends")
```
before installing the package using `remotes`:
```r
remotes::install_github("giabaio/BCEA", ref="dev")
```

Under Linux or MacOS, it is sufficient to install the package via `remotes`:
```r
install.packages("remotes")
remotes::install_github("giabaio/BCEA", ref="dev")
```

## Licence
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Contributing
Please submit contributions through `Pull Requests`. To report issues and/or seek support, please file a new ticket in the [issue](https://github.com/giabaio/bmhe_utils/issues) tracker.
