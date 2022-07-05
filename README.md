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

## Use
Load the package into the `R` workspace as usual
```
library(bmhe)
```
and use all the available functions. Roughly speaking, these can be divided into "plotting", "printing" and "utility".

### Plotting

- `betaplot`	Trial-and-error Beta plot (using `[manipulate](https://cran.r-project.org/web/packages/manipulate/index.html)`)
- `coefplot`	"Coefplot" for the parameters in the model (using `[tidyverse](https://www.tidyverse.org/)`)
- `diagplot`	Specialised diagnostic plots to check convergence and autocorrelation of the MCMC run
- `gammaplot`	Trial-and-error Gamma plot (using `[manipulate](https://cran.r-project.org/web/packages/manipulate/index.html)`)
- `posteriorplot`	Various plots for the posteriors in a 'bugs' or 'jags' object
- `traceplot`	Makes a traceplot (eg to visualise MCMC simulations from multiple chains, using `[tidyverse](https://www.tidyverse.org/)`)

### Printing

- `print.bugs`	Modifies the built-in print method for the `R2OpenBUGS` package to provide a few more options and standardisation
- `print.rjags`	Modifies the built-in print method for the `R2jags` package to provide a few more options and standardisation
- `stats`	Computes and prints summary statistics for a vector or matrix of simulated values

### Utility
- `betaPar`	Computes the parameters of a Beta distribution so that the mean and standard dev are the input (m,s)
- `betaPar2`	Compute the parameters of a Beta distribution, given a prior guess for key parameters. Based on "Bayesian ideas and data analysis", page 100. Optimisation method to identify the values of a,b that give required conditions on the Beta distribution
- `ilogit`	Computes the inverse logit of a number between -infinity and +infinity
- `logit`	Computes the logit of a number
- `lognPar`	Computes mean and variance of a logNormal distribution so that the parameters on the natural scale are mu and sigma
- `odds2probs`	Maps from odds to probabilities
- `OR`	Computes the odds ratio between two probabilities

## Licence
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Contributing
Please submit contributions through `Pull Requests`. To report issues and/or seek support, please file a new ticket in the [issue](https://github.com/giabaio/bmhe_utils/issues) tracker.
