% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{diagplot}
\alias{diagplot}
\title{Specialised diagnostic plots}
\usage{
diagplot(x, what = "Rhat", ...)
}
\arguments{
\item{x}{an object of class `bugs', see \code{\link{bugs}}, or of class
'jags', see \code{\link{jags}} for details}

\item{what}{A string indicating what diagnostic measure should be plotted.
Options are 'Rhat' (default), indicating the PSR statistic, or 'n.eff',
indicating the effective sample size}

\item{...}{Additional options}
}
\description{
Creates a plot showing the output of convergence indicators, such as
the Potential Scale Reduction and the effective sample size
}
\examples{
\dontrun{
}
}
\seealso{
\code{\link{bugs}}, \code{\link{jags}}
}
\author{
Gianluca Baio
}
\keyword{Diagnostic}
\keyword{plots}
