\name{SimLongi-package}
\alias{SimLongi-package}
\docType{package}
\title{
A package for simultaneous inference in longitudinal data settings
}
\description{
Simultaneous p-values and confidence intervals for comparing several groups per time point and/or several time points per group using multiple contrasts.
}
\details{
\tabular{ll}{
Package: \tab SimLongi\cr
Type: \tab Package\cr
Version: \tab 1.1-2\cr
Date: \tab 2014-04-04\cr
License: \tab GPL (>= 2)\cr
}
Multiple contrast tests using parameter and covariance estimates from a general linear model (with \code{SimLongi}), or a linear mixed-effects model (with \code{SimLongiMix}), or from combining multiple marginal models using the method of Pipper et al. (2012) (with \code{SimLongiMMM}). A convenient plotting function for simultaneous confidence intervals is also available (\code{PlotCI}), as are several example datasets (\code{brady}, \code{heart1}, \code{vitaminE}).
}
\author{
Philip Pallmann <pallmann@biostat.uni-hannover.de>
}
\references{
Pipper, C. B., Ritz, C., Bisgaard, H. (2012) A versatile method for confirmatory evaluation of the effects of a covariate in multiple models. Journal of the Royal Statistical Society, Series C: Applied Statistics, 61(2), 315--326.
}
\keyword{
package
}
\seealso{
\code{\link{SimLongi}}, \code{\link{SimLongiMix}}, \code{\link{SimLongiMMM}}
}
\examples{
## SimLongi()
## SimLongiMix()
## SimLongiMMM()
}
