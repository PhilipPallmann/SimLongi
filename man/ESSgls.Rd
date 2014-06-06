\name{ESSgls}
\alias{ESSgls}
\title{
Effective Sample Size for \code{gls} Parameters
}
\description{
Compute the effective sample sizes associated with the parameter estimates of a general linear model fitted via \code{gls}.
}
\usage{
ESSgls(model, IDcolumn)
}
\arguments{
  \item{model}{A general linear model fit of class \code{gls} from \code{nlme}.}
  \item{IDcolumn}{An integer giving the number of the column in the \code{data.frame} on which the model fit is based that contains the ID variable (cluster identifier).}
}
\details{
The \code{data.frame} used for model fitting must be sorted along the ID variable, otherwise you might get mock results!
}
\value{
A vector giving the effective sample sizes associated with the estimated parameters of the \code{gls} model.
}
\references{
Faes, C., Molenberghs, G., Aerts, M., Verbeke, G., Kenward, M. G. (2009) The effective sample size and an alternative small-sample degrees-of-freedom method. The American Statistician, 63(4), 389--399.
}
\author{
Philip Pallmann \email{pallmann@biostat.uni-hannover.de}
}
\note{
Models with multiple hierarchical levels are currently not supported.
}
\seealso{
\code{\link{ESSlme}}
}
\examples{
### Using the example from the gls help page
library(nlme)
data(Ovary)
gls1 <- gls(follicles ~ sin(2 * pi * Time) + cos(2 * pi * Time), Ovary,
            correlation=corAR1(form=~1|Mare))
#ESSgls(model=gls1, IDcolumn=1)
}
\keyword{htest}