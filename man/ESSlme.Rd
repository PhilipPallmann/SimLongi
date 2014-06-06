\name{ESSlme}
\alias{ESSlme}
\title{
Effective Sample Size for \code{lme} Parameters
}
\description{
Compute the effective sample sizes associated with the parameter estimates of a linear mixed-effects model fitted via \code{lme}.
}
\usage{
ESSlme(model, IDcolumn)
}
\arguments{
  \item{model}{A linear mixed-effects model fit of class \code{lme} from \code{nlme}.}
  \item{IDcolumn}{An integer giving the number of the column in the \code{data.frame} on which the model fit is based that contains the ID variable (cluster identifier).}
}
\details{
The \code{data.frame} used for model fitting must be sorted along the ID variable, otherwise you might get mock results!
}
\value{
A vector giving the effective sample sizes associated with the estimated parameters of the \code{lme} model.
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
\code{\link{ESSgls}}
}
\examples{
### Using the example from the lme help page
library(nlme)
data(Orthodont)
lme1 <- lme(distance ~ age + Sex, Orthodont, random=~1)
#ESSlme(model=lme1, IDcolumn=3)
}
\keyword{htest}