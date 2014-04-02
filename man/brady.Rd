\name{brady}
\alias{brady}
\docType{data}
\title{
Bradykinin receptor antagonism
}
\description{
xxx.
}
\usage{data(brady)}
\format{
  A data frame with 560 observations on the following 4 variables.
  \describe{
    \item{\code{Drug}}{A factor with 3 levels giving the treatment arm.}
    \item{\code{ID}}{An integer vector running from 1 to 112 giving the patient ID.}
    \item{\code{Time}}{A factor with 5 levels giving the occasion of measurement.}
    \item{\code{logConc}}{A numeric vector giving the measured log-concentration of D-dimer.}
  }
}
\details{

SIMULATED DATA

}
\source{
For the original trial: \url{http://clinicaltrials.gov/ct2/show/NCT00223704}
}
\references{
Balaguer, J. M., Yu, C., Byrne, J. G., Ball, S. K., Petracek, M. R., Brown, N. J., Pretorius, M. (2013) Contribution of endogenous bradykinin to fibrinolysis, inflammation, and blood product transfusion following cardiac surgery: a randomized clinical trial. Clinical Pharmacology & Therapeutics, 93(4), 326--334.
}
\examples{
data(brady)
str(brady)
}
\keyword{datasets}