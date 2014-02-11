\name{satterthwaite}
\alias{satterthwaite}
\title{
Satterthwaite Approximation
}
\description{
Computes Satterthwaite's approximation to the degrees of freedom.
}
\usage{
satterthwaite(var, n, df=n-1)
}
\arguments{
  \item{var}{A vector of numeric values giving the independent sample variances.}
  \item{n}{A vector of integers giving the sample sizes.}
  \item{df}{A vector of numeric values giving the degrees of freedom for the independent sample variances. Default is \code{n-1}.}
}
\value{
A numeric value giving the pooled degrees of freedom according to the Satterthwaite approximation.
}
\references{
Satterthwaite, F. E. (1941) Synthesis of variance. Psychometrika, 6(5), 309--316.

Satterthwaite, F. E. (1946) An approximate distribution of estimates of variance components. Biometrics Bulletin, 2(6), 110--114.
}
\author{
Philip Pallmann \email{pallmann@biostat.uni-hannover.de}
}
\examples{
# Perfectly homoscedastic
satterthwaite(rep(2, 3), rep(10, 3))

# Mildly heteroscedastic
satterthwaite(c(1.9, 2, 2.1), rep(10, 3))

# Strongly heteroscedastic
satterthwaite(1:3, rep(10, 3))
}
\keyword{htest}