\name{vitaminE}
\alias{vitaminE}
\docType{data}
\title{
Vitamin E Supplementation
}
\description{
Longitudinal study on the effect of vitamin E supplements on weight gain in guinea pigs.
}
\usage{data(vitaminE)}
\format{
  A data frame with 90 observations on the following 4 variables.
  \describe{
    \item{\code{dose}}{A factor with 3 levels giving the dose level.}
    \item{\code{animal}}{A factor with 15 levels giving the ID of the animal.}
    \item{\code{weight}}{A numeric vector giving the animal's body weight.}
    \item{\code{week}}{A numeric vector giving the week by the end of which the measurement was taken.}
  }
}
\details{
15 guinea pigs were treated with a growth inhibitor during the first week of a dietary study, then their body weight (in grams) was measured at the end of weeks 1, 3, 4, 5, 6, and 7. After week 4 the animals were randomly assigned to a dose group (\code{zero}, \code{low}, or \code{high}), and vitamin E therapy started at the beginning of week 5.
}
\source{
Crowder, M. J. and Hand, D. J. (1990) Analysis of Repeated Measures. Chapman & Hall, London, UK, pp. 27--29.
}
\references{
Hand, D. J., Daly, F., Lunn, A. D., McConway, K. J., Ostrowski, E. (1994) A Handbook of Small Data Sets. Chapman & Hall, London, UK, p. 325.
}
\examples{
data(vitaminE)
str(vitaminE)
}
\keyword{datasets}