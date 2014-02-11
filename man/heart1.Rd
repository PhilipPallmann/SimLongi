\name{heart}
\alias{heart}
\docType{data}
\title{
Heart Rates
}
\description{
Heart rates dataset from Milliken & Johnson
}
\usage{data(heart)}
\format{
  A data frame with 96 observations on the following 4 variables.
  \describe{
    \item{\code{drug}}{a numeric vector}
    \item{\code{time}}{a numeric vector}
    \item{\code{person}}{a numeric vector}
    \item{\code{heartrate}}{a numeric vector}
  }
}
\details{
%%  If necessary, more details than the __description__ above
}
\source{
%%  reference to a publication or URL from which the data were obtained
}
\references{
Milliken, G. A. and Johnson, D. E. (1992) Analysis of Messy Data. Volume I: Designed Experiments. Chapman & Hall, London, UK, p.
}
\examples{
data(heart)
str(heart)
}
\keyword{datasets}