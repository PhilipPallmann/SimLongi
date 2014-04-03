\name{SimLongiMMM}
\alias{SimLongiMMM}
\title{
Simultaneous Inference for Longitudinal Data Using Multiple Marginal Models
}
\description{
Simultaneous inference for a set of contrasts (linear combinations) of means in longitudinal scenarios. Computes multiplicity-adjusted p-values and simultaneous confidence intervals for comparing groups at multiple time points by combining time-point-specific marginal models as described by Pipper et al. (2012).
}
\usage{
SimLongiMMM(data, response, group, time, id, covariates=NULL,
            contrasts=NULL, type="Dunnett", base=1,
            alternative="two.sided", level=0.95, refdist="normal")
}
\arguments{
  \item{data}{A data frame.}
  \item{response}{A character string giving the name of the response variable in \code{data}.}
  \item{group}{A character string giving the name of the (treatment) group variable in \code{data}.}
  \item{time}{A character string giving the name of the time variable in \code{data}.}
  \item{id}{A character string giving the name of the subject variable in \code{data}.}
  \item{covariates}{ccc}
  \item{contrasts}{. Default to \code{NULL}. See examples.}
  \item{type}{A character string defining the type of contrast matrix (i.e., the set of comparisons); ignored unless \code{contrasts=NULL}. Two widespread and often useful choices are \code{Dunnett} (many-to-one i.e., compare each group against a common reference) and \code{Tukey} (all-pairs i.e., pairwise comparisons of all groups). Further options are \code{Sequen}, \code{AVE}, \code{Changepoint}, \code{Williams}, \code{Marcus}, \code{McDermott}, \code{UmbrellaWilliams}, and \code{GrandMean}.}
  \item{base}{An integer specifying the reference group with many-to-one comparisons; ignored unless \code{contrasts=NULL} and \code{type="Dunnett"}.}
  \item{alternative}{The direction of the alternative to be tested against. Default is \code{two.sided}. Options for one-sided testing are \code{greater} or \code{less}.}
  \item{level}{A numeric value giving the simultaneous confidence level (1 - alpha).}
  \item{refdist}{.}
}
\details{
%%  If necessary, more details than the description above
}
\value{
A list of class \code{silo} with elements
  \item{Results}{A table listing comparisonwise the estimated difference with standard error, lower and upper simultaneous confidence bounds, value of the test statistic, and multiplicity-adjusted p-value.}
  \item{CovStat}{The covariance matrix of test statistics.}
  \item{CritValue}{The critical value (equicoordinate quantile from a multivariate t-distribution).}
  \item{Alternative}{The direction of the alternative.}
  \item{ConfLevel}{The confidence level as specified via \code{level}.}
  \item{RefDist}{The reference distribution (multivariate normal or t).}
  \item{DF}{The degrees of freedom used for the multivariate t distribution (zero if multivariate normal).}
  \item{Corr}{The estimated correlation matrix of time points and test statistics.}
  \item{ContMat}{The contrast matrix applied to each point in time.}
}
\references{
Hothorn, T., Bretz, F., Westfall, P. (2008) Simultaneous inference in general parametric models. Biometrical Journal, 50(3), 346--363.

Pipper, C. B., Ritz, C., Bisgaard, H. (2012) A versatile method for confirmatory evaluation of the effects of a covariate in multiple models. Journal of the Royal Statistical Society, Series C: Applied Statistics, 61(2), 315--326.
}
\author{
Philip Pallmann \email{pallmann@biostat.uni-hannover.de}
}
\note{
bla
}
\seealso{
\code{\link{SimLongi}}, \code{\link{SimLongiMix}}
}
\examples{
data(heart)

# Many-to-one comparisons of groups per time point
# taking the third group ("control") as reference

SimLongiMMM(data=heart, response="heartrate", group="drug",
            time="time", id="person", type="Dunnett",
            base=3)$Results

# For comparison:

SimLongi(data=heart, response="heartrate", group="drug",
         time="time", id="person", direction="gpt",
         type="Dunnett", base=3, df="normal")$Results

# With multivariate t as reference distribution:

SimLongiMMM(data=heart, response="heartrate", group="drug",
            time="time", id="person", type="Dunnett",
            base=3, refdist="t")$Results

% WARUM SIND DIE SE'S JETZT IMMER GLEICH FUER GLEICHE ZEITPUNKTE???
}
\keyword{htest}