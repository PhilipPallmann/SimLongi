\name{SimLongiMix}
\alias{SimLongiMix}
\title{
Simultaneous Inference for Longitudinal Data Using Linear Mixed-Effects Modeling
}
\description{
Simultaneous inference for a set of contrasts (linear combinations) of means in longitudinal scenarios. Computes multiplicity-adjusted p-values and simultaneous confidence intervals for comparing groups at multiple time points, or comparing time points in multiple groups.
}
\usage{
SimLongiMix(data, response, group, time, id, covariates=NULL,
            rand=list("1|id", "time|id", "group|id", "timegroup|id"),
            contrasts=NULL, type="Dunnett", base=1, direction="gpt",
            alternative="two.sided", level=0.95, df="kr")
}
\arguments{
  \item{data}{A data frame.}
  \item{response}{A character string giving the name of the response variable in \code{data}.}
  \item{group}{A character string giving the name of the (treatment) group variable in \code{data}.}
  \item{time}{A character string giving the name of the time variable in \code{data}.}
  \item{id}{A character string giving the name of the subject variable in \code{data}.}
  \item{covariates}{ccc}
  \item{rand}{A list containing the random effects structures to be employed. There are four options of increasing complexity: \code{1|id} for random subject effects; \code{time|id} for time-point-specific random subject effects; \code{group|id} for treatment-specific random subject effects; and \code{timegroup|id} for time-point- and treatment-specific random subject effects. Note that for the latter option the correlation between time-point- and group-specific random effects is (modelled as) zero. Also note that complex models can cause instabilities when being fitted.}
  \item{contrasts}{An optional matrix of appropriate dimensions defining the contrasts to be applied. Default to \code{NULL}. See examples.}
  \item{type}{A character string defining the type of contrast matrix (i.e., the set of comparisons); ignored unless \code{contrasts=NULL}. Two widespread and often useful choices are \code{Dunnett} (many-to-one i.e., compare each group against a common reference) and \code{Tukey} (all-pairs i.e., pairwise comparisons of all groups). Further options are \code{Sequen}, \code{AVE}, \code{Changepoint}, \code{Williams}, \code{Marcus}, \code{McDermott}, \code{UmbrellaWilliams}, and \code{GrandMean}.}
  \item{base}{An integer specifying the reference group with many-to-one comparisons; ignored unless \code{contrasts=NULL} and \code{type="Dunnett"}.}
  \item{direction}{Defines which factor's levels are to be compared at each level of the other factor; ignored unless \code{contrasts=NULL}. \code{gpt} ("groups per time") invokes comparisons among groups separately and simultaneously for each time point. By contrast, \code{tpg} ("times per group") compares time points separately and simultaneously for each group.}
  \item{alternative}{The direction of the alternative to be tested against. Default is \code{two.sided}. Options for one-sided testing are \code{greater} or \code{less}.}
  \item{level}{A numeric value giving the simultaneous confidence level (1 - alpha).}
  \item{df}{A character string specifying the approximation to the degrees of freedom for the multivariate t-distribution. Must be one of \code{kr} (Kenward-Roger), \code{pb} (Pinheiro-Bates), \code{ess} (effective sample size), \code{adj} (adjusted), \code{naive} (naive), \code{res} (residual), or \code{normal} (multivariate normal instead of t). See details.}
}
\details{
The function performs time-point-wise comparisons of treatment groups, or treatment-group-wise comparisons of points in time, using multiple contrast as described by Hothorn et al. (2008). Test statistics are built with fixed-effects and covariance estimates from an appropriately parameterized linear mixed-effects model (e.g., Verbeke & Molenberghs 2000). If \code{rand} contains more than one element, AICc model selection (Burnham & Anderson 2002) is employed for selecting a "best-fitting" model to base further inferences on. Both multiplicity-adjusted p-values and simultaneous confidence intervals are provided.

Several approximations to the degrees of freedom for the multivariate t-distribution can be chosen. \code{kr} computes the approximation of Kenward & Roger (1997) as implemented in package \code{pbkrtest}. \code{pb} invokes the containment degrees of freedom as described by Pinheiro & Bates (2000, p. 91) and implemented in their \code{nlme} package. \code{naive} calculates the number of independent units minus the number of cell means; it is prone to make results conservative (i.e., not exploit their type I error level). In contrast, \code{residual} compute the total number of observations minus the number of cell means and is therefore likely cause anticonservatism. \code{normal} uses a critical point from a multivariate normal distribution (i.e., a multivariate t-distribution at "infinite" degrees of freedom).
}
\value{
A list of class \code{silo} with elements
  \item{Results}{A table listing comparisonwise the estimated difference with standard error, lower and upper simultaneous confidence bounds, value of the test statistic, and multiplicity-adjusted p-value.}
  \item{CovStat}{The covariance matrix of test statistics.}
  \item{CritValue}{The critical value (equicoordinate quantile from a multivariate t-distribution).}
  \item{Alternative}{The direction of the alternative.}
  \item{ConfLevel}{The confidence level as specified via \code{level}.}
  \item{DFMethod}{The approximation to the degrees of freedom.}
  \item{DF}{The degrees of freedom used for the multivariate t-distribution (zero if multivariate normal).}
  \item{ContMat}{The contrast matrix.}
  \item{BestMod}{The formula (in \code{R} syntax) of the AICc-chosen model.}
  \item{ModSelTab}{A model selection table.}
  \item{AWBest}{The Akaike weight of the AICc-chosen model (with respect to the set of models considered).}
  \item{CovBest}{The estimated covariance matrix of the AICc-chosen model.}
  \item{Model}{The fit of the AICc-chosen model.}
}
\references{
Burnham, K. P. & Anderson, D. R. (2002) Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach. Second Edition. Springer, New York, NY.

Hothorn, T., Bretz, F., Westfall, P. (2008) Simultaneous inference in general parametric models. Biometrical Journal, 50(3), 346--363.

Kenward, M. G. & Roger, J. H. (1997) Small sample inference for fixed effects from restricted maximum likelihood. Biometrics, 53(3), 983--997.

Pinheiro, J. C. & Bates, D. M. (2000) Mixed-Effects Models in S and S-PLUS. Springer, New York, NY.

Verbeke, G. & Molenberghs, G. (2000) Linear Mixed Models for Longitudinal Data. Springer, New York, NY.
}
\author{
Philip Pallmann \email{pallmann@biostat.uni-hannover.de}
}
\seealso{
\code{\link{SimLongi}}, \code{\link{SimLongiMMM}}
}
\examples{
data(heart)

# Many-to-one comparisons of groups per time point
# taking the third group ("control") as reference

Mix <- SimLongiMix(data=heart, response="heartrate", group="drug",
                   time="time", id="person", rand=list("1|id", "time|id"), 
                   direction="gpt", type="Dunnett", base=3)

Mix$Results

# The simplest model was chosen:
Mix$BestMod

# Kenward-Roger-approximated denominator degrees of freedom:
Mix$DF

# A graphical display of simultaneous confidence intervals:
PlotCI(Mix)
}
\keyword{htest}