\name{SimLongi}
\alias{SimLongi}
\title{
Simultaneous Inference for Longitudinal Data
}
\description{
Simultaneous inference for a set of contrasts (linear combinations) of means in longitudinal scenarios. Computes multiplicity-adjusted p-values and simultaneous confidence intervals for comparing groups at multiple time points, or comparing time points in multiple groups, or both.
}
\usage{
SimLongi(data, response, group, time, id, covariates=NULL,
         var=list("hom", "het", "hett"), cor=list("CS", "AR1", "UN"),
         contrasts=NULL, type="Dunnett", base=1, direction="gpt",
         alternative="two.sided", level=0.95, df="ess")
}
\arguments{
  \item{data}{A data frame.}
  \item{response}{A character string giving the name of the response variable in \code{data}.}
  \item{group}{A character string giving the name of the (treatment) group variable in \code{data}.}
  \item{time}{A character string giving the name of the time variable in \code{data}.}
  \item{id}{A character string giving the name of the subject variable in \code{data}.}
  \item{covariates}{ccc}
  \item{var}{A list of character strings itemizing the patterns of heteroscedasticity to be included in model selection. The alternatives are \code{hom} (equal variances across groups and time points), \code{het} (different variances for all groups and time points), \code{hett} (variances equal across groups but different for time points), and \code{hetg} (variances equal across time points but different for groups).}
  \item{cor}{A list of character strings itemizing the correlation structures to be included in model selection. The alternatives are \code{CS} (compound symmetry), \code{AR1} (first-order autoregressive), \code{CAR1} (continuous first-order autoregressive), \code{AR2} (second-order autoregressive), \code{MA1} (first-order moving average), \code{MA2} (second-order moving average), \code{ARMA11} (first-order autoregressive moving average), and \code{UN} (unstructured).}
  \item{contrasts}{. Default to \code{NULL}. See examples.}
  \item{type}{A character string defining the type of contrast matrix (i.e., the set of comparisons); ignored unless \code{contrasts=NULL}. Two widespread and often useful choices are \code{Dunnett} (many-to-one i.e., compare each group against a common reference) and \code{Tukey} (all-pairs i.e., pairwise comparisons of all groups). Further options are \code{Sequen}, \code{AVE}, \code{Changepoint}, \code{Williams}, \code{Marcus}, \code{McDermott}, \code{UmbrellaWilliams}, and \code{GrandMean}. For \code{direction="both"} a vector of two \code{character} strings may be inserted: the first one specifies the \code{gpt} comparisons and the second one the \code{tpg} comparisons. If only a single \code{character} string is inserted when \code{direction="both"}, this set of contrasts is applied to \code{gpt} and \code{tpg} comparisons.}
  \item{base}{An integer specifying the reference group with many-to-one comparisons; ignored unless \code{contrasts=NULL} and \code{type="Dunnett"}. For \code{direction="both"} and \code{type="Dunnett"} for both \code{gpt} and \code{tpg} comparisons, a vector of two \code{integer}s may be inserted to define reference groups for both factors (group and time). See examples.}
  \item{direction}{Defines which factor's levels are to be compared at each level of the other factor; ignored unless \code{contrasts=NULL}. \code{gpt} ("groups per time") invokes comparisons among groups separately and simultaneously for each time point. By contrast, \code{tpg} ("times per group") compares time points separately and simultaneously for each group. \code{both} carries out jointly both types of comparisons.}
  \item{alternative}{ccc}
  \item{level}{A numeric value giving the simultaneous confidence level (1 - alpha).}
  \item{df}{A character string specifying the approximation to the degrees of freedom for the multivariate t-distribution. Must be one of \code{ess} (effective sample size), \code{adj} (adjusted), \code{pb} (Pinheiro-Bates), \code{satt} (Satterthwaite), \code{kr} (Kenward-Roger), \code{con} (containment), \code{naive} (naive), \code{res} (residual), or \code{normal} (multivariate normal i.e., "infinite" degrees of freedom). See details.}
}
\details{

AICc model selection

four variance structures

specified corration structures


\code{ess} uses as degrees of freedom the effective sample size computed under an AR(1)-assumption (Faes et al. 2009).

Computation of
in the book of Pinheiro & Bates (2000, p. 91)


Kenward & Roger (1997, 2009)

\code{kr} uses the Kenward-Roger implementation of Halekoh & Hojsgaard (2013 ZITAT)


\code{naive} calculates the number of independent units minus the number of cell means; it is prone to make results conservative (i.e., not exploit their type I error level). In contrast, \code{residual} calculates the total number of observations minus the number of cell means and is therefore likely cause anticonservatism.

override model selection procedure by just inserting one var and one cor pattern

}
\value{
An list with elements
  \item{Results}{A table listing comparisonwise the estimated difference with standard error, lower and upper simultaneous confidence bounds, values of the test statistics, and multiplicity-adjusted p-values.}
  \item{CovStat}{The covariance matrix of test statistics.}
  \item{CritValue}{.}
  \item{Alternative}{.}
  \item{ConfLevel}{The confidence level as specified via \code{level}.}
  \item{DFMethod}{.}
  \item{DF}{.}
  \item{BestMod}{Formula of the AICc-best model.}
  \item{ModSelTab}{Model selection table.}
  \item{AWBest}{Akaike weight of the AICc-best model (with respect to the set of models considered).}
  \item{CovBest}{The estimated covariance matrix for the AICc-best model.}
  \item{Model}{The fit of the AICc-best model.}
}
\references{
Burnham, K. P. & Anderson, D. R. (2002) Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach. Second Edition. Springer, New York, NY.

Faes, C., Molenberghs, G., Aerts, M., Verbeke, G., Kenward, M. G. (2009) The effective sample size and an alternative small-sample degrees-of-freedom method. The American Statistician, 63(4), 389--399.

Hothorn, T., Bretz, F., Westfall, P. (2008) Simultaneous inference in general parametric models. Biometrical Journal, 50(3), 346--363.

Kenward, M. G. & Roger, J. H. (1997) Small sample inference for fixed effects from restricted maximum likelihood. Biometrics, 53(3), 983--997.

Kenward, M. G. & Roger, J. H. (2009) An improved approximation to the precision of fixed effects from restricted maximum likelihood. Computational Statistics and Data Analysis, 53(7), 2583--2595.

Pinheiro, J. C. & Bates, D. M. (2000) Mixed-Effects Models in S and S-PLUS. Springer, New York, NY.

Satterthwaite, F. E. (1941) Synthesis of variance. Psychometrika, 6(5), 309--316.
}
\author{
Philip Pallmann \email{pallmann@biostat.uni-hannover.de}
}
\note{
Vorsicht bei CAR1
}
\seealso{
\code{\link{SimLongiMMM}}
}
\examples{
data(heart)

# Many-to-one comparisons of groups per time point
# taking the third group ("control") as reference

SimLongi(data=heart, response="heartrate", group="drug",
         time="time", id="person", direction="gpt",
         type="Dunnett", base=3)$Results

# The same test using a self-built contrast matrix

library(multcomp)
cmat <- kronecker(contrMat(numeric(3), "Dunnett", base=3), diag(4))
SimLongi(data=heart, response="heartrate", group="drug",
         time="time", id="person", contrasts=cmat)$Results

# All-pairs comparisons of groups per time point,
# and many-to-one comparisons of time points per group
# taking (by default) the first time point ("T1") as reference

SimLongi(data=heart, response="heartrate", group="drug",
         time="time", id="person", direction="both",
         type=c("Tukey", "Dunnett"))$Results

# Many-to-one comparisons of groups per time point
# taking the third group ("control") as reference,
# and many-to-one comparisons of time points per group
# taking the last time point ("T4") as reference

SimLongi(data=heart, response="heartrate", group="drug",
         time="time", id="person", direction="both",
         type="Dunnett", base=c(3, 4))$Results

# DFS HIER NOCH EINBAUEN! (UND EVTL. WEITERE NEUERUNGEN)

}
\keyword{htest}