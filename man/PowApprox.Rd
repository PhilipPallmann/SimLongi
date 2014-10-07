\name{PowApprox}
\alias{PowApprox}
\title{
Bla
}
\description{
Bla
}
\usage{
PowApprox(mu, cov, n, df, cmat, rhs = 0, alpha = 0.05, alternative = c("two.sided", "less", "greater"), power = c("global", "anypair", "allpairs"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{cov}{
%%     ~~Describe \code{cov} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{df}{
%%     ~~Describe \code{df} here~~
}
  \item{cmat}{
%%     ~~Describe \code{cmat} here~~
}
  \item{rhs}{
%%     ~~Describe \code{rhs} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{alternative}{
%%     ~~Describe \code{alternative} here~~
}
  \item{power}{
%%     ~~Describe \code{power} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mu, cov, n, df, cmat, rhs = 0, alpha = 0.05, alternative = c("two.sided", 
    "less", "greater"), power = c("global", "anypair", "allpairs")) 
{
    MU <- matrix(mu, ncol = 1)
    COV <- cmat \%*\% (cov/n) \%*\% t(cmat)
    R <- cov2cor(COV)
    M <- nrow(cmat)
    RHS <- rep(rhs, M)
    ExpT <- (as.numeric(cmat \%*\% MU) - RHS)/sqrt(diag(COV))
    switch(alternative, two.sided = {
        crit <- qmvt(p = 1 - alpha, tail = "both.tails", df = df, 
            corr = R)[["quantile"]]
    }, less = {
        crit <- qmvt(p = 1 - alpha, tail = "upper.tail", df = df, 
            corr = R)[["quantile"]]
    }, greater = {
        crit <- qmvt(p = 1 - alpha, tail = "lower.tail", df = df, 
            corr = R)[["quantile"]]
    })
    switch(power, global = {
        switch(alternative, two.sided = {
            beta <- pmvt(lower = rep(-crit, M), upper = rep(crit, 
                M), delta = ExpT, df = df, corr = R)
            whichHA <- which(ExpT != 0)
        }, less = {
            beta <- pmvt(lower = rep(crit, M), upper = rep(Inf, 
                M), delta = ExpT, df = df, corr = R)
            whichHA <- which(ExpT < 0)
        }, greater = {
            beta <- pmvt(lower = rep(-Inf, M), upper = rep(crit, 
                M), delta = ExpT, df = df, corr = R)
            whichHA <- which(ExpT > 0)
        })
    }, anypair = {
        switch(alternative, two.sided = {
            whichHA <- which(ExpT != 0)
            MHA <- length(whichHA)
            if (MHA < 1) {
                warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated.")
                beta <- 1 - alpha
            } else {
                beta <- pmvt(lower = rep(-crit, MHA), upper = rep(crit, 
                  MHA), delta = ExpT[whichHA], df = df, corr = R[whichHA, 
                  whichHA])
            }
        }, less = {
            whichHA <- which(ExpT < 0)
            MHA <- length(whichHA)
            if (MHA < 1) {
                warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated.")
                beta <- 1 - alpha
            } else {
                beta <- pmvt(lower = rep(crit, MHA), upper = rep(Inf, 
                  MHA), delta = ExpT[whichHA], df = df, corr = R[whichHA, 
                  whichHA])
            }
        }, greater = {
            whichHA <- which(ExpT > 0)
            MHA <- length(whichHA)
            if (MHA < 1) {
                warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated.")
                beta <- 1 - alpha
            } else {
                beta <- pmvt(lower = rep(-Inf, MHA), upper = rep(crit, 
                  MHA), delta = ExpT[whichHA], df = df, corr = R[whichHA, 
                  whichHA])
            }
        })
    }, allpairs = {
        switch(alternative, two.sided = {
            whichHA <- which(ExpT != 0)
            MHA <- length(whichHA)
            if (MHA < 1) {
                warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated.")
                beta <- 1 - alpha
            } else {
                nsim <- 1e+05
                RT <- rmvt(n = nsim, delta = ExpT[whichHA], df = df, 
                  sigma = as.matrix(R[whichHA, whichHA]), method = "svd")
                nreject <- sum(apply(RT, 1, function(x) {
                  min(abs(x))
                }) > abs(crit))
                beta <- 1 - (nreject/nsim)
            }
        }, less = {
            whichHA <- which(ExpT < 0)
            MHA <- length(whichHA)
            if (MHA < 1) {
                warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated.")
                beta <- 1 - alpha
            } else {
                beta <- 1 - pmvt(lower = rep(-Inf, MHA), upper = rep(crit, 
                  MHA), delta = ExpT[whichHA], df = df, corr = R[whichHA, 
                  whichHA])
            }
        }, greater = {
            whichHA <- which(ExpT > 0)
            MHA <- length(whichHA)
            if (MHA < 1) {
                warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated.")
                beta <- 1 - alpha
            } else {
                beta <- 1 - pmvt(lower = rep(crit, MHA), upper = rep(Inf, 
                  MHA), delta = ExpT[whichHA], df = df, corr = R[whichHA, 
                  whichHA])
            }
        })
    })
    pow <- round(1 - beta[[1]], 3)
    return(pow)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
