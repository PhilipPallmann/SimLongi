PowApprox <- function(mu, cov, n, df, cmat, rhs=0, alpha=0.05, alternative=c("two.sided", "less", "greater"),
                      power=c("global", "anypair", "allpairs")){
  
  MU <- matrix(mu, ncol=1)
  COV <- cmat %*% (cov/n) %*% t(cmat)
  R <- cov2cor(COV)
  M <- nrow(cmat)
  RHS <- rep(rhs, M)
  ExpT <- (as.numeric(cmat %*% MU) - RHS)/sqrt(diag(COV))
  
  switch(alternative,
         two.sided={
           crit <- qmvt(p=1-alpha, tail="both.tails", df=df, corr=R)[["quantile"]]
         },
         less={
           crit <- qmvt(p=1-alpha, tail="upper.tail", df=df, corr=R)[["quantile"]]
         },
         greater={
           crit <- qmvt(p=1-alpha, tail="lower.tail", df=df, corr=R)[["quantile"]]
         })
  
  switch(power,
         global={
           switch(alternative,
                  two.sided={
                    beta <- pmvt(lower=rep(-crit, M), upper=rep(crit, M), delta=ExpT, df=df, corr=R)
                    whichHA <- which(ExpT != 0)
                  },
                  less={
                    beta <- pmvt(lower=rep(crit, M), upper=rep(Inf, M), delta=ExpT, df=df, corr=R)
                    whichHA <- which(ExpT < 0)
                  },
                  greater={
                    beta <- pmvt(lower=rep(-Inf, M), upper=rep(crit, M), delta=ExpT, df=df, corr=R)
                    whichHA <- which(ExpT > 0)
                  })
         },
         anypair={
           switch(alternative,
                  two.sided={
                    whichHA <- which(ExpT != 0)
                    MHA <- length(whichHA)
                    if(MHA < 1){
                      warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated.")
                      beta <- 1 - alpha
                    }else{
                      beta <- pmvt(lower=rep(-crit, MHA), upper=rep(crit, MHA), delta=ExpT[whichHA], df=df, corr=R[whichHA, whichHA])
                    }
                  },
                  less={
                    whichHA <- which(ExpT < 0)
                    MHA <- length(whichHA)
                    if(MHA < 1){
                      warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated.")
                      beta <- 1 - alpha
                    }else{
                      beta <- pmvt(lower=rep(crit, MHA), upper=rep(Inf, MHA), delta=ExpT[whichHA], df=df, corr=R[whichHA, whichHA])
                    }
                  },
                  greater={
                    whichHA <- which(ExpT > 0)
                    MHA <- length(whichHA)
                    if(MHA < 1){
                      warning("All contrasts are under their corresponding null hypothesis, anypair power can not be calculated.")
                      beta <- 1 - alpha
                    }else{
                      beta <- pmvt(lower=rep(-Inf, MHA), upper=rep(crit, MHA), delta=ExpT[whichHA], df=df, corr=R[whichHA, whichHA])
                    }
                  })
         },
         allpairs={
           switch(alternative,
                  two.sided={
                    whichHA <- which(ExpT != 0)
                    MHA <- length(whichHA)
                    if(MHA < 1){
                      warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated.")
                      beta <- 1 - alpha
                    }else{
                      #beta <- pmvt(lower=rep(-crit, MHA), upper=rep(crit, MHA), delta=ExpT[whichHA], df=df, corr=R[whichHA, whichHA])
                      # geht nur so einfach, wenn GENAU eine Hypothese unter H_A ist!!!
                      nsim <- 100000
                      RT <- rmvt(n=nsim, delta=ExpT[whichHA], df=df, sigma=as.matrix(R[whichHA, whichHA]), method="svd")
                      nreject <- sum(apply(RT, 1, function(x){min(abs(x))}) > abs(crit))
                      beta <- 1 - (nreject/nsim)
                      #simerror <- sqrt(beta * (1 - beta)/nsim)
                      #attr(beta, which = "simerror") <- simerror
                    }
                  },
                  less={
                    whichHA <- which(ExpT < 0)
                    MHA <- length(whichHA)
                    if(MHA < 1){
                      warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated.")
                      beta <- 1 - alpha
                    }else{
                      beta <- 1 - pmvt(lower=rep(-Inf, MHA), upper=rep(crit, MHA), delta=ExpT[whichHA], df=df, corr=R[whichHA, whichHA])
                    }
                  },
                  greater={
                    whichHA <- which(ExpT > 0)
                    MHA <- length(whichHA)
                    if(MHA < 1){
                      warning("All contrasts are under the corresponding null hypotheses, all pair power can not be calculated.")
                      beta <- 1 - alpha
                    }else{
                      beta <- 1 - pmvt(lower=rep(crit, MHA), upper=rep(Inf, MHA), delta=ExpT[whichHA], df=df, corr=R[whichHA, whichHA])
                    }
                  })
         })
  
  pow <- round(1 - beta[[1]], 3)
  #HAtrue <- numeric(M)
  #HAtrue[whichHA] <- 1
  #settings <- data.frame(cmat, expContrast=ExpL, rhs = RHS, 
  #                       ExpTstat = ExpTeststat, underHA = HAtrue)
  #out <- list(power = pow, mu = mu, n = n, conexp = settings, 
  #            crit = crit, alternative = alternative, ptype = ptype, 
  #            alpha = alpha)
  #return(out)
  return(pow)
  
}