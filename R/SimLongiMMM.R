SimLongiMMM <- function(data, response, group, time, id, covariates=NULL,
                        contrasts=NULL, type="Dunnett", base=1,
                        alternative="two.sided", level=0.95, df="essdf"){
  
  df <- match.arg(df, c("ess", "essdf", "adj", "pb", "satt", "kr", "con", "naive", "res", "normal"))
  
  if(is.null(contrasts)==T){
    
    typeset <- c("Dunnett", "Tukey", "Sequen", "AVE", "Changepoint", "Williams",
                 "Marcus", "McDermott", "UmbrellaWilliams", "GrandMean")
    
    type <- match.arg(type, typeset)
    
  }else{
    
    if(dim(contrasts)[2]!=nlevels(as.factor(data[, group])) * nlevels(as.factor(data[, time]))){
      stop("The number of columns in your contrast matrix must equal the
           number of groups multiplied by the number of time points!")
    }
    
    if(any(rowSums(contrasts)!=0)){
      stop("The sum of elements in each row of your contrast matrix must be zero!")
    }
    
    }
  
  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
  
  dat <- data.frame(response=data[, response], group=as.factor(data[, group]),
                    time=as.factor(data[, time]), id=as.factor(data[, id]))
  
  dat$tg <- dat[, "time"]:dat[, "group"]
  dat$gt <- dat[, "group"]:dat[, "time"]
  
  groups <- nlevels(dat[, "group"])
  times <- nlevels(dat[, "time"])
  ids <- nlevels(dat[, "id"])
  tgs <- nlevels(dat[, "tg"])
  
  ngroup <- summary(dat[, "group"])
  ntime <- summary(dat[, "time"])
  nid <- summary(dat[, "id"])
  ntg <- summary(dat[, "tg"])
  
  ########## Marginal Models ##########
  
  dat$time <- as.factor(dat$time)
  
  datlist <- split(dat, dat$time)
  
  modlist <- list()
  
  for(z in 1:length(datlist)){
    modlist[[z]] <- lm(response ~ group, datlist[[z]])
  }
  
  names(modlist) <- levels(dat$time)
  
  ########## Contrast Matrices ##########
  
  if(is.null(contrasts)==T){
    cmat <- contrMat(ngroup, type, base=base)
  }else{
    cmat <- contrasts
  }
  
  ########## Degrees of Freedom ##########
  
  if(df=="ess"){
    essmod <- gls(response ~ gt - 1, dat, correlation=corAR1(form=~1|id))
    phi <- cov2cor(vcov(essmod))[1, 2]
    def <- floor(ids * ((times - (times - 2) * phi) / (1 + phi))) # - times * groups)
  }
  
  if(df=="essdf"){
    essmod <- gls(response ~ gt - 1, dat, correlation=corAR1(form=~1|id))
    phi <- cov2cor(vcov(essmod))[1, 2]
    def <- floor(ids * ((times - (times - 2) * phi) / (1 + phi)) - times * groups)
  }
  
  if(df=="adj"){
    def <- floor(ids + (times^2 * groups^2) / ids - tgs)
  }
  
  if(df=="pb"){
    def <- dim(dat)[1] - ids - tgs + 1
  }
  
  if(df=="satt"){
    vvv <- as.vector(tapply(dat[, "response"], dat[, "tg"], var))
    def <- floor(satterthwaite(vvv, ntg))
  }
  
  if(df=="kr"){
    keromod <- lmer(response ~ tg + (1|id), dat)
    keromod0 <- lmer(response ~ 1 + (1|id), dat)
    def <- floor(KRmodcomp(keromod, keromod0)$stats$ddf)
  }
  
  if(df=="con"){
    conmod <- lmer(response ~ tg - 1 + (1|id), dat)
    def <- dim(dat)[1] - rankMatrix(cbind(getME(conmod, "X"), as.matrix(getME(conmod, "Z"))))[1]
  }
  
  if(df=="naive"){
    def <- ids - (groups * times)
  }
  
  if(df=="res"){
    def <- dim(dat)[1] - (groups * times)
  }
  
  if(df=="normal"){
    def <- 0
  }
  
  if(df!="normal"){
    def <- max(def, 2)
  }
  
  ########## Tests and SCIs ##########
  
  # geht leider nicht:
  # test <- glht(mmm(modlist), mlf(mcp(group="Dunnett")))
  
  ### deshalb ganz primitiv:
  if(nlevels(dat$time)==2){
    test <- glht(mmm(T1=modlist[[1]], T2=modlist[[2]]), mlf(mcp(group=cmat)),
                 alternative=alternative, df=def)
  }
  if(nlevels(dat$time)==3){
    test <- glht(mmm(T1=modlist[[1]], T2=modlist[[2]], T3=modlist[[3]]), mlf(mcp(group=cmat)),
                 alternative=alternative, df=def)
  }
  if(nlevels(dat$time)==4){
    test <- glht(mmm(T1=modlist[[1]], T2=modlist[[2]], T3=modlist[[3]], T4=modlist[[4]]), mlf(mcp(group=cmat)),
                 alternative=alternative, df=def)
  }
  if(nlevels(dat$time)==5){
    test <- glht(mmm(T1=modlist[[1]], T2=modlist[[2]], T3=modlist[[3]], T4=modlist[[4]],
                     T5=modlist[[5]]), mlf(mcp(group=cmat)), alternative=alternative, df=def)
  }
  if(nlevels(dat$time)==6){
    test <- glht(mmm(T1=modlist[[1]], T2=modlist[[2]], T3=modlist[[3]], T4=modlist[[4]],
                     T5=modlist[[5]], T6=modlist[[6]]), mlf(mcp(group=cmat)), alternative=alternative, df=def)
  }
  if(nlevels(dat$time)==7){
    test <- glht(mmm(T1=modlist[[1]], T2=modlist[[2]], T3=modlist[[3]], T4=modlist[[4]],
                     T5=modlist[[5]], T6=modlist[[6]], T7=modlist[[7]]), mlf(mcp(group=cmat)),
                 alternative=alternative, df=def)
  }
  if(nlevels(dat$time)==8){
    test <- glht(mmm(T1=modlist[[1]], T2=modlist[[2]], T3=modlist[[3]], T4=modlist[[4]],
                     T5=modlist[[5]], T6=modlist[[6]], T7=modlist[[7]], T8=modlist[[8]]),
                 mlf(mcp(group=cmat)), alternative=alternative, df=def)
  }
  
  tests <- summary(test)
  testci <- confint(test, level=level)
  
  est <- coef(test)
  covm <- vcov(test)
  corm <- cov2cor(vcov(test))
  se <- tests$test$sigma[1:dim(covm)[1]]
  stat <- tests$test$tstat[1:dim(covm)[1]]
  padj <- tests$test$pvalues[1:dim(covm)[1]]
  low <- testci$confint[, "lwr"]
  up <- testci$confint[, "upr"]
  crit <- attr(testci$confint, "calpha")
  
  tab <- data.frame(Estimate=est, SE=se, Lower=low, Upper=up, Statistic=stat, P=padj)
  tab <- round(tab, 4)
  
  comprownames <- matrix(unlist(strsplit(rownames(tab), ": ")), byrow=T, ncol=2)[, 2]
  rownames(tab) <- paste(rep(levels(dat$time), each=nrow(cmat)), ":", comprownames, sep="")
  
  ########## Output ##########
  
  out <- list()
  
  out[["Results"]] <- tab
  out[["CovStat"]] <- covm
  out[["CritValue"]] <- crit
  out[["Alternative"]] <- alternative
  out[["ConfLevel"]] <- level
  out[["DFMethod"]] <- df
  out[["DF"]] <- def
  out[["Corr"]] <- corm
  out[["ContMat"]] <- cmat
  class(out) <- "silo"
  
  return(out)
  
}