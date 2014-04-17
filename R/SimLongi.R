SimLongi <- function(data, response, group, time, id, covariates=NULL,
                     var=list("hom", "het", "hett"), cor=list("CS", "AR1", "UN"),
                     contrasts=NULL, type="Dunnett", base=1, direction="gpt",
                     alternative="two.sided", level=0.95, df="ess"){

  if(min(var %in% list("hom", "het", "hett", "hetg"))==0){
    stop("Your heteroscedasticity patterns must be among the following:
         'hom', 'het', 'hett', 'hetg'.")
  }
  
  if(min(cor %in% list("CS", "AR1", "CAR1", "AR2", "MA1", "MA2", "ARMA11", "UN"))==0){
    stop("Your correlation structures must be among the following:
         'CS', 'AR1', 'AR2', 'MA1', 'MA2', 'ARMA11', 'CAR1', 'UN'.")
  }
  
  df <- match.arg(df, c("ess", "essdf", "adj", "pb", "satt", "kr", "con", "naive", "res", "normal"))
  
  if(is.null(contrasts)==T){

    direction <- match.arg(direction, c("gpt", "tpg", "both"))
      
    if(direction=="both" & length(type)==1){
      type <- rep(type, 2)
    }
    
    if(direction=="both" & length(base)==1){
      base <- rep(base, 2)
    }
    
    typeset <- c("Dunnett", "Tukey", "Sequen", "AVE", "Changepoint", "Williams",
                 "Marcus", "McDermott", "UmbrellaWilliams", "GrandMean")
    
    if(direction!="both"){
      type <- match.arg(type, typeset)
    }else{
      type1 <- type[1]
      type2 <- type[2]
      type1 <- match.arg(type1, typeset)
      type2 <- match.arg(type2, typeset)
    }
    
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
  
  ########## Contrast Matrices ##########
  
  if(is.null(contrasts)==T){
    
    if(direction=="gpt"){
      mat1 <- contrMat(ngroup, type[1], base=base)
      cmat <- diag(times) %x% mat1
      rownames(cmat) <- paste(rep(levels(dat[, "time"]), each=dim(mat1)[1]), rownames(mat1), sep=":")
    }
    
    if(direction=="tpg"){
      mat2 <- contrMat(ntime, type[1], base=base)
      cmat <- mat2 %x% diag(groups)
      rownames(cmat) <- paste(rep(levels(dat[, "group"]), dim(mat2)[1]), rep(rownames(mat2), each=nlevels(dat[, "group"])), sep=":")
    }
    
    if(direction=="both"){
      mat1 <- contrMat(ngroup, type[1], base=base[1])
      cmat1 <- diag(times) %x% mat1
      rownames(cmat1) <- paste(rep(levels(dat[, "time"]), each=dim(mat1)[1]), rownames(mat1), sep=":")
      mat2 <- contrMat(ntime, type[2], base=base[2])
      cmat2 <- mat2 %x% diag(groups)
      rownames(cmat2) <- paste(rep(levels(dat[, "group"]), dim(mat2)[1]), rep(rownames(mat2), each=nlevels(dat[, "group"])), sep=":")
      cmat <- rbind(cmat1, cmat2)
    }
    
  }else{
    
    cmat <- contrasts
    
  }
  
  ########## Covariance Model Selection ##########
  
  mlist <- list()
  count <- 1
  
  if("hom" %in% var){
    if("CS" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            correlation=corCompSymm(form=~1|id))
      count <- count + 1
    }
    if("AR1" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            correlation=corAR1(form=~1|id))
      count <- count + 1
    }
    if("AR2" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            correlation=corARMA(form=~1|id, p=2))
      count <- count + 1
    }
    if("MA1" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            correlation=corARMA(form=~1|id, q=1))
      count <- count + 1
    }
    if("MA2" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            correlation=corARMA(form=~1|id, q=2))
      count <- count + 1
    }
    if("ARMA11" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            correlation=corARMA(form=~1|id, p=1, q=1))
      count <- count + 1
    }
    if("CAR1" %in% cor){
      dat$timeC <- as.numeric(dat[, "time"])
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            correlation=corCAR1(form=~timeC|id))
      count <- count + 1
    }
    if("UN" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            correlation=corSymm(form=~1|id))
      count <- count + 1
    }
  }
  
  if("het" %in% var){
    if("CS" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|tg),
                            correlation=corCompSymm(form=~1|id))
      count <- count + 1
    }
    if("AR1" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|tg),
                            correlation=corAR1(form=~1|id))
      count <- count + 1
    }
    if("AR2" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|tg),
                            correlation=corARMA(form=~1|id, p=2))
      count <- count + 1
    }
    if("MA1" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|tg),
                            correlation=corARMA(form=~1|id, q=1))
      count <- count + 1
    }
    if("MA2" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|tg),
                            correlation=corARMA(form=~1|id, q=2))
      count <- count + 1
    }
    if("ARMA11" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|tg),
                            correlation=corARMA(form=~1|id, p=1, q=1))
      count <- count + 1
    }
    if("CAR1" %in% cor){
      dat$timeC <- as.numeric(dat[, "time"])
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|tg),
                            correlation=corCAR1(form=~timeC|id))
      count <- count + 1
    }
    if("UN" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|tg),
                            correlation=corSymm(form=~1|id))
      count <- count + 1
    }
  }
  
  if("hett" %in% var){
    if("CS" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|time),
                            correlation=corCompSymm(form=~1|id))
      count <- count + 1
    }
    if("AR1" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|time),
                            correlation=corAR1(form=~1|id))
      count <- count + 1
    }    
    if("AR2" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|time),
                            correlation=corARMA(form=~1|id, p=2))
      count <- count + 1
    }
    if("MA1" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|time),
                            correlation=corARMA(form=~1|id, q=1))
      count <- count + 1
    }
    if("MA2" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|time),
                            correlation=corARMA(form=~1|id, q=2))
      count <- count + 1
    }
    if("ARMA11" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|time),
                            correlation=corARMA(form=~1|id, p=1, q=1))
      count <- count + 1
    }
    if("CAR1" %in% cor){
      dat$timeC <- as.numeric(dat[, "time"])
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|time),
                            correlation=corCAR1(form=~timeC|id))
      count <- count + 1
    } 
    if("UN" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|time),
                            correlation=corSymm(form=~1|id))
      count <- count + 1
    }
  }
  
  if("hetg" %in% var){
    if("CS" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|group),
                            correlation=corCompSymm(form=~1|id))
      count <- count + 1
    }
    if("AR1" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|group),
                            correlation=corAR1(form=~1|id))
      count <- count + 1
    }
    if("AR2" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|group),
                            correlation=corARMA(form=~1|id, p=2))
      count <- count + 1
    }
    if("MA1" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|group),
                            correlation=corARMA(form=~1|id, q=1))
      count <- count + 1
    }
    if("MA2" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|group),
                            correlation=corARMA(form=~1|id, q=2))
      count <- count + 1
    }
    if("ARMA11" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|group),
                            correlation=corARMA(form=~1|id, p=1, q=1))
      count <- count + 1
    }
    if("CAR1" %in% cor){
      dat$timeC <- as.numeric(dat[, "time"])
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|group),
                            correlation=corCAR1(form=~timeC|id))
      count <- count + 1
    }
    if("UN" %in% cor){
      mlist[[count]] <- gls(response ~ tg - 1, dat,
                            weights = varIdent(form=~1|group),
                            correlation=corSymm(form=~1|id))
      count <- count + 1
    }
  }
  
  mstab <- mod.sel(mlist)
  
  bestweight <- mstab[1, "weight"]
  
  best <- as.numeric(attr(mstab, "order")[1])
  
  ########## Degrees of Freedom ##########
  
  if(df=="ess"){
    essmod <- gls(response ~ gt - 1, dat, correlation=corAR1(form=~1|id))
    phi <- cov2cor(vcov(essmod))[1, 2]
    def <- floor(ids * ((times - (times - 2) * phi) / (1 + phi)))
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
  
  test <- glht(mlist[[best]], cmat, alternative=alternative, df=def)
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
    
  ########## Output ##########
  
  out <- list()
  
  out[["Results"]] <- tab
  out[["CovStat"]] <- covm
  out[["CritValue"]] <- crit
  out[["Alternative"]] <- alternative
  out[["ConfLevel"]] <- level
  out[["DFMethod"]] <- df
  out[["DF"]] <- def
  out[["ContMat"]] <- cmat
  out[["BestMod"]] <- mlist[[best]]$call
  out[["ModSelTab"]] <- mstab
  out[["AWBest"]] <- bestweight
  out[["CovBest"]] <- vcov(mlist[[best]])
  out[["Model"]] <- mlist[[best]]
  class(out) <- "silo"
  
  return(out)
  
}