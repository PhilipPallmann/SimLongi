ESSlme <- function(model, IDcolumn){
  
  if(class(model)!="lme"){
    stop("Your model must be of class 'lme'!")
  }
  
  N <- mix2$dims$ngrps[1]
  p <- length(model$coefficients$fixed)
  
  ni <- table(model$groups)
  
  Vilist <- Wilist <- Viprodlist <- Wiprodlist <- list()
  
  for(i in 1:N){
    Vilist[[i]] <- getVarCov(model, individuals=i, type="marginal")[[1]]
    if(length(diag(Vilist[[i]])) > 1){
      Wilist[[i]] <- diag(diag(Vilist[[i]]))
    }else{
      Wilist[[i]] <- diag(Vilist[[i]])
    }
  }
  
  X <- matrix((model.matrix(model, getData(model))), ncol=p)
  
  Xilist <- split(X, getData(model)[, IDcolumn])
  makemat <- function(x) matrix(x, ncol=p)
  Xilist <- lapply(Xilist, makemat)
  
  for(i in 1:N){
    Viprodlist[[i]] <- t(Xilist[[i]]) %*% solve(Vilist[[i]]) %*% Xilist[[i]]
    Wiprodlist[[i]] <- t(Xilist[[i]]) %*% solve(Wilist[[i]]) %*% Xilist[[i]]
  }
  
  Vprod <- Reduce("+", Viprodlist)
  Wprod <- Reduce("+", Wiprodlist)
  
  Varhat <- solve(Vprod)
  Vartilde <- solve(Wprod)
  
  w <- Vartilde/Varhat
  
  Ntilde <- diag(w) * sum(ni)
  
  return(Ntilde)
  
}