ESSgls <- function(model, IDcolumn=1){
  
  if(class(model)!="gls"){
    stop("Your model must be of class 'gls'!")
  }
  
  N <- attr(model$modelStruct$corStruct, "Dim")$M
  p <- model$dims$p
  
  ni <- attr(model$modelStruct$corStruct, "Dim")$len
  
  Vilist <- Wilist <- Viprodlist <- Wiprodlist <- list()
  
  for(i in 1:N){
    Vilist[[i]] <- getVarCov(model, individual=i)
    if(length(diag(Vilist[[i]])) > 1){
      Wilist[[i]] <- diag(diag(Vilist[[i]]))
    }else{
      Wilist[[i]] <- diag(Vilist[[i]])
    }
  }
  
  X <- matrix(model.matrix(formula(model), data=getData(model)), ncol=p)
  
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