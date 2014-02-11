satterthwaite <- function(var, n, df=n-1){
  m <- (1/n) / sum(1/n)
  v <- sum(m * var)
  df <- v^2 / sum(m^2 * var^2/df)
  return(df)
}