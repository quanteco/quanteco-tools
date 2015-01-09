#' Simulate Autocorrelated Data
#' 
#' Simulates autocorrelated data given a vector (\strong{X}), \eqn{\rho}, and a neighborhood matrix (\strong{W})
#' 
#' @param X a vector of values
#' @param rho stength of correlation
#' @param W a neighbohood matrix
#' @export
#' 


sim.autocorrelated<-function(X,rho,W){
  n<-length(X)
  e<-eigen(W,only.values=TRUE)$values
  if (is.complex(e))
    feasible <- 1/(range(Re(e)))
  else feasible <- 1/(range(e))
  if (rho <= feasible[1] || rho >= feasible[2])
    stop(paste("Rho outside feasible range:", feasible))
  mat<-diag(n)-rho * W
  res<-chol2inv(chol(mat))
  return(res%*%X)
}