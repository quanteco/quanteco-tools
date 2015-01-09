#' Calculate Moran's I
#' 
#' Calculates Moran's I with a given vector (\strong{X}) and a neighborhood matrix (\strong{W})
#' 
#' @param X a vector of values
#' @param W a neighborhood matrix
#' @return Moran's I
#' @export 
#' @details Moran's I is a value that represents the degree of autocorrelation in a vector (\strong{x}) given a weight matrix (\strong{W}). It is similar to calculating the correlation between two variables and is typically restricted to be between -1 and 1; however, this depends on whether the weight matrix is weighted by row. Negative values of Moran's I represent negative autocorrelation and positive values represent positive autocorrelation. Although these values can be compared to the spatial autocorrelation of other vectors, there is little utility of the Moran's I value itself. 
#' @seealso \code{\link{weightmatrix}}, \code{\link{morans.test}}
#' @examples
#' ############
#' #Not run   #
#' ############
#' library(matrix)
#' n<-100
#' y<-rnorm(n)
#' D<-array(dim=rep(length(y),2))
#' diag(D)<-0
#' D[lower.tri(D,diag=F)]<-runif(length(D[lower.tri(D,diag=F)]),0,500)
#' library(Matrix)
#' D<-forceSymmetric(D,uplo="L")

morans.I<-function(X,W){
  n<-length(X)
  I<-diag(n)
  J<-rep(1,n)
  H<-I-J%*%t(J)/n
  morans<-(n/(t(J)%*%W%*%J))%*%((t(X)%*%H%*%W%*%H%*%X)/(t(X)%*%H%*%X))
  return(as.vector(morans))
}
