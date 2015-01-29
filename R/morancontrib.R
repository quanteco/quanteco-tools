#' Calculates Positive or Negative Contributions to Moran's I
#' 
#' Calculates either postive or negative autocorrelation contributions to Moran's I using the method described by Dray (2011)
#' 
#' @param X a vector of values
#' @param vectors an optional matrix of eigenvectors if no neighborhood matrix is supplied
#' @param values an optional vector Moran's I coefficients corresponding to each eigenvector supplied
#' @param W a neighborhood matrix
#' @param stat a character sting specifying the whether "positve" (default) or "negative" contributions should be calculated.
#' @export
#' @details When calculating Moran's I for a vector, the value observed represents a combination of both negative and positive autocorrelation forces. This can be problematic when negative autocorrelation forces are as stong as the positive autcorrelation forces. When this is the case, Moran's I might be low despite the presence of both positive and autocorrelation. Dray (2000) demonstrates that the positive and negative contributions to Moran's I can be calculated seperately. 
morancontrib<-function(X,vectors=NULL,values=NULL,W=NULL,stat=c("positive","negative")){
  if(is.null(W)==F){
    vectors=mem(W)$vectors
    values=mem(W)$values
  }
  else{
    if (is.null(vectors) | is.null(values)) stop("Either vectors and values or W needs to be non-NULL")
  }
  Z<-scale(X)
  ER.MC= -(1/(length(Z)-1))
  tmp<-rep(0,length(values))
  if(length(stat)>1)stat="positive"
  if(stat=="positive"|stat=="pos"|stat=="p"){
    tmp[values>ER.MC]<-values[values>ER.MC]*((t(Z)%*%vectors[,values>ER.MC]*t(Z)%*%vectors[,values>ER.MC])/length(Z)^2)
    return(tmp)
  }
  else if(stat=="negative"|stat=="neg"|stat=="n"){
    tmp[values<ER.MC]<-values[values<ER.MC]*((t(Z)%*%vectors[,values<ER.MC]*t(Z)%*%vectors[,values<ER.MC])/length(Z)^2)
    return(tmp)
  }
  else(stop("stat must be 'positive' or 'negative'"))
}