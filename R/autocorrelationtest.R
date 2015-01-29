#' Monte Carlo Test for Autocorrelation
#' 
#' Performs a Monte Carlo test for autocorrelation using the method described by Dray (2011)
#' #' 
#' @param X a vector of values
#' @param vectors an optional matrix of eigenvectors if no neighborhood matrix is supplied
#' @param values an optional vector Moran's I coefficients corresponding to each eigenvector supplied
#' @param W a neighborhood matrix
#' @param test a character sting specifying the alternative hypothesis, must be one of "positve" (default), "negative", or "two-sided".
#' @param N the number of permutations or randomizations; default is 999
#' @export
#' @details When calculating Moran's I for a vector, the value observed represents a combination of both negative and positive autocorrelation forces. This can be problematic when negative autocorrelation forces are as stong as the positive autcorrelation forces. When this is the case, Moran's I might be low despite the presence of both positive and autocorrelation. Dray (2000) demonstrates that the positive and negative contributions to Moran's I can be calculated seperately, and this calculation can be used to have a more robust test for positive and/or negative autocorrelation. 
#' @seealso \code{\link{weightmatrix}}, \code{\link{morans.i}}
autocorrelation.test<-function(X,vectors=NULL,values=NULL,W=NULL,test=c("positive","negative","two-sided"),N=999){
  if(length(test)>1)test="positive"
  if(test=="positive"|test=="pos"|test=="p"){
    obs<-sum(morancontrib(X=X,vectors=vectors,values=values,W=W,stat="positive"))
    return((sum(pbapply(replicate(N,sample(X)),2,function(x)sum(morancontrib(x,vectors=vectors,values=values,W=W,stat="positive")))>=obs)+1)/(N+1))
  }
  else if(test=="negative"|test=="neg"|test=="n"){
    obs<-sum(morancontrib(X=X,vectors=vectors,values=values,stat="negative"))
    return((sum(pbapply(replicate(N,sample(X)),2,function(x)sum(morancontrib(x,vectors=vectors,values=values,W=W,stat="negative")))<=obs)+1)/(N+1))
  }
  else if(test=="two-sided"|test=="two"|test=="t"){
    obs<-sum(morancontrib(X=X,vectors=vectors,values=values,stat="positive"))
    p.pos<-(sum(pbapply(replicate(N,sample(X)),2,function(x)sum(morancontrib(x,vectors=vectors,values=values,W=W,stat="positive")))>=obs)+1)/(N+1)
    obs<-sum(morancontrib(X=X,vectors=vectors,values=values,stat="negative"))
    p.neg<-(sum(pbapply(replicate(N,sample(X)),2,function(x)sum(morancontrib(x,vectors=vectors,values=values,W=W,stat="negative")))<=obs)+1)/(N+1)
    return(2*min(p.pos,p.neg))
  }
  else(stop("test must be 'positive','negative', or 'two-sided'"))
}