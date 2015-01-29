#' Monte-Carlo Test for Positive or Negative Contributions to Moran's I
#' 
#' Performs a randomization test for either postive or negative autocorrelation contributions to Moran's I using the method described by Dray (2011)
#' #' 
#' @param X a vector of values
#' @param vectors an optional matrix of eigenvectors if no neighborhood matrix is supplied
#' @param values an optional vector Moran's I coefficients corresponding to each eigenvector supplied
#' @param W a neighborhood matrix
#' @param stat a character sting specifying the whether "positve" (default) or "negative" contributions should be calculated.
#' @param N the number of permutations (default=999)
#' @export
#' @details When calculating Moran's I for a vector, the value observed represents a combination of both negative and positive autocorrelation forces. This can be problematic when negative autocorrelation forces are as stong as the positive autcorrelation forces. When this is the case, Moran's I might be low despite the presence of both positive and autocorrelation. Dray (2000) demonstrates that the positive and negative contributions to Moran's I can be calculated seperately. Using the formulation, we can test for signficant contribuitons to either positive or negative autocorrelation. 

morancontrib.test<-function(X,vectors=NULL,values=NULL,W=NULL,stat=c("positive","negative"),N=999){
  if(length(stat)>1)stat="positive"
  obs<-morancontrib(X=X,vectors=vectors,values=values,W=W,stat=stat)
  if(stat=="positive"|stat=="pos"|stat=="p"){
    return((rowSums(pbapply(replicate(N,sample(X)),2,function(x)morancontrib(x,vectors=vectors,values=values,W=W,stat=stat))>=as.vector(obs))+1)/(N+1))
  }
  else if(stat=="negative"|stat=="neg"|stat=="n"){
    return((rowSums(pbapply(replicate(N,sample(X)),2,function(x)morancontrib(x,vectors=vectors,values=values,W=W,stat=stat))<=as.vector(obs))+1)/(N+1))
  }
  else(stop("stat must be 'positive' or 'negative'"))
}