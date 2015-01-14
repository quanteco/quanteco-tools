#' Monte Carlo Test for Moran's I
#' 
#' Performs a Monte Carlo test for Moran's I, given a vector (\strong{X}) and a neighborhood matrix (\strong{W})
#' 
#' @param X a vector of values
#' @param W a neighborhood matrix
#' @param N the number of permutations or randomizations; default is 999
#' @param test a character sting specifying the alternative hypothesis, must be one of "positve" (default), "negative", or "two-sided".
#' @param graph a logical statement specifying if a graph of the test should be returned. 
#' @param print.results a logical statement specifying if results should be printed (\code{TRUE}; default) or not (\code{FALSE})
#' @export
#' @details While calculating Moran's I is helpful, we need some way to determine if the value is higher (or lower) than expected by chance. There are a variety of ways to perform such a test. A fairly robust approach that doesn't require the assumption of normality is to perform a randomization or permutation test. Basically, we can rearrange the elements of x and calculate Moran's I for each rearrangement. Depending on the number of elements in x, we can either perform all possible rearrangements (i.e., a permutation test) or a large number of combinations (i.e., a randomization test). 
#' @seealso \code{\link{weightmatrix}}, \code{\link{morans.i}}
#' @examples
#' library(matrix)
#' set.seed(200)
#' n<-100
#' y<-rnorm(n)
#' D<-array(dim=rep(length(y),2))
#' diag(D)<-0
#' D[lower.tri(D,diag=F)]<-runif(length(D[lower.tri(D,diag=F)]),0,500)
#' library(Matrix)
#' D<-forceSymmetric(D,uplo="L")
#' W<-weightmatrix(D)
#' y<-sim.autocorrelated(y,0.1,W)
#' morans.test(y,W,graph=T)

morans.test<-function(X,W,N=999,test=c("positive","negative","two-sided"),graph=F,print.results=T){  
  writeLines("Moran's I test for Autocorrelation")
  require(pbapply)
  observed<-morans.I(X,W)
  options("pbapply.pb"="txt")
  if(length(X)<8){
    if(N>gamma(length(X)+1))writeLines(paste("Only",gamma(length(X)+1),"permutations were used due to small sample size"))
    require(combinat)
    store<-pbapply(array(unlist(permn(X)),dim=c(length(X),gamma(length(X)+1)))[,1:(ifelse(N<=gamma(length(X)+1),N,gamma(length(X)+1)))],2,function(x)morans.I(x,W))
  }
  else store<-pbapply(replicate(N,sample(X)),2,function(x)morans.I(x,W))
  if(length(test)>1){test=test[1]}
  if(test=="positive"){
    p.value<-(sum(ifelse(store>=observed,1,0))+1)/(length(store)+1)
    expected=(-1/(length(X)-1))  
  }
  else if(test=="negative"){
    p.value<-(sum(ifelse(store<=observed,1,0))+1)/(length(store)+1)
    expected=(-1/(length(X)-1))
  }
  else if(test=="two-sided"){
    observed<-abs(observed-mean(store))
    store<-abs(store-mean(store))
    expected=abs(-1/(length(X)-1))
    p.value<-(sum(ifelse(store>=observed,1,0))+1)/(length(store)+1)
  }
  else{stop("test must be at 'positive', 'negative' or 'two-sided'")}
  if(graph==T){
    require(ggplot2)
    tmp.dat<-data.frame(store=store,observed=observed)
    export.graph<-ggplot(tmp.dat,aes(x=store))+
      scale_y_sqrt()+geom_density()+
      geom_vline(aes(xintercept=observed),color="red",size=1)+
      xlab("Moran's I Coefficient") +
      ylab("Empirical Density")+theme_bw(base_family="serif")
    print(export.graph)
  }
  if(print.results==T){
    print(list(Observed=observed,Expected=expected,p.value=p.value))
  }
  z<-list(Observed=observed,Expected=expected,p.value=p.value,Values=store)
}