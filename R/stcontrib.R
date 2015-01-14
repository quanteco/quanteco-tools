#' Tests for Spatiotemporal Contribution
#' 
#' Tests for significant temporal (TCV) and spatial (SCV) contributions to variance
#' 
#' @param Y a time by site matrix of abudances
#' @param trans a string determining the transformation used on Y; currently implemented options are: 'chord' (default) and 'hellinger'
#' @param stat a string determining whether to test for temporal (TCV; default) or spatial contributions (SCV) to variance
#' @param N the number of randomizations to perform
#' @export
#' @return a vector of emperical p-values for either temporal (TCV) or spatial (SCV) contributions to variance
#' @details This function is used to test for signficant temporal (TCV) or spatial (SCV) contributions to spatiotemporal variance in abundance or counts of a species.
#' @seealso \code{\link{stvariance}}


stcontrib.test<-function(Y,trans=c("chord","hellinger"),stat=c('TCV','SCV'),N=999,...){
  if(length(trans)>1) trans="chord" 
  if(length(stat)>1) stat='TCV'
  obs<-stvariance(Y,trans,print.results=F,...)
  if(stat=='TCV'){
    tmp<-matrix((rowSums(apply(apply(replicate(N,apply(Y,2,
                                                       sample)),3,function(x) stvariance(x,trans=trans,print.results=F,...)$TCV),2,
                               function(x) x>obs$TCV))+1)/(N+1),nrow=1)
    colnames(tmp)<-rownames(Y)
    
  }
  else if(stat=='SCV'){
    tmp<-matrix((rowSums(apply(apply(replicate(N,t(apply(Y,1,
                                                         sample))),3,function(x) stvariance(x,trans=trans,print.results=F,...)$SCV),2,
                               function(x) x>obs$SCV))+1)/(N+1),nrow=1)
    colnames(tmp)<-colnames(Y)
  }
  else{stop("stat must be 'TCV' or 'SCV'")}
  return(tmp)
}