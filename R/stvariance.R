#' Calculate Spatiotemporal Variation
#' 
#' Calculates and partitions spatiotemporal variation given a n (time) x p (site) matrix
#' 
#' @param Y a time by site matrix of abudances
#' @param trans a string determining the transformation used on Y; currently implemented options are: 'chord' (default) and 'hellinger'
#' @param print.results a logical result to print results (default) or be silent
#' @export
#' @return a list including (1) total spatiotemporal variance, (2) temporal contributions to variance, and (3) spatial contributions to variance
#' @details This function is used to partition spatiotemporal variance in abundance into: (1) total spatiotemporal variance, (2) temporal contributions to variance, and (3) spatial contributions to variance.
#' @examples
#' Create Data
#' Y<-array(c(0,2,1,5,0,0,0,3,1,2,4,2,0,0,2,2,0,2,2,0,0,2,3,1,5,8,4,2,1,2,2,2,4,3,1,2,4,4,3,2),dim=c(8,5))
#' rownames(Y)<-c(2006:2012,2014)
#' colnames(Y)<-1:5

stvariance<-function(Y,trans=c("chord","hellinger"),print.results=T,...){
  if(length(trans)>1){trans="chord"}
  if(trans=="chord"|trans=="c"){
    X<-Y/sqrt(rowSums(Y^2,...))
  }
  else if(trans=="hellinger"|trans=="h"){
    X<-sqrt(Y/rowSums(Y,...))
  }
  else{stop("trans must be 'chord' or 'hellinger'")}
  VarY<-sum(colSums((X-rowMeans(X,...))^2,...))/(nrow(X)-1)
  TCV<-matrix(rowSums(((X-rowMeans(X,...))^2)/sum(colSums((X-rowMeans(X,...))^2,...),...),...),nrow=1)
  colnames(TCV)<-rownames(Y)
  SCV<-matrix(colSums(((X-rowMeans(X,...))^2)/sum(colSums((X-rowMeans(X,...))^2,...),...),...),nrow=1)
  colnames(SCV)<-colnames(Y)
  if(print.results==T){
    writeLines(paste("Total Spatiotemporal Variance:",round(VarY,2)))
    writeLines("")
    writeLines("Temporal Contributions to Variance:")
    print(TCV)
    writeLines("")
    writeLines("Spatial Contributions to Variance:")
    print(SCV)
  }
  z<-list('VarY'=VarY,'TCV'=TCV,'SCV'=SCV)
}