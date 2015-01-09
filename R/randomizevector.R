#' Randomize vector for Monte Carlo tests.
#' 
#' Creates a list of N randomizations of a vector (X)
#' 
#' @param N the number of randomizations to be made
#' @param X the vector to randomize
#' @return A list of N vectors that represent randomizations of X

randomize_vector<-function(X,N){
  lst<-list()
  for(i in 1:N){
    lst[[i]]<-sample(X,length(X))
  }
  return(lst)
}