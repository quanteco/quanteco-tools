#' Create a weight matrix
#' 
#' Creates a weighted neighborhood matrix from a distance matrix using different rule sets.
#' 
#'  @param D a distance matrix or n x n matrix of distances to be used
#'  @param method the method to be used to create the wegith matrix; must be one of "dbmem" (default), "linear", "concave-up", "concave-down", or "connectivity".
#'  @param thresh a threshold value to be used in which sites are no longer connected; the default is the maximum distance on a minimum spanning tree.
#'  @param alpha the power used for the "concave-down" method.
#'  @param beta the power used for the "concave-up" method.
#'  @param diag.mat a 0 or 1 to be used as the diagonal weight in the output matrix.
#'  @details There are a variety of rules that one can use to define a neighborhood matrix, especially with lattice data (i.e., data located in regions). However, often we deal with data within a spatial (or temporal) domain that does not explicitly fit into regions, but rather they are points in space, time, or from an  evolutionary spectrum. For these cases, we need some way of defining whether a site, time, or species is a neighbor or not with another and a way of expressing the difficulty of exchange between other sites, times, or species.  This is why we create a weighted neighborhood matrix.
#'  
#'  Using a distance matrix, we can use a variety of rules to calculate a weighted neighborhood matrix.  This function uses four different methods that describe the increase in difficultly or exchange between points or species.  These functions are the ones described by Drey et al. (2006) and include: distance-based Moran's Eigenvector Maps (dbmem), Linear Fuction, Concave-Down, Concave-Up, and Connectivity.
#'  
#'  For the first four methods we need to define both (1) a connectivity matrix (B) and (2) an edge weight matrix (A); for the latter we only need a connectivity matrix.   To create a connectivity matrix, we need a thresdhold in which to say sites are no longer connected. The weightmatrix function allows for the user to provide a value; however, there is a standard theshold used for dbmem.  This is usually difined as the maximum distance on a minimum span treee.   This span tree is simply a tree that connects all sites together with edges, and the minimum span tree is the tree with the lowest cost. With the connectivity matrix defined, an edge matrix can be calculated using the different rule sets.
#'  
#'  To compute the weighted neighborhood we just need to take the Hadamard product of B and A.  However, the diagonal will be ones representing a connection.  This is good if one is going to try to derive more complex relationships with other weight matrices (e.g., spatiotemporal designs), but we might want to have zeros on the diagonal instead (default of function).  This is controled by the diag.mat argument. 
#'  @return A weighted neighboorhod matrix. 
#'  @examples
#'  #Not run
#'  D<-as.matrix(dist(cbind(runif(100,0,100),runif(100,0,100))))
#'  #Distance-based Moran's Eigenvector Maps
#'  W<-weightmatrix(D,method='dbmem')
#'  W
#'  
#'  #Concave-Up
#'  W<-weightmatrix(D,method='concave-up',beta=2)
#'  W
#'  @export 
#'  
#'  
weightmatrix<-function (D, method = c("dbmem", "linear", "concave-up", "concave-down","connectivity"), thresh = NULL, alpha = NULL, beta = NULL, diag.mat = 0) {
  if (class(D) == "dist")
    D = as.matrix(D)
  require(vegan)
  n <- dim(D)[1]
  if (is.null(thresh))
    thresh = max(spantree(D)$dist)
  out = A = (B <- array(dim = dim(D)))
  if (length(method) > 1)
    method = "dbmem"
  B<-apply(D,2,function(x)ifelse(x>thresh,0,1))
  if(method=="dbmem") A<-apply(D,2,function(x)1 - ((x/(4 * thresh))^2))
  else if(method=="linear")A<-apply(D,2,function(x)1 - (x/max(D)))
  else if (method == "concave-down") {
    if (is.null(alpha)) stop("alpha must be specified for the method 'concave-down'")
    A<-apply(D,2,function(x)1 - (x/max(D))^alpha)
  }
  else if (method == "concave-up") {
    if (is.null(beta)) stop("beta must be specified for the method 'concave-up")
    A<-apply(D,2,function(x)1/x^beta)
  }
  else if (method == "connectviity")A<-B
  else stop("method must be 'dbmem', 'linear', 'concave-up' or 'concave-down'")
  out <- A * B
  if (diag.mat == 0) {
    diag(out) <- 0
  }
  return(out)
}