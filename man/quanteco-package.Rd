\name{quanteco-package}
\alias{quanteco-package}
\alias{quanteco}
\docType{package}
\title{
Quantitative Functions and Tests for Ecology
}
\description{
Functions and tests for spatioautocorrelation, model averaging, data visualization, and data analysis.
}
\details{
\tabular{ll}{
Package: \tab quanteco\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-12-19\cr
License: \tab GPL \cr
}
}
\author{
Gregory B. Anderson

Maintainer: Gregory B. Anderson <quanteco@vt.edu>
}
\references{
Dray, S., P. Legendre, and P. R. Peres-Neto. 2006. Spatial modeling: a comprehensive framework for principal coordinate analysis of neighborhood matrices (PCNM). Ecological Modeling 196:483-493.

Griffith, D. A., and P. R. Peres-Neto. 2006. Spatial modeling in ecology: the flexibility of eigenfunction spatial analyses. Ecology 87:2603-2613.
}
\keyword{ package }
\examples{
#Create Weighted Neighborhood Matrix
library(matrix)
set.seed(200)
n<-100
y<-rnorm(n)
D<-array(dim=rep(length(y),2))
diag(D)<-0
D[lower.tri(D,diag=F)]<-runif(length(D[lower.tri(D,diag=F)]),0,500)
library(Matrix)
D<-forceSymmetric(D,uplo="L")
W<-weighmatrix(D)
W

#Test for autocorrelation using Monte Carolo tests of Moran's I
y<-sim.autocorrelated(y,0.1,W)
morans.test(y,W,graph=T)
}
