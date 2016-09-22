\name{BMEA-package}
\alias{BMEA-package}
\alias{BMEA}
\docType{package}
\title{Bayesian Modelling for Exon Arrays}
\description{
  Designed to operate on the aroma.affymetrix framework. Builds on the probe-level modelling approaches for 3' Affymetrix Arrays by incorporating a term for the proportion of exons present in a cellType. Uses a Bayesian (MCMC) approach.
}
\details{
\tabular{ll}{
Package: \tab BMEA\cr
Type: \tab Package\cr
Version: \tab 0.0.1\cr
Date: \tab 2012-01-04\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
The BMEA package is designed for analysis of Exon Array data, and uses a
computationally intensive MCMC approach with batches of genes being analysed in parallel. 

Users are also advised to be aware that analysis of an entire dataset can take more than 24hrs, depending on the computing framework & size of the dataset.
}
\author{
Steve Pederson

Maintainer: Steve Pederson <stephen.pederson@adelaide.edu.au>

}
\keyword{ BMEA }
\keyword{ Bayesian }
\keyword{ Microarray }
