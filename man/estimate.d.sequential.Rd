% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/d_estimation_fnxs.R
\name{estimate.d.sequential}
\alias{estimate.d.sequential}
\title{Estimate d using sequential method}
\usage{
estimate.d.sequential(geno.matrix, fit.matrix, d.hat.MLE, d.hat.RDB, d.range,
  d.max.adj = 1.1)
}
\arguments{
\item{geno.matrix}{Genotype matrix generated in
\code{\link{generate.geno.matrix}} or read in.}

\item{fit.matrix}{Fitness matrix generated in
\code{\link{sim.stick.data}} or read in.}

\item{d.hat.MLE}{The MLE estimate of d from \code{\link{estimate.d.MLE}}}

\item{d.hat.RDB}{The estimate of d from \code{\link{estimate.d.RDB}}}

\item{d.range}{Interval of d where MLE was searched. Defines when valid MLE exists (see details).}

\item{d.max.adj}{When forced to use the maximum estimator, the estimate is adjusted upwards
by this factor (see details). Default = 1.1 (inflate observation 10\%).}
}
\value{
Estimate of d. Name indicates methods used.
}
\description{
Estimate d using sequential method
}
\details{
If a valid MLE exists, function returns it. When MLE is at the boundary
as defined by \code{d.range} or less than the observed distance from wild type fitness to the maximum fitness,
this is not considered a valid estimate. If MLE is not valid, but RDB
estimate exists, it returns RDB estimate. If neither exists, returns the Max estimated.
The Max estimate is based on the largest observed fitness times a factor \code{d.max.adj}.
Name of the returned object indicates the method estimate is based on (MLE, RDB or Max).
}
\examples{
n.muts <- length(Khan.data[1,])-1
geno.matrix <- Khan.data[,seq(1, n.muts)]
fit.matrix <- as.matrix(Khan.data[,(n.muts+1)])
d.hat.MLE <- estimate.d.MLE(geno.matrix, fit.matrix,c(0.1, 10),0.001,c(2,1))
d.hat.RDB <- estimate.d.RDB(geno.matrix, fit.matrix,-100)$d.hat.RDB
estimate.d.sequential(geno.matrix, fit.matrix, d.hat.MLE, d.hat.RDB, c(0.1, 10), 1.1)
}

