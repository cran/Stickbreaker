% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/batch_analysis_fnxs.R
\name{analyze.mult.add.data.batch}
\alias{analyze.mult.add.data.batch}
\title{Analyze batch data generated by \code{\link{sim.fit.mult.add.data.batch}}}
\usage{
analyze.mult.add.data.batch(mut.vals, coe.vals, sig.vals, inpath, relative,
  coes)
}
\arguments{
\item{mut.vals}{Vector of number of mutations to simulate}

\item{coe.vals}{Vector of stickbreaking coefficients to analyze}

\item{sig.vals}{Vector of sigma values to analyze}

\item{inpath}{The path to read input files from. Generally same as \code{outpath} from
\code{\link{sim.fit.stick.data.batch}}.}

\item{relative}{\code{TRUE/FALSE}. Should bias and rMSE be divided by parameter value to make them relative values?}

\item{coes}{Vector of stickbreaking coefficients for each mutation}
}
\value{
\code{error.table} A data frame with the parametric conditions and the failure rate, bias and rMSE
for the parameters d, the stickbreaking coefficients, and sigma. \cr
}
\description{
Analyze batch data generated by \code{\link{sim.fit.mult.add.data.batch}}
}
\details{
This function assesses bias and relative root mean squared error (rMSE) from the data generated by
\code{\link{sim.fit.mult.add.data.batch}}.
}

