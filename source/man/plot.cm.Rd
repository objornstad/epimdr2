% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network.R
\name{plot.cm}
\alias{plot.cm}
\title{Function to plot an object of class CM}
\usage{
\method{plot}{cm}(x, ...)
}
\arguments{
\item{x}{an object of class cm}

\item{...}{other arguments}
}
\value{
A plot of the contract matrix
}
\description{
Function to plot an object of class CM
}
\examples{
cm=ringlattice(N=20,K=4)
\dontrun{plot(cm)}
}
