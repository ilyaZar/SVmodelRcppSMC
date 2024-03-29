% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_pmmh_r.R
\name{bpfPMMH}
\alias{bpfPMMH}
\title{Boostrap Particle Filter for the toy SV model}
\usage{
bpfPMMH(NN, yt, phiX, sigmaX, betaY)
}
\arguments{
\item{NN}{number of particles}

\item{yt}{measurements/data of the unobservered latent state process}

\item{phiX}{phi parameter (see model equations for details)}

\item{sigmaX}{sigma parameter (see model equations for details)}

\item{betaY}{beta parameter (see model equations for details)}
}
\value{
a list of one element (was previously returniong more) giving the 
  log-likelihood estimate
}
\description{
Computes, besides other things, an estimate of the log-likelihood.
}
