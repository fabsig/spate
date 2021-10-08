\name{ffbs.spectral}
\alias{ffbs.spectral}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Forward Filtering Backward Sampling algorithm in the spectral space of
  the SPDE.
}
\description{
  Forward Filtering Backward Sampling algorithm for sampling from the
  joint full conditional of the coefficients \eqn{\alpha} and for
  evaluation of the log-likelihood.
}
\usage{
ffbs.spectral(w=NULL,wFT=NULL,spec=NULL,Gvec=NULL,tau2=NULL,par=NULL,n,T,lglk=FALSE,
              BwSp=TRUE,NF=n*n,indCos=(1:((n*n-4)/2)*2+3),ns=4,nu=1,dt=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{w}{Observed data or latent process w (depending on which data
    model is used) in an T x n*n matrix with columns and rows (points on a
    grid stacked into a vector) corresponding to time and space, respectively.
}
  \item{wFT}{
Vector of length T*n*n containing the real Fourier transform of 'w'.
}
  \item{spec}{
Spectrum of the innovations \eqn{\hat{\epsilon}} in a vector of length
n*n. If 'spec' is not given, it is constructed based on 'par'.
}
\item{Gvec}{
  The propagator matrix G in vector format obtained from 'get.G.vec'. If 'Gvec' is not given, it is constructed based on 'par'.
}
\item{tau2}{
  Measurement error variance tau2. If 'NULL'; tau2=par[9].
}
  \item{par}{ Vector of parameters for the SPDE in the following order: rho_0, sigma^2, zeta,
    rho_1, gamma, alpha, mu_x, mu_y, tau^2. If 'spec' and 'Gvec' are
    given, 'par' will not be used.
  }
\item{n}{Number of grid points on each axis. n*n is the total number of spatial points.
}
  \item{T}{
Number of points in time.
}
  \item{lglk}{
Logical; if 'TRUE' the value of the log-likelihood is returned as well.
}
  \item{BwSp}{
Logical; if 'TRUE' a sample from the full conditional of \eqn{\alpha} is returned.
}
  \item{NF}{
Number of Fourier functions used.
}
\item{indCos}{Vector of integers indicating the position cosine terms in
  the 1:NF real Fourier functions. The first 'ns' cosine wavenumbers in 'wave' are not included in 'indCos'.
}
  \item{ns}{
Number of real Fourier functions that have only a cosine and no sine
term. 'ns' is maximal 4.
}
 \item{nu}{ Smoothness parameter of the Matern covariance function for the innovations. By
default this equals 1 corresponding to the Whittle covariance function.}
  \item{dt}{ Temporal lag between two time points. By default, this equals 1.}
}
\value{
A list with entries (depending on whether 'lglk' are 'BwSp' are
'TRUE' or 'FALSE'):
\item{simAlpha}{A T x n*n matrix with a sample from the full conditional
  of latent process \eqn{\alpha},}
\item{ll}{The evaluated log-likelihood,}
}

\author{
Fabio Sigrist
}
