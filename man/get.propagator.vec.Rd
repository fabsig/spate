\name{get.propagator.vec}
\alias{get.propagator.vec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Propagator matrix G in vector form.
}
\description{
Function for obtaining the spectral propagator matrix G of the vector autoregressive model
for the Fourier coefficients in vector form.
}
\usage{
get.propagator.vec(wave, indCos, zeta, rho1, gamma, alpha, muX, muY, dt = 1,ns=4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave}{
    Spatial wavenumbers.
  }
   \item{indCos}{
    Vector of integers indicating the position of columns in 'wave' of wavenumbers of cosine terms.
} 
\item{zeta}{
  Damping parameter
}
  \item{rho1}{
Range parameter of the diffusion term
}
  \item{gamma}{
Parameter that determines the amount of anisotropy in the diffusion term
}
  \item{alpha}{
Parameter that determines the direction of anisotropy in the diffusion term
}
  \item{muX}{
X component of the drift vector.
}
  \item{muY}{
Y component of the drift vector.
}
\item{dt}{
  Temporal lag between two time points. By default, this equals 1.
}
\item{ns}{Number of real Fourier functions that have only a cosine and no sine
term. 'ns' is maximal 4.
}
}

\value{
A list with three elements 'G11C', 'G11', and 'G12'. The first element contains a vector of
length 'ns' which corresponds to the diagonal propagator of the
cosin-only terms. The second element contains the remaining diagonal
entries of G, i.e., the diagonal entries of the propagator for the
cosine / sine pairs. Note that for each pair, only one value is taken
since the diagonal elements for both the cosin and sine terms are equal. The third element is a vector with the off-diagonal
terms of the propagator for the cosine / sine pairs.
}

\author{
Fabio Sigrist
}

\examples{
##For illustration, four grid points on each axis
n <- 4
wave <- wave.numbers(n)
G <- get.propagator(wave=wave$wave,indCos=wave$indCos,zeta=0.5, rho1=0.1,
           gamma=2,alpha=pi/4, muX=0.2, muY=-0.15,dt=1,ns=4)
diag(G)[1:4]
diag(G[wave$indCos,wave$indCos])
diag(G[wave$indCos,wave$indCos+1])
get.propagator.vec(wave=wave$wave,indCos=wave$indCos,zeta=0.5, rho1=0.1,
gamma=2,alpha=pi/4, muX=0.2, muY=-0.15,dt=1,ns=4)
}
