\name{innov.spec}
\alias{innov.spec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spectrum of the innovation term epsilon.
}
\description{
  Spectrum of the innovation term epsilon.
}
\usage{
innov.spec(wave,n,ns=4,rho0,sigma2,zeta,rho1,alpha,gamma,nu=1,dt=1,norm=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave}{
    Spatial wavenumbers.
  }
   \item{n}{
 Number of grid points on each axis. n x n is the total number of spatial points.
    }
\item{ns}{Number of real Fourier functions that have only a cosine and no sine
term. 'ns' is maximal 4.
}
\item{rho0}{
  Range of the Matern covariance funtion for the innovation term epsilon
}
  \item{sigma2}{
Marginal variance of the Matern covariance funtion for the innovation term epsilon
}
  \item{zeta}{
Damping parameter
}
  \item{rho1}{
Range parameter of the diffusion term
}
  \item{alpha}{
Parameter that determines the direction of anisotropy in the diffusion term
}
\item{gamma}{
  Parameter that determines the amount of anisotropy in the
  diffusion term
}
  \item{nu}{
Smoothness parameter of the Matern covariance function for the innovations. By
default, this equals 1 corresponding to the Whittle covariance function.
}  
\item{dt}{
Temporal lag between two time points. By default, this equals 1.
}
\item{norm}{logical; if 'TRUE' the spectrum is multiplied by n*n so that after
  applying the real Fourier transform 'real.FFT' one has the correct normalization.
}
}

\value{
  Vector with the spectrum of the integrated innovation term epsilon hat.
}

\author{
Fabio Sigrist
}

\examples{
n <- 100
spec <- innov.spec(wave=spate.init(n=n,T=1)$wave,n=n,rho0=0.05,sigma2=0.5,zeta=0.5,
                   rho1=0.05,alpha=pi/4,gamma=2,norm=TRUE)
sim <- real.fft(sqrt(spec)*rnorm(n*n),n=n,inv=FALSE)
image(1:n,1:n,matrix(sim,nrow=n),main="Sample from the integrated
stochastic innovation",xlab="",ylab="",col=cols())
}
