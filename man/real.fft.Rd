\name{real.fft}
\alias{real.fft}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Fast calculation of the two-dimensional real Fourier transform.
}
\description{
  Fast calculation of the real Fourier transform. This is a wrapper
  function for a C function which uses the complex FFT function from the 'fftw3' library.
}
\usage{
real.fft(w,n,inv=TRUE,indFFT=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{w}{
A spatial field in a stacked vector of length N=n^2.
}
\item{n}{
  Number of grid points on each axis. n x n is the total number of spatial points.
}
  \item{inv}{
Indicates whether the inverse Fourier transform should be calculated or not.
}
  \item{indFFT}{
A list of containing vectors of natural numbers representing indices
used to transform between the real and the complex Fourier transform.
  }
}

\value{
A vector of length n*n containing the real (inverse) Fourier transformation of 'w'.
}

\author{
Fabio Sigrist
}

\examples{
n <- 100
spec <- matern.spec(wave=spate.init(n=n,T=1)$wave,n=n,rho0=0.05,sigma2=1,norm=TRUE)
sim <- real.fft(sqrt(spec)*rnorm(n*n),n=n,inv=FALSE)
image(1:n,1:n,matrix(sim,nrow=n),main="Sample from Matern field",xlab="",ylab="")
}
