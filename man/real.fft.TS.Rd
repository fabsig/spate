\name{real.fft.TS}
\alias{real.fft.TS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Fast calculation of the two-dimensional real Fourier transform of a space-time field. For
 each time point, the spatial field is transformed.
}
\description{
  This function calculates the two-dimensional real Fourier transform of
  a space-time field. This is a wrapper
  function for a C function which uses the complex FFT function from the
  'fftw3' library. In contrast to using T times the function 'real.FFT',
  R needs to communicate with C only once and not T times which saves
  computational time.
}
\usage{
real.fft.TS(w,n,T,inv=TRUE,indFFT=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{w}{
    Spatio-temporal field in a stacked vector of length T x N. Stacking is done
    first over space and then time. E.g., the first N=n^2 entries contain
    the spatial field at time t=1. Note that the spatial field itself is
    stacked as well, i.e., each spatial field is in a vector of length N=n^2.
}
\item{n}{
  Number of grid points on each axis. n x n is the total number of spatial points.
}
\item{T}{
  Number of time points.
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
A vector of length T x N containing the (inverse) Fourier transform of 'w'.
}

\author{
Fabio Sigrist
}

\examples{
n <- 100
T <- 4
spec <- matern.spec(wave=spate.init(n=n,T=1)$wave,n=n,rho0=0.05,sigma2=1,norm=TRUE)
specsim <- matrix(0,nrow=T,ncol=n*n)
for(t in 1:T) specsim[t,] <- rnorm(n*n)*sqrt(spec)
maternsim <- vect.to.TSmat(real.fft.TS(TSmat.to.vect(specsim),n=n,T=T,inv=FALSE),T=T)
opar <- par(no.readonly = TRUE)
par(mfrow=c(2,2))
for(t in 1:T) image(1:n,1:n,matrix(maternsim[t,],nrow=n),
                    main="Sample from Matern field",xlab="",ylab="")
par(opar) # Reset par() settings
}
