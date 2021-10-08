\name{matern.spec}
\alias{matern.spec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spectrum of the Matern covariance function.
}
\description{
  Spectrum of the Matern covariance function. Note that the spectrum is
  renormalized, by dividing with the sum over all frequencies so that
  they sum to one, so that
  \eqn{\sigma^2} is the marginal variance no matter how many
  wavenumbers are included.
}
\usage{
matern.spec(wave, n, ns=4, rho0, sigma2, nu = 1, norm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wave}{
    Spatial wavenumbers.
  }
  \item{n}{
 Number of grid points on each axis. n x n is the total number of spatial points.
    }
\item{ns}{Integer indicating the number of cosine-only terms. Maximally
  this is 4.}
  \item{rho0}{
Range parameter.
}
  \item{sigma2}{
Marginal variance parameter.
}
\item{nu}{
  Smoothness parameter of the Matern covariance function. By default this equals 1 corresponding to the Whittle covariance function.
}
\item{norm}{
  logical; if 'TRUE' the spectrum is multiplied by n*n so that after
  applying the real Fourier transform 'real.FFT' one has the correct normalization.
}
}
\details{
  The Matern covariance function is of the form
  \deqn{\sigma^2 2^(1-\nu) \Gamma(\nu)^{-1} (d/\rho_0)^{\nu} K_{\nu}(d/\rho_0)}
  with 'd' being the  Euclidean distance between two points and K_nu(.)
  a modified Bessel function. Its spectrum is given by
  \deqn{2^{\nu-1} \nu ((1/\rho_0)^(2\nu)) (\pi*((1/\rho_0)^2 + w)^(\nu + 1))^{-1}}  
  where 'w' is a spatial wavenumber.
}
\value{
Vector with the spectrum of the Matern covariance function.
}

\author{
Fabio Sigrist
}

\examples{
n <- 100
spec <- matern.spec(wave=spate.init(n=n,T=1)$wave,n=n,rho0=0.05,sigma2=1,norm=TRUE)
sim <- real.fft(sqrt(spec)*rnorm(n*n),n=n,inv=FALSE)
image(1:n,1:n,matrix(sim,nrow=n),main="Sample from a Gaussian process
with Matern covariance function",xlab="",ylab="",col=cols())
}

