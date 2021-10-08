\name{sample.four.coef}
\alias{sample.four.coef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Sample from the full conditional of the Fourier coefficients.
}
\description{
Sample from the full conditional of the Fourier coefficients.
}
\usage{
sample.four.coef(w=NULL,wFT=NULL,spec=NULL,Gvec=NULL,tau2=NULL,par=NULL,n,T,
                 NF=n*n,indCos=(1:((n*n-4)/2)*2+3),ns=4,nu=1,dt=1)
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
  A T x n*n matrix with a sample from the full conditional
  of latent process \eqn{\alpha}.
}

\examples{
##Specifications for simulated example
n <- 50
T <- 4
par <- c(rho0=0.1,sigma2=0.2,zeta=0.5,rho1=0.1,gamma=2,alpha=pi/4,muX=0.2,muY=-0.2,tau2=0.01)
spateSim <- spate.sim(par=par,n=n,T=T,seed=4)
w <- spateSim$w
##Sample from full conditional
Nmc <- 50
alphaS <- array(0,c(T,n*n,Nmc))
wFT <- real.fft.TS(w,n=n,T=T)
for(i in 1:Nmc){
  alphaS[,,i] <- sample.four.coef(wFT=wFT,par=par,n=n,T=T,NF=n*n)
}
##Mean from full conditional
alphaMean <- apply(alphaS,c(1,2),mean)
xiMean <- real.fft.TS(alphaMean,n=n,T=T,inv=FALSE)

opar <- par(no.readonly = TRUE)
par(mfrow=c(2,4),mar=c(1,1,1,1))
for(t in 1:4) image(1:n,1:n,matrix(w[t,],nrow=n),xlab="",ylab="",col=cols(),
                    main=paste("w(",t,")",sep=""),xaxt='n',yaxt='n')
for(t in 1:4) image(1:n,1:n,matrix(xiMean[t,],nrow=n),xlab="",ylab="",col=cols(),
                    main=paste("xiPost(",t,")",sep=""),xaxt='n',yaxt='n')
par(opar) # Reset par() settings
}

\author{
Fabio Sigrist
}
