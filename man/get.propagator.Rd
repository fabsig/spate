\name{get.propagator}
\alias{get.propagator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Propagator matrix G.
}
\description{
Function for obtaining the spectral propagator matrix G of the vector autoregressive model
for the Fourier coefficients.
}
\usage{
get.propagator(wave, indCos, zeta, rho1, gamma, alpha, muX, muY, dt = 1, ns=4)
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
Propagator matrix G.
}

\author{
Fabio Sigrist
}

\examples{
##For illustration, four grid points on each axis
n <- 4
wave <- wave.numbers(n)
G <- get.propagator(wave=wave$wave,indCos=wave$indCos,zeta=0.5, rho1=0.1, gamma=2,
           alpha=pi/4, muX=0.2, muY=-0.15,dt=1,ns=4)
round(G,digits=2)
## View(round(G,digits=2))

##An example
n <- 50
spec <- matern.spec(wave=spate.init(n=n,T=1)$wave,n=n,rho0=0.05,sigma2=1,norm=TRUE)
alphat <- sqrt(spec)*rnorm(n*n)
##Propagate initial state
wave <- wave.numbers(n)
G <- get.propagator(wave=wave$wave,indCos=wave$indCos,zeta=0.5, rho1=0.02, gamma=2,
            alpha=pi/4, muX=0.2, muY=0.2,dt=1,ns=4)
alphat1 <- G\%*\%alphat

opar <- par(no.readonly = TRUE)
par(mfrow=c(1,2))
image(1:n,1:n,matrix(real.fft(alphat,n=n,inv=FALSE),nrow=n),main="Whittle
field",xlab="",ylab="",col=cols())
image(1:n,1:n,matrix(real.fft(alphat1,n=n,inv=FALSE),nrow=n),main="Propagated
field",xlab="",ylab="",col=cols())
par(opar) # Reset par() settings
}
