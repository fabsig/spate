\name{propagate.spectral}
\alias{propagate.spectral}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that propagates a state (spectral coefficients). 
}
\description{
Function that propagates the vector 'alphat'. This is
equivalent to multiplying 'alphat' with the propagator matrix G. It is a
lot faster though, due to the block-diagonal structure of G. This is a
wrapper function of a C function.
}
\usage{
propagate.spectral(alphat,spateFT=NULL,n=NULL,Gvec=NULL,par=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alphat}{
    A vector of spectral coefficients.
  }
   \item{spateFT}{
   A 'spateFT' obtained from 'spate.init'. Either this or 'n' needs to
   be given.
 }
   \item{n}{
 Number of points on each axis. n x n is the total number of spatial
 points. Either this or 'spateFT' needs to be given.
}
\item{Gvec}{
  The propagator matrix G in vector format obtained from
  'get.propagator.vec'. If 'Gvec' is not given, it is constructed based on 'par'.
}
  \item{par}{
 Parameters for the SPDE in the following order: rho_0, sigma^2, zeta,
 rho_1, gamma, alpha, mu_x, mu_y, tau^2. If 'Gvec' is not given, 'par'
 needs to be given.
}
}


\value{
  A vector of propagated coefficients G*alphat.
}

\author{
Fabio Sigrist
}

\examples{
n <- 50
spec <- matern.spec(wave=spate.init(n=n,T=1)$wave,n=n,rho0=0.05,sigma2=1,norm=TRUE)
alphat <- sqrt(spec)*rnorm(n*n)
##Propagate initial state
wave <- wave.numbers(n)
Gvec <- get.propagator.vec(wave=wave$wave,indCos=wave$indCos,zeta=0.1,rho1=0.02,gamma=2,
                  alpha=pi/4,muX=0.2,muY=0.2,dt=1,ns=4)
alphat1 <- propagate.spectral(alphat,n=n,Gvec=Gvec)

opar <- par(no.readonly = TRUE)
par(mfrow=c(1,2))
image(1:n,1:n,matrix(real.fft(alphat,n=n,inv=FALSE),nrow=n),main="Whittle
field",xlab="",ylab="",col=cols())
image(1:n,1:n,matrix(real.fft(alphat1,n=n,inv=FALSE),nrow=n),main="Propagated
field",xlab="",ylab="",col=cols())
par(opar) # Reset par() settings
}
