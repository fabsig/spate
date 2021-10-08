\name{print.spateSim}
\alias{print.spateSim}
\title{
Print function for 'spateSim' objects.
}
\description{
Print function for 'spateSim' objects.  
}
\usage{
\method{print}{spateSim}(x,...)
}
\arguments{
  \item{x}{
    'spateSim' object obtained from 'spate.sim'.
  }
  \item{...}{Arguments to be passed to 'spate.plot'}
}

\author{
Fabio Sigrist
}
\examples{
spateSim <-spate.sim(par=c(rho0=0.1,sigma2=0.2,zeta=0.5,rho1=0.1,gamma=2,
                     alpha=pi/4,muX=0.2,muY=-0.2,tau2=0.01),n=50,T=9)
spateSim
}
