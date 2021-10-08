\name{summary.spateSim}
\alias{summary.spateSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary function for 'spateSim' objects.
}
\description{
Summary function for 'spateSim' objects.  
}
\usage{
\method{summary}{spateSim}(object,...)
}
\arguments{
  \item{object}{
    'spateSim' object obtained from 'spate.sim()'.
  }
  \item{...}{not used.}
}
\author{
Fabio Sigrist
}

\examples{
spateSim <- spate.sim(par=c(rho0=0.1,sigma2=0.2,zeta=0.5,rho1=0.1,gamma=2,
                      alpha=pi/4,muX=0.2,muY=-0.2,tau2=0.01),n=50,T=9)
summary(spateSim)
}
