\name{plot.spateSim}
\alias{plot.spateSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting function for 'spateSim' objects.
}
\description{
This is the plotting function for 'spateSim' objects. It calles the
function 'spate.plot()'.
}
\usage{
\method{plot}{spateSim}(x,..., plotXi =TRUE,plotW = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
'spateSim' object obtained from 'spate.sim'.
}
\item{...}{Arguments to be passed to 'spate.plot'}
  \item{plotXi}{
Logical; if 'TRUE' the latent process 'xi' is plotted.
}
  \item{plotW}{
Logical; if 'TRUE' the observed process 'w' is plotted.
  }
}
\value{
Plots illustrating the simulated space-time field.
}

\author{
Fabio Sigrist
}

\examples{
spateSim <-spate.sim(par=c(rho0=0.1,sigma2=0.2,zeta=0.5,rho1=0.1,gamma=2,
                     alpha=pi/4,muX=0.2,muY=-0.2,tau2=0.01),n=50,T=9)
plot(spateSim)
}
