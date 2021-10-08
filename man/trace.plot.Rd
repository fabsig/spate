\name{trace.plot}
\alias{trace.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Trace plots for MCMC output analysis.
}
\description{
  Trace plots for MCMC output analysis.
}
\usage{
trace.plot(data, true = NULL, BurnIn = NULL,BurnInAdaptive=NULL)
}
\arguments{
  \item{data}{
A p x Nmc data.frame of matrix where p denotes the number of parameters
and Nmc the number of Monte Carlo samples.
}
  \item{true}{
The true value of the parameters (for simulation only).
}
  \item{BurnIn}{
 The number of samples used as burn-in if the burn-in has not yet been
 removed from the sample.
}
  \item{BurnInAdaptive}{
 The number of samples used as burn-in before starting the adaptive
 estimation of Metropolis-Hastings
 proposal covariance for the hyper-parameters.
}
}
\value{
Trace plots.
}
\author{
Fabio Sigrist
}
\examples{
data <- matrix(rnorm(1200),nrow=6)
opar <- par(no.readonly = TRUE)
par(mfrow=c(2,3))
trace.plot(data,true=rep(0,6))
par(opar) # Reset par() settings
}
