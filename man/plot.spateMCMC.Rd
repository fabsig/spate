\name{plot.spateMCMC}
\alias{plot.spateMCMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot fitted spateMCMC objects.
}
\description{
Plots trace plots, pair plots, the posterior of the hyperparameters and the posterior of the latent spatio-temporal process. 
}
\usage{
\method{plot}{spateMCMC}(x,..., trace = TRUE, hist = TRUE,
              medianHist=TRUE, pairs = FALSE,ask = TRUE, ToFile = FALSE,
              path = NULL,file = NULL,true=NULL,BurnInAdaptive=NULL,
              postProcess = FALSE)
}

\arguments{
  \item{x}{
A 'spateMCMC' object obtained from 'spate.mcmc'.
}
\item{...}{Arguments to be passed to 'spate.plot' in case
  'postProcess=TRUE' is selected.
  }
  \item{trace}{
logical; if 'TRUE' trace plots are made
}
\item{hist}{
  logical; if 'TRUE' histograms of the posterior distributions for the hyper-parameters are plotted
}  \item{medianHist}{
  logical; if 'TRUE' medians are added to the histograms.
}
  \item{pairs}{
logical; if 'TRUE' trace plots are made
}
\item{ask}{logical;  if 'TRUE' (and the R session is interactive) the user is asked for input, before a new figure is drawn. 
}
  \item{ToFile}{
logical; if 'TRUE' the plots are save to a file.
}
  \item{path}{
The path.
}
  \item{file}{
The file name.
}
  \item{true}{
The true value of the parameters (for simulation only).
}
  \item{BurnInAdaptive}{
 The number of samples used as burn-in before starting the adaptive
 estimation of Metropolis-Hastings
 proposal covariance for the hyper-parameters.
}
  \item{postProcess}{
logical; if 'TRUE' the posterior of the spatio-temporal process xi is plotted as well.
}
}
\value{
Plots illustrating a fitted model saved in a 'spateMCMC' object.
}
\author{
Fabio Sigrist
}
\examples{
data("spateMCMC")
plot(spateMCMC,medianHist=TRUE,pairs=TRUE)
}
