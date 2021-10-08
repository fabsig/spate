\name{mcmc.summary}
\alias{mcmc.summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary function for MCMC output.
}
\description{
Auxilary function for summarizing MCMC output and illustrating the posterior distributions.
}
\usage{
mcmc.summary(data, probs = c(0.025, 0.5, 0.975), mean = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Matrix of size p x Nmc where p denotes the number of parameters and Nmc
the number of MCMC samples.
}
  \item{probs}{
 Vector of quantiles that should be computed for each parameter.
}
  \item{mean}{
logical; if 'TRUE' the mean of the posterior distributions is computed
as well.
}
}

\value{
Matrix with quantiles and the mean of the posterior distributions.
}

\author{
Fabio Sigrist
}

\examples{
data("spateMCMC")
mcmc.summary(spateMCMC$Post,mean=TRUE)
}
