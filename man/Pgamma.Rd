\name{Pgamma}
\alias{Pgamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prior for amount of anisotropy in diffusion parameter gamma.
}
\description{
Default prior for amount of anisotropy in diffusion parameter gamma. A uniform
prior on log(gamma) over the interval [1/100,100] is used.
}
\usage{
Pgamma(gamma, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gamma}{
A quantile
}
\item{log}{
  Indicates whether the logarithm should be calculated or not
}
}

\value{
  Value of (log) prior at 'gamma'.
}

\author{
Fabio Sigrist
}

