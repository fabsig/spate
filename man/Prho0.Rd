\name{Prho0}
\alias{Prho0}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
  Prior for range parameter rho0 of innovation epsilon.
}
\description{
  Default prior for range parameter rho0 of stochastic source-sink term epsilon. A uniform prior on [0,100] is used.
}

\usage{
Prho0(rho0, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rho0}{
A quantile
}
\item{log}{
  Indicates whether the logarithm should be calculated or not
}
}

\value{
Value of (log) prior at 'rho0'.
}

\author{
Fabio Sigrist
}
