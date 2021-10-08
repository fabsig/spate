\name{Prho1}
\alias{Prho1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Prior for range parameter rho1 of diffusion.
}
\description{
  Default prior for range parameter rho1 of diffusive term. A uniform prior on
  [0,100] is used.
}
\usage{
Prho1(rho1, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rho1}{
A quantile.
}
\item{log}{
  Indicates whether the logarithm should be calculated or not
}
}

\value{
Value of (log) prior at 'rho1'.
}

\author{
Fabio Sigrist
}
