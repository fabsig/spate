\name{Psigma2}
\alias{Psigma2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Prior for for variance parameter sigma2 of innovation epsilon.

  hyperparameter.
}
\description{
  Default prior for marginal variance parameter sigma2 (=sigma^2) of the stochastic source-sink term epsilon. A uniform, improper prior on sigma (P[sigma] propto 1 or P[sigma2] propto 1/tau) is used.
}
\usage{
Psigma2(sigma2, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sigma2}{
A quantile
}
  \item{log}{
Indicates whether the logarithm should be calculated or not

}
}

\value{
  Value of (log) prior at 'sigma2'.
}

\author{
Fabio Sigrist
}
