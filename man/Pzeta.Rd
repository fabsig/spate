\name{Pzeta}
\alias{Pzeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Prior for damping parameter zeta.
}
\description{
Default prior for damping parameter zeta. A uniform, improper prior on
the positive real line (P[zeta] propto 1) is used.
}
\usage{
Pzeta(zeta, log = FALSE)
}

\arguments{
  \item{zeta}{
A quantile
}
  \item{log}{
Indicates whether the logarithm should be calculated or not
}
}

\value{
Value of (log) prior at 'zeta'.
}

\author{
Fabio Sigrist
}
