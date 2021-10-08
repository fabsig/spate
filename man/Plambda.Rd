\name{Plambda}
\alias{Plambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prior for transformation parameter of the Tobit model.
}
\description{
  Default prior for transformation parameter of the Tobit model. A locally constant, improper prior on
the positive real line is used.
}
\usage{
Plambda(lambda, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{
A quantile
}
\item{log}{
  Indicates whether the logarithm should be calculated or not
}
}

\value{
  Value of (log) prior at 'lambda'.
}

\author{
  Fabio Sigrist
}
