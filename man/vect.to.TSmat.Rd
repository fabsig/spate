\name{vect.to.TSmat}
\alias{vect.to.TSmat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts a stacked vector into matrix.
}
\description{
Converts a stacked N*T vector into a time-space matrix with columns and rows
corresponding to time and space, respectively.
}
\usage{
vect.to.TSmat(vect, T = 1)
}
\arguments{
  \item{vect}{
A vector of stacked values. Stacking is done first over space and then time.
}
  \item{T}{
Number of time points.
}
}

\value{
A T x N matrix with columns and rows corresponding to time and space, respectively.
}

\author{
Fabio Sigrist
}

\examples{
vect <- 1:12
vect
vect.to.TSmat(vect,T=3)
}
