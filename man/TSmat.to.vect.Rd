\name{TSmat.to.vect}
\alias{TSmat.to.vect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Converts a matrix stacked vector.
}
\description{
Converts a time-space matrix with columns and rows
corresponding to time and space into a stacked N*T vector.
}
\usage{
TSmat.to.vect(mat)
}
\arguments{
  \item{mat}{
A T x N matrix with columns and rows corresponding to time and space, respectively.
}
}

\value{
A vector of stacked values. Stacking is done first over space and then time.
}

\author{
Fabio Sigrist
}
\examples{
vect <- 1:12
mat <- vect.to.TSmat(vect,T=3)##Convert vector to matrix
TSmat.to.vect(mat)##Convert matrix to vector.
}
