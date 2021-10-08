\name{wave.numbers}
\alias{wave.numbers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wave numbers.
}
\description{
Returns wave numbers used in real Fourier transform.
}
\usage{
wave.numbers(n)
}
\arguments{
  \item{n}{
 Number of grid points on each axis. n x n is the total number of spatial points.
}
}

\value{
  Returns a list with
  \item{wave}{A 2 x n^2 matrix with wavenumbers used in the real Fourier
  transform. The first four columns contain the wavenumbers that are only used by cosine terms and not by sine terms. Subsequent columns alternate between wavenumbers of cosine and sine terms.}
  \item{indCos}{Vector of integers indicating the position of columns in
    'wave' of wavenumbers of cosine terms. The first four cosine wavenumbers in 'wave' are not included in 'indCos'.}
}

\author{
Fabio Sigrist
}
