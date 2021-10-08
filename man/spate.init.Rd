\name{spate.init}
\alias{spate.init}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Constructor for 'spateFT' object which are used for the two-dimensional
 Fourier transform.
}
\description{
Auxilary function for constructing 'spateFT' objects which are used for the two-dimensional
 Fourier transform.
}
\usage{
spate.init(n,T,NF=n*n)
}
\arguments{
  \item{n}{
 Number of points on each axis. n x n is the total number of spatial points.
}
  \item{T}{
Number of temporal points.
}
  \item{NF}{
This integer specifies the number of Fourier functions. If NF<n*n,
 dimension reduction is obtained. In this case, Fourier functions with
 wavenumbers closest to the origin (0,0) are first included. If a given
 'NF' implies a basis with anisotropy, i.e., higher frequencies in one
 direction than in another, this is automatically corrected by using the
 next higher integer NF' such that the basis has the same resolution in
 all directions.
}

}
\value{
  A 'spateFT' object. This is a list with
  \item{wave}{a matrix containing the wavenumbers}
  \item{indCos}{a vector indicating the position of the cosine terms (excluding the frist 'ns')}
  \item{ns}{an integer indicating the number of cosine-only terms}
  \item{indFFT}{a list of indices used for the conversion between the complex FFT and the real Fourier transform.}
  \item{n}{number of points on each axis}
  \item{T}{number of points in time}
}
\author{
Fabio Sigrist
}
