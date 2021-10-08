\name{ffbs}
\alias{ffbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Forward Filtering Backward Sampling algorithm.
}
\description{
Forward Filtering Backward Sampling algorithm for sampling from the
joint full conditional of the hidden state of a linear, Gaussian state
space model. To be more specific, one samples from \eqn{P[\alpha|.]} where
\eqn{\alpha} is specified through
\deqn{y_t = lp_t + H xi_t + nu_t,  \nu_t ~ N(0,\Omega)}
and
\deqn{\alpha_t = G \alpha_{t-1} + \epsilon_t,  \epsilon_t ~
  N(0,\Sigma).}
}
\details{In the context of the SPDE, \eqn{\alpha} are the Fourier coefficients.}
\usage{
ffbs(y, lp, G, Sigma, H, Omega, N = dim(y)[2],T = dim(y)[1],
      NF = dim(G)[1], lglk = FALSE, BwSp = TRUE, filt = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Observed data in an T x N matrix with columns and rows corresponding to time and space, respectively.
}
\item{lp}{
  Mean (linear predictor) in an T x N matrix with columns and rows corresponding to time and space, respectively.
}
  \item{G}{
Propagator matrix of the latent process \eqn{\alpha}.
}
  \item{Sigma}{
Innovation covariance matrix of the latent process \eqn{\alpha}.
}
  \item{H}{
Observation matrix relating y to \eqn{\alpha}.
}
  \item{Omega}{
Covariance matrix of the observation error \eqn{\nu}.
}
  \item{N}{
Number of points in space.
}
  \item{T}{
Number of points in time.
}
  \item{NF}{
Dimension of the latent process \eqn{\alpha}.
}
  \item{lglk}{
Logical; if 'TRUE' the value of the log-likelihood is returned as well.
}
  \item{BwSp}{
Logical; if 'TRUE' a sample from the full conditional of \eqn{\alpha} is returned.
}
  \item{filt}{
Logical; if 'TRUE' the filtered values for \eqn{\alpha} are returned.
}
}

\value{
A list with entries (depending on whether 'lglk', 'BwSp', 'filt' are
'TRUE' or 'FALSE'):
\item{simAlpha}{A T x N matrix with a sample from the full conditional
  of latent process \eqn{\alpha},}
\item{ll}{The evaluated log-likelihood,}
\item{mtt}{A T x N matrix with the mean of the full conditional of latent process \eqn{\alpha}.}
}

\author{
Fabio Sigrist
}
