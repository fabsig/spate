\name{spate-package}
\alias{spate-package}
\alias{spate}
\docType{package}
\title{
  Spatio-temporal modeling of large data with the spectral SPDE approach
}
\description{This is an R package for spatio-temporal modeling of large data sets. It provides tools for modeling of Gaussian processes in space and time defined through a stochastic partial differential equation (SPDE). The SPDE is solved in the spectral space, and after discretizing in time and space, a linear Gaussian state space model is obtained. When doing inference, the main computational difficulty consists in evaluating the likelihood and in sampling from the full conditional of the spectral coefficients, or equivalently, the latent space-time process. In comparison to the traditional approach of using a spatio-temporal covariance function, the spectral SPDE approach is computationally advantageous. This package aims at providing tools for two different modeling approaches. First, the SPDE based spatio-temporal model can be used as a component in a customized hierarchical Bayesian model (HBM). The functions of the package then provide parametrizations of the process part of the model as well as computationally efficient algorithms needed for doing inference with the HBM. Alternatively, the adaptive MCMC algorithm implemented in the package can be used as an algorithm for doing inference without any additional modeling. The MCMC algorithm supports data that follow a Gaussian or a censored distribution with point mass at zero. Covariates can be included in the model through a regression term.
}
\author{
Fabio Sigrist, Hans R. Kuensch, Werner A. Stahel

Maintainer: Fabio Sigrist <sigrist@stat.math.ethz.ch> 
}
\references{
  Fabio Sigrist, Hans R. K\"unsch, and Werner A. Stahel, "Stochastic Partial Differential Equation Based Modeling of Large Space-Time Data Sets", Journal of the Royal Statistical Society: Series B, Volume 77, Issue 1, 2015, pages 3-33

  Fabio Sigrist, Hans R. Kuensch, Werner A. Stahel, "spate: An R Package for Spatio-Temporal Modeling with a Stochastic Advection-Diffusion Process.", Journal of Statistical Software, Volume 63, Number 14, 2015, pages 1-23, URL http://www.jstatsoft.org/v63/i14/
}

\keyword{ package }

