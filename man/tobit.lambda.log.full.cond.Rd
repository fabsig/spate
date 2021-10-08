\name{tobit.lambda.log.full.cond}
\alias{tobit.lambda.log.full.cond}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Full conditional for transformation parameter lambda.
}
\description{
 Full conditional for transformation parameter lambda of the Tobit
 model. This is used in the Metropolis-Hastings step of the MCMC algorithm.
}
\usage{
tobit.lambda.log.full.cond(y, z, tau2, lambda)
}
\arguments{
  \item{y}{
 Observed data.
}
  \item{z}{
 Latent Gaussian variable.
}
  \item{tau2}{
 Value of variance (corresponds to nugget effect).
}
  \item{lambda}{
Value of transformation parameter lambda.
}
}
\value{
Value of evaluated full conditional for transformation parameter lambda.
}
\author{
Fabio Sigrist
}

