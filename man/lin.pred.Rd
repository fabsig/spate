\name{lin.pred}
\alias{lin.pred}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear predictor.
}
\description{
Calculates the linear predictor.
}
\usage{
lin.pred(x, beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Covariates in an array of dimensions p x T X N, where p denotes the
 number of covariates, T the number of time points, and N the number of
 spatial points.
}
  \item{beta}{
Coefficients of covariates in a vector of length p.
}
}

\value{
Matrix of dimension T x N with linear predictors.
}

\author{
Fabio Sigrist
}
