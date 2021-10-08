\name{Ptau2}
\alias{Ptau2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prior for nugget effect parameter tau2.
}
\description{
Default prior for measurment error or small scale variation tau2 (nugget effect). A uniform, improper prior on tau (P[tau] propto 1 or P[tau2] propto 1/tau) is used.
}
\usage{
Ptau2(tau2, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tau2}{
A quantile
}
  \item{log}{
Indicates whether the logarithm should be calculated or not
}
}

\value{
  Value of (log) prior at tau2.
}

\author{
Fabio Sigrist
}

