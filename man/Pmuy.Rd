\name{Pmuy}
\alias{Pmuy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prior for y-component of drift.
}
\description{
  Default prior for y-component of drift vector mu. A uniform prior on the
  interval [-0.5,0.5] is used.
}
\usage{
Pmuy(muy, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{muy}{
A quantile
}
\item{log}{
  Indicates whether the logarithm should be calculated or not
}
}

\value{
Value of (log) prior at 'muy'.
}

\author{
Fabio Sigrist
}

