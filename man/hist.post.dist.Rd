\name{post.dist.hist}
\alias{post.dist.hist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Histogram of posterior distributions.
}
\description{
Histogram of posterior distributions.
}
\usage{
post.dist.hist(data, true=NULL, breaks = 20, mean = FALSE, median = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Matrix of size p x Nmc where p denotes the number of
    parameters and Nmc the number of MCMC samples.
  }
  \item{true}{
The true value of the parameters (for simulation only).
}
  \item{breaks}{
Parameter for 'hist()' function.
}
  \item{mean}{
logical; if 'TRUE' the mean is added to the histogram.
}
\item{median}{
  logical; if 'TRUE' the median is added to the histogram.
}
}

\value{
Histograms illustrating posterior distributions.
}

\author{
Fabio Sigrist
}
