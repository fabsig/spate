\name{Palpha}
\alias{Palpha}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prior for direction of anisotropy in diffusion parameter alpha.
}
\description{
Default prior for direction of anisotropy in diffusion parameter alpha. A uniform
prior on [0,pi/4] is used.
}
\usage{
Palpha(alpha, log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{
A quantile
}
\item{log}{
    Indicates whether the logarithm should be calculated or not
}
}

\value{
  Value of (log) prior at 'alpha'.

}

\author{
Fabio Sigrist
}
