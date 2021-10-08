\name{map.obs.to.grid}
\alias{map.obs.to.grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Maps non-gridded data to a grid.
}
\description{
Maps non-gridded data to a grid based on the coordinates supplied. Cells with no data are NA. For cells with more than one data point, the average is taken.
}
\usage{
map.obs.to.grid(n,y.non.grid,coord,lengthx=NULL,lengthy=NULL)
}
\arguments{
  \item{y.non.grid}{
Observed data in an T x N matrix with columns and rows corresponding to
time and space, respectively. The coordinates of each observation point need to be specified in 'coord'. 
}  
\item{coord}{
Matrix of dimension N x 2 with coordinates of the N observation points. Based on to these coordinates, each observation location is then mapped to a grid cell.
}
  \item{lengthx}{
Use together with 'coord' to specify the length of the x-axis. This is
usefull if the observations lie in a rectangular area instead of a
square. The length needs to be at least as large as the largest
x-distance in 'coord.
}
\item{lengthy}{
  Use together with 'coord' to specify the length of the y-axis. This is
usefull if the observations lie in a rectangular area instead of a
square. The length needs to be at least as large as the largest
y-distance in 'coord.
}
  \item{n}{
Number of point per axis of the square into which the points are
mapped. In total, the process is modeled on a grid of size n*n.
  }
}

\value{
  The function returns data in an T x n^2 matrix with columns and rows corresponding to time and space, respectively. Cells with no data are NA. For cells with more than one data point, the average is taken.
}

\author{
Fabio Sigrist
}

\examples{
## See code of 'spate.mcmc'.

}
