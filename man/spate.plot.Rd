\name{spate.plot}
\alias{spate.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a spatio-temporal field.
}
\description{
Generates a figure or an animation of a spatio-temporal field.
}
\usage{
spate.plot(xi,nx=NULL,whichT=NULL,format="ImgTogether",ToFile=FALSE,path=NULL,
           file=NULL,indScale=FALSE,main=NULL,mfrow=NULL,
           imagesize=c(1000,1000),zlim=NULL,breaks=NULL,...)
}
\arguments{
  \item{xi}{
    A spatio-temporal field stored in an T x N matrix with columns and rows corresponding to time and space, respectively.
  }
  \item{nx}{Integer specifying the number of points on the x-axis. If 'NULL',
    a quadratic grid is assumed.}
  \item{whichT}{Vector of integers specifying the time points that
    are plotted. If 'NULL', all time points are plotted.}
  \item{format}{
 A string specifying how the spatio-temporal field should be ploted. "ImgTogether" produces one single plot containing all spatial fields at all time
 points. With "ImgSeparate", the spatial fields at each
 time point are plotted in separate plots.
}  
\item{ToFile}{
 Indicates whether the output should be saved to a file.
}
  \item{path}{
Path indicating where to save the file.
}
  \item{file}{
File name.
}
  \item{indScale}{
 Indicates whether the color scale for the spatial plots is the same
 for all time points (indScale=FALSE) or separate for each time point (indScale=TRUE).
}  
\item{main}{Titles for the plots. Can be either be NULL or a character vector of length equal to the number of time points or 1.
}
\item{mfrow}{See 'par'. Can be either NULL or an integer vector of length two. If
it is NULL, the functions determines mfrow automatically.}
  \item{imagesize}{
 The size of the .jpeg image if ToFile=TRUE.
}
\item{zlim}{Graphical parameter to be passed to 'image'. Determines the scale on the z-axis
  of the plots. If 'indScale=FALSE' one can specify the common scale on
  the z-axis of the plots through this argument.}
\item{breaks}{Graphical parameter to be passed to 'image'.}
\item{...}{
Other graphical parameters that are passed to 'image' and 'par'.
  }
}
\value{
Plots illustrating a space-time field.
}
\author{
Fabio Sigrist
}
\examples{
spateSim <- spate.sim(par=c(rho0=0.1,sigma2=0.2,zeta=0.5,rho1=0.1,gamma=2,
                      alpha=pi/4,muX=0.2,muY=-0.2,tau2=0.01),n=50,T=9)
spate.plot(spateSim$xi)
}
