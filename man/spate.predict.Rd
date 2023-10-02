\name{spate.predict}
\alias{spate.predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Obtain samples from predictive distribution in space and time.
}
\description{
 Obtain samples from predictive distribution in space and time given the
 posterior of the hyperparameters.
}
\usage{
spate.predict(y,tPred,sPred=NULL,xPred=NULL,yPred=NULL,spateMCMC,Nsim=200,
              BurnIn=5,coord=NULL,lengthx=NULL,lengthy=NULL,Sind=NULL,
              n=NULL,IncidenceMat=FALSE,x=NULL,DataModel="Normal",
              DimRed=FALSE,NFour=NULL,seed=NULL,nu =1,trace=FALSE)
}
\arguments{
  \item{y}{
    Observed data in an T x N matrix with columns and rows corresponding to time and space, respectively.
}
\item{x}{
   Covariates in an array of dimensions p x T X N, where p denotes the
 number of covariates, T the number of time points, and N the number of
 spatial points.
}
  \item{tPred}{
Time points where predictions are made.This needs to be a vector if
predictions are made at multiple times.  For instance, if T is the number
of time points in the data 'y', then tPred=c(T+1, T+2) means that
predictions are made at time 'T+1' and 'T+2'.  If 'xPred' and 'yPred' are empty, then
predictions are made at all spatial points for each time point in
'tPred'. Otherwise 'xPred' and 'yPred', or 'sPred, need to have the same length as 'tPred', and
predictions are made at the points (tPred,xPred,yPred), ore (tPred,
sPred), respectively.
}
\item{sPred}{
  Vector of indices of grid cells (positions of locations in the stacked
  spatial vector) where predictions are made. This is an alternative to
  specifying the coordinates 'xPred' and 'yPred'.
}
\item{xPred}{
  Vector of x-coordinates of spatial points where predictions are made. This is an alternative to
  specifying the grid cell in 'sPred'.
}
\item{yPred}{
  Vector of y-coordinates of spatial points where predictions are made. This is an alternative to
  specifying the grid cell in 'sPred'.
}
  \item{spateMCMC}{
'spateMCMC' object obtained from 'spate.mcmc' containing the posterior of the hyper-parameters and
information on the model used.
}
  \item{Nsim}{
Number of samples used to characterize the predictive distribution.
}
  \item{BurnIn}{
Length of burn-in period.
}
\item{coord}{
If specified, this needs to be a matrix of dimension N x 2 with coordinates of the N observation points. Observations in 'y' can either be on a square grid or not. If not, the coordinates of each observation point need to be specified in 'coord'. According to these coordinates, each observation location is then mapped to a grid cell. If 'coord' is not specified, the observations in 'y' are assumed to lie on a square grid with each axis scaled so that it has unit length.
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
\item{Sind}{
  Vector of indices of grid cells where observations are made, in case,
  the observation are not made at every grid cell. Alternatively, the
  coordinates of the observation locations can be specfied in 'coord'.
}
 \item{n}{
Number of point per axis of the square into which the points are
mapped. In total, the process is modeled on a grid of size n*n.
  }
  \item{IncidenceMat}{
Logical; if 'TRUE' and incidence matrix relating the latent process to
observation locations. This is only recommended to use when the
observations are relatively low-dimensional and when the latent process
is modeled in a reduced dimensional spaceas well.
}
\item{DataModel}{
  Specifies the data model. "Normal" or "SkewTobit".
}
\item{DimRed}{
Logical; if 'TRUE' dimension reduction is applied. This means that not
the full number (n*n) of Fourier functions is used but rather only a
reduced dimensional basis of dimension 'NFour'.
}
  \item{NFour}{
If 'DimRed' is 'TRUE', this specifies the number of Fourier functions.
}
\item{seed}{
  Seed for random generator.
}
  \item{nu}{
Smoothness parameter of the Matern covariance function for the innovations. By
default this equals 1 corresponding to the Whittle covariance function.
} 
\item{trace}{
Logical; if 'TRUE' tracing information on the progress of the MCMC algorithm is produced.
}
}
\value{
  Depending on whether 'xPred' and 'yPred' are empty or not, either
  \item{(i)}{an array of size t x s x Nsim where the first index is for time, the second
    for space, and the third for the number of samples 'Nsim'}
  or
  \item{(ii)}{a matrix of size length(tPred) x Nsim}
}
\author{
Fabio Sigrist
}
\examples{
par <- c(rho0=0.1,sigma2=0.2,zeta=0.5,rho1=0.1,gamma=2,alpha=pi/4,muX=0.2,muY=-0.2,tau2=0.01)
##Simulate data
spateSim <- spate.sim(par=par,n=20,T=20,seed=4)
w <- spateSim$w
data("spateMCMC")
##Make predictions. Takes a couple of seconds
predict <- spate.predict(y=w,tPred=(17:25),spateMCMC=spateMCMC,Nsim =200,
                      BurnIn=10,DataModel="Normal")

Pmean <- apply(predict,c(1,2),mean)
Psd <- apply(predict,c(1,2),sd)

opar <- par(no.readonly = TRUE)
par(mfrow=c(2,2))
image(1:20,1:20,matrix(w[19,],nrow=20),main="Observed field at t=19",xlab="x",ylab="y")
image(1:20,1:20,matrix(Pmean[3,],nrow=20),main="Fitted field at t=19",xlab="x",ylab="y")
image(1:20,1:20,matrix(w[20,],nrow=20),main="Observed field at t=20",xlab="x",ylab="y")
image(1:20,1:20,matrix(Pmean[4,],nrow=20),main="Fitted field at t=20",xlab="x",ylab="y")

par(mfrow=c(3,3))
zlim=c(min(Pmean),max(Pmean))
for(i in 1:9){
  image(1:20,1:20,matrix(Pmean[i,],nrow=20),zlim=zlim,
        main=paste("Mean t=",i+16,sep=""),xlab="x",ylab="y")
}

par(mfrow=c(3,3))
zlim=c(min(Psd),max(Psd))
for(i in 1:9){
  image(1:20,1:20,matrix(Psd[i,],nrow=20),zlim=zlim,
        main=paste("Std.dev. t=",i+16,sep=""),xlab="x",ylab="y")
par(opar) # Reset par() settings
}
}
