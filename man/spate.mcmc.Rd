\name{spate.mcmc}
\alias{spate.mcmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 MCMC algorithm for fitting the model.
}
\description{
MCMC algorithm for fitting the model.
}
\usage{
spate.mcmc(y,coord=NULL,lengthx=NULL,lengthy=NULL,Sind=NULL,n=NULL,
          IncidenceMat=FALSE,x=NULL,SV=c(rho0=0.2,sigma2=0.1,
          zeta=0.25,rho1=0.2,gamma=1,alpha=0.3,muX=0,muY=0,tau2=0.005),
          betaSV=rep(0,dim(x)[1]),RWCov=NULL,parh=NULL,tPred=NULL,
          sPred=NULL,P.rho0=Prho0,P.sigma2=Psigma2,P.zeta=Pzeta,P.rho1=Prho1,
          P.gamma=Pgamma,P.alpha=Palpha,P.mux=Pmux,P.muy=Pmuy,P.tau2=Ptau2,
          lambdaSV=1,sdlambda=0.01,P.lambda=Plambda,DataModel="Normal",
          DimRed=FALSE,NFour=NULL,indEst=1:9,Nmc=10000,BurnIn =1000,
          path=NULL,file=NULL,SaveToFile=FALSE,PlotToFile=FALSE,
          FixEffMetrop=TRUE,saveProcess=FALSE,Nsave=200,seed=NULL,
          Padding=FALSE,adaptive=TRUE,NCovEst=500,BurnInCovEst=500,
          MultCov=0.5,printRWCov=FALSE,MultStdDevLambda=0.75,
          Separable=FALSE,Drift=!Separable,Diffusion=!Separable,
          logInd=c(1,2,3,4,5,9),nu=1,plotTrace=TRUE,
          plotHist=FALSE,plotPairs=FALSE,trueVal=NULL,
          plotObsLocations=FALSE,trace=TRUE,monitorProcess=FALSE,
          tProcess=NULL,sProcess=NULL)
}
\arguments{
  \item{y}{
Observed data in an T x N matrix with columns and rows corresponding to
time and space (observations on a grid stacked into a vector),
respectively. By default, at each time point, the observations are
assumed to lie on a square grid with each axis scaled so that it has unit length.
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
Logical; if 'TRUE' an incidence matrix relating the latent process to
observation locations is used. This is only recommended to use when the
observations are relatively low-dimensional and when the latent process
is modeled in a reduced dimensional space as well.
}
  \item{x}{
 Covariates in an array of dimensions p x T X N, where p denotes the
 number of covariates, T the number of time points, and N the number of
 spatial points.
}
  \item{SV}{
Starting values for parameters. Parameters for the SPDE in the following
order: rho_0, sigma^2, zeta, rho_1, gamma, alpha, mu_x, mu_y, tau^2. rho_0 and sigma^2 are the range
 and marginal variance of the Matern covariance funtion for the
 innovation term epsilon. zeta is the damping parameter. rho_1, gamma,
 and alpha parametrize the diffusion matrix with rho_1 being a range
 parameter, gamma and alpha determining the amount and the direction,
 respectively, of anisotropy. mu_x and mu_y are the two components of
 the drift vector. tau^2 denotes the nugget effect or measurment error.
}
  \item{betaSV}{
Starting values for regression coefficients.
}
\item{RWCov}{
 Covariance matrix of the proposal distribution used in the random walk
 Metropolis-Hastings step for the hyper-parameters.
}
\item{parh}{Only used in prediction mode. If 'parh'
  is not 'NULL', this indicates that 'spate.mcmc' is used for making
  predictions at locations (tPred,sPred) instead of applying the traditional MCMC
  algorithm. In case 'parh' is not 'NULL', it is a Npar x Nsim
  matrix containing Nsim samples from the posterior of the Npar
  parameters. This argument is used by the wrapper function 'spate.predict'.
 }
  \item{tPred}{
Time points where predictions are made.This needs to be a vector if
predictions are made at multiple times.  For instance, if T is the number
of time points in the data 'y', then tPred=c(T+1, T+2) means that
predictions are made at time 'T+1' and 'T+2'. This argument is used by
the wrapper function 'spate.predict'.
}
\item{sPred}{
  Vector of indices of grid cells (positions of locations in the stacked
  spatial vector) where predictions are made. This argument is used by
the wrapper function 'spate.predict'.
}
\item{P.rho0}{
Function specifying the prior for rho0.
}
\item{P.sigma2}{
Function specifying the prior for sigma2.
}
\item{P.zeta}{
Function specifying the prior for zeta.
}
\item{P.rho1}{
Function specifying the prior for rho1.
} 
\item{P.gamma}{
Function specifying the prior for gamma.
}
\item{P.alpha}{
Function specifying the prior for alpha.
}
\item{P.mux}{
Function specifying the prior for mux.
}
\item{P.muy}{
Function specifying the prior for muy.
}
\item{P.tau2}{
Function specifying the prior for tau2.
}
\item{lambdaSV}{
Starting value for transformation parameter lambda in the Tobit model.
}
\item{sdlambda}{
 Standard deviation of the proposal distribution used in the random walk
 Metropolis-Hastings step for lambda.
}
\item{P.lambda}{
Function specifying the prior for lambda.
}
\item{DataModel}{
Specifies the data model. "Normal" or "SkewTobit" are available options.
}  
\item{DimRed}{
Logical; if 'TRUE' dimension reduction is applied. This means that not
the full number (n*n) of Fourier functions is used but rather only a
reduced dimensional basis of dimension 'NFour'.
}
  \item{NFour}{
If 'DimRed' is 'TRUE', this specifies the number of Fourier functions.
}
  \item{indEst}{
 A vector of numbers specifying which for which parameters the posterior
 should be computed and which should be held fix (at their starting
 value). If the corresponding to the index of rho_0, sigma^2, zeta,
 rho_1, gamma, alpha, mu_x, mu_y, tau^2 is present in the vector, the
 parameter will be estimated otherwise not. Default is indEst=1:9 which
 means that one samples from the posterior for all parameters.
}
  \item{Nmc}{
Number of MCMC samples.
}
  \item{BurnIn}{
Length of the burn-in period.
}
  \item{path}{
Path, in case plots and / or the spateMCMC object should be save in a file.
}
\item{file}{
  File name, in case plots and / or the spateMCMC object should be save in a file.
}
  \item{SaveToFile}{
Indicates whether the spateMCMC object should be save in a file.
}
\item{PlotToFile}{
  Indicates whether the MCMC output analysis plots should be save in a file.
}
  \item{FixEffMetrop}{
The fixed effects, i.e., the regression coefficients, can either be
sampled in a Gibbs step or updated together with the hyperparameters in
the Metropolis-Hastings step. The latter is the default and recommended
option since correlations between fixed effects and the random process
can result in slow mixing.
}
  \item{saveProcess}{
Logical; if 'TRUE' samples from the posterior of the latent
spatio-temporal process xi are saved.
}
  \item{Nsave}{
Number of samples from the posterior of the latent
spatio-temporal process xi that should be save.
}
\item{seed}{
   Seed for random generator.
}
  \item{Padding}{
Indicates whether padding is applied or not. If the range parameters are
large relative to the domain, this is recommended since otherwise
spurious periodicity can occur.
}
  \item{adaptive}{
Indicates whether an adaptive Metropolis-Hastings algorithm is used or
not. If yes, the proposal covariance matrix 'RWCov' is adaptively
estimated during the algorithm and tuning does not need to be done by hand.
}
  \item{NCovEst}{
 Minimal number of samples to be used for estimating the proposal matrix.
}
  \item{BurnInCovEst}{
Burn-in period for estimating the proposal matrix.
}
\item{MultCov}{Numeric used as multiplier for the adaptively estimated
  proposal cocariance matrix 'RWCov' of the hyper-parameters. I.e., the estimated covariance
  matrix is multiplied by 'MultCov'.}
\item{printRWCov}{Logical, if 'TRUE' the estimated
  proposal cocariance matrix is printed each time.}
\item{MultStdDevLambda}{Numeric used as multiplier for the adaptively
  estimated proposal
  standard deviation of the Tobit transformation parameter lambda. I.e.,
  the estimated standard deviation is multiplied by 'MultStdDevLambda'.}
\item{Separable}{
Indicates whether a separable model, i.e., no transport / drift and no
diffusion, should be estimated.
}
  \item{Drift}{
Indicates whether a drift term should be included.
}
  \item{Diffusion}{
Indicates whether a diffusion term should be included.
}
 \item{logInd}{
Indicates which parameters are sampled on the log-scale. Default is
logInd=c(1, 2, 3, 4, 5, 9) corresponding to rho_0, sigma2,  zeta,
rho_1, gamma, and tau^2.
}
  \item{nu}{
Smoothness parameter of the Matern covariance function for the innovations. By
default this equals 1 corresponding to the Whittle covariance function.
} 
\item{plotTrace}{
Indicates whether trace plots are made.
}
  \item{plotHist}{
Indicates whether histograms of the posterior distributions are made.
}
  \item{plotPairs}{
Indicates whether scatter plots of the hyper-parameters and the
regression coefficients are made.
}
  \item{trueVal}{
In simulations, true values can be supplied for comparison with the MCMC output.
}  
\item{plotObsLocations}{
Logical; if 'TRUE' the observations locations are ploted together with
the grid cells.
}
\item{trace}{
Logical; if 'TRUE' tracing information on the progress of the MCMC algorithm is produced.
}
\item{monitorProcess}{
  Logical; if 'TRUE' in addition to the trace plots of the hyper-parameters, the mixing
  properties of the latent process xi=Phi*alpha is monitored.  This is
  done by plotting the current sample of the process. More specifically,
  the time series at locations 'sProcess' and the spatial fieldd at time points 'tProcess'.
  }
  \item{tProcess}{
    To be secified if 'monitorProcess=TRUE'. Time points at which
    spatial fields of the sampled process should be plotted.
    }
  \item{sProcess}{
    To be secified if 'monitorProcess=TRUE'. Locations at which time
    series of the sampled process should be plotted. 
    }
}

\value{
  The function returns a 'spateMCMC' object with, amongst others, the following entries
  \item{Post}{Matrix containing samples from the posterior of the hyper-parameters and the regression coefficient}
  \item{xiPost}{Array with samples from the posterior of the spatio-temporal process}
  \item{RWCov}{(Estimated) proposal covariance matrix}
}

\author{
Fabio Sigrist
}

\examples{
##Specify hyper-parameters
par <- c(rho0=0.1,sigma2=0.2,zeta=0.5,rho1=0.1,gamma=2,alpha=pi/4,muX=0.2,muY=-0.2,tau2=0.01)
##Simulate data
spateSim <- spate.sim(par=par,n=20,T=20,seed=4)
w <- spateSim$w

##Below is an example to illustrate the use of the MCMC algorithm.
##In practice, more samples are needed for a sufficiently large effective sample size.

##The following takes a couple of minutes.
##Load the precomputed object some lines below to save time.
##spateMCMC <- spate.mcmc(y=w,x=NULL,SV=c(rho0=0.2,sigma2=0.1,
##              zeta=0.25,rho1=0.2,gamma=1,alpha=0.3,muX=0,muY=0,tau2=0.005),
##              RWCov=diag(c(0.005,0.005,0.05,0.005,0.005,0.001,0.0002,0.0002,0.0002)),
##              Nmc=10000,BurnIn=2000,seed=4,Padding=FALSE,plotTrace=TRUE,NCovEst=500,
##              BurnInCovEst=500,trueVal=par,saveProcess=TRUE)
##spateMCMC
##plot(spateMCMC.fit,true=par,postProcess=TRUE)

##Instead of waiting, you can also use this precomputed object
data("spateMCMC")
spateMCMC
plot(spateMCMC,true=par,medianHist=FALSE)

}
