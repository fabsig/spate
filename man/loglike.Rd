\name{loglike}
\alias{loglike}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log-likelihood of the hyperparameters.
}
\description{
Evaluates the log-likelihood of the hyperparameters given the data (Gaussian case) or
given the latent variable w (in the Tobit case).
}
\usage{
loglike(par=NULL,w=NULL,wFT=NULL,x=NULL,spec=NULL,Gvec=NULL,tau2=NULL,n,T,
         NF=n*n,indCos=(1:((n*n-4)/2)*2+3),ns=4,nu=1,dt=1,logScale=FALSE,
         logInd=c(1,2,3,4,5,9),negative=FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{ Vector of parameters for the SPDE in the following order: rho_0, sigma^2, zeta,
    rho_1, gamma, alpha, mu_x, mu_y, tau^2, regression coefficients beta. rho_0 and sigma^2 are the range
    and marginal variance of the Whittle covariance funtion for the
    innovation term epsilon. zeta is the damping parameter. rho_1, gamma,
    and alpha parametrize the diffusion matrix with rho_1 being a range
    parameter, gamma and alpha determining the amount and the direction,
    respectively, of anisotropy. mu_x and mu_y are the two components of
    the drift vector. tau^2 denotes the variance of nugget effect or
    measurment error. Subsequently in par are the
    regression coefficients beta, if there are covariates.
  }
  \item{w}{
    Matrix of size T x N, where T and N denote the number of points in time
    and space. In the case of a Gaussian data model, w contains the
    observed values, with the Tobit model, w denotes the latent normal variable.
  }
  \item{wFT}{ A vector with the (discrete) Fourier transform of the observed or latent w, depending
    on which data model is used. Note that, in contrast to w, this needs to
    be in stacked vector format. Use 'TSmat.to.vect'.
  }
  \item{x}{ Covariates in an array of dimensions p x T X N, where p denotes the
    number of covariates, T the number of time points, and n the number of
    spatial points.
  }
  \item{spec}{
    A vector containing the spectrum of the innovation term epsilon. If 'spec' is not given, it is constructed based on 'par'.
  }
  \item{Gvec}{
    The propagator matrix G in vector format obtained from 'get.G.vec'. If 'Gvec' is not given, it is constructed based on 'par'.
  }
  \item{tau2}{
    Measurement error variance tau2. If 'NULL'; tau2=par[9].
  }
  \item{n}{
    Number of grid points on each axis. n x n is the total number of spatial points.
  }
  \item{T}{
    Number of points in time.
  }
  \item{NF}{Number of Fourier functions.}

  \item{indCos}{ Vector of integers indicating the position of wavenumbers of cosine-only terms. }

  \item{ns}{ Number of real Fourier functions that have only a cosine and no sine
term. 'ns' is maximal 4.}
  \item{nu}{ Smoothness parameter of the Matern covariance function for the innovations. By
default this equals 1 corresponding to the Whittle covariance function.}
  \item{dt}{ Temporal lag between two time points. By default, this equals 1.}

  \item{logScale}{
    logical; if 'TRUE' the parameters specified in 'logInd' are on the
    logarithmic scale. This is used for constraining parameters to be positive.
  }
  \item{logInd}{
    Vector of integers indicating which parameters are on the log-scale.
  }
  \item{negative}{
    logical; if 'TRUE' the negative log-likelihood is returned otherwise the
    positive log-likelihood is returned.
  }
}



\value{
Value of the log-likelihood evaluated at 'par'.
}

\author{
Fabio Sigrist
}

\examples{
n <- 20
T <- 20
##Specify hyper-parameters
par <- c(rho0=0.1,sigma2=0.2,zeta=0.5,rho1=0.1,gamma=2,alpha=pi/4,muX=0.2,muY=-0.2,tau2=0.01)
##Simulate data
spateSim <- spate.sim(par=par,n=n,T=T,seed=4)
w <- spateSim$w

##Initial values for optim. This takes a couple of seconds.
parI <- c(rho0=0.2,sigma2=0.1,zeta=0.25,rho1=0.01,gamma=1,alpha=0.3,muX=0,muY=0,tau2=0.005)
logInd=c(1,2,3,4,5,9)
##Transform to log-scale
parI[logInd] <- log(parI[logInd])

##Fourier transform needs to be done only once
wFT <- real.fft.TS(w,n=n,T=T)
##ML estimation using optim, takes a couple of seconds
##Load the precomputed object a line below to save time
##spateMLE <- optim(par=parI,loglike,control=list(trace=TRUE,maxit=1000),wFT=wFT,method="L-BFGS-B",
##     lower=c(-10,-10,-10,-10,-10,0,-0.5,-0.5,-10),
##     upper=c(10,10,10,10,10,pi/2,0.5,0.5,10),negative=TRUE,
##     logScale=TRUE,hessian=TRUE,n=n,T=T)
data("spateMLE")

mle <- spateMLE$par
mle[logInd] <- exp(mle[logInd])
sd=sqrt(diag(solve(spateMLE$hessian)))

MleConfInt <- data.frame(array(0,c(4,9)))
colnames(MleConfInt) <- names(par)
rownames(MleConfInt) <- c("True","Estimate","Lower","Upper")
MleConfInt[1,] <- par
MleConfInt[2,] <- mle
MleConfInt[3,] <- spateMLE$par-2*sd
MleConfInt[4,] <- spateMLE$par+2*sd
MleConfInt[c(3,4),logInd] <- exp(MleConfInt[c(3,4),logInd])
cat("\n")
round(MleConfInt,digits=4)

}
