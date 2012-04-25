\name{SemiParSampleSel}
\alias{SemiParSampleSel}

\title{Semiparametric Sample Selection Modelling with Continuous Response}

\description{ 
\code{SemiParSampleSel} can be used to fit continuous response sample selection models where the linear predictors are flexibly specified using parametric and 
regression spline components. During the model fitting process, the possible presence of correlated error equations is accounted for. Regression 
spline bases are extracted from the package \code{mgcv}. Multi-dimensional smooths are available 
via the use of penalized thin plate regression splines (isotropic). The current implementation does not support scale invariant tensor 
product smooths. 

}



\usage{
SemiParSampleSel(formula.eq1, formula.eq2, data=list(),  
                 iterlimSP=50, pr.tol=1e-6, 
                 gamma=1, aut.sp=TRUE, fp=FALSE, start.v=NULL, 
                 rinit=1, rmax=100, fterm=sqrt(.Machine$double.eps), 
                 mterm=sqrt(.Machine$double.eps),   
                 control=list(maxit=50,tol=1e-6,step.half=25,
                              rank.tol=sqrt(.Machine$double.eps)))
}


\arguments{ 

\item{formula.eq1}{A GAM formula for equation 1. \code{s} terms are used to specify smooth smooth functions of 
                  predictors. \code{SemiParSampleSel} supports the use shrinkage smoothers for variable
                  selection purposes. See the examples below and the documentation of \code{mgcv} for further 
                  details on GAM formula specifications. Note that the formula MUST refer to the selection equation.} 

\item{formula.eq2}{A GAM formula for equation 2.} 

\item{data}{An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the
            variables are taken from \code{environment(formula)}, typically the environment from which \code{SemiParSampleSel} is called.} 

\item{iterlimSP}{A positive integer specifying the maximum number of loops to be performed before the smoothing parameter estimation step is terminated.} 

\item{pr.tol}{Tolerance to use in judging convergence of the algorithm when automatic smoothing parameter selection is used.} 

\item{gamma}{It is is an inflation factor for the model degrees of freedom in the UBRE score. Smoother models can be obtained setting 
             this parameter to a value greater than 1. Typically \code{gamma=1.4} achieves this.}

\item{aut.sp}{If \code{TRUE}, then automatic multiple smoothing parameter selection is carried out. If \code{FALSE}, then smoothing parameters 
              are set to the values obtained from the univariate fits.}

\item{fp}{If \code{TRUE}, then a fully parametric model with regression splines if fitted. This only makes sense 
          if used jointly with \code{aut.sp=FALSE}. See the example below.}  

\item{start.v}{Although strictly not recommended, starting values for all model parameters can be provided here. Otherwise these are obtained 
               using the Heckman sample selection correction approach.}  

\item{rinit}{Starting trust region radius. The trust region radius is adjusted as the algorithm proceeds. Bad initial 
             values waste a few steps while the radius is adjusted, but do not prevent the algorithm from working properly. See the documentation
             of \code{trust} for further details.}

\item{rmax}{Maximum allowed trust region radius. This may be set very large. If set small, the algorithm traces a steepest 
            descent path.}        
        
\item{fterm}{Positive scalar giving the tolerance at which the difference in objective function values in a step is considered close 
             enough to zero to terminate the algorithm.}

\item{mterm}{Positive scalar giving the tolerance at which the two-term Taylor-series approximation to the difference in objective 
             function values in a step is considered close enough to zero to terminate the algorithm.}

\item{control}{It is a list containing iteration control constants with the following elements: \code{maxit}: maximum number of iterations of the 
               \code{magic} algorithm; \code{tol}: tolerance to use in judging convergence; \code{step.half}: if a trial step fails then 
               the method tries halving it up to a maximum of \code{step.half} times; \code{rank.tol}: constant used to test for numerical rank 
               deficiency of the problem. See the documentation of \code{magic} in \code{mgcv} for further details.}

}


\details{

The association between the responses is modelled by the correlation parameter \eqn{\rho}{\rho} of a bivariate normal 
distribution. In a semiparametric bivariate sample selection model the linear predictors are flexibly specified using 
parametric components and smooth functions of covariates. Replacing the smooth components with their regression spline expressions yields a fully parametric bivariate sample selection model. In principle, classic 
maximum likelihood estimation can be employed. However, to avoid overfitting, penalized likelihood maximization has to be employed instead. Here the use of 
penalty matrices allows for the suppression of that part of smooth term complexity which has no support from the data. The tradeoff between smoothness 
and fitness is controlled by smoothing parameters associated with the penalty matrices. Smoothing parameters are chosen to 
minimize the approximate Un-Biased Risk Estimator (UBRE) score.

The optimization problem is solved by Newton-Raphson's method. Automatic smoothing parameter selection is integrated 
using a performance-oriented iteration approach (Gu, 1992; Wood, 2004) combined with a `leapfrog' algorithm (Smith, 1996). Roughly 
speaking, at each iteration, (i) the penalized weighted least squares 
problem is solved, then (ii) the smoothing parameters of that problem estimated by approximate UBRE. Steps (i) and (ii) are iterated until 
convergence. Details of the underlying fitting methods are given in Marra and Radice (submitted). 

}


\value{
The function returns an object of class \code{SemiParSampleSel} as described in \code{SemiParSampleSelObject}.
}


\author{

Maintainer: Giampiero Marra \email{giampiero@stats.ucl.ac.uk}

}



\references{

Gu C. (1992), Cross validating non-Gaussian data. \emph{Journal of Computational and Graphical Statistics}, 1(2), 169-179.

Marra G. and Radice R. (submitted), Semiparametric Continuous Response Models for Controlling Sample Selection Bias.

Smith G.K. (1996), Partitioned algorithms for maximum likelihood and other non-linear estimation. \emph{Statistics and Computing}, 6(3), 201-216.

Wood S.N. (2004), Stable and efficient multiple smoothing parameter estimation for generalized additive models. \emph{Journal of the American Statistical Association}, 99(467), 673-686.

}


\section{WARNINGS}{

Any automatic smoothing parameter selection procedure is not likely to work well when the data have low information content. In the current context,
convergence failure may especially occur when \eqn{\rho}{\rho} is high and the total number and selected number of observations is low. If this happens, then one might either (i) lower 
the total number of parameters to estimate by reducing the dimension of the regression spline 
bases, (ii) set the smoothing parameters to the values obtained from the univariate fits (\code{aut.sp=FALSE}), or (iii) set the smoothing parameters to the values obtained from
the non-converged algorithm. The default option is (iii).

Fully parametric modelling is allowed for. However, it is not possible to specify one linear predictor as a function of parametric and smooth
components, and the other as a function of parametric terms only. If continuous covariates are available, then we should let the data determine which
effects are linear or non-linear and for which equations.


} 

\seealso{
\code{\link{InfCr}}, \code{\link{plot.SemiParSampleSel}}, \code{\link{SemiParSampleSel-package}}, \code{\link{SemiParSampleSelObject}}, \code{\link{summary.SemiParSampleSel}} 
}

\examples{

library(SemiParSampleSel)

############
## Generate data
## Correlation between the two equations and covariate correlation 0.5 
## Sample size 2000 

set.seed(0)

n <- 2000

rhC <- rhU <- 0.5      

SigmaU <- matrix(c(1,rhU,rhU,1),2,2)
U      <- rmvnorm(n,rep(0,2),SigmaU)

SigmaC <- matrix( c(1,rhC,rhC,
                    rhC,1,rhC,
                    rhC,rhC,1), 3 , 3)

cov    <- rmvnorm(n,rep(0,3),SigmaC, method="svd")
cov    <- pnorm(cov)

bi <- round(cov[,1]); x1 <- cov[,2]; x2 <- cov[,3]
  
f11 <- function(x) -0.7*(4*x + 2.5*x^2 + 0.7*sin(5*x) + cos(7.5*x))
f12 <- function(x) -0.4*( -0.3 - 1.6*x + sin(5*x))  
f21 <- function(x) 0.6*(exp(x) + sin(2.9*x)) 

ys <-  0.58 + 2.5*bi + f11(x1) + f12(x2) + U[, 1] > 0
y  <- -0.68 - 1.5*bi + f21(x1) +         + U[, 2]
yo <- y*(ys > 0)
  

dataSim <- data.frame(ys,yo,bi,x1,x2)


## CLASSIC SAMPLE SELECTION MODEL
## the first equation must be the selection equation

out <- SemiParSampleSel(ys ~ bi + x1 + x2, 
                        yo ~ bi + x1, 
                        data=dataSim)

summary(out)

InfCr(out)
InfCr(out,cr="BIC")

## SEMIPARAMETRIC SAMPLE SELECTION MODEL
## the first equation MUST be the selection equation

## "cr" cubic regression spline basis      - "cs" shrinkage version of "cr"
## "tp" thin plate regression spline basis - "ts" shrinkage version of "tp"
## for smooths of one variable, "cr/cs" and "tp/ts" achieve similar results 
## k is the basis dimension - default is 10
## m is the order of the penalty for the specific term - default is 2

out <- SemiParSampleSel(ys ~ bi + s(x1,bs="cr",k=10,m=2) + s(x2,bs="cr",k=10), 
                        yo ~ bi + s(x1,bs="cr",k=10), 
                        data=dataSim)

InfCr(out)

## compare the two summary outputs
## the second output produces a summary of the results obtained when only 
## the outcome equation is fitted, i.e. selection bias is not accounted for

summary(out)
summary(out$gam2)


## estimated smooth function plots
## the red line is the true curve
## the blue line is the naive curve not accounting for selection bias

x1 <- sort(x1)
f21.x1 <- f21(x1)[order(x1)]-mean(f21(x1))

plot(out, eq=2, select=1, ylim=c(-1,0.8)); lines(x1, f21.x1, col="red")
par(new=TRUE)
plot(out$gam2, select=1, se=FALSE, col="blue",ylim=c(-1,0.8),ylab="",rug=FALSE)

#
#

}


\keyword{sample selection}
\keyword{semiparametric sample selection modelling} 
\keyword{sample selection model} 
\keyword{smooth} 
\keyword{regression spline}
\keyword{shrinkage smoother}
\keyword{variable selection}







