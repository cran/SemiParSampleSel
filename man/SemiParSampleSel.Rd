\name{SemiParSampleSel}
\alias{SemiParSampleSel}

\title{Semiparametric Sample Selection Modelling with Continuous Response}

\description{ 
\code{SemiParSampleSel} can be used to fit continuous response sample selection models where the linear predictors 
are flexibly specified using parametric and 
regression spline components. The depedence between the selection and outcome equations is modelled through the use of copulas. Regression 
spline bases are extracted from the package \code{mgcv}. Multi-dimensional smooths are available 
via the use of penalized thin plate regression splines (isotropic). The current implementation does not support scale invariant tensor 
product smooths. 

}



\usage{
SemiParSampleSel(formula.eq1, formula.eq2, data=list(), weights=NULL,  
                 BivD="N", margins=c("N","N"),
                 start.theta=NULL, start.v=NULL, gamma=1, aut.sp=TRUE, fp=FALSE, 
                 rinit=1, rmax=100, fterm=sqrt(.Machine$double.eps), 
                 mterm=sqrt(.Machine$double.eps),  
                 iterlimsp=50, pr.tolsp=1e-6, 
                 control=list(maxit=50,tol=1e-6,step.half=25,
                              rank.tol=sqrt(.Machine$double.eps)))
}


\arguments{ 

\item{formula.eq1}{A formula for equation 1, i.e. the selection equation. \code{s} terms are used to specify smooth functions of 
                  predictors. \code{SemiParSampleSel} supports the use shrinkage smoothers for variable
                  selection purposes. See the examples below and the documentation of \code{mgcv} for further 
                  details on formula specifications. Note that \code{formula.eq1} MUST refer to the selection equation.} 

\item{formula.eq2}{A formula for equation 2, i.e. the outcome equation.} 

\item{data}{An optional data frame, list or environment containing the variables in the model. If not found in \code{data}, the
            variables are taken from \code{environment(formula)}, typically the environment from which \code{SemiParSampleSel} is called.} 

\item{weights}{Optional vector of prior weights to be used in fitting.} 


\item{BivD}{Type of bivariate copula employed. Possible choices are "N", "C", "J", "FGM", "F", "AMH" and "G" which stand for bivariate normal, Clayton, 
	   Joe, Farlie-Gumbel-Morgenstern, Frank, Ali-Mikhail-Haq and Gumbel.} 

\item{margins}{A two-dimensional vector which specifies the marginal distributions of the selection and outcome equations. Possible choices 
               are c("N","N") for normal margins and c("N","G") for normal and gamma margins.}

\item{start.theta}{A starting value for the association parameter of the copula given in BivD.}

\item{start.v}{Starting values for all model parameters can be provided here; these have to be given in the following order: parametric and 
               regression spline (if present) coefficients for the first equation, parametric and regression spline
                    coefficients for the second equation, dispersion parameter and association coefficient. Otherwise, these are obtained 
               using an adaptation of the two-stage Heckman sample selection correction approach.}  

\item{gamma}{Inflation factor for the model degrees of freedom in the UBRE score. Smoother models can be obtained setting 
             this parameter to a value greater than 1. Typically \code{gamma=1.4} achieves this.}

\item{aut.sp}{If \code{TRUE}, then automatic multiple smoothing parameter selection is carried out. If \code{FALSE}, then smoothing parameters 
              are set to the values obtained from the univariate fits.}

\item{fp}{If \code{TRUE}, then a fully parametric model with regression splines if fitted.}  


\item{rinit}{Starting trust region radius. The trust region radius is adjusted as the algorithm proceeds. See the documentation
             of \code{trust} for further details.}

\item{rmax}{Maximum allowed trust region radius. This may be set very large. If set small, the algorithm traces a steepest 
            descent path.}        
        
\item{fterm}{Positive scalar giving the tolerance at which the difference in objective function values in a step is considered close 
             enough to zero to terminate the algorithm.}

\item{mterm}{Positive scalar giving the tolerance at which the two-term Taylor-series approximation to the difference in objective 
             function values in a step is considered close enough to zero to terminate the algorithm.}

\item{iterlimsp}{A positive integer specifying the maximum number of loops to be performed before the smoothing parameter estimation step is terminated.} 

\item{pr.tolsp}{Tolerance to use in judging convergence of the algorithm when automatic smoothing parameter selection is used.} 

\item{control}{List containing iteration control constants for the automatic smoothing parameter selection method. These are: \code{maxit}: maximum 
               number of iterations of the 
               \code{magic} algorithm; \code{tol}: tolerance to use in judging convergence; \code{step.half}: if a trial step fails then 
               the method tries halving it up to a maximum of \code{step.half} times; \code{rank.tol}: constant used to test for numerical rank 
               deficiency of the problem. See the documentation of \code{magic} in \code{mgcv} for further details.}

}


\details{

The association between the responses is modelled by parameter \eqn{\theta}{\theta}. In a semiparametric bivariate sample selection model 
the linear predictors are flexibly specified using 
parametric components and smooth functions of covariates. Replacing the smooth components with their regression spline expressions yields a 
fully parametric bivariate sample selection model. In principle, classic 
maximum likelihood estimation can be employed. However, to avoid overfitting, penalized likelihood maximization is used instead. Here the use of 
penalty matrices allows for the suppression of that part of smooth term complexity which has no support from the data. The tradeoff between smoothness 
and fitness is controlled by smoothing parameters associated with the penalty matrices. Smoothing parameters are chosen to 
minimize the approximate Un-Biased Risk Estimator (UBRE) score.

The optimization problem is solved by Newton-Raphson's method. Automatic smoothing parameter selection is integrated 
using a performance-oriented iteration approach (Gu, 1992; Wood, 2004) combined with a `leapfrog' algorithm (Smith, 1996). Roughly 
speaking, at each iteration, (i) the penalized weighted least squares 
problem is solved, and (ii) the smoothing parameters of that problem estimated by approximate UBRE. Steps (i) and (ii) are iterated until 
convergence. Details of the underlying fitting methods are given in Marra and Radice (in press) and Wojtys et. al (submitted). 

}


\value{
The function returns an object of class \code{SemiParSampleSel} as described in \code{SemiParSampleSelObject}.
}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\references{

Gu C. (1992), Cross validating non-Gaussian data. \emph{Journal of Computational and Graphical Statistics}, 1(2), 169-179.

Marra G. and Radice R. (in press), Estimation of a Regression Spline Sample Selection Model. \emph{Computational Statistics and Data Analysis}.

Smith G.K. (1996), Partitioned algorithms for maximum likelihood and other non-linear estimation. \emph{Statistics and Computing}, 6(3), 201-216.

Wojtys M., Marra G. and Radice R. (submitted), Copula Regression Spline Sample Selection Models: The R Package SemiParSampleSel.

Wood S.N. (2004), Stable and efficient multiple smoothing parameter estimation for generalized additive models. \emph{Journal of the American Statistical Association}, 99(467), 673-686.

}


\section{WARNINGS}{

In the current context, convergence failure may occur when \eqn{\theta}{\theta} is very high, the total number and selected number of 
observations is low, and the data have low information content. 


} 

\seealso{
\code{\link{InfCr}}, \code{\link{plot.SemiParSampleSel}}, \code{\link{SemiParSampleSel-package}}, \code{\link{SemiParSampleSelObject}}, \code{\link{ss.checks}}, \code{\link{predict.SemiParSampleSel}}, \code{\link{summary.SemiParSampleSel}} 
}

\examples{

library(SemiParSampleSel)

######################################################################
## Generate data
## Correlation between the two equations and covariate correlation 0.5 
## Sample size 2000 
######################################################################

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
ss.checks(out)
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

x1.s <- sort(x1[dataSim$ys>0])
f21.x1 <- f21(x1.s)[order(x1.s)]-mean(f21(x1.s))

plot(out, eq=2, ylim=c(-1,0.8)); lines(x1.s, f21.x1, col="red")
par(new=TRUE)
plot(out$gam2, se=FALSE, col="blue", ylim=c(-1,0.8), ylab="",rug=FALSE)

#
#

###################################################
## example using Clayton copula with normal margins
###################################################

\dontrun{ 

set.seed(0)

teta <- 10
sig  <- 1.5

myCop <- archmCopula(family = "clayton", dim = 2, param = teta)

# other copula options are for instance: "amh", "frank", "gumbel", "joe"
# for FGM use the following code:
# myCop <- fgmCopula(teta, dim=2)

bivg  <- mvdc(copula = myCop, c("norm","norm"), list(list(mean=0,sd=1), list(mean=0,sd=sig)))
er    <- rMvdc(n,bivg)

ys <-  0.58 + 2.5*bi + f11(x1) + f12(x2) + er[, 1] > 0
y  <- -0.68 - 1.5*bi + f21(x1) +         + er[, 2]
yo <- y*(ys > 0)
  
dataSim <- data.frame(ys,yo,bi,x1,x2)

out <- SemiParSampleSel(ys ~ bi + s(x1) + s(x2), 
                        yo ~ bi + s(x1), 
                        data=dataSim, BivD="C")
ss.checks(out)
summary(out)

x1.s <- sort(x1[dataSim$ys>0])
f21.x1 <- f21(x1.s)[order(x1.s)]-mean(f21(x1.s))

plot(out, eq=2, ylim=c(-1.1,1.6)); lines(x1.s, f21.x1, col="red")
par(new=TRUE)
plot(out$gam2, se=FALSE, col="blue", ylim=c(-1.1,1.6), ylab="",rug=FALSE)

#
#

}

}

\keyword{copula}
\keyword{sample selection}
\keyword{semiparametric sample selection modelling} 
\keyword{sample selection model} 
\keyword{smooth} 
\keyword{regression spline}
\keyword{shrinkage smoother}
\keyword{variable selection}







