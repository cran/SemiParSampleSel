\name{SemiParSampleSel}
\alias{SemiParSampleSel}

\title{Semiparametric Sample Selection Modelling with Continuous or Discrete Response}

\description{ 
\code{SemiParSampleSel} can be used to fit continuous or discrete response sample selection models where the linear predictors 
are flexibly specified using parametric and 
regression spline components. The depedence between the selection and outcome equations is modelled through the use of copulas. Regression 
spline bases are extracted from the package \code{mgcv}. Multi-dimensional smooths are available 
via the use of penalized thin plate regression splines. If it makes sense, the dependence parameter of 
the chosen bivariate distribution as well as the shape and dispersion parameters of the outcome distribution can be specified as 
functions of semiparametric predictors.

}


\usage{

SemiParSampleSel(formula, data = list(), weights = NULL, subset = NULL, 
                  start.v = NULL, start.theta = NULL,
                  BivD = "N", margins = c("probit","N"), fp = FALSE, infl.fac = 1,  
                  rinit = 1, rmax = 100, iterlimsp = 50, pr.tolsp = 1e-6, bd = NULL,
                  parscale)
}


\arguments{ 

\item{formula}{A list of two formulas, one for selection equation and the other for the outcome equation. \code{s} terms 
               are used to specify smooth smooth functions of 
                  predictors. \code{SemiParSampleSel} supports the use shrinkage smoothers for variable
                  selection purposes. See the examples below and the documentation of \code{mgcv} for further 
                  details on formula specifications. Note that the first formula MUST refer to the selection equation.
                  Furthermore,
                  if it makes sense, more equations can be specified for the other model parameters (see Example 1 below).} 

\item{data}{An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the
            variables are taken from \code{environment(formula)}, typically the environment from which \code{SemiParSampleSel} is called.} 

\item{weights}{Optional vector of prior weights to be used in fitting.} 

\item{subset}{Optional vector specifying a subset of observations to be used in the fitting process.}

\item{start.v}{Starting values for all model parameters can be provided here. Otherwise, these are obtained using an 
               adaptation of the two-stage Heckman sample selection correction approach.}  

\item{start.theta}{A starting value for the association parameter of the copula given in BivD.}

\item{BivD}{Type of bivariate error distribution employed. Possible choices are "N", "C0", "C90", "C180", "C270", "J0", "J90", "J180", "J270", 
            "G0", "G90", "G180", "G270", "F", "FGM" and "AMH" which stand for bivariate normal, Clayton, rotated Clayton (90 degrees), 
             survival Clayton,
             rotated Clayton (270 degrees), Joe, rotated Joe (90 degrees), survival Joe, rotated Joe (270 degrees),
             Gumbel, rotated Gumbel (90 degrees), survival Gumbel, rotated Gumbel (270 degrees), Frank, 
             Farlie-Gumbel-Morgenstern, and Ali-Mikhail-Haq.}

\item{margins}{A two-dimensional vector which specifies the marginal distributions of the selection and outcome 
               equations. The first margin currently admits only "probit" or equivalently "N".
               The second margin can be "N", "GA", "P", "NB", "D", "PIG", "S", "BB", "BI", "GEOM", "LG", 
               "NBII", "WARING", "YULE", "ZIBB", "ZABB", "ZABI", "ZIBI", "ZALG", "ZANBI", "ZINBI", "ZAP", 
               "ZIP", "ZIP2", "ZIPIG" which stand for normal, gamma, Poisson, negative binomial type I, Delaporte, 
               Poisson inverse Gaussian, Sichel, beta binomial, binomial, geometric, logarithmic, negative binomial
               type II, Waring, Yule, zero inflated beta binomial, zero altered beta binomial, zero altered binomial,
               zero inflated binomial, zero altered logarithmic, zero altered negative binomial type I, zero inflated
               negative binomial type I, zero altered Poisson, zero inflated Poisson, zero inflated Poisson type II
               and zero inflated Poisson inverse Gaussian.}
               
\item{fp}{If \code{TRUE}, then a fully parametric model with regression splines if fitted.}  

\item{infl.fac}{Inflation factor for the model degrees of freedom in the UBRE score. Smoother models can be obtained setting 
             this parameter to a value greater than 1.}

\item{rinit}{Starting trust region radius. The trust region radius is adjusted as the algorithm proceeds. See the documentation
             of \code{trust} for further details.}

\item{rmax}{Maximum allowed trust region radius. This may be set very large. If set small, the algorithm traces a steepest 
            descent path.}        
       
\item{iterlimsp}{A positive integer specifying the maximum number of loops to be performed before the smoothing parameter estimation step 
                 is terminated.} 

\item{pr.tolsp}{Tolerance to use in judging convergence of the algorithm when automatic smoothing parameter estimation is used.} 

\item{bd}{Binomial denominator. To be used in the case of "BB", "BI", "ZIBB", "ZABB", "ZABI", "ZIBI".}

\item{parscale}{The algorithm will operate as if optimizing objfun(x / parscale, ...). If missing then no rescaling is done. See the 
documentation of \code{trust} for more details.}   

}


\details{

The association between the responses is modelled by parameter \eqn{\rho}{\rho} or \eqn{\theta}{\theta}. In a semiparametric 
bivariate sample selection model 
the linear predictors are flexibly specified using 
parametric components and smooth functions of covariates. Replacing the smooth components with their regression spline expressions yields a 
fully parametric bivariate sample selection model. In principle, classic 
maximum likelihood estimation can be employed. However, to avoid overfitting, penalized likelihood maximization is used instead. Here the use of 
penalty matrices allows for the suppression of that part of smooth term complexity which has no support from the data. The tradeoff 
between smoothness 
and fitness is controlled by smoothing parameters associated with the penalty matrices. Smoothing parameters are chosen to 
minimize the approximate Un-Biased Risk Estimator (UBRE) score, which can also be viewed as an approximate AIC.

The optimization problem is solved by a trust region algorithm. Automatic smoothing parameter selection is integrated 
using a performance-oriented iteration approach (Gu, 1992; Wood, 2004). Roughly 
speaking, at each iteration, (i) the penalized weighted least squares 
problem is solved, and (ii) the smoothing parameters of that problem estimated by approximate UBRE. Steps (i) and (ii) are 
iterated until convergence. Details of the underlying fitting methods are given in Marra and Radice (2013) and Wojtys et. al (in press). 


}


\value{
The function returns an object of class \code{SemiParSampleSel} as described in \code{SemiParSampleSelObject}.
}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\references{

Gu C. (1992), Cross validating non-Gaussian data. \emph{Journal of Computational and Graphical Statistics}, 1(2), 169-179.

Marra G. and Radice R. (2013), Estimation of a Regression Spline Sample Selection Model. \emph{Computational Statistics and Data Analysis}, 61, 158-173.

Wojtys M., Marra G. and Radice R. (in press), Copula Regression Spline Sample Selection Models: The R Package SemiParSampleSel. \emph{Journal of Statistical Software}.

Wood S.N. (2004), Stable and efficient multiple smoothing parameter estimation for generalized additive models. \emph{Journal of the American Statistical Association}, 99(467), 673-686.

}


\section{WARNINGS}{

Convergence failure may occur when \eqn{\rho}{\rho} or \eqn{\theta}{\theta} is very high, and/or the total number and selected number of 
observations are low, and/or there are important mistakes in the model specification (i.e., using C90 when the model equations 
are positively associated), and/or there are many smooth components in the model as compared to the number of observations. Convergence 
failure may also mean that an infinite cycling between steps (i) and (ii) occurs. In this case, the smoothing parameters are set to the values 
obtained from the non-converged algorithm (\code{\link{conv.check}} will give a warning). In such cases, we 
recommend re-specifying the model, and/or using some rescaling (see \code{parscale}).

In the context of non-random sample selection, it would not make much sense to specify the dependence
parameter as function of covariates. This is because the assumption is that 
the dependence parameter models the 
association between the unobserved confounders in the two equations. However, this option does make sense when it is believed that the  
association coefficient varies across geographical areas, for instance. 


} 

\seealso{
\code{\link{aver}}, \code{\link{plot.SemiParSampleSel}}, \code{\link{SemiParSampleSel-package}}, \code{\link{SemiParSampleSelObject}}, \code{\link{conv.check}}, \code{\link{predict.SemiParSampleSel}}, \code{\link{summary.SemiParSampleSel}} 
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

SigmaU <- matrix(c(1, rhU, rhU, 1), 2, 2)
U      <- rmvnorm(n, rep(0,2), SigmaU)

SigmaC <- matrix(rhC, 3, 3); diag(SigmaC) <- 1

cov    <- rmvnorm(n, rep(0,3), SigmaC, method = "svd")
cov    <- pnorm(cov)

bi <- round(cov[,1]); x1 <- cov[,2]; x2 <- cov[,3]
  
f11 <- function(x) -0.7*(4*x + 2.5*x^2 + 0.7*sin(5*x) + cos(7.5*x))
f12 <- function(x) -0.4*( -0.3 - 1.6*x + sin(5*x))  
f21 <- function(x) 0.6*(exp(x) + sin(2.9*x)) 

ys <-  0.58 + 2.5*bi + f11(x1) + f12(x2) + U[, 1] > 0
y  <- -0.68 - 1.5*bi + f21(x1) +         + U[, 2]
yo <- y*(ys > 0)
  

dataSim <- data.frame(ys, yo, bi, x1, x2)

## CLASSIC SAMPLE SELECTION MODEL
## the first equation MUST be the selection equation

out <- SemiParSampleSel(list(ys ~ bi + x1 + x2, 
                             yo ~ bi + x1), 
                        data = dataSim)
conv.check(out)
summary(out)

AIC(out)
BIC(out)
aver(out)


## SEMIPARAMETRIC SAMPLE SELECTION MODEL

## "cr" cubic regression spline basis      - "cs" shrinkage version of "cr"
## "tp" thin plate regression spline basis - "ts" shrinkage version of "tp"
## for smooths of one variable, "cr/cs" and "tp/ts" achieve similar results 
## k is the basis dimension - default is 10
## m is the order of the penalty for the specific term - default is 2

out <- SemiParSampleSel(list(ys ~ bi + s(x1, bs = "tp", k = 10, m = 2) + s(x2), 
                             yo ~ bi + s(x1)), 
                        data = dataSim)
conv.check(out)                        
AIC(out)
aver(out)

## compare the two summary outputs
## the second output produces a summary of the results obtained when only 
## the outcome equation is fitted, i.e. selection bias is not accounted for

summary(out)
summary(out$gam2)

## estimated smooth function plots
## the red line is the true curve
## the blue line is the naive curve not accounting for selection bias

x1.s <- sort(x1[dataSim$ys>0])
f21.x1 <- f21(x1.s)[order(x1.s)] - mean(f21(x1.s))

plot(out, eq = 2, ylim = c(-1, 0.8)); lines(x1.s, f21.x1, col = "red")
par(new = TRUE)
plot(out$gam2, se = FALSE, col = "blue", ylim = c(-1, 0.8), ylab = "", rug = FALSE)




\dontrun{ 

## SEMIPARAMETRIC SAMPLE SELECTION MODEL with association and dispersion parameters 
## depending on covariates as well

out <- SemiParSampleSel(list(ys ~ bi + s(x1) + s(x2), 
                             yo ~ bi + s(x1),
                                ~ bi, 
                                ~ bi + x1), 
                        data = dataSim)
conv.check(out)                        
summary(out)
out$sigma
out$theta

#
#

###################################################
## example using Clayton copula with normal margins
###################################################

set.seed(0)

theta <- 5
sig  <- 1.5

myCop <- archmCopula(family = "clayton", dim = 2, param = theta)

# other copula options are for instance: "amh", "frank", "gumbel", "joe"
# for FGM use the following code:
# myCop <- fgmCopula(theta, dim=2)

bivg  <- mvdc(copula = myCop, c("norm", "norm"), 
              list(list(mean = 0, sd = 1), 
              list(mean = 0, sd = sig)))
er    <- rMvdc(n, bivg)

ys <-  0.58 + 2.5*bi + f11(x1) + f12(x2) + er[, 1] > 0
y  <- -0.68 - 1.5*bi + f21(x1) +         + er[, 2]
yo <- y*(ys > 0)
  
dataSim <- data.frame(ys, yo, bi, x1, x2)

out <- SemiParSampleSel(list(ys ~ bi + s(x1) + s(x2), 
                             yo ~ bi + s(x1)), 
                        data = dataSim, BivD = "C0")
conv.check(out)
summary(out)
aver(out)

x1.s <- sort(x1[dataSim$ys>0])
f21.x1 <- f21(x1.s)[order(x1.s)] - mean(f21(x1.s))

plot(out, eq = 2, ylim = c(-1.1, 1.6)); lines(x1.s, f21.x1, col = "red")
par(new = TRUE)
plot(out$gam2, se = FALSE, col = "blue", ylim = c(-1.1, 1.6), ylab = "", rug = FALSE)

#
#

########################################################
## example using Gumbel copula with normal-gamma margins
########################################################

set.seed(0)

k <- 2                                # shape of gamma distribution
miu  <- exp(-0.68 - 1.5*bi + f21(x1)) # mean values of y's (log m = Xb)
lambda <- k/miu	                      # rate of gamma distribution

theta <- 6

# Two-dimensional Gumbel copula with unif margins
gumbel.cop <- onacopula("Gumbel", C(theta, 1:2))

# Random sample from two-dimensional Gumbel copula with uniform margins
U <- rnacopula(n = n, gumbel.cop)		  

# Margins: normal and gamma
er <- cbind(qnorm(U[,1], 0, 1), qgamma(U[, 2], shape = k, rate = lambda))

ys <- 0.58 + 2.5*bi + f11(x1) + f12(x2) + er[, 1] > 0
y  <- er[, 2]
yo <- y*(ys > 0)

dataSim <- data.frame(ys, yo, bi, x1, x2)

out <- SemiParSampleSel(list(ys ~ bi + s(x1) + s(x2), 
                             yo ~ bi + s(x1)), 
                        data = dataSim, BivD = "G0", margins = c("N", "G"))
conv.check(out)
summary(out)
aver(out)

x1.s <- sort(x1[dataSim$ys>0])
f21.x1 <- f21(x1.s)[order(x1.s)] - mean(f21(x1.s))

plot(out, eq = 2, ylim = c(-1.1, 1)); lines(x1.s, f21.x1, col = "red")
par(new = TRUE)
plot(out$gam2, se = FALSE, col = "blue", ylim = c(-1.1, 1), ylab = "", rug = FALSE)

#
#


########################################################
## Example for discrete margins and normal copula
########################################################


# Creating simulation function
bcds <- function(n, s.tau = 0.2, s.sigma = 1, s.nu = 0.5, 
                rhC = 0.2, outcome.margin = "PO", copula = "FGM")  {
  
# Generating covariates     
  SigmaC <- matrix( c(1,rhC,rhC,rhC,rhC,1,rhC,rhC,rhC,rhC,1,rhC,rhC,rhC,rhC,1), 4 , 4)
  covariates    <- rmvnorm(n,rep(0,4),SigmaC, method="svd")
  covariates    <- pnorm(covariates)
  x1 <- covariates[,1]
  x2 <- covariates[,2]
  x3 <- round(covariates[,3])
  x4 <- round(covariates[,4])
  
  # Establishing copula object
  if (copula == "FGM") {
      Cop <- fgmCopula(dim = 2, param = iTau(fgmCopula(), s.tau))
    } else if (copula == "N") {
      Cop <- ellipCopula(family = "normal", dim = 2, param = iTau(normalCopula(), s.tau))
    } else if (copula == "AMH") {
      Cop <- archmCopula(family = "amh", dim = 2, param = iTau(amhCopula(), s.tau))
    } else if (copula == "C0") {
      Cop <- archmCopula(family = "clayton", dim = 2, param = iTau(claytonCopula(), s.tau))
    } else if (copula == "F") {
      Cop <- archmCopula(family = "frank", dim = 2, param = iTau(frankCopula(), s.tau))
    } else if (copula == "G0") {
      Cop <- archmCopula(family = "gumbel", dim = 2, param = iTau(gumbelCopula(), s.tau))
    } else if (copula == "J0") {
      Cop <- archmCopula(family = "joe", dim = 2, param = iTau(joeCopula(), s.tau))
    }  
  
  # Setting up equations
  f1 <- function(x) 0.4*(-4 - (5.5*x-2.9) + 3*(4.5*x-2.3)^2 - (4.5*x-2.3)^3)
  f2 <- function(x) x*sin(8*x)
  mu_s  <- 1.0 + f1(x1) - 2.0*x2 + 3.1*x3 - 2.2*x4
  mu_o  <-  exp(1.3 + f2(x1) - 1.9*x2 + 2.4*x3 - 0.1*x4)
  
  # Creating margin dependent object
  if (outcome.margin == "P") {
        speclist <- list(mu = mu_o)
        outcome.margin2 <- "PO"
     } else if (outcome.margin == "NB") {
        speclist  <- list(mu = mu_o, sigma = s.sigma)
        outcome.margin2 <- "NBI"
     } else if (outcome.margin == "D") {
        speclist  <- list(mu = mu_o, sigma = s.sigma, nu = s.nu)
        outcome.margin2 <- "DEL"
     } else if (outcome.margin == "PIG") {
        speclist  <- list(mu = mu_o, sigma = s.sigma)
        outcome.margin2 <- "PIG"
     } else if (outcome.margin == "S") {
        speclist  <- list(mu = mu_o, sigma = s.sigma, nu = s.nu)
        outcome.margin2 <- "SICHEL"
   }
  spec  <- mvdc(copula = Cop, c("norm", outcome.margin2), 
                list(list(mean = mu_s, sd = 1), speclist))  
  
  # Simulating data
  simGen <- rMvdc(n, spec)
  y <- ifelse(simGen[,1]>0, simGen[,2], -99)
  
  dataSim <- data.frame(y, x1, x2, x3, x4)  
  dataSim  
}



# Creating plots of the true functional form of x1 in both equations
xt1  <- seq(0, 1, length.out=200)
xt2  <- seq(0,1, length.out=200)
f1t <- function(x) 0.4*(-4 - (5.5*x-2.9) + 3*(4.5*x-2.3)^2 - (4.5*x-2.3)^3)
f2t <- function(x) x*sin(8*x)
plot(xt1, f1t(xt1))
plot(xt2, f2t(xt2))




# Simulating 1000 deviates

set.seed(0)

dataSim<- bcds(1000, s.tau = 0.6, s.sigma = 0.1, s.nu = 0.5, 
                rhC = 0.5, outcome.margin = "NB", copula = "N")
dataSim$y.probit<-ifelse(dataSim$y >= 0, 1, 0)


# Estimating SemiParSampleSel

out1 <- SemiParSampleSel(list(y.probit ~ s(x1) + x2 + x3 + x4, y ~ s(x1) + x2 + x3 + x4), 
                          data = dataSim, BivD = "N", margins = c("N", "P")) 
conv.check(out1)

out2 <- SemiParSampleSel(list(y.probit ~ s(x1) + x2 + x3 + x4, y ~ s(x1) + x2 + x3 + x4), 
                         data = dataSim, BivD = "N", margins = c("N", "NB")) 
conv.check(out2)


# Model comparison

AIC(out1)
AIC(out2)
VuongClarke(out1, out2)


# Model diagnostics

summary(out2, cm.plot = TRUE)
plot(out2, eq = 1)
plot(out2, eq = 2)
aver(out2, univariate = TRUE)
aver(out2, univariate = FALSE)

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







