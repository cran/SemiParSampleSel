\name{SemiParSampleSel-package}
\alias{SemiParSampleSel-package}
\docType{package}

\title{Semiparametric Sample Selection Modelling with Continuous or Discrete Response}

\description{\code{SemiParSampleSel} provides a function for fitting continuous and discrete response (copula) sample selection models 
with parametric and nonparametric predictor effects. Several bivariate copula distributions are supported. The dependence parameter of 
the copula distribution as well as the shape and dispersion parameters of the outcome distribution can be specified as 
functions of semiparametric predictors as well. Smoothness selection is
achieved automatically and interval calculations are based on a Bayesian approach.
}


\details{\code{SemiParSampleSel} provides a function for flexible sample selection modelling with continuous or discrete response. The underlying 
representation and 
estimation of the model is based on a penalized regression spline approach, with automatic smoothness selection. The numerical routine carries 
out function minimization using a trust region algorithm from the package \code{trust} in combination 
with an adaptation of a low level smoothness selection fitting procedure from the package \code{mgcv}. 
 
\code{\link{SemiParSampleSel}} supports the use of many smoothers as extracted from \code{mgcv}. Scale invariant tensor product smooths
are not currently supported. Estimation is by penalized maximum likelihood with automatic smoothness selection by approximate Un-Biased 
Risk Estimator (UBRE) score, which can also be viewed as an approximate AIC. 
The depedence between the selection and outcome equations is modelled through the use of copulas.

Confidence intervals for smooth components and nonlinear functions of the model
parameters are derived using a Bayesian approach. Approximate p-values for testing 
individual smooth terms for equality to the zero function are also provided and based on the approach
implemented in \code{mgcv}. Functions \code{\link{plot.SemiParSampleSel}} and 
\code{\link{summary.SemiParSampleSel}} extract such information from a fitted \code{\link{SemiParSampleSel}} object. Model/variable 
selection is also possible via the use of shrinakge smoothers or information criteria. 

Function \code{\link{aver}} calculates the average outcome corrected for non-random sample selection.

If it makes sense, the dependence parameter of 
the copula function as well as the shape and dispersion parameters of the outcome distribution can be specified as 
functions of semiparametric predictors.

}


\author{

Giampiero Marra (University College London, Department of Statistical Science), Rosalba Radice (Birkbeck, University of London, Department of Economics, Mathematics and Statistics), Malgorzata Wojtys (University of Plymouth, School of Computing and Mathematics), Karol Wyszynski (University College London, Department of Statistical Science)

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}

\references{

Marra G. and Radice R. (2013), Estimation of a Regression Spline Sample Selection Model. \emph{Computational Statistics and Data Analysis}, 61, 158-173.

Wojtys M., Marra G. and Radice R. (in press), Copula Regression Spline Sample Selection Models: The R Package SemiParSampleSel. \emph{Journal of Statistical Software}.

}


\seealso{
\code{\link{SemiParSampleSel}}
}


\keyword{package}
\keyword{sample selection} 
\keyword{semiparametric sample selection modelling} 
\keyword{smooth} 
\keyword{regression}
\keyword{variable selection}

