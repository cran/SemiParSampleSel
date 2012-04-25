\name{SemiParSampleSel-package}
\alias{SemiParSampleSel-package}
\docType{package}

\title{Semiparametric Sample Selection Modelling with Continuous Response}

\description{\code{SemiParSampleSel} provides a function for fitting continuous response sample selection models with semiparametric predictors, including linear and nonlinear effects.
}

\details{\code{SemiParSampleSel} provides a function for flexible sample selection modelling with continuous response. The underlying representation and 
estimation of the model is based on a penalized regression spline approach, with automatic smoothness selection. The numerical routine carries 
out function minimization using a trust region algorithm from the package \code{trust} in combination 
with an adaptation of a low level smoothness selection fitting procedure from the package \code{mgcv} combined with a `leapfrog' algorithm. 
 
\code{\link{SemiParSampleSel}} supports the use of many smoothers as extracted from \code{mgcv}. Scale invariant tensor product smooths
are not currently supported. Estimation is by penalized maximum likelihood with automatic smoothness selection by approximate Un-Biased Risk Estimator (UBRE) score. 

Confidence intervals for smooth components are derived using a Bayesian approach. Approximate p-values for testing 
individual smooth terms for equality to the zero function are also provided. Functions \code{\link{plot.SemiParSampleSel}} and 
\code{\link{summary.SemiParSampleSel}} extract such information from a fitted \code{\link{SemiParSampleSelObject}}. Model/variable selection is also 
possible via the use of shrinakge smoothers or information criteria. 

}


\author{

Giampiero Marra (University College London) and Rosalba Radice (London School of Hygiene and Tropical Medicine)

Maintainer: Giampiero Marra \email{giampiero@stats.ucl.ac.uk}
}

\references{

Marra G. and Radice R. (submitted), Semiparametric Continuous Response Models for Controlling Sample Selection Bias.

}


\seealso{
\code{\link{SemiParSampleSel}}
}


\keyword{package}
\keyword{sample selection} 
\keyword{semiparametric sample selection modelling} 
\keyword{smooth} 
\keyword{regression}
\keyword{model selection}

