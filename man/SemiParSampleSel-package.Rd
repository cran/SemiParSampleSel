\name{SemiParSampleSel-package}
\alias{SemiParSampleSel-package}
\docType{package}

\title{Semiparametric Sample Selection Modelling with Continuous Response}

\description{\code{SemiParSampleSel} provides a function for fitting continuous response (copula) sample selection models 
with parametric and nonparametric predictor effects.
}

\details{\code{SemiParSampleSel} provides a function for flexible sample selection modelling with continuous (normal or gamma) response. The underlying 
representation and 
estimation of the model is based on a penalized regression spline approach, with automatic smoothness selection. The numerical routine carries 
out function minimization using a trust region algorithm from the package \code{trust} in combination 
with an adaptation of a low level smoothness selection fitting procedure from the package \code{mgcv}, via a `leapfrog' algorithm. 
 
\code{\link{SemiParSampleSel}} supports the use of many smoothers as extracted from \code{mgcv}. Scale invariant tensor product smooths
are not currently supported. Estimation is by penalized maximum likelihood with automatic smoothness selection by approximate Un-Biased 
Risk Estimator (UBRE) score. 
The depedence between the selection and outcome equations is modelled through the use of copulas.

Confidence intervals for smooth components are derived using a Bayesian approach. P-values for testing 
individual smooth terms for equality to the zero function are also provided. Functions \code{\link{plot.SemiParSampleSel}} and 
\code{\link{summary.SemiParSampleSel}} extract such information from a fitted \code{\link{SemiParSampleSelObject}}. Variable selection is also 
possible via the use of shrinakge smoothers or information criteria. 

}


\author{

Giampiero Marra (University College London, Department of Statistical Science), Rosalba Radice (Birkbeck, University of London, Department of Economics, Mathematics and Statistics) and Malgorzata Wojtys (University College London, Department of Statistical Science)

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}

\references{

Marra G. and Radice R. (in press), Estimation of a Regression Spline Sample Selection Model. \emph{Computational Statistics and Data Analysis}.

Wojtys M., Marra G. and Radice R. (submitted), Copula Regression Spline Sample Selection Models: The R Package SemiParSampleSel.

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

