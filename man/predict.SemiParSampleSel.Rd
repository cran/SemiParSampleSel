\name{predict.SemiParSampleSel}
\alias{predict.SemiParSampleSel}

\title{Prediction from fitted SemiParSampleSel model}

\description{

It takes a fitted \code{SemiParSampleSel} object produced by \code{SemiParSampleSel()} and produces predictions 
for a new set of values of the model covariates or the original values used for the model fit. 
Standard errors of predictions can be produced and are based on the posterior distribution of the model coefficients.

This function is a wrapper for \code{predict.gam()} in \code{mgcv}. Please see the documentation of \code{mgcv} for full details. 
}

\usage{

\method{predict}{SemiParSampleSel}(object, eq, ...)


}


\arguments{ 

\item{object}{A fitted \code{SemiParSampleSel} object as produced by \code{SemiParSampleSel()}.}

\item{eq}{The equation to be considered for prediction.}

\item{...}{Other arguments as in \code{predict.gam()} in \code{mgcv}.}

}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}


\seealso{
\code{\link{SemiParSampleSel}}, \code{\link{aver}}, \code{\link{plot.SemiParSampleSel}}, \code{\link{summary.SemiParSampleSel}}
}




\keyword{model}
\keyword{prediction}







