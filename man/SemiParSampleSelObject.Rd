\name{SemiParSampleSelObject}
\alias{SemiParSampleSelObject}

\title{Fitted SemiParSampleSel object}

\description{A fitted semiparametric bivariate sample selection object returned by function \code{SemiParSampleSel} and of class.}


\value{ 


\item{fit}{A list of values and diagnostics extracted from the output of the algorithm. For instance, \code{fit$argument} and \code{fit$S.h} return
           the estimated parameters and overall penalty matrix scaled by its smoothing parameters, for the bivariate sample selection model. See the 
           documentation of \code{trust} for diagnostics.}

\item{gam1}{A univariate GAM object for equation 1. See the documentation of \code{mgcv} for full details.}

\item{gam2}{A univariate GAM object for equation 2.}

\item{gam2.1}{A univariate GAM object for equation 2, estimated using the Heckman sample selection correction procedure.}

\item{sp}{Estimated smoothing parameters of the smooth components for the fitted bivariate sample selection model.}

\item{iter.sp}{Number of iterations performed for the smoothing parameter estimation step.}

\item{start.v}{Starting values for all model parameters of the semiparametric bivariate sample selection algorithm. These are obtained using the Heckman sample selection 
               correction approach when starting values are not provided.}

\item{sigma}{Estimated standard deviation of the outcome equation error term.}

\item{rho}{Estimated correlation coefficient between the two equations.}

\item{n}{Sample size.}

\item{n.sel}{Number of selected observations.}

\item{X1}{It returns the design matrix associated with the first linear predictor.}

\item{X2}{It returns the design matrix associated with the second linear predictor.} 

\item{X1.d2}{Number of columns of the design matrix for equation 1. This is used for internal calculations.}

\item{X2.d2}{Number of columns of the design matrix for equation 2.} 

\item{l.sp1}{Number of smooth components in equation 1.}

\item{l.sp2}{Number of smooth components in equation 2.}

\item{He}{Penalized hessian.}

\item{HeSh}{Unpenalized hessian.}

\item{Vb}{Inverse of the penalized hessian. This corresponds to the Bayesian variance-covariance matrix used for `confidence' interval calculations.}

\item{F}{This is given by \code{Vb*HeSh}.} 

\item{t.edf}{Total degrees of freedom of the estimated bivariate sample selection model. It is calculated as \code{sum(diag(F))}.}

\item{bs.mgfit}{A list of values and diagnostics extracted from \code{magic}.}

\item{conv.sp}{If \code{TRUE} then the smoothing parameter selection algorithm converged.}
                
\item{wor.c}{It contains the working model quantities given by the square root of 
             the weight matrix times the pseudo-data vector and design matrix, \code{rW.Z} and \code{rW.X}.}                

\item{eta1,eta2}{The estimated linear predictors for the two equations.}

\item{dat}{It returns the full design matrix associated with the response variables and two linear predictors.}


}


\author{

Maintainer: Giampiero Marra \email{giampiero@stats.ucl.ac.uk}

}






\seealso{
\code{\link{SemiParSampleSel}}, \code{\link{plot.SemiParSampleSel}}, \code{\link{summary.SemiParSampleSel}}
}











