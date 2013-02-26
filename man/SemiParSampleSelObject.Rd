\name{SemiParSampleSelObject}
\alias{SemiParSampleSelObject}

\title{Fitted SemiParSampleSel object}

\description{A fitted semiparametric sample selection object returned by function \code{SemiParSampleSel} and of class "SemiParSampleSel".}


\value{ 


\item{fit}{List of values and diagnostics extracted from the output of the algorithm. For instance, \code{fit$gradient} and \code{fit$S.h} return
           the gradient vector and overall penalty matrix scaled by its smoothing parameters, for the bivariate probit model. See the 
           documentation of \code{trust} for details on the diagnostics provided.}

\item{gam1}{Univariate fit for equation 1. See the documentation of \code{mgcv} for full details.}

\item{gam2}{Univariate fit for equation 2.}

\item{gam2.1}{Univariate fit for equation 2, estimated using an adaptation of the Heckman sample selection correction procedure.}

\item{coefficients}{The coefficients of the fitted semiparametric sample selection model. They are given in the following 
                    order: parametric and regression spline (if present)
                    coefficients for the first equation, parametric and regression spline
                    coefficients for the second equation, logarithm of standard deviation of the output (in bivariate normal case) or logarithm of 
                    scale parameter (in normal-gamma case) and transformed association coefficient.}

\item{weights}{Prior weights used during model fitting.}

\item{sp}{Estimated smoothing parameters of the smooth components for the fitted sample selection model.}

\item{iter.if}{Number of iterations performed in the initial step of the Fisher scoring algorithm.}

\item{iter.sp}{Number of iterations performed for the smoothing parameter estimation step.}

\item{iter.fi}{Final number of iterations performed after smoothing parameter/leap frog estimation.}

\item{start.v}{Starting values for all model parameters of the semiparametric sample selection algorithm. These are obtained using the  
              Heckman sample selection 
               correction approach when starting values are not provided.}

\item{phi}{Estimated dispersion for the response of the outcome equation. In the normal bivariate case, this corresponds to the variance. }

\item{sigma}{Estimated standard deviation for the response of the outcome equation in the case of normal marginal distribution of the outcome.}

\item{shape}{Estimated shape parameter for the response of the outcome equation in the case of gamma marginal distribution of the outcome.}

\item{theta}{Estimated coefficient linking the two equations. In the normal bivariate case, this corresponds to the correlation coefficient.}

\item{tau}{Estimated Kendall's tau association coefficient between the two equations.}

\item{n}{Sample size.}

\item{n.sel}{Number of selected observations.}

\item{X1}{Design matrix associated with the first linear predictor.}

\item{X2}{Design matrix associated with the second linear predictor.} 

\item{X1.d2}{Number of columns of the design matrix for equation 1. This is used for internal calculations.}

\item{X2.d2}{Number of columns of the design matrix for equation 2.} 

\item{l.sp1}{Number of smooth components in equation 1.}

\item{l.sp2}{Number of smooth components in equation 2.}

\item{He}{Penalized hessian.}

\item{HeSh}{Unpenalized hessian.}

\item{Vb}{Inverse of the penalized hessian. This corresponds to the Bayesian variance-covariance matrix used for `confidence' interval calculations.}

\item{F}{This is given by \code{Vb*HeSh}.} 

\item{BivD}{Type of bivariate copula distribution employed.} 

\item{margins}{Margins used in the bivariate copula specification.}

\item{t.edf}{Total degrees of freedom of the estimated sample selection model. It is calculated as \code{sum(diag(F))}.}

\item{bs.mgfit}{A list of values and diagnostics extracted from \code{magic} in \code{mgcv}.}

\item{conv.sp}{If \code{TRUE} then the smoothing parameter selection algorithm converged.}
                
\item{wor.c}{Working model quantities given by the square root of 
             the weight matrix times the pseudo-data vector and design matrix, \code{rW.Z} and \code{rW.X}.}                

\item{eta1,eta2}{Estimated linear predictors for the two equations.}

\item{y1}{Binary outcome of the selection equation.}

\item{y2}{Dependent variable of the outcome equation.}

\item{logL}{Value of the (unpenalized) log-likelihood evaluated at the (penalized or unpenalized) parameter 
            estimates.}

\item{fp}{If \code{TRUE}, then a fully parametric model was fitted.} 


}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}






\seealso{
\code{\link{SemiParSampleSel}}, \code{\link{plot.SemiParSampleSel}}, \code{\link{predict.SemiParSampleSel}}, \code{\link{summary.SemiParSampleSel}}
}











