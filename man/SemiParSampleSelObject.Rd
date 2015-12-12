\name{SemiParSampleSelObject}
\alias{SemiParSampleSelObject}

\title{Fitted SemiParSampleSel object}

\description{A fitted semiparametric sample selection object returned by function \code{SemiParSampleSel} and of class "SemiParSampleSel".}


\value{ 


\item{fit}{List of values and diagnostics extracted from the output of the algorithm. For instance, \code{fit$gradient} and \code{fit$S.h} return
           the gradient vector and overall penalty matrix scaled by its smoothing parameters, for the bivariate model. See the 
           documentation of \code{trust} for details on the diagnostics provided.}

\item{gam1}{Univariate fit for selection equation. See the documentation of \code{mgcv} for full details.}

\item{gam2,gam3,gam4,gam5}{Univariate fit for the outcome equation and equations 3, 4 and 5 when present.}

\item{gam2.1}{Univariate fit for equation 2, estimated using an adaptation of the Heckman sample selection correction procedure.}

\item{coefficients}{The coefficients of the fitted semiparametric sample selection model.}

\item{weights}{Prior weights used during model fitting.}

\item{sp}{Estimated smoothing parameters of the smooth components for the fitted sample selection model.}

\item{iter.sp}{Number of iterations performed for the smoothing parameter estimation step.}

\item{iter.if}{Number of iterations performed in the initial step of the algorithm.}

\item{iter.inner}{Number of iterations performed inside smoothing parameter estimation step.}

\item{start.v}{Starting values for all model parameters of the semiparametric sample selection algorithm. These are obtained using the  
               Heckman sample selection correction approach when starting values are not provided and the dependence parameter is not
               specified as a function of a linear predictor.}

\item{phi}{Estimated dispersion for the response of the outcome equation. In the normal bivariate case, this corresponds to the variance. }

\item{sigma}{Estimated standard deviation for the response of the outcome equation in the case of normal marginal distribution of the outcome.}

\item{shape}{Estimated shape parameter for the response of the outcome equation in the case of gamma marginal distribution of the outcome.}

\item{nu}{Estimated shape parameter for the response of the outcome equation in the case of a discrete distribution.}

\item{theta}{Estimated coefficient linking the two equations. In the normal bivariate case, this corresponds to the correlation coefficient.}

\item{n}{Sample size.}

\item{n.sel}{Number of selected observations.}

\item{X1,X2,X3,X4,X5}{Design matrices associated with the linear predictors.}

\item{X1.d2,X2.d2,X3.d2,X4.d2,X5.d2}{Number of columns of the design matrices.}

\item{l.sp1,l.sp2,l.sp3,l.sp4,l.sp5}{Number of smooth components in the equations.}

\item{He}{Penalized hessian.}

\item{HeSh}{Unpenalized hessian.}

\item{Vb}{Inverse of the penalized hessian. This corresponds to the Bayesian variance-covariance matrix used for `confidence' interval calculations.}

\item{F}{This is given by \code{Vb*HeSh}.} 

\item{BivD}{Type of bivariate copula distribution employed.} 

\item{margins}{Margins used in the bivariate copula specification.}

\item{t.edf}{Total degrees of freedom of the estimated sample selection model. It is calculated as \code{sum(diag(F))}.}

\item{bs.mgfit}{A list of values and diagnostics extracted from \code{magic} in \code{mgcv}.}

\item{conv.sp}{If \code{TRUE} then the smoothing parameter selection algorithm converged.}
                
\item{wor.c}{Working model quantities.}                

\item{eta1,eta2}{Estimated linear predictors for the two equations.}

\item{y1}{Binary outcome of the selection equation.}

\item{y2}{Dependent variable of the outcome equation.}

\item{logLik}{Value of the (unpenalized) log-likelihood evaluated at the (penalized or unpenalized) parameter 
            estimates.}

\item{fp}{If \code{TRUE}, then a fully parametric model was fitted.} 

\item{X2s}{Full design matrix of outcome equation.}


}



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}






\seealso{
\code{\link{aver}}, \code{\link{SemiParSampleSel}}, \code{\link{plot.SemiParSampleSel}}, \code{\link{predict.SemiParSampleSel}}, \code{\link{summary.SemiParSampleSel}}
}











