\name{summary.SemiParSampleSel}
\alias{summary.SemiParSampleSel}

\title{SemiParSampleSel summary}

\description{
It takes a fitted \code{SemiParSampleSel} object produced by \code{SemiParSampleSel()} and produces some summaries from it.
}

\usage{

\method{summary}{SemiParSampleSel}(object, n.sim=1000, s.meth="svd", prob.lev=0.05, 
                  cm.plot = FALSE, xlim = c(-3, 3), ylab = "Outcome margin",
                  xlab = "Selection margin", ...)

}


\arguments{ 

\item{object}{A fitted \code{SemiParSampleSel} object as produced by \code{SemiParSampleSel()}.}

\item{n.sim}{The number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
to calculate `confidence' intervals for \eqn{\theta}{\theta} and \eqn{\phi}{\phi}.} 

\item{s.meth}{Matrix decomposition used to determine the matrix root of the covariance matrix. See the documentation of \code{mvtnorm} for further details.} 

\item{prob.lev}{Probability of the left and right tails of the posterior distribution used for interval calculations.} 

\item{cm.plot}{If \code{TRUE} display contour plot of the model based on average parameter values.}

\item{xlim}{Maximum and minimum values of the selection margin to be displayed by \code{cm.plot}.}

\item{ylab}{Label of the outcome margin axis.}

\item{xlab}{Label of the selection margin axis.}

\item{...}{Other arguments.}

}


\details{ 

Using a low level function in \code{mgcv}, based on the results of Marra and Wood (2012), `Bayesian p-values' are returned for the smooth terms. These have 
better frequentist performance than their frequentist counterpart. Let \eqn{\hat{\bf f}}
and \eqn{{\bf V}_f}{V_f} denote the vector of values of a smooth term evaluated at the original covariate values and the
corresponding Bayesian covariance matrix, and let \eqn{{\bf V}_f^{r-}}{V_f^{r-}} denote 
the rank \eqn{r}{r} pseudoinverse of \eqn{{\bf V}_f}{V_f}. The statistic used 
is \eqn{T=\hat{\bf f}^\prime {\bf V}_f^{r-} \hat{\bf f}}. This is 
compared to a chi-squared distribution with degrees of freedom given by \eqn{r}{r}, which is obtained by 
biased rounding of the estimated degrees of freedom.

Covariate selection can also be achieved using a single penalty shrinkage approach as shown in Marra and Wood (2011). 

See Wojtys et al. (in press) for further details.

}



\value{

\item{tableP1}{Table containing parametric estimates, their standard errors, z-values and p-values for equation 1.}
\item{tableP2,tableP3,tableP4,tableP5}{As above but for equation 2, and equations 3, 4 and 5 if present.}
\item{tableNP1}{Table of nonparametric summaries for each smooth component including estimated degrees of freedom, estimated rank, 
                approximate Wald statistic for testing the null hypothesis that the smooth term is zero and corresponding p-value, for equation 1.}
\item{tableNP2,tableNP3,tableNP4,tableNP5}{As above but for equation 2 and equations 3, 4 and 5 if present.}
\item{n}{Sample size.}
\item{n.sel}{Number of selected observations.}
\item{sigma}{Estimated standard deviation for the response of the outcome equation in the case of normal marginal distribution of the outcome.}
\item{shape}{Estimated shape parameter for the response of the outcome equation in the case of gamma marginal distribution of the outcome.}
\item{phi}{Estimated dispersion for the response of the outcome equation.} 
\item{theta}{Estimated coefficient linking the two equations.} 
\item{nu}{Estimated coefficient for the response of the outcome equation when the Delaporte and Sichel distributions are employed.} 
\item{formula1,formula2,formula3,formula4,formula5}{Formulas used for equations 1, 2, 3, 4 and 5.}
\item{l.sp1,l.sp2,l.sp3,l.sp4,l.sp5}{Number of smooth components in equations 1, 2, 3, 4 and 5.}
\item{t.edf}{Total degrees of freedom of the estimated sample selection model.}
\item{CIsig}{`Confidence' interval for \eqn{\sigma}{\sigma} in the case of normal marginal distribution of the outcome.}
\item{CIshape}{`Confidence' interval for the shape parameter in the case of gamma distribution of the outcome.}
\item{CInu}{`Confidence' interval for the shape parameter in the case of a discrete distribution of the outcome.}
\item{CIth}{`Confidence' intervals for \eqn{\theta}{\theta}.}
\item{BivD}{Selected copula function.}
\item{margins}{Margins used in the bivariate copula specification.}
\item{n.sel}{Number of selected observations.}

}

\references{

Marra G. and Wood S.N. (2011), Practical Variable Selection for Generalized Additive Models. \emph{Computational Statistics and Data Analysis}, 55(7), 2372-2387.

Marra G. and Wood S.N. (2012), Coverage Properties of Confidence Intervals for Generalized Additive Model Components. \emph{Scandinavian Journal of Statistics}, 39(1), 53-74.

Wojtys M., Marra G. and Radice R. (in press), Copula Regression Spline Sample Selection Models: The R Package SemiParSampleSel. \emph{Journal of Statistical Software}.

}



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{SemiParSampleSelObject}}, \code{\link{plot.SemiParSampleSel}}, \code{\link{predict.SemiParSampleSel}}
}

\examples{
## see examples for SemiParSampleSel
}


\keyword{model}
\keyword{information criteria}
\keyword{smooth}
\keyword{regression}
\keyword{shrinkage smoother}
\keyword{variable selection}






