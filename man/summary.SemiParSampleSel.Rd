\name{summary.SemiParSampleSel}
\alias{summary.SemiParSampleSel}

\title{SemiParSampleSel summary}

\description{
It takes a fitted \code{SemiParSampleSel} object produced by \code{SemiParSampleSel()} and produces some summaries from it.
}

\usage{

\method{summary}{SemiParSampleSel}(object, n.sim = 1000, s.meth = "svd", prob.lev = 0.05, ...)


}


\arguments{ 

\item{object}{A fitted \code{SemiParSampleSel} object as produced by \code{SemiParSampleSel()}.}

\item{n.sim}{The number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
to calculate `confidence' intervals for \eqn{\theta}{\theta} and \eqn{\phi}{\phi}.} 

\item{s.meth}{Matrix decomposition used to determine the matrix root of the covariance matrix. See the documentation of \code{mvtnorm} for further details.} 

\item{prob.lev}{Probability of the left and right tails of the posterior distribution used for interval calculations.} 

\item{...}{Other arguments.}

}


\details{ 

As in the package \code{mgcv}, based on the results of Marra and Wood (2012), `Bayesian p-values' are returned for the smooth terms. These have 
better frequentist performance than their frequentist counterpart. Let \eqn{\hat{\bf f}}
and \eqn{{\bf V}_f}{V_f} denote the vector of values of a smooth term evaluated at the original covariate values and the
corresponding Bayesian covariance matrix, and let \eqn{{\bf V}_f^{r-}}{V_f^{r-}} denote 
the rank \eqn{r}{r} pseudoinverse of \eqn{{\bf V}_f}{V_f}. The statistic used 
is \eqn{T=\hat{\bf f}^\prime {\bf V}_f^{r-} \hat{\bf f}}. This is 
compared to a chi-squared distribution with degrees of freedom given by \eqn{r}{r}, which is obtained by 
biased rounding of the estimated degrees of freedom.

Covariate selection can also be achieved using a single penalty shrinkage approach as shown in Marra and Wood (2011). 

See Wojtys et al. (submitted) for further details.

}



\value{

\item{tableP1}{Table containing parametric estimates, their standard errors, z-values and p-values for equation 1.}
\item{tableP2}{As above but for equation 2.}
\item{tableNP1}{Table of nonparametric summaries for each smooth component including estimated degrees of freedom, estimated rank, 
                approximate Wald statistic for testing the null hypothesis that the smooth term is zero and corresponding p-value, for equation 1.}
\item{tableNP2}{As above but for equation 2.}
\item{n}{Sample size.}
\item{n.sel}{Number of selected observations.}
\item{sigma}{Estimated standard deviation for the response of the outcome equation in the case of normal marginal distribution of the outcome.}
\item{shape}{Estimated shape parameter for the response of the outcome equation in the case of gamma marginal distribution of the outcome.}
\item{phi}{Estimated dispersion for the response of the outcome equation.} 
\item{theta}{Estimated coefficient linking the two equations.} 
\item{tau}{Estimated Kendall's tau association coefficient between the two equations.}
\item{formula1}{Formula used for equation1.}
\item{formula2}{Formula used for equation2.}
\item{l.sp1}{Number of smooth components in equation 1.}
\item{l.sp2}{Number of smooth components in equation 2.}
\item{t.edf}{Total degrees of freedom of the estimated sample selection model.}
\item{CIsig}{`Confidence' interval for \eqn{\sigma}{\sigma} in the case of normal marginal distribution of the outcome.}
\item{CIshape}{`Confidence' interval for the shape parameter in the case of gamma distribution of the outcome.}
\item{CIth}{`Confidence' intervals for \eqn{\theta}{\theta}.}
\item{CIkt}{`Confidence' intervals for \eqn{\tau}{\tau}.}
\item{BivD}{Selected copula function.}
\item{margins}{Margins used in the bivariate copula specification.}
\item{n.sel}{Number of selected observations.}

}

\references{

Marra G. and Wood S.N. (2011), Practical Variable Selection for Generalized Additive Models. \emph{Computational Statistics and Data Analysis}, 55(7), 2372-2387.

Marra G. and Wood S.N. (2012), Coverage Properties of Confidence Intervals for Generalized Additive Model Components. \emph{Scandinavian Journal of Statistics}, 39(1), 53-74.

Wojtys M., Marra G. and Radice R. (submitted), Copula Regression Spline Sample Selection Models: The R Package SemiParSampleSel.

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






