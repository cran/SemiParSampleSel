\name{aver}
\alias{aver}

\title{Estimated overall average from sample selection model}

\description{ 
\code{aver} can be used to calculate the overall estimated average from a sample selection model, with corresponding confidence interval
obtained using the delta method.}



\usage{

aver(x, sig.lev = 0.05, sw = NULL, univariate = FALSE, delta = TRUE, n.sim = 100)

}

\arguments{ 

\item{x}{A fitted \code{SemiParSampleSel} object as produced by \code{SemiParSampleSel()}.}
 
\item{sig.lev}{Significance level.} 

\item{sw}{Survey weights.} 

\item{univariate}{It indicates whether the prevalence is calculated using a (naive/classic) univariate equation model 
             or sample selection model. This option has been introduced to compared adjusted and unadjusted estimates.} 

\item{delta}{If TRUE then the delta method is used for confidence interval calculations, otherwise Bayesian posterior 
             simulation is employed.} 

\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used when delta =              FALSE. It may be increased if more precision is required.} 

}


\details{

\code{aver} estimates the overall average of an outcome of interest when there are missing values that are not at random. 

}



\value{

\item{res}{It returns three values: lower confidence interval limit, estimated average and upper confidence interval limit.}
\item{sig.lev}{Significance level used.}

}




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}


\seealso{
\code{\link{SemiParSampleSel-package}}, \code{\link{SemiParSampleSel}}, \code{\link{summary.SemiParSampleSel}}  
}

\examples{
## see examples for SemiParSampleSel
}

\keyword{semiparametric sample selection modelling} 
\keyword{average outcome} 









