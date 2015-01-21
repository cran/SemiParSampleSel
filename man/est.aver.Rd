\name{est.aver}
\alias{est.aver}

\title{Estimated overall average from sample selection model}

\description{ 
\code{est.aver} can be used to calculate the overall estimated average from a sample selection model, with corresponding confidence interval
obtained using the delta method.}



\usage{

est.aver(x, sig.lev = 0.05, sw = NULL, naive = FALSE)

}

\arguments{ 

\item{x}{A fitted \code{SemiParSampleSel} object as produced by \code{SemiParSampleSel()}.}
 
\item{sig.lev}{Significance level.} 

\item{sw}{Survey weights.} 

\item{naive}{It indicates whether the prevalence is calculated using a (naive/classic) univariate equation model 
             or sample selection model. This option has been introduced to compared adjusted and unadjusted estimates.} 


}


\details{

est.aver estimates the overall average of an outcome of interest when there are missing values that are not at random. 

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









