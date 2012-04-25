\name{InfCr}
\alias{InfCr}

\title{Information Criteria}

\description{ 
\code{InfCr} can be used to calculate the Akaike information criterion (AIC), corrected AIC and Bayesian information criterion.}


\usage{

InfCr(x,cr="AIC")

}

\arguments{ 

\item{x}{A fitted \code{SemiParSampleSel} object as produced by \code{SemiParSampleSel()}.}
 
\item{cr}{The criterion to be used. Default is \code{AIC}; alternatives are \code{AICc} and \code{BIC}.} 

}


\details{
When comparing fitted objects, the smaller the chosen criterion, the better the fit.
}


\value{

It returns a numeric value corresponding to the selected criterion.
}



\author{

Maintainer: Giampiero Marra \email{giampiero@stats.ucl.ac.uk}

}




\seealso{
\code{\link{SemiParSampleSel-package}}, \code{\link{SemiParSampleSel}}, \code{\link{summary.SemiParSampleSel}}
}

\examples{
## see examples for SemiParBIVProbit
}

\keyword{semiparametric sample selection modelling} 
\keyword{sample selection model} 
\keyword{information criterion} 








