\name{logLik.SemiParSampleSel}
\alias{logLik.SemiParSampleSel}

\title{Extract the log likelihood for a fitted SemiParSampleSel}


\description{It extracts the log-likelihood for a fitted \code{SemiParSampleSel} model.

}

\usage{
\method{logLik}{SemiParSampleSel}(object, ...)
}


\arguments{ 
  
\item{object}{A fitted \code{SemiParSampleSel} object as produced by \code{SemiParSampleSel()}.}

\item{...}{Un-used for this function.}


}


\details{

Modification of the classic \code{logLik} which accounts for the estimated degrees of freedom used in \code{SemiParSampleSel} objects.
This function is provided so that information criteria work correctly with \code{SemiParSampleSel} objects by using the correct degrees 
of freedom. 

}



\value{

Standard \code{\link{logLik}} object.

}



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}




\seealso{
\code{\link{AIC}}, \code{\link{BIC}}
}



\examples{
## see examples for SemiParSampleSel
}

\keyword{logLik} 
\keyword{AIC} 
\keyword{BIC}






