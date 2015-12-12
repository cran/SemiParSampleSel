\name{post.check}
\alias{post.check}

\title{Post-estimation response diagnostics}

\description{ It takes the response vector and produces a QQ-plot based on the estimates obtained 
from the sample selection model. For discrete responses normalized and randomized QQ-plots are produced.

}


\usage{

post.check(x, bd = 1)

}



\arguments{ 

\item{x}{A fitted \code{SemiParSampleSel} object as produced by \code{SemiParSampleSel()}.} 

\item{bd}{Binomial denominator - maximum value the response can take. This argument is only relevant in the case of BI, BB, ZIBI, ZABI, ZIBB and ZABB.} 
 

}





\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}




