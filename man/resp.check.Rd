\name{resp.check}
\alias{resp.check}

\title{Preliminary response diagnostics}

\description{ It takes takes the response vector and produces a histogram and QQ-plot. For discrete responses normalized and randomized QQ-plots
are produced.


}


\usage{

resp.check(y, margin = "N", bd = 1)

}



\arguments{ 

\item{y}{reponse vector} 

\item{margin}{margin of interest. Is the response distributed according to the margin? The default is normal outcome.} 

\item{bd}{Binomial denominator - maximum value the response can take. This argument is only relevant in the case of BI, BB, ZIBI, ZABI, ZIBB and ZABB.} 
 

}





\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}


\examples{

library(SemiParSampleSel)

# Generating normal response with n = 2000
ys <- rnorm(2000, mean=3, sd=1)
resp.check(ys, margin = "N")

# Generating gamma response with n = 2000
ys <- rgamma(2000, shape=2, rate=3)
resp.check(ys, margin = "G")

# Generating Poisson response with n = 2000
ys <- rPO(2000, mu=3)
resp.check(ys, margin = "P")

# Generating negative binomial response with n = 2000
ys <- rNBI(2000, mu=3, sigma=1)
resp.check(ys, margin = "NB")



}






