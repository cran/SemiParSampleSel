theta.tau <- function( BivD, theta.star ) {

  if(BivD=="N")   { theta <- tanh(theta.star);	KendTau <- tau(normalCopula(theta))  }
  else if(BivD=="C")   { theta <- exp(theta.star);	KendTau <- tau(claytonCopula(theta)) }
  else if(BivD=="rC")  { theta <- exp(theta.star);	KendTau <- -tau(claytonCopula(theta))}
  else if(BivD=="J")   { theta <- 1+exp(theta.star);	KendTau <- tau(joeCopula(theta))     }
  else if(BivD=="rJ")  { theta <- 1+exp(theta.star);	KendTau <- -tau(joeCopula(theta))    }
  else if(BivD=="FGM") { theta <- tanh(theta.star);	KendTau <- tau(fgmCopula(theta))     }
  else if(BivD=="F")   { theta <- theta.star;           KendTau <- tau(frankCopula(theta))   }
  else if(BivD=="AMH") { theta <- tanh(theta.star);	KendTau <- tau(amhCopula(theta))     }
  else if(BivD=="G")   { theta <- 1+exp(theta.star);	KendTau <- tau(gumbelCopula(theta))  }
  else if(BivD=="rG")  { theta <- 1+exp(theta.star);	KendTau <- -tau(gumbelCopula(theta)) }
  
  c(theta,KendTau)

}