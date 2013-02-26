theta.tau <- function( BivD, theta.star ) {

  if(BivD=="N")   { theta <- tanh(theta.star);	KendTau <- tau(normalCopula(theta))  }
  if(BivD=="C")   { theta <- exp(theta.star);	KendTau <- tau(claytonCopula(theta)) }
  if(BivD=="J")   { theta <- 1+exp(theta.star);	KendTau <- tau(joeCopula(theta))     }
  if(BivD=="FGM") { theta <- tanh(theta.star);	KendTau <- tau(fgmCopula(theta))     }
  if(BivD=="F")   { theta <- theta.star;        KendTau <- tau(frankCopula(theta))   }
  if(BivD=="AMH") { theta <- tanh(theta.star);	KendTau <- tau(amhCopula(theta))     }
  if(BivD=="G")   { theta <- 1+exp(theta.star);	KendTau <- tau(gumbelCopula(theta))  }
  
  c(theta,KendTau)

}