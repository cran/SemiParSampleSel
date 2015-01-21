theta.tau <- function( BivD, theta.star ) {

  if(BivD=="N")                   { theta <- tanh(theta.star) }#;	KendTau <- tau(normalCopula(theta))  }
  if(BivD %in% c("C0", "C180") )  { theta <- exp(theta.star)}#;	KendTau <- tau(claytonCopula(theta)) }
  if(BivD %in% c("C90","C270") )  { theta <- exp(theta.star)}#;	KendTau <- -tau(claytonCopula(theta))}
  if(BivD %in% c("J0", "J180"))   { theta <- 1+exp(theta.star)}#;	KendTau <- tau(joeCopula(theta))     }
  if(BivD %in% c("J90", "J270"))  { theta <- 1+exp(theta.star)}#;	KendTau <- -tau(joeCopula(theta))    }
  if(BivD=="FGM")                 { theta <- tanh(theta.star)}#;	KendTau <- tau(fgmCopula(theta))     }
  if(BivD=="F")                   { theta <- theta.star}#;        KendTau <- tau(frankCopula(theta))   }
  if(BivD=="AMH")                 { theta <- tanh(theta.star)}#;	KendTau <- tau(amhCopula(theta))     }
  if(BivD %in% c("G0", "G180"))   { theta <- 1+exp(theta.star)}#;	KendTau <- tau(gumbelCopula(theta))  }
  if(BivD %in% c("G90", "G270"))  { theta <- 1+exp(theta.star)}#;	KendTau <- -tau(gumbelCopula(theta)) }
  
  KendTau <- NULL
  
  c(theta, KendTau)

}