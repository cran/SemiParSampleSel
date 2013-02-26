print.SemiParSampleSel <- function(x,...){

  if(x$margins[2]=="N")      nn <- "with normal margins"
  else if(x$margins[2]=="G") nn <- "with normal-gamma margins"

  if(x$BivD=="N")   cop <- paste("Bivariate Normal",nn)
  if(x$BivD=="C")   cop <- paste("Clayton Copula",nn)
  if(x$BivD=="J")   cop <- paste("Joe Copula",nn) 
  if(x$BivD=="FGM") cop <- paste("FGM Copula",nn) 
  if(x$BivD=="F")   cop <- paste("Frank Copula",nn)
  if(x$BivD=="AMH") cop <- paste("AMH Copula",nn) 
  if(x$BivD=="G")   cop <- paste("Gumbel Copula",nn) 

  cat("\nFamily: SAMPLE SELECTION",cop,"\n\nSELECTION EQ.: ")
  print(x$gam1$formula)

  cat("OUTCOME   EQ.: ")
  print(x$gam2$formula)
  cat("\n")


  if(x$BivD=="N") cp <- ("rho = ") else cp <- ("theta = ")

if(x$margins[2]=="N") cat("n = ",x$n,"  n.sel = ",x$n.sel,"  sigma = ",round(x$sigma,3),"\n",cp,round(x$theta,3),"  Kendall's Tau = ",round(x$tau,3),"\ntotal edf = ",round(x$t.edf,3),"\n\n",sep="")
  
if(x$margins[2]=="G") cat("n = ",x$n,"  n.sel = ",x$n.sel,"  shape = ",round(x$shape,3),"\n",cp,round(x$theta,3),"  Kendall's Tau = ",round(x$tau,3),"\ntotal edf = ",round(x$t.edf,3),"\n\n",sep="")


  invisible(x)

}
