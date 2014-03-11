print.summary.SemiParSampleSel <- function(x,digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){


  if(x$margins[2]=="N") nn <- "with normal margins"
  if(x$margins[2]=="G") nn <- "with normal-gamma margins"

  if(x$BivD=="N")   cop <- paste("Bivariate Normal",nn)
  else if(x$BivD=="C")   cop <- paste("Clayton Copula",nn)
  else if(x$BivD=="rC")  cop <- paste("90 degrees rotated Clayton Copula",nn)
  else if(x$BivD=="J")   cop <- paste("Joe Copula",nn) 
  else if(x$BivD=="rJ")  cop <- paste("90 degrees rotated Joe Copula",nn) 
  else if(x$BivD=="FGM") cop <- paste("FGM Copula",nn) 
  else if(x$BivD=="F")   cop <- paste("Frank Copula",nn)
  else if(x$BivD=="AMH") cop <- paste("AMH Copula",nn) 
  else if(x$BivD=="G")   cop <- paste("Gumbel Copula",nn) 
  else if(x$BivD=="rG")  cop <- paste("90 degrees rotated Gumbel Copula",nn) 

  cat("\nFamily: SAMPLE SELECTION",cop,"\n\nSELECTION EQ.: ")
  print(x$formula1)
  cat("\n") 
  printCoefmat(x$tableP1,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

  if(x$l.sp1!=0){
      cat("Smooth components' approximate significance:\n")
      printCoefmat(x$tableNP1,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
      cat("\n")
  }

  cat("\nOUTCOME EQ.: ")
  print(x$formula2)
  cat("\n")
  printCoefmat(x$tableP2,digits = digits, signif.stars = signif.stars,na.print = "NA",)
  cat("\n")

  if(x$l.sp2!=0){
      cat("Smooth components' approximate significance:\n")
      printCoefmat(x$tableNP2,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
      cat("\n")
  }

  if(x$BivD=="N") cp <- ("rho = ") else cp <- ("theta = ")
                  
                  
if(x$margins[2]=="N"){                     
                   dis <-  "  sigma = "
                   diss <- x$sigma
                   CIdiss <- c(x$CIsig[1],x$CIsig[2]) 
                     }else{
                   dis <-  "  shape = "
		   diss <- x$shape
                   CIdiss <- c(x$CIshape[1],x$CIshape[2]) 
                   }

cat("n = ",x$n,"  n.sel = ",x$n.sel,dis,round(diss,3),"(",round(CIdiss[1],3),",",round(CIdiss[2],3),")","\n",cp,round(x$theta,3),"(",round(x$CIth[1],3),",",round(x$CIth[2],3),")","  Kendall's Tau = ",round(x$tau,3),"(",round(x$CIkt[1],3),",",round(x$CIkt[2],3),")","\n","total edf = ",round(x$t.edf,3),"\n\n",sep="") 
  
  invisible(x)

}


















