print.summary.SemiParSampleSel <- function(x,digits = max(3, getOption("digits") - 3),
                                             signif.stars = getOption("show.signif.stars"),...){

  cat("\nFamily: SAMPLE SELECTION\n\nSELECTION EQ.: ")
  print(x$formula1)
  cat("\n") 
  printCoefmat(x$tableP1,digits = digits, signif.stars = signif.stars,na.print = "NA",...)
  cat("\n")

    if(x$l.sp1!=0 && x$l.sp2!=0){
    cat("Smooth terms' approximate significance:\n")
    printCoefmat(x$tableNP1,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }

  cat("\nOUTCOME EQ.: ")
  print(x$formula2)
  cat("\n")
  printCoefmat(x$tableP2,digits = digits, signif.stars = signif.stars,na.print = "NA",)
  cat("\n")

    if(x$l.sp1!=0 && x$l.sp2!=0){
    cat("Smooth terms' approximate significance:\n")
    printCoefmat(x$tableNP2,digits = digits, signif.stars = signif.stars,has.Pvalue = TRUE,na.print = "NA",cs.ind = 1,...)
    cat("\n")
    }

    if(x$l.sp1!=0 && x$l.sp2!=0){ cat("n = ",x$n,"  n.sel = ",x$n.sel,"  sigma = ",round(x$sigma,3),"(",round(x$CIsi[1],3),",",round(x$CIsi[2],3),")","\n","rho = ",round(x$rho,3),"(",round(x$CIrs[1],3),",",round(x$CIrs[2],3),")","  total edf = ",round(x$t.edf,3),"\n\n", sep="") }
    else{ cat("n = ",x$n,"  n.sel = ",x$n.sel,"  sigma = ",round(x$sigma,3),"(",round(x$CIsi[1],3),",",round(x$CIsi[2],3),")","\n","rho = ",round(x$rho,3),"(",round(x$CIrs[1],3),",",round(x$CIrs[2],3),")","  total edf = ",round(x$t.edf,3),"\n\n",sep="")  }

invisible(x)

}

