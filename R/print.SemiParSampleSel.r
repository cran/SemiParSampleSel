print.SemiParSampleSel <- function(x,...){

  cat("\nFamily: SAMPLE SELECTION\n\nSELECTION EQ.: ")
  print(x$gam1$formula)

  cat("OUTCOME   EQ.: ")
  print(x$gam2$formula)
  cat("\n")

  cat("n = ",x$n,"  n.sel = ",x$n.sel,"  sigma = ",round(x$sigma,3),"  rho = ",round(x$rho,3),"  total edf = ",round(x$t.edf,3),"\n\n")

invisible(x)

}
