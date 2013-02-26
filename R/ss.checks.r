ss.checks <- function(x){

  e.v <- eigen(x$He)$values

  cat("\nLargest absolute gradient value:",max(abs(x$fit$gradient)))

  if (min(e.v) > 0) cat("\nInformation matrix is positive definite\n")
  else cat("\nEigenvalue range of the information matrix: [",min(e.v),",",max(e.v),"]\n", sep = "")

  if((x$l.sp1!=0 || x$l.sp2!=0) && x$fp==FALSE){
      cat("\nTrust region iterations before smoothing parameter estimation:",x$iter.if)
      cat("\nSmoothing parameter/leapfrog loops:",x$iter.sp)
      cat("\nTrust region iterations after smoothing parameter/leapfrog step:",x$iter.fi,"\n\n")
      

  }
  else cat("\nTrust region iterations:",x$iter.if,"\n\n")


}
