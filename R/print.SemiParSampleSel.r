print.SemiParSampleSel <- function(x,...){

  if(x$margins[2]=="N")      nn <- "with normal margins"
  if(x$margins[2]=="G")      nn <- "with normal-gamma margins"
  if(x$margins[2]=="P")      nn <- "with normal-poisson margins"
  if(x$margins[2]=="NB")     nn <- "with normal-negative binomial margins"
  if(x$margins[2]=="D")      nn <- "with normal-delaporte margins"
  if(x$margins[2]=="PIG")    nn <- "with normal-poisson inverse gaussian margins"
  if(x$margins[2]=="S")      nn <- "with normal-sichel margins"
  
  if(x$BivD=="N")   cop <- paste("Bivariate Normal",nn)
  if(x$BivD=="C0")  cop <- paste("Clayton Copula",nn)
  if(x$BivD=="C90") cop <- paste("Rotated Clayton Copula (90 degrees)",nn)
  if(x$BivD=="C180")cop <- paste("Survival Clayton Copula",nn)
  if(x$BivD=="C270")cop <- paste("Rotated Clayton Copula (270 degrees)",nn) 
  if(x$BivD=="J0")  cop <- paste("Joe Copula",nn)
  if(x$BivD=="J90") cop <- paste("Rotated Joe Copula (90 degrees)",nn)
  if(x$BivD=="J180")cop <- paste("Survival Joe Copula",nn)
  if(x$BivD=="J270")cop <- paste("Rotated Joe Copula (270 degrees)",nn) 
  if(x$BivD=="FGM") cop <- paste("FGM Copula",nn) 
  if(x$BivD=="F")   cop <- paste("Frank Copula",nn)
  if(x$BivD=="AMH") cop <- paste("AMH Copula",nn) 
  if(x$BivD=="G0")  cop <- paste("Gumbel Copula",nn)
  if(x$BivD=="G90") cop <- paste("Rotated Gumbel Copula (90 degrees)",nn)
  if(x$BivD=="G180")cop <- paste("Survival Gumbel Copula",nn)
  if(x$BivD=="G270")cop <- paste("Rotated Gumbel Copula (270 degrees)",nn) 


  cat("\nFamily: ",cop,"\n\nSELECTION EQ.: ")
  print(x$gam1$formula)

  cat("OUTCOME   EQ.: ")
  print(x$gam2$formula)
  cat("\n")


  if(x$BivD=="N") cp <- ("rho = ") else cp <- ("theta = ")


if(x$margins[2] %in% c("N", "NB", "PIG", "D", "S")){dis <-  "  sigma = "; diss <- x$sigma}
if(x$margins[2]=="G"){dis <-  "  shape = "; diss <- x$shape}
if(x$margins[2] %in% c("D", "S")){dis1 <-  "  nu = "; diss1 <- x$nu}


#
#if(x$margins[2] %in% c("G", "N", "NB", "PIG")){
#cat("n = ",x$n,"  n.sel = ",x$n.sel,dis,round(diss,3),"\n",cp,round(x$theta,3),"  Kendall's Tau = ",round(x$tau,3),"\ntotal edf = ",round(x$t.edf,3),"\n\n",sep="")
#} 
#
#if(x$margins[2] =="P"){
#cat("n = ",x$n,"  n.sel = ",x$n.sel,"\n",cp,round(x$theta,3),"  Kendall's Tau = ",round(x$tau,3),"\ntotal edf = ",round(x$t.edf,3),"\n\n",sep="")
#}
#
#
#if(x$margins[2] %in% c("D", "S")){
#cat("n = ",x$n,"  n.sel = ",x$n.sel,dis,round(diss,3),dis1,round(diss1,3),"\n",cp,round(x$theta,3),"  Kendall's Tau = ",round(x$tau,3),"\ntotal edf = ",round(x$t.edf,3),"\n\n",sep="")
#}  
#
#

if(x$margins[2] %in% c("G", "N", "NB", "PIG")){
cat("n = ",x$n,"  n.sel = ",x$n.sel,dis,round(diss,3),"\n",cp,round(x$theta,3),"  total edf = ",round(x$t.edf,3),"\n\n",sep="")
} 

if(x$margins[2] =="P"){
cat("n = ",x$n,"  n.sel = ",x$n.sel,"\n",cp,round(x$theta,3),"  total edf = ",round(x$t.edf,3),"\n\n",sep="")
}


if(x$margins[2] %in% c("D", "S")){
cat("n = ",x$n,"  n.sel = ",x$n.sel,dis,round(diss,3),dis1,round(diss1,3),"\n",cp,round(x$theta,3),"  total edf = ",round(x$t.edf,3),"\n\n",sep="")
}  




  invisible(x)

}
















