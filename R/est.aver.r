est.aver <- function(x, sig.lev = 0.05, sw = NULL, naive = FALSE){

X2sg <- x$X2s

if(naive==TRUE)  eta2 <- X2sg%*%coef(x$gam2) 
if(naive==FALSE) eta2 <- x$eta2

if(is.null(sw)) sw <- rep(1,length(eta2)) 

if(x$margins[2] == "N") core <- apply( c(eta2)*X2sg,      2, weighted.mean,  w = sw) 
if(x$margins[2] != "N") core <- apply( c(exp(eta2))*X2sg, 2, weighted.mean,  w = sw) 

if(naive==FALSE){

  if(x$margins[2] %in% c("N","G","NB","PIG") )  zs <- c(0,0) 
  if(x$margins[2] == "P")                       zs <- c(0)
  if(x$margins[2] %in% c("D", "S"))             zs <- c(0,0,0) 
       
}


if(naive==FALSE) G <- c( rep(0, x$X1.d2), core , zs)  
if(naive==TRUE)  G <- c( core ) 
  
if(x$margins[2] == "N") wm <- weighted.mean(eta2, w=sw)
if(x$margins[2] != "N") wm <- weighted.mean(exp(eta2), w=sw)
  
if(naive==FALSE) Vv <- x$Vb
if(naive==TRUE)  Vv <- x$gam2$Vp  


  sv <- sqrt( t(G)%*%Vv%*%G ) 

  qz <- qnorm(sig.lev/2, lower.tail = FALSE)
  lb <- wm - qz*sv 
  ub <- wm + qz*sv 

  res <- c(lb, wm, ub)

  out <- list(res = res, sig.lev = sig.lev)
 
  class(out) <- "est.aver"

  out

}


# importFrom(gamlss.dist,dNBI,pNBI,dDEL,pDEL,dPIG,pPIG,dSI,pSI,tofyDEL2)
