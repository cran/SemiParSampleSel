SemiParSampleSel <- function(formula.eq1, formula.eq2, data=list(), weights=NULL, 
                             BivD="N", margins=c("N","N"), start.theta=NULL, start.v=NULL,
                             gamma=1, aut.sp=TRUE, fp=FALSE, rinit=1, rmax=100, 
                             fterm=sqrt(.Machine$double.eps), mterm=sqrt(.Machine$double.eps),
                             iterlimsp=50, pr.tolsp=1e-6,  
        		     control=list(maxit=50,tol=1e-6,step.half=25,rank.tol=sqrt(.Machine$double.eps)) ){


  if(margins[1]!="N") stop("Error in margin names.")
  if(!(margins[2] %in% c("N","G")) ) stop("Error in margin names.")
  if(margins[2]=="G") stop("Gamma case not finalised yet. Check the next release.")
  if(!(BivD %in% c("N", "C", "J", "FGM", "F", "AMH", "G"))) stop("Error in parameter BivD value. It should be one of: N, C, J, FGM, F, AMH, G.")

  qu.mag <- sp <- NULL
  
  ########################################
  # Naive models needed for starting point
  ########################################

  gam1 <- eval(substitute(gam(formula.eq1, binomial(link="probit"), gamma=gamma, weights=weights, data=data),list(weights=weights)))
  y1 <- gam1$y; inde <- y1 > 0
  if(margins[2]=="N")      gam2 <- eval(substitute(gam(formula.eq2, gamma=gamma, data=data, weights=weights, subset=inde),list(weights=weights,inde=inde)))
  else if(margins[2]=="G") gam2 <- eval(substitute(gam(formula.eq2, gamma=gamma, data=data, weights=weights, subset=inde, family=Gamma(link = "log")),list(weights=weights,inde=inde))) else stop("Error in margin names.")

  ########################################
  # Data Objects
  ########################################

  n <- length(inde)
  X1 <- model.matrix(gam1);  X1.d2 <- dim(X1)[2];  X2.d2 <- length(coef(gam2))
  X2 <- matrix(0,n,X2.d2,dimnames = list(c(1:n),c(names(coef(gam2)))) )
  X2[inde, ] <- model.matrix(gam2)
  y2 <- rep(0,n); y2[inde] <- gam2$y
  l.sp1 <- length(gam1$smooth); l.sp2 <- length(gam2$smooth)
  if(is.null(weights)) weights <- rep(1,length(gam1$y)) else weights <- model.weights(gam1$model) 
  if( (l.sp1!=0 || l.sp2!=0) && fp==FALSE) qu.mag <- S.m(gam1,gam2,l.sp1,l.sp2)

  ############################################
  # Heckmann's procedure for linear predictors
  ############################################

  p.g1 <- predict(gam1)
  imr <- data$imr <- dnorm(p.g1)/pnorm(p.g1)

  formula.eq2.1 <- update.formula(formula.eq2, ~. + imr)
  if(margins[2]=="N")      gam2.1 <- eval(substitute(gam(formula.eq2.1, gamma=gamma, data=data, weights=weights, subset=inde),list(weights=weights,inde=inde)))
  else if(margins[2]=="G") gam2.1 <- eval(substitute(gam(formula.eq2.1, gamma=gamma, data=data, weights=weights, subset=inde, family=Gamma(link = "log")),list(weights=weights,inde=inde)))

  sigma <- sqrt(mean(resid(gam2.1)^2)+mean(imr[inde]*(imr[inde]+p.g1[inde]))*gam2.1$coef["imr"]^2)[[1]]; names(sigma) <- "sigma"
  co  <- (gam2.1$coef["imr"]/sigma)[[1]] 

  if(margins[2]=="G") { k <- (summary(gam2)$dispersion)^(-1); names(k) <- "shape" }

  	if(l.sp1!=0 && l.sp2!=0) sp <- c(gam1$sp,gam2.1$sp)
  	if(l.sp1==0 && l.sp2!=0) sp <- c(gam2.1$sp)
        if(l.sp1!=0 && l.sp2==0) sp <- c(gam1$sp)

  ########################################
  # Starting value for theta
  ########################################

  a.theta <- st.theta.star(start.theta, co, BivD);  names(a.theta) <- "theta.star"


  if(margins[2]=="N") { 
        if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"],log(sigma),a.theta)
        names(start.v)[length(start.v)-1] <- "log.sigma"
  }
  
  if(margins[2]=="G") {
        if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"],log(k),a.theta)
        names(start.v)[length(start.v)-1] <- "log.shape"
  }



  ########################################
  # Fitting procedure
  ########################################

  fit <- fit.SemiParSampleSel( dat=cbind(y1,y2), X1=X1, X2=X2, qu.mag=qu.mag, sp=sp, X1.d2=X1.d2, X2.d2=X2.d2, gp1=gam1$nsdf, gp2=gam2$nsdf, n=n,
                          l.sp1=l.sp1, l.sp2=l.sp2, weights=weights, iterlimSP=iterlimsp, pr.tol=pr.tolsp, BivD=BivD, margins=margins,
                          aut.sp=aut.sp, gamma=gamma, fp=fp, start.v=start.v, rinit=rinit, rmax=rmax, fterm=fterm, mterm=mterm, control=control )

  He <- fit$fit$hessian
  He.eig <- eigen(He,symmetric=TRUE)
  k.e <- sum(as.numeric(He.eig$val<sqrt(.Machine$double.eps)))
   
  if(k.e!=0){
      ind.e <- (length(He.eig$val)-(k.e-1)):length(He.eig$val)
      min.e <- min(He.eig$val[1:(ind.e[1]-1)])
      for(i in 1:k.e) He.eig$val[ind.e[i]] <- min.e/10^i  
      Vb <- He.eig$vec%*%diag(1/He.eig$val)%*%t(He.eig$vec)      
  }
  else {
      Vb <- He.eig$vec%*%diag(1/He.eig$val)%*%t(He.eig$vec) 
  }
         
  if((l.sp1!=0 || l.sp2!=0) && fp==FALSE){ HeSh <- He - fit$fit$S.h; F <- Vb%*%HeSh } else { HeSh <- He; F <- diag(rep(1,dim(Vb)[1])) }      
  t.edf <- sum(diag(F))


  ###########################################################################
  # Transforming theta back to its original scale and computing Kendall's tau
  ###########################################################################

  theta <- theta.tau(BivD=BivD, theta.star=fit$fit$argument["theta.star"])
  KendTau <- theta[2]
  theta   <- theta[1]

  names(theta) <- names(KendTau) <- NULL

  if(margins[2]=="N") { sigma <- exp(fit$fit$argument["log.sigma"]); names(sigma) <- k <- NULL; phi <- sigma^2 }  # dispersion: variance for normal case
  if(margins[2]=="G") { k     <- exp(fit$fit$argument["log.shape"]); names(k) <- sigma <- NULL; phi <- k^{-1} }   # dispersion: 1/k for gamma case


  ###########################################################################
  # Final object
  ###########################################################################

       L <- list(fit=fit$fit, gam1=gam1, gam2=gam2, gam2.1=gam2.1, coefficients=fit$fit$argument, weights=weights, sp=fit$sp, 
            iter.if=fit$iter.if, iter.sp=fit$iter.sp, iter.fi=fit$iter.fi, start.v=start.v, phi=phi, sigma=sigma, shape=k, theta=theta, tau=KendTau, 
            n=n, n.sel=length(gam2$y), X1=X1, X2=X2, X1.d2=X1.d2, X2.d2=X2.d2, 
            l.sp1=l.sp1, l.sp2=l.sp2, He=He, HeSh=HeSh, Vb=Vb, F=F, BivD=BivD, margins=margins,
            t.edf=t.edf, bs.mgfit=fit$bs.mgfit, conv.sp=fit$conv.sp, wor.c=fit$wor.c, eta1=fit$fit$eta1, eta2=fit$fit$eta2, 
            y1=y1, y2=y2, logL=fit$fit$l, fp=fp )

  class(L) <- "SemiParSampleSel"

  L

}
