SemiParSampleSel <- function(formula, data = list(), weights = NULL, subset = NULL, 
                             method = "trwlF", start.v = NULL, start.theta = NULL, 
                             BivD = "N", margins = c("N","N"), fp = FALSE, infl.fac = 1, pPen1 = NULL, pPen2 = NULL, 
                             rinit = 1, rmax = 100, iterlimsp = 50, pr.tolsp = 1e-6, 
                             parscale){ # pPen1 = NULL, pPen2 = NULL,
 
 if(margins[2] == "G") stop("For tested version of Gamma case, check next release.")
 if(margins[2] %in% c("P", "NB", "D", "PIG", "S")) stop("Check next release for final tested version of discrete case.")
 
  marg2 <- c("N", "G", "P", "NB", "D", "PIG", "S")
  copul <- c("N", "FGM", "F", "AMH", 
             "C0", "C90", "C180", "C270", 
             "J0", "J90", "J180", "J270", 
             "G0", "G90", "G180", "G270")
  copul1 <- c("C180", "C270", "J180", "J270", "G180", "G270")           
 
  if(margins[1]!="N") stop("Error in selection equation's margin name. It can currently only be N.")
  if(!(margins[2] %in% marg2) ) stop("Error in outcome equation's margin name. It should be one of: N, G, P, NB, D, PIG, S.")
  if(!(BivD %in% copul)) stop("Error in parameter BivD value. It should be one of: N, FGM, F, AMH, C0, C90, C180, C270, J0, J90, J180, J270, G0, G90, G180, G270.")
   
  if(margins[2] %in% c("N","G") && BivD %in% copul1) stop("Rotations by 180 and 270 for normal and gamma outcome margins not currently supported.")
   
  qu.mag <- sp <- nu <- NULL
  
  ############
  # Data check 
  ############  
  
  ig <- interpret.gam(formula)
  mf <- match.call(expand.dots = FALSE)
  
  pred.n <- union(ig[[1]]$pred.names,c(ig[[2]]$pred.names,ig[[2]]$response))
  fake.formula <- paste(ig[[1]]$response, "~", paste(pred.n, collapse = " + ")) 
  environment(fake.formula) <- environment(ig$fake.formula)
  mf$formula <- fake.formula 
  mf$start.v <- mf$start.theta <- mf$BivD <- mf$margins <- mf$infl.fac <- mf$fp <- mf$rinit <-  mf$pPen1 <- mf$pPen2 <- mf$rmax <- mf$iterlimsp <- mf$pr.tolsp <- mf$parscale <- NULL
  mf$drop.unused.levels <- TRUE 
  mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())

  indS <- as.logical(data[,ig[[1]]$response])==FALSE 
  indS <- ifelse( is.na(indS), FALSE, indS) 
  data[indS, ig[[2]]$response] <- ifelse( is.na(data[indS, ig[[2]]$response]), 0, data[indS, ig[[2]]$response]) 
  data <- na.omit(data)
   
  if(is.null(weights)) weights <- rep(1,dim(data)[1]) else weights <- data[,"(weights)"]    
    
  formula.eq1 <- formula[[1]]
  formula.eq2 <- formula[[2]] 
    

  #########################################
  # Naive models needed for starting values 
  #########################################

  gam1 <- eval(substitute(gam(formula.eq1, binomial(link="probit"), gamma=infl.fac, weights=weights, data=data, paraPen=pPen1),list(weights=weights)))
  y1 <- gam1$y
  inde <- y1 > 0
  
  if(margins[2]=="N")               gam2 <- eval(substitute(gam(formula.eq2, gamma=infl.fac, data=data, paraPen=pPen2, weights=weights, subset=inde),                               list(weights=weights,inde=inde)))
  if(margins[2]=="G")               gam2 <- eval(substitute(gam(formula.eq2, gamma=infl.fac, data=data, paraPen=pPen2, weights=weights, subset=inde, family=Gamma(link = "log")),   list(weights=weights,inde=inde)))
  if(!(margins[2] %in% c("N","G"))) gam2 <- eval(substitute(gam(formula.eq2, gamma=infl.fac, data=data, paraPen=pPen2, weights=weights, subset=inde, family=poisson(link = "log")), list(weights=weights,inde=inde)))

  ##############
  # Data Objects
  ##############

  n <- length(y1)
  X1 <- model.matrix(gam1); X1.d2 <- dim(X1)[2]; X2.d2 <- length(coef(gam2))
  X2 <- matrix(0, n, X2.d2, dimnames = list(c(1:n),c(names(coef(gam2)))) )
  X2[inde, ] <- model.matrix(gam2)
  y2 <- rep(0,n); y2[inde] <- gam2$y
  l.sp1 <- length(gam1$smooth)
  l.sp2 <- length(gam2$smooth)
  gp1 <- gam1$nsdf
  gp2 <- gam2$nsdf
  
  if( (l.sp1!=0 || l.sp2!=0) && fp==FALSE) qu.mag <- S.m(gam1, gam2, l.sp1, l.sp2)

  ########################################################
  # Approximate Heckmann's procedure for linear predictors
  ########################################################

  p.g1 <- predict(gam1)
  imr <- data$imr <- dnorm(p.g1)/pnorm(p.g1)
  formula.eq2.1 <- update.formula(formula.eq2, ~. + imr)
  
  if(margins[2]=="N")               gam2.1 <- eval(substitute(gam(formula.eq2.1, gamma=infl.fac, data=data, weights=weights, subset=inde),                               list(weights=weights,inde=inde)))
  if(margins[2]=="G")               gam2.1 <- eval(substitute(gam(formula.eq2.1, gamma=infl.fac, data=data, weights=weights, subset=inde, family=Gamma(link = "log")),   list(weights=weights,inde=inde)))
  if(!(margins[2] %in% c("N","G"))) gam2.1 <- eval(substitute(gam(formula.eq2.1, gamma=infl.fac, data=data, weights=weights, subset=inde, family=poisson(link = "log")), list(weights=weights,inde=inde)))


  sigma <- sqrt(mean(residuals(gam2.1, type = "deviance")^2)+mean(imr[inde]*(imr[inde]+p.g1[inde]))*gam2.1$coef["imr"]^2)[[1]]
  names(sigma) <- "sigma"
  co  <- (gam2.1$coef["imr"]/sigma)[[1]] # this is an approximate correlation

  if(margins[2]=="G") { k <- (summary(gam2)$dispersion)^(-1); names(k) <- "shape" }
  
  #####################################
  # Smoothing parameter starting values
  #####################################

  if(l.sp1!=0 && l.sp2!=0) sp <- c(gam1$sp, gam2.1$sp)
  if(l.sp1==0 && l.sp2!=0) sp <- c(gam2.1$sp)
  if(l.sp1!=0 && l.sp2==0) sp <- c(gam1$sp)

  #################
  # Starting values
  #################

  a.theta <- st.theta.star(start.theta, co, BivD); names(a.theta) <- "theta.star"

  if(margins[2] %in% c("N", "NB", "PIG")) { 
        if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"],log(sigma),a.theta)
        names(start.v)[length(start.v)-1] <- "log.sigma"
  }
  
  if(margins[2]=="G") {
        if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"],log(k),a.theta)
        names(start.v)[length(start.v)-1] <- "log.shape"
  }
  
  
  if(margins[2] == "P") {
          if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"],a.theta)
  }
  
  
  if(margins[2] == "D") { 
        if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"], log(sigma), qlogis(0.51) , a.theta)
        names(start.v)[length(start.v)-2] <- "log.sigma"
        names(start.v)[length(start.v)-1] <- "logi.nu" 
  }
  
    if(margins[2] == "S") { 
          if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"], log(sigma), 0.5 , a.theta)
          names(start.v)[length(start.v)-2] <- "log.sigma"
          names(start.v)[length(start.v)-1] <- "nu"  
    }
  

  ###################
  # Fitting procedure
  ###################
  
    if(missing(parscale)) parscale <- rep(1, length(start.v))
  
    if(margins[2] %in% c("N", "G") ) funcop <- ghss else funcop <- ghssD  

    respvec <- cbind(y1, y2)
    
    VC <- list(X1 = X1, 
               X2 = X2, 
               X1.d2 = X1.d2, 
               X2.d2 = X2.d2,
               gp1 = gp1, 
               gp2 = gp2,
               pPen1 = pPen1,
               pPen2 = pPen2,
               n = n,
               l.sp1 = l.sp1, 
               l.sp2 = l.sp2,
               infl.fac = infl.fac,
               weights = weights,
               fp = fp,
               BivD = BivD, 
               margins = margins, parscale = parscale)
  
  
    if(method == "trwl")  fit.func <- fit.SemiParSampleSel
    if(method == "trwlF") fit.func <- fit.SemiParSampleSel1
  
    fit <- fit.func(funcop = funcop, start.v = start.v, dat = respvec, VC = VC, qu.mag = qu.mag, sp = sp, 
                    iterlimsp = iterlimsp, pr.tolsp = pr.tolsp, 
                    rinit = rinit, rmax = rmax )

  ############################
  # Post-estimation quantities
  ############################

  epsilon <- sqrt(.Machine$double.eps)
  He <- fit$fit$hessian
  He.eig <- eigen(He, symmetric=TRUE)
  if(min(He.eig$values) < epsilon) He.eig$values[which(He.eig$values < epsilon)] <- 0.0000001
  Vb <- He.eig$vectors%*%tcrossprod(diag(1/He.eig$values),He.eig$vectors)  
  
  if((l.sp1!=0 || l.sp2!=0) && fp==FALSE){ HeSh <- He - fit$fit$S.h
                                           F <- Vb%*%HeSh
                                          }else{ HeSh <- He; F <- diag(rep(1,dim(Vb)[1])) } 
  t.edf <- sum(diag(F))
  
  resp  <- rep(1, n)
  fs    <- as.formula( paste("resp","~",formula.eq2[3], sep="") ) 
  X2s   <- gam(fs, data=data, fit = FALSE)$X # no good indicator here, so we never drop obs.
  param <- fit$fit$argument[-c(1:length(gam1$coef))]
  param <- param[1:length(gam2$coef)]
  eta2  <- X2s%*%param
  
  ###########################################################################
  # Transforming theta back to its original scale and computing Kendall's tau
  ###########################################################################

  theta <- theta.tau(BivD = BivD, theta.star = fit$fit$argument["theta.star"])
  KendTau <- theta[2] # we probably do not want kendall's tau in output
  theta   <- theta[1]

  names(theta) <- names(KendTau) <- NULL

  # dispersion: variance for normal case 
  # dispersion: 1/k for gamma case

  if(margins[2] %in% c("N", "NB", "PIG", "D", "S") ) { sigma <- exp(fit$fit$argument["log.sigma"]); names(sigma) <- k <- NULL; phi <- sigma^2 } 

  if(margins[2] == "P") { sigma <- k <- phi <- nu <- NULL} 

  if(margins[2]=="G") { k <- exp(fit$fit$argument["log.shape"]); names(k) <- sigma <- NULL; phi <- k^{-1} }   

  if(margins[2] == "S"){nu <- fit$fit$argument["nu"]; names(nu) <- NULL}
  
  if(margins[2] == "D"){nu <- plogis(fit$fit$argument["logi.nu"]); names(nu) <- NULL }


  ###########################################################################
  # Final object
  ###########################################################################

       L <- list(fit = fit$fit, gam1 = gam1, gam2 = gam2, gam2.1 = gam2.1, 
                 coefficients = fit$fit$argument, weights = weights, 
                 sp = fit$sp, 
                 iter.if = fit$iter.if, iter.sp = fit$iter.sp, iter.inner = fit$iter.inner, 
                 start.v = start.v, phi = phi, sigma = sigma, shape = k, nu = nu, 
                 theta = theta, tau = KendTau, 
                 n = n, n.sel = length(gam2$y), 
                 X1 = X1, X2 = X2, X1.d2 = X1.d2, X2.d2 = X2.d2, 
                 l.sp1 = l.sp1, l.sp2 = l.sp2, 
                 He = He, HeSh = HeSh, Vb = Vb, F = F, 
                 BivD = BivD, margins = margins,
                 t.edf = t.edf, bs.mgfit = fit$bs.mgfit, 
                 conv.sp = fit$conv.sp, wor.c = fit$wor.c, pPen1 = pPen1, pPen2 = pPen2,  
                 eta1 = fit$fit$eta1, eta2 = eta2, 
                 y1 = y1, y2 = y2, logLik = -fit$fit$l, fp = fp,
                 gp1 = gp1, gp2 = gp2, X2s = X2s, magpp = fit$magpp)

  class(L) <- "SemiParSampleSel"

  L

}
