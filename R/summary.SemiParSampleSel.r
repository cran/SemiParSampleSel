summary.SemiParSampleSel <- function(object,n.sim=1000,s.meth="svd",sig.lev=0.05,...){

  testStat <- function (p, X, V, rank = NULL) {
      qrx <- qr(X)
      R <- qr.R(qrx)
      V <- R %*% tcrossprod(V[qrx$pivot, qrx$pivot],R)
      V <- (V + t(V))/2
      ed <- eigen(V, symmetric = TRUE)
      k <- max(0, floor(rank))
      nu <- abs(rank - k)

          if (rank > k + 0.05 || k == 0) 
              k <- k + 1
          nu <- 0
          rank <- k
      
      if (nu > 0) 
          k1 <- k + 1
      else k1 <- k
      r.est <- sum(ed$values > max(ed$values) * .Machine$double.eps^0.9)
      if (r.est < k1) {
          k1 <- k <- r.est
          nu <- 0
          rank <- r.est
      }
      vec <- ed$vectors
      if (k1 < ncol(vec)) 
          vec <- vec[, 1:k1, drop = FALSE]
      if (k == 0) {
          vec <- t(t(vec) * sqrt(nu/ed$val[1]))
      }
      if (nu > 0 && k > 0) {
          if (k > 1) 
              vec[, 1:(k - 1)] <- t(t(vec[, 1:(k - 1)])/sqrt(ed$val[1:(k - 
                  1)]))
          b12 <- 0.5 * nu * (1 - nu)
          if (b12 < 0) 
              b12 <- 0
          b12 <- sqrt(b12)
          B <- matrix(c(1, b12, b12, nu), 2, 2)
          ev <- diag(ed$values[k:k1]^-0.5)
          B <- ev %*% B %*% ev
          eb <- eigen(B, symmetric = TRUE)
          rB <- eb$vectors %*% tcrossprod(diag(sqrt(eb$values)),eb$vectors)
          vec[, k:k1] <- t(tcrossprod(rB,vec[, k:k1]))
      }
      else {
          vec <- t(t(vec)/sqrt(ed$val[1:k]))
      }
      d <- crossprod(vec,R%*%p)
      d <- sum(d^2)
      attr(d, "rank") <- rank
      d
  }

  F  <- object$F
  Vr <- object$Vb
          
  SE <- sqrt(diag(Vr))
  n  <- object$n 

  bs <- rmvnorm(n.sim, mean = coef(object), sigma=object$Vb, method=s.meth)
  d.sig <- dim(object$Vb)[1]-1
  d.rho <- dim(object$Vb)[1]
  est.SIGb <- est.THETAb <- est.KeTb <- rep(NA,n.sim)

  for(i in 1:n.sim){ 
      est.SIGb[i] <- exp(bs[i,d.sig])
        if(object$BivD=="N")   { est.THETAb[i] <- tanh(bs[i,d.rho]);	est.KeTb[i] <- tau(normalCopula(est.THETAb[i]))	}
        if(object$BivD=="C")   { est.THETAb[i] <- exp(bs[i,d.rho]);	est.KeTb[i] <- tau(claytonCopula(est.THETAb[i]))}
        if(object$BivD=="J")   { est.THETAb[i] <- 1+exp(bs[i,d.rho]);	est.KeTb[i] <- tau(joeCopula(est.THETAb[i]))	}  
	if(object$BivD=="FGM") { est.THETAb[i] <- tanh(bs[i,d.rho]);	est.KeTb[i] <- tau(fgmCopula(est.THETAb[i]))	}
	if(object$BivD=="F")   { est.THETAb[i] <- bs[i,d.rho];		est.KeTb[i] <- tau(frankCopula(est.THETAb[i]))	}
	if(object$BivD=="AMH") { est.THETAb[i] <- tanh(bs[i,d.rho]);	est.KeTb[i] <- tau(amhCopula(est.THETAb[i]))	}
        if(object$BivD=="G")   { est.THETAb[i] <- 1+exp(bs[i,d.rho]);	est.KeTb[i] <- tau(gumbelCopula(est.THETAb[i]))	}  
  }

  CIphi <- as.numeric(quantile(est.SIGb,c(sig.lev/2,1-sig.lev/2),na.rm=TRUE))
  CIth  <- as.numeric(quantile(est.THETAb,c(sig.lev/2,1-sig.lev/2),na.rm=TRUE))
  CIkt  <- as.numeric(quantile(est.KeTb,c(sig.lev/2,1-sig.lev/2),na.rm=TRUE))

  tableN <- list(NULL,NULL)
  table  <- list()
  
  ind <- list(ind1=1:1:object$gam1$nsdf, ind2=object$X1.d2+(1:object$gam2$nsdf) )

  for(i in 1:2){
       estimate <- coef(object)[ind[[i]]]
       se       <- SE[ind[[i]]]
       ratio    <- estimate/se
       pv       <- 2*pnorm(abs(ratio), lower.tail = FALSE)
       table[[i]] <- cbind(estimate,se,ratio,pv)
       dimnames(table[[i]])[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  }


  if((object$l.sp1!=0 || object$l.sp2!=0)){
  
    	pTerms.df <- pTerms.chi.sq <- pTerms.pv <- edf <- tableN <- list(0,0)
          
             for(i in 1:2){

             if(i==1) {mm <- object$l.sp1; if(mm==0) next}
             if(i==2) {mm <- object$l.sp2; if(mm==0) break}
    
  		for(k in 1:mm){
  
                        if(i==1){gam <- object$gam1; ind <- (gam$smooth[[k]]$first.para):(gam$smooth[[k]]$last.para)} 
                            else{gam <- object$gam2; ind <- (gam$smooth[[k]]$first.para:gam$smooth[[k]]$last.para)+object$X1.d2}
  			edf[[i]][k] <- sum(diag(F)[ind])
  			names(edf[[i]])[k] <- gam$smooth[[k]]$label 
  			b  <- coef(object)[ind]
  			V  <- Vr[ind,ind]
  			if(i==1) Xt <- object$X1[, 1:length(ind)+gam$nsdf] else Xt <- object$X2[, 1:length(ind)+gam$nsdf]
  			pTerms.df[[i]][k] <- min(ncol(Xt), edf[[i]][k])
  			pTerms.chi.sq[[i]][k] <- Tp <- testStat(b, Xt, V, pTerms.df[[i]][k])
  			pTerms.df[[i]][k] <- attr(Tp, "rank")
                        pTerms.pv[[i]][k] <- pchisq(pTerms.chi.sq[[i]][k], df = pTerms.df[[i]][k], lower.tail = FALSE)
  			                 
                }
                tableN[[i]] <- cbind(edf[[i]], pTerms.df[[i]], pTerms.chi.sq[[i]], pTerms.pv[[i]])
                dimnames(tableN[[i]])[[2]] <- c("edf", "Est.rank", "Chi.sq", "p-value")
             }
  
  }
  

     res <- list(tableP1=table[[1]], tableP2=table[[2]], 
                 tableNP1=tableN[[1]], tableNP2=tableN[[2]], 
                 n=n, phi=object$phi, sigma=object$sigma, shape=object$shape, theta=object$theta, tau=object$tau, 
                 formula1=object$gam1$formula, formula2=object$gam2$formula, 
                 l.sp1=object$l.sp1, l.sp2=object$l.sp2, 
                 t.edf=object$t.edf, CIsig=CIphi, CIshape=CIphi, CIth=CIth, CIkt=CIkt, BivD=object$BivD, margins=object$margins, n.sel=object$n.sel)
  

  class(res) <- "summary.SemiParSampleSel"

  res


}

