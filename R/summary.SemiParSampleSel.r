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
  Vf <- object$Vb
          
  SE <- sqrt(diag(Vf))
  n  <- object$n 

  bs <- rmvnorm(n.sim, mean = object$fit$argument, sigma=object$Vb, method=s.meth)
  d.sig <- dim(object$Vb)[1]-1
  d.rho <- dim(object$Vb)[1]
  est.SIGb <- est.RHOb <- rep(NA,n.sim)
  for(i in 1:n.sim){ est.SIGb[i] <- exp(bs[i,d.sig]); est.RHOb[i] <- tanh(bs[i,d.rho]) }
  CIsi <- as.numeric(quantile(est.SIGb,c(sig.lev/2,1-sig.lev/2)))
  CIrs <- as.numeric(quantile(est.RHOb,c(sig.lev/2,1-sig.lev/2)))

  estimate1 <- object$fit$argument[1:object$gam1$nsdf]
  se1       <- SE[1:object$gam1$nsdf]
  ratio1    <- estimate1/se1
  pv1       <- 2*pnorm(abs(ratio1), lower.tail = FALSE)
  table1    <- cbind(estimate1,se1,ratio1,pv1)

  estimate2 <- object$fit$argument[object$X1.d2+(1:object$gam2$nsdf)]
  se2       <- SE[object$X1.d2+(1:object$gam2$nsdf)]
  ratio2    <- estimate2/se2
  pv2       <- 2*pnorm(abs(ratio2), lower.tail = FALSE)
  table2    <- cbind(estimate2,se2,ratio2,pv2)

  dimnames(table1)[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  dimnames(table2)[[2]] <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")

  if(object$l.sp1!=0 && object$l.sp2!=0){
  	pTerms.df1 <- pTerms.chi.sq1 <- pTerms.pv1 <- edf1 <- NA
  	pTerms.df2 <- pTerms.chi.sq2 <- pTerms.pv2 <- edf2 <- NA
		for(k in 1:length(object$gam1$smooth)){
			ind <- object$gam1$smooth[[k]]$first.para:object$gam1$smooth[[k]]$last.para
			edf1[k] <- sum(diag(F)[ind])
			names(edf1)[k] <- object$gam1$smooth[[k]]$label 
			b  <- object$fit$argument[ind]
			V  <- Vf[ind,ind]
			Xt <- object$X1[, 1:length(ind)+object$gam1$nsdf]
			pTerms.df1[k] <- min(ncol(Xt), edf1[k])
			pTerms.chi.sq1[k] <- Tp <- testStat(b, Xt, V, pTerms.df1[k])
			pTerms.df1[k] <- attr(Tp, "rank")
                        pTerms.pv1[k] <- pchisq(pTerms.chi.sq1[k], df = pTerms.df1[k], lower.tail = FALSE)			
            }
  	table1.1 <- cbind(edf1, pTerms.df1, pTerms.chi.sq1, pTerms.pv1)

		for(k in 1:length(object$gam2$smooth)){
			ind <- (object$gam2$smooth[[k]]$first.para:object$gam2$smooth[[k]]$last.para)+object$X1.d2
			edf2[k] <- sum(diag(F)[ind])
			names(edf2)[k] <- object$gam2$smooth[[k]]$label 
			b  <- object$fit$argument[ind]
			V  <- Vf[ind,ind]
			Xt <- object$X2[, 1:length(ind)+object$gam2$nsdf]
			pTerms.df2[k] <- min(ncol(Xt), edf2[k])
			pTerms.chi.sq2[k] <- Tp <- testStat(b, Xt, V, pTerms.df2[k])   
			pTerms.df2[k] <- attr(Tp, "rank")
			pTerms.pv2[k] <- pchisq(pTerms.chi.sq2[k], df = pTerms.df2[k], lower.tail = FALSE)
		}
	table2.2 <- cbind(edf2, pTerms.df2, pTerms.chi.sq2, pTerms.pv2)
  dimnames(table1.1)[[2]] <- c("edf", "Est.rank", "Chi.sq", "p-value")
  dimnames(table2.2)[[2]] <- c("edf", "Est.rank", "Chi.sq", "p-value")
  }


  if(object$l.sp1!=0 && object$l.sp2!=0){res <- list(tableP1=table1, tableP2=table2, 
                                           tableNP1=table1.1, tableNP2=table2.2, 
                                           n=n, sigma=object$sigma, rho=object$rho, 
                                           formula1=object$gam1$formula, formula2=object$gam2$formula, 
                                           l.sp1=object$l.sp1, l.sp2=object$l.sp2, 
                                           t.edf=object$t.edf, CIsi=CIsi, CIrs=CIrs, n.sel=object$n.sel)
                               class(res) <- "summary.SemiParSampleSel"
                               res
  }
  else{res <- list(tableP1=table1,tableP2=table2, 
                   n=n, sigma=object$sigma, rho=object$rho, 
                   formula1=object$gam1$formula, formula2=object$gam2$formula, 
                   l.sp1=0, l.sp2=0, 
                   t.edf=object$t.edf, CIsi=CIsi, CIrs=CIrs, n.sel=object$n.sel)
       class(res) <- "summary.SemiParSampleSel"
       res
  }

}

