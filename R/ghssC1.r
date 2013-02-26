ghssC1 <- function(params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, sp=NULL, l.sp1, l.sp2, gp1, gp2, exl1, exl2, qu.mag=NULL){

  eta1 <- dat1%*%params[1:X1.d2]  
  eta2 <- dat2%*%params[(X1.d2+1):(X1.d2+X2.d2)]   
  sqv.st <- exl1
  teta.st <- exl2
  sqv <- exp(sqv.st)
  teta <- exp(teta.st)

  i0 <- 1-dat[, 1]
  i1 <- dat[, 1]

  e2 <- (dat[,2]-eta2)/sqv
  F1 <- pnorm(-eta1);
  F1 <- ifelse(F1>0.0000000001,F1,0.0000000001)
  F2 <- pnorm(e2);
  F2 <- ifelse(F2>0.0000000001,F2,0.0000000001)
  f2 <- dnorm(e2)/sqv;
  f2 <- ifelse(f2>0.0000000001,f2,0.0000000001)
  lnF1 <- log(F1)
  teta1 <- teta+1
  u <- F1^(-teta)+F2^(-teta)-1
  z <- F2^(-teta1)*u^(-teta1/teta);
  z <- ifelse(z<0.9999999999,z,0.9999999999)
  zt <- z/(1-z)
  P <- F1^(-teta1)*u^(-1)*zt
  K <- f2*F2^(-teta1)*(F2^teta-u^(-1))
  ph <- dnorm(-eta1)
  d <- P*ph*i1
  A <- teta1*d*(K*teta1/(1-z)-teta*f2*F2^(-teta1)*u^(-1))
  E <- K*teta1/(1-z)+e2/sqv+f2*F2^(-teta1)*(F2^teta-teta/u)
  h <- -zt*teta1*K

  settings <- cbind(e2,F1,F2,f2,lnF1,u,z,zt,P,K,ph,d,A,E,h)
  if(any(abs(settings)=='Inf') | any(settings=='NaN') | any(settings=='NA')) stop("Ill-conditioned task.")


# Likelihood:

  l.par <- weights*(i0*lnF1 + i1*(log(1-z)+log(f2))) 

 
# Gradient:

  dl.dbe1 <- weights*(-i0*F1^(-1)+teta1*i1*P)*ph
  dl.dbe2 <- weights*i1*(h+e2/sqv)


# Hessian:

  d2l.be1.be1 <- -weights*( (ph*(-i0*(ph/F1-eta1)/F1+
		teta1*i1*P*(ph*(teta1-F1^(-teta)*(teta+teta1/(1-z))/u)/F1-eta1))) )
  d2l.be1.be2 <- -weights*A
  d2l.be2.be2 <- -weights*i1*(h*E-sqv^(-2))
  
  be1.be1 <- crossprod(dat1*c(d2l.be1.be1),dat1)
  be1.be2 <- crossprod(dat1*c(d2l.be1.be2),dat2) 
  be2.be2 <- crossprod(dat2*c(d2l.be2.be2),dat2) 
  

  H <- rbind( cbind( be1.be1    , be1.be2 ), 
              cbind( t(be1.be2) , be2.be2 )   ) 


  res <- -sum(l.par)
  G   <- c( -colSums( c(dl.dbe1)*dat1 ) ,
                   -colSums( c(dl.dbe2)*dat2 )   )


      S <- mapply("*", qu.mag$Ss, sp, SIMPLIFY=FALSE)
      S <- do.call(adiag, lapply(S, unlist))

      if(l.sp1!=0 && l.sp2!=0) S.h <- adiag(matrix(0,gp1,gp1),
                                            S[1:(X1.d2-gp1),1:(X1.d2-gp1)],
                                            matrix(0,gp2,gp2),
                                            S[(X1.d2-(gp1-1)):dim(S)[2],(X1.d2-(gp1-1)):dim(S)[2]])

      if(l.sp1==0 && l.sp2!=0) S.h <- adiag(matrix(0,gp1,gp1), matrix(0,gp2,gp2), S)
      if(l.sp1!=0 && l.sp2==0) S.h <- adiag(matrix(0,gp1,gp1), S, matrix(0,gp2,gp2))

  S.h1 <- 0.5*crossprod(params,S.h)%*%params
  S.h2 <- S.h%*%params

  S.res <- res
  res <- S.res + S.h1
  G   <- G + S.h2
  H   <- H + S.h  


 list(value=res, gradient=G, hessian=H, S.h=S.h, l=S.res, eta1=eta1, eta2=eta2, dat1=dat1, dat2=dat2, 
      dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, 
      d2l.be1.be1=d2l.be1.be1, d2l.be1.be2=d2l.be1.be2, d2l.be2.be2=d2l.be2.be2)
    

}


