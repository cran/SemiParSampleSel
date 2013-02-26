ghssJ1 <- function(params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, sp=NULL, l.sp1, l.sp2, gp1, gp2, exl1, exl2, qu.mag=NULL){

  eta1 <- dat1%*%params[1:X1.d2]  
  eta2 <- dat2%*%params[(X1.d2+1):(X1.d2+X2.d2)]   
  sqv.st <- exl1
  teta.st <- exl2
  sqv <- exp(sqv.st)
  teta <- 1+exp(teta.st)
  
# Settings:

  i0 <- 1-dat[, 1]
  i1 <- dat[, 1]

  e2 <- (dat[,2]-eta2)/sqv
  F1 <- pnorm(-eta1)		
  F1b <- 1-F1
  F1b <- ifelse(F1b>0.0000000001,F1b,0.0000000001)
  F2b <- 1-pnorm(e2)
  F2b <- ifelse(F2b>0.0000000001,F2b,0.0000000001)
  f2 <- dnorm(e2)/sqv
  f2 <- ifelse(f2>0.0000000001,f2,0.0000000001)
  lnF1b <- log(F1b)
  lnF2b <- log(F2b)
  teta1 <- teta-1
  u <- F1b^teta+F2b^teta-(F1b*F2b)^teta
  ut <- u^(1/teta-1)
  z <- (1-F1b^teta)*(F2b^teta1)*ut
  z <- ifelse(z<0.9999999999,z,0.9999999999)
  zt <- (z-1)^(-1)
  p <- F2b^teta1*ut*(u+teta1)*teta^(-1)*u^(-1)
  d <- teta1*F1b^teta*u^(-1)
  h <- z*zt*f2*F2b^(-1)*d
  E <- f2*F2b^(-1)*( teta*u^(-1)*F1b^teta-teta-1-d*zt ) + e2/sqv
  ph <- dnorm(-eta1)
  A <- -i1*zt*teta*teta1*F1b^teta1*F2b^(-1)*f2*(
	p*(1-z*zt*F1b^teta*u^(-1))-F2b^teta*(z/u)*(2*teta+u-1)/(teta*u)
	)*ph


  settings <- cbind(e2,F1b,F2b,f2,lnF1b,lnF2b,u,ut,z,zt,p,ph,d,E,h,A)
  if(any(abs(settings)=='Inf') | any(settings=='NaN') | any(settings=='NA')) stop("Ill-conditioned task.")


# Likelihood:

  l.par <- weights*(i0*log(F1) + i1*(log(1-z)+log(f2))) 

 
# Gradient:

  dl.dbe1 <- weights*(-i0*F1^(-1) - i1*teta*zt*F1b^teta1*p)*ph
  dl.dbe2 <- weights*i1*(h+e2/sqv)


# Hessian:

  d2l.be1.be1 <- -weights*( (ph*( -i0*(ph/F1-eta1)/F1 -
		i1*teta*zt*(F1b^(teta-2))*(
			ph*( p*(teta1+teta*p*F1b^teta*zt) - teta1*u^(-1)*(u-F2b^teta)*(2*p+z*(teta*u)^(-1)*(1-F2b^teta)) )
			-F1b*p*eta1
		))) )
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


