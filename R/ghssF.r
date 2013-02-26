ghssF <- function(params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, sp=NULL, l.sp1, l.sp2, gp1, gp2, fp, qu.mag=NULL){


  eta1 <- dat1%*%params[1:X1.d2]  
  eta2 <- dat2%*%params[(X1.d2+1):(X1.d2+X2.d2)]   
  sqv.st <- params[(X1.d2+X2.d2+1)]
  teta <- params[(X1.d2+X2.d2+2)]
  sqv <- exp(sqv.st)

  
# Settings:

  i0 <- 1-dat[, 1]
  i1 <- dat[, 1]

  e2 <- (dat[,2]-eta2)/sqv
  F1 <- pnorm(-eta1);
  F1 <- ifelse(F1>0.0000000001,F1,0.0000000001)
  F2 <- pnorm(e2);
  F2 <- ifelse(F2>0.0000000001,F2,0.0000000001)
  f2 <- dnorm(e2)/sqv;
  f2 <- ifelse(f2>0.0000000001,f2,0.0000000001)
  u <- exp(teta*(F1+F2))-exp(teta*(1+F2))
  z <- u/(u-exp(teta*(1+F1))+exp(teta))
  p <- z*exp(teta*(1+F1+F2))*u^(-2)
  b <- (z^2)*(u^(-1))*((1-z)*z^(-1)*(u*(F1+F2-1)+(F1-1)*exp(teta*(1+F2))) + F1*exp(teta*(1+F1)))
  h <- z*u^(-1)*f2*teta*exp(teta)*(exp(teta*F1)-1)
  ph <- dnorm(-eta1)
  A <- i1*(1-exp(teta))*(teta^2)*f2*(p^2)*u*(exp(-teta*F1)-exp(-teta))*ph
  E <- teta*z*f2+e2/sqv
  B <- i1*f2*exp(teta)*z*u^(-1)*( (1+teta+teta*F1)*exp(teta*F1)-teta-1-
		teta*(exp(teta*F1)-1)*z*u^(-1)*(u*(F1+F2)+(F1-1)*exp(teta*(1+F2))-(1+F1)*exp(teta*(1+F1))+exp(teta) ) )


  settings <- cbind(e2,F1,F2,f2,u,z,p,ph,A,b,B,E,h)
  if(any(abs(settings)=='Inf') | any(settings=='NaN') | any(settings=='NA')) stop("Ill-conditioned task.")


# Likelihood:

  l.par <- weights*(i0*log(F1) + i1*(log(z)+log(f2))) 


# Gradient:

  dl.dbe1 <- weights*(-i0*F1^(-1)-i1*(1-exp(teta))*teta*p)*ph
  dl.dbe2 <- weights*i1*(h+e2/sqv)
  dl.dsqv.st <- weights*i1*(h*e2*sqv+e2^2-1)
  dl.dteta <- weights*i1*b*z^(-1)


# Hessian:

  d2l.be1.be1 <- -weights*( ph*(-i0*(ph/F1-eta1)/F1-
	teta*(1-exp(teta))*i1*p*(ph*teta*p*(exp(teta*F1)*(exp(teta*(F2-1))-1)-exp(teta*(1-F1))*(exp(teta*F2)-1))-eta1)) )
  d2l.be1.be2 <- -weights*A
  d2l.be1.sqv.st <- -weights*A*e2*sqv
  d2l.be1.teta <- -weights*( i1*p*( (1+teta)*exp(teta)-1+teta*(1-exp(teta))*p*(
		(F1+F2-1)*exp(teta*(F1+F2-1))-2*F2*exp(teta*F2)-F1*exp(teta*F1)+
		(1-F1+F2)*exp(teta*(1-F1+F2))-(1-F1)*exp(teta*(1-F1))+exp(teta)
	) )*ph )

  d2l.be2.be2 <- -weights*i1*(h*E-sqv^(-2))
  d2l.be2.sqv.st <- -weights*i1*(h*(E*e2-sqv^(-1))-2*e2/sqv^2)*sqv
  d2l.be2.teta <- -weights*B

  d2l.sqv.st.sqv.st <- -weights*i1*(sqv^2)*e2*(h*(E*e2-sqv^(-1))-2*e2*(sqv^(-2)))
  d2l.sqv.st.teta <- -weights*e2*B*sqv

  d2l.teta.teta <- -weights*( i1*( b*z^(-1)*( b*z^(-1)-u^(-1)*(u*(F1+F2)+exp(teta*(1+F2))*(F1-1)) )
		+z*u^(-1)*F1*(1+F1)*exp(teta*(1+F1))+(F1+F2-1)*((1-z)*(F1+F2)-b*z^(-1))
		+u^(-1)*exp(teta*(1+F2))*(F1-1)*( (1-z)*(F1+2*F2)-b*z^(-1) )
		) )


  be1.be1 <- crossprod(dat1*c(d2l.be1.be1),dat1)
  be1.be2 <- crossprod(dat1*c(d2l.be1.be2),dat2) 
  be1.sqv.st <- t(t(rowSums(t(dat1*c(d2l.be1.sqv.st)))))
  be1.teta <- t(t(rowSums(t(dat1*c(d2l.be1.teta)))))

  be2.be2 <- crossprod(dat2*c(d2l.be2.be2),dat2) 
  be2.sqv.st <- t(t(rowSums(t(dat2*c(d2l.be2.sqv.st)))))
  be2.teta <- t(t(rowSums(t(dat2*c(d2l.be2.teta)))))


  H <- rbind( cbind( be1.be1    , be1.be2    , be1.sqv.st  ,  be1.teta ), 
              cbind( t(be1.be2) , be2.be2    , be2.sqv.st  ,  be2.teta ),
              cbind( t(be1.sqv.st) , t(be2.sqv.st) , sum(d2l.sqv.st.sqv.st), sum(d2l.sqv.st.teta) ) ,
              cbind( t(be1.teta) , t(be2.teta) , sum(d2l.sqv.st.teta), sum(d2l.teta.teta) )
            )



  res <- -sum(l.par)
  G   <- c( -colSums( c(dl.dbe1)*dat1 ) ,
            -colSums( c(dl.dbe2)*dat2 )    ,
            -sum( dl.dsqv.st ) ,  
            -sum( dl.dteta )   )


  if( ( l.sp1==0 && l.sp2==0 ) || fp==TRUE) S.h <- S.h1 <- S.h2 <- 0
  else{

      S <- mapply("*", qu.mag$Ss, sp, SIMPLIFY=FALSE)
      S <- do.call(adiag, lapply(S, unlist))

      if(l.sp1!=0 && l.sp2!=0) S.h <- adiag(matrix(0,gp1,gp1),
                                            S[1:(X1.d2-gp1),1:(X1.d2-gp1)],
                                            matrix(0,gp2,gp2),
                                            S[(X1.d2-(gp1-1)):dim(S)[2],(X1.d2-(gp1-1)):dim(S)[2]],
                                            0,0)

      if(l.sp1==0 && l.sp2!=0) S.h <- adiag(matrix(0,gp1,gp1), matrix(0,gp2,gp2), S, 0, 0)
      if(l.sp1!=0 && l.sp2==0) S.h <- adiag(matrix(0,gp1,gp1), S, matrix(0,gp2,gp2), 0, 0)
      

     S.h1 <- 0.5*crossprod(params,S.h)%*%params
     S.h2 <- S.h%*%params
  }

     
    
  S.res <- res
  res <- S.res + S.h1
  G   <- G + S.h2
  H   <- H + S.h  


  list(value=res, gradient=G, hessian=H, S.h=S.h, l=S.res, eta1=eta1, eta2=eta2, dat1=dat1, dat2=dat2,
       dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2,
       dl.dsqv.st=dl.dsqv.st,
       dl.dteta=dl.dteta, 
       d2l.be1.be1=d2l.be1.be1, d2l.be1.be2=d2l.be1.be2, d2l.be2.be2=d2l.be2.be2,
       d2l.be1.sqv.st=d2l.be1.sqv.st,
       d2l.be1.cor.st=d2l.be1.teta,
       d2l.be2.sqv.st=d2l.be2.sqv.st, 
       d2l.be2.cor.st=d2l.be2.teta,
       d2l.sqv.st.sqv.st=d2l.sqv.st.sqv.st,
       d2l.sqv.st.cor.st=d2l.sqv.st.teta,    
       d2l.cor.st.cor.st=d2l.teta.teta)

}

