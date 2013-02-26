ghssF2 <- function(params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, params.f, sp=NULL, qu.mag=NULL){

  eta1 <- dat1%*%params.f[1:X1.d2]  
  eta2 <- dat2%*%params.f[(X1.d2+1):(X1.d2+X2.d2)]   
  sqv.st <- params[1]
  teta <- params[2]
  sqv <- exp(sqv.st)

  
# Settings:

  i0 <- 1-dat[, 1]
  i1 <- dat[, 1]

  e2 <- (dat[,2]-eta2)/sqv
  F1 <- pnorm(-eta1)
  F1 <- ifelse(F1>0.0000000001,F1,0.0000000001)
  F2 <- pnorm(e2)
  F2 <- ifelse(F2>0.0000000001,F2,0.0000000001)
  f2 <- dnorm(e2)/sqv
  f2 <- ifelse(f2>0.0000000001,f2,0.0000000001)
  u <- exp(teta*(F1+F2))-exp(teta*(1+F2))
  z <- u/(u-exp(teta*(1+F1))+exp(teta))
  b <- (z^2)*(u^(-1))*((1-z)*z^(-1)*(u*(F1+F2-1)+(F1-1)*exp(teta*(1+F2))) + F1*exp(teta*(1+F1)))
  h <- z*u^(-1)*f2*teta*exp(teta)*(exp(teta*F1)-1)
  ph <- dnorm(-eta1)
  E <- teta*z*f2+e2/sqv
  B <- i1*f2*exp(teta)*z*u^(-1)*( (1+teta+teta*F1)*exp(teta*F1)-teta-1-
		teta*(exp(teta*F1)-1)*z*u^(-1)*(u*(F1+F2)+(F1-1)*exp(teta*(1+F2))-(1+F1)*exp(teta*(1+F1))+exp(teta) ) )


  settings <- cbind(e2,F1,F2,f2,u,z,ph,b,B,E,h)
  if(any(abs(settings)=='Inf') | any(settings=='NaN') | any(settings=='NA')) stop("Ill-conditioned task.")


# Likelihood:

  l.par <- weights*(i0*log(F1) + i1*(log(z)+log(f2))) 
  

# Gradient:

  dl.dsqv.st <- weights*i1*(h*e2*sqv+e2^2-1)
  dl.dteta <- weights*i1*b*z^(-1)


# Hessian:

  d2l.sqv.st.sqv.st <- -weights*i1*(sqv^2)*e2*(h*(E*e2-sqv^(-1))-2*e2*(sqv^(-2)))
  d2l.sqv.st.teta <- -weights*e2*B*sqv

  d2l.teta.teta <- -weights*( i1*( b*z^(-1)*( b*z^(-1)-u^(-1)*(u*(F1+F2)+exp(teta*(1+F2))*(F1-1)) )
		+z*u^(-1)*F1*(1+F1)*exp(teta*(1+F1))+(F1+F2-1)*((1-z)*(F1+F2)-b*z^(-1))
		+u^(-1)*exp(teta*(1+F2))*(F1-1)*( (1-z)*(F1+2*F2)-b*z^(-1) )
		) )


  H <- rbind( 
              cbind( sum(d2l.sqv.st.sqv.st), sum(d2l.sqv.st.teta) ) ,
              cbind( sum(d2l.sqv.st.teta), sum(d2l.teta.teta) )
            ) 

  
  res <- -sum(l.par)
  G   <- c( -sum( dl.dsqv.st ) ,  
            -sum( dl.dteta )   )
  H   <- H

  list(value=res, gradient=G, hessian=H, l=res, dat1=dat1, dat2=dat2) 
     
}


