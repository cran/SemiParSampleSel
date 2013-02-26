ghssFGM2 <- function(params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, params.f, sp=NULL, qu.mag=NULL){

  eta1 <- dat1%*%params.f[1:X1.d2]  
  eta2 <- dat2%*%params.f[(X1.d2+1):(X1.d2+X2.d2)]   
  sqv.st <- params[1]
  teta.st <- params[2]
  sqv <- exp(sqv.st)
  teta <- tanh(teta.st)
  

# Settings:

  i0 <- 1-dat[, 1]
  i1 <- dat[, 1]

  e2 <- (dat[,2]-eta2)/sqv
  F1 <- pnorm(-eta1)
  F2 <- pnorm(e2)
  f2 <- dnorm(e2)/sqv
  u <- 1-teta*F1*(1-2*F2)
  h <- 2*teta*F1*f2*u^(-1)


# Likelihood:

  l.par <- weights*(i0*log(F1) + i1*(log(1-F1)+log(1-teta*F1*(1-2*F2))+log(f2))) 
  

# Gradient:

  dl.dsqv.st <- weights*i1*(-h*e2*sqv+e2^2-1)
  dl.dteta.st <- weights*i1*(1-u^(-1))*(1-teta^2)*teta^(-1)


# Hessian:

  d2l.sqv.st.sqv.st <- -weights*i1*(-h*u^(-1)*((e2^2-1)*u*sqv^(-1) + h*u*e2)*sqv^2-2*e2)*e2
  d2l.sqv.st.teta.st <- weights*i1*2*F1*f2*u^(-2)*e2*sqv*(1-teta^2)

  d2l.teta.st.teta.st <- -weights*i1*(u-1)*u^(-1)*((1-teta^2)*teta^(-2)*u^(-1)-(1+teta^2)*teta^(-2))*(1-teta^2)


  H <- rbind( 
              cbind( sum(d2l.sqv.st.sqv.st), sum(d2l.sqv.st.teta.st) ) ,
              cbind( sum(d2l.sqv.st.teta.st), sum(d2l.teta.st.teta.st) )
            ) 

  
  res <- -sum(l.par)
  G   <- c( -sum( dl.dsqv.st ) ,  
            -sum( dl.dteta.st )   )


  list(value=res, gradient=G, hessian=H, l=res, dat1=dat1, dat2=dat2) 
     
}


