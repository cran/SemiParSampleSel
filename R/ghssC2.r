ghssC2 <- function(params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, params.f, sp=NULL, qu.mag=NULL){


  eta1 <- dat1%*%params.f[1:X1.d2]  
  eta2 <- dat2%*%params.f[(X1.d2+1):(X1.d2+X2.d2)]   
  sqv.st <- params[1]
  teta.st <- params[2]
  sqv <- exp(sqv.st)
  teta <- exp(teta.st)


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
  lnF1 <- log(F1)
  lnF2 <- log(F2)
  teta1 <- teta+1
  u <- F1^(-teta)+F2^(-teta)-1
  lnu <- log(u)
  z <- F2^(-teta1)*u^(-teta1/teta);
  z <- ifelse(z<0.9999999999,z,0.9999999999)
  zt <- z/(1-z)
  K <- f2*F2^(-teta1)*(F2^teta-u^(-1))
  C <- F1^(-teta)*lnF1+F2^(-teta)*lnF2
  b <- (lnF2-lnu/teta^2-C*teta1/(teta*u))/(z-1)
  B <- -zt*i1*(K+teta1*(K*b+f2*F2^(-teta1)*u^(-2)*(u*lnF2-C)))
  E <- K*teta1/(1-z)+e2/sqv+f2*F2^(-teta1)*(F2^teta-teta/u)
  h <- -zt*teta1*K

  settings <- cbind(e2,F1,F2,f2,lnF1,lnF2,u,lnu,z,zt,K,C,b,B,E,h)
  if(any(abs(settings)=='Inf') | any(settings=='NaN') | any(settings=='NA')) stop("Ill-conditioned task.")


# Likelihood:

  l.par <- weights*(i0*lnF1 + i1*(log(1-z)+log(f2))) 
  

# Gradient:

  dl.dsqv.st <- weights*i1*(h*e2*sqv+e2^2-1)
  dl.dteta.st <- weights*i1*zt*(lnF2*teta-teta^(-1)*lnu-C*teta1/u)


# Hessian:

  d2l.sqv.st.sqv.st <- -weights*i1*(h*(E*e2-sqv^(-1))-2*e2/sqv^2)*e2*sqv^2
  d2l.sqv.st.teta.st <- -weights*e2*B*sqv*teta

  d2l.teta.st.teta.st <- -weights*( i1*zt*((z-1)*teta^2*b^2+teta*lnF2+lnu/teta+C*(1-teta)/u
		-teta*teta1*(C^2/u^2-(F1^(-teta)*lnF1^2+F2^(-teta)*lnF2^2)/u)) )


  H <- rbind( 
              cbind( sum(d2l.sqv.st.sqv.st), sum(d2l.sqv.st.teta.st) ) ,
              cbind( sum(d2l.sqv.st.teta.st), sum(d2l.teta.st.teta.st) )
            ) 

  
  res <- -sum(l.par)
  G   <- c( -sum( dl.dsqv.st ) ,  
            -sum( dl.dteta.st )   )
  H   <- H

 list(value=res, gradient=G, hessian=H, l=res, dat1=dat1, dat2=dat2) 
     
}


