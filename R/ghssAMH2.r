ghssAMH2 <- function(params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, params.f, sp=NULL, qu.mag=NULL){

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
  F1 <- pnorm(-eta1);
  F1 <- ifelse(F1>0.0000000001,F1,0.0000000001)
  F2 <- pnorm(e2);
  F2 <- ifelse(F2>0.0000000001,F2,0.0000000001)
  f2 <- dnorm(e2)/sqv;
  f2 <- ifelse(f2>0.0000000001,f2,0.0000000001)
  u <- 1-teta*(1-F1)*(1-F2)
  z <- F1*(1-teta+teta*F1)*u^(-2);
  z <- ifelse(z<0.9999999999,z,0.9999999999)
  b <- (1-F1)*(2*z*(1-F2)*u^(-1)-F1*u^(-2))
  h <- 2*teta*(z/(z-1))*(1-F1)*f2*u^(-1)
  ph <- dnorm(-eta1)
  E <- teta*u^(-1)*(1-F1)*f2*(z-3)/(z-1)+e2/sqv
  B <- i1*2*(1-teta^2)*(1-F1)*f2*u^(-1)*(z-1)^(-1)*(z*u^(-1)-teta*b*(z-1)^(-1))

  settings <- cbind(e2,F1,F2,f2,u,z,b,B,E,h)
  if(any(abs(settings)=='Inf') | any(settings=='NaN') | any(settings=='NA')) stop("Ill-conditioned task.")


# Likelihood:

  l.par <- weights*(i0*log(F1) + i1*(log(1-z)+log(f2))) 
  

# Gradient:

  dl.dsqv.st  <- weights*i1*(h*e2*sqv+e2^2-1)
  dl.dteta.st <- weights*i1*(1-teta^2)*b*(z-1 )^(-1)


# Hessian:

  d2l.sqv.st.sqv.st   <- -weights*i1*(sqv^2)*e2*(h*(E*e2-sqv^(-1))-2*e2*(sqv^(-2)))
  d2l.sqv.st.teta.st  <- -weights*e2*B*sqv

  d2l.teta.st.teta.st <- -weights*( i1*(1-teta^2)*(z-1)^(-1)*(
		(1-teta^2)*(
			-b^2*(z-1)^(-1)+2*(1-F1)*(1-F2)*u^(-1)*(b+(1-F1)*(z*u^(-1)*(1-F2)-F1*u^(-2)))
			) - 2*teta*b 
		) )


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


