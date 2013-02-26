ghssJ2 <- function(params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, params.f, sp=NULL, qu.mag=NULL){

  eta1 <- dat1%*%params.f[1:X1.d2]  
  eta2 <- dat2%*%params.f[(X1.d2+1):(X1.d2+X2.d2)]   
  sqv.st <- params[1]
  teta.st <- params[2]
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
  d <- teta1*F1b^teta*u^(-1)
  C <- (F1b^teta)*lnF1b+(F2b^teta)*lnF2b-((F1b*F2b)^teta)*(lnF1b+lnF2b)
  b <- lnF2b+(1/teta-1)*(C/u)-log(u)/(teta^2)
  B <- z*b-(F2b^teta1)*ut*(F1b^teta)*lnF1b
  Ct <- (F1b^teta)*(lnF1b^2)+(F2b^teta)*(lnF2b^2)-((F1b*F2b)^teta)*((lnF1b+lnF2b)^2)
  h <- z*zt*f2*F2b^(-1)*d
  E <- f2*F2b^(-1)*( teta*u^(-1)*F1b^teta-teta-1-d*zt ) + e2/sqv
  Bt <- i1*h*(1-teta1*(C/u-lnF1b+B*z^(-1)*(z-1)^(-1)))


  settings <- cbind(e2,F1b,F2b,f2,lnF1b,lnF2b,u,ut,z,zt,C,d,b,B,Bt,E,h)
  if(any(abs(settings)=='Inf') | any(settings=='NaN') | any(settings=='NA')) stop("Ill-conditioned task.")


# Likelihood:

  l.par <- weights*(i0*log(F1) + i1*(log(1-z)+log(f2))) 
  

# Gradient:

  dl.dsqv.st <- weights*i1*(h*e2*sqv+e2^2-1)
  dl.dteta.st <- weights*i1*zt*B*teta1


# Hessian:

 d2l.sqv.st.sqv.st <- -weights*sqv*e2*i1*( h*(E*sqv*e2-1)-2*e2/sqv )
  d2l.sqv.st.teta.st <- -weights*e2*Bt*sqv

  d2l.teta.st.teta.st <- -weights*( i1*zt*teta1*( B-teta1*(
			zt*B^2-B*b+z*(teta1*teta^(-1)*(Ct/u-C^2/u^2) + 2*teta^(-2)*(C/u-log(u)/teta)) + (z*b-B)*(lnF1b+b)
			) ) )



  H <- rbind( 
              cbind( sum(d2l.sqv.st.sqv.st), sum(d2l.sqv.st.teta.st) ) ,
              cbind( sum(d2l.sqv.st.teta.st), sum(d2l.teta.st.teta.st) )
            ) 

  
  res <- -sum(l.par)
  G   <- c( -sum( dl.dsqv.st ) ,  
            -sum( dl.dteta.st )   )


  list(value=res, gradient=G, hessian=H, l=res, dat1=dat1, dat2=dat2) 
     
}


