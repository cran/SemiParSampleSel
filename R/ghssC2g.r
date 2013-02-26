ghssC2g <- function(params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, params.f, sp=NULL, qu.mag=NULL){

  eta1 <- dat1%*%params.f[1:X1.d2]  
  eta2 <- dat2%*%params.f[(X1.d2+1):(X1.d2+X2.d2)]   
  k.st <- params[1]
  teta.st <- params[2]
  k <- exp(k.st)
  teta <- exp(teta.st)
  

# Settings:

  i0 <- 1-dat[, 1]
  i1 <- dat[, 1]
  i2 <- dat[, 2]

  F1 <- pnorm(-eta1);					
  F1 <- ifelse(F1>0.0000000001,F1,0.0000000001)
  F1 <- ifelse(F1<0.9999999999,F1,0.9999999999)

  F2 <- pgamma(i2,shape=k,rate=k*exp(-eta2));		
  F2 <- ifelse(F2>0.0000000001,F2,0.0000000001)
  F2 <- ifelse(F2<0.9999999999,F2,0.9999999999)

  f2 <- dgamma(i2,shape=k,rate=k*exp(-eta2));		
  f2[i2==0] <- 0
  f2 <- ifelse(f2>0.0000000001,f2,0.0000000001)

  lnF1 <- log(F1)
  lnF2 <- log(F2)
  teta1 <- teta+1
  u <- F1^(-teta)+F2^(-teta)-1
  z <- F2^(-teta1)*u^(-teta1/teta);
  z <- ifelse(z<0.9999999999,z,0.9999999999)
  z <- ifelse(z>0.0000000001,z,0.0000000001)
  zt <- z/(1-z)
  C <- F1^(-teta)*lnF1+F2^(-teta)*lnF2
  h <- -zt*teta1*F2^(-teta1)*(F2^teta-u^(-1))
  b <- (lnF2-log(u)/teta^2-C*teta1*teta^(-1)*u^(-1))/(z-1)
  B <- -zt*i1*F2^(-teta1)*( (1+teta1*b)*(F2^teta-u^(-1))+teta1*(lnF2/u-C*u^(-2)))
  E <- F2^(-1)-F2^(-teta1)*teta*u^(-1)-h/z
  ye <- i2*exp(-eta2)
  G <- gamma(k)
  psi <- psigamma(k)
  psiprim <- psigamma(k,deriv=1)
  fun1 <- function(t,k){t^(k-1)*exp(-t)*log(t)}
  Gkf <- function(k, y) { sapply(y, function(y) integrate(fun1,y,Inf,k)$value ) }
  Gk <- Gkf(k,k*ye)
  fun2 <- function(t,k){t^(k-1)*exp(-t)*(log(t))^2}
  G2kf <- function(k, y) { sapply(y, function(y) integrate(fun2,y,Inf,k)$value ) }
  G2k <- G2kf(k,k*ye)
  dF2k <- i2*f2*k^(-1)+psi*(1-F2)-Gk/G
  df2kbyf2 <- -psi+log(k*i2)+1-eta2-ye;
  df2kbyf2[i1==0]=0;
  d2F2k <-i2*k^(-1)*f2*(2*df2kbyf2+ye-1-1/k) + (1-F2)*(psiprim-psi^2) + 2*psi*Gk*G^(-1)-G2k/G 			


  settings <- cbind(F1,F2,f2,lnF1,lnF2,u,z,zt,C,h,b,B,E,ye,Gk,G2k,dF2k,df2kbyf2,d2F2k)
  if(any(abs(settings)=='Inf') | any(settings=='NaN') | any(settings=='NA')) stop("Ill-conditioned task.")


# Likelihood:

  l.par <- weights*(i0*lnF1 + i1*(log(1-z)+log(f2))) 
  

# Gradient:

  dl.dk.st    <- weights*i1*( -dF2k*h + df2kbyf2 )*k
  dl.dteta.st <- -weights*i1*teta*z*b


# (Minus) Hessian:

  d2l.k.st.k.st       <- -weights*i1*( -dF2k*h + df2kbyf2 - h*k*(d2F2k-dF2k^2*E) + 1 - k*psiprim)*k
  d2l.k.st.teta.st    <- weights*dF2k*B*teta*k
  d2l.teta.st.teta.st <- -weights*i1*zt*( teta^2*b^2*(z-1) + teta*lnF2 + log(u)/teta + (1-teta)*C*u^(-1) +
	teta*teta1*((F1^(-teta)*lnF1^2+F2^(-teta)*lnF2^2)*u^(-1)-C^2*u^(-2)) )


  H <- rbind( 
              cbind( sum(d2l.k.st.k.st), sum(d2l.k.st.teta.st) ) ,
              cbind( sum(d2l.k.st.teta.st), sum(d2l.teta.st.teta.st) )
            ) 

  
  res <- -sum(l.par)
  G   <- c( -sum( dl.dk.st ) ,  
            -sum( dl.dteta.st )   )
  H   <- H

  list(value=res, gradient=G, hessian=H, l=res, dat1=dat1, dat2=dat2) 
     
}


