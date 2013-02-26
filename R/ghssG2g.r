ghssG2g <- function(params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, params.f, sp=NULL, qu.mag=NULL){

  eta1 <- dat1%*%params.f[1:X1.d2]  
  eta2 <- dat2%*%params.f[(X1.d2+1):(X1.d2+X2.d2)]   
  k.st <- params[1]
  teta.st <- params[2]
  k <- exp(k.st)
  teta <- 1+exp(teta.st)

  
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

  lnF1 <- -log(F1)
  lnF2 <- -log(F2)
  u <- lnF1^teta + lnF2^teta
  ut <- u^(1/teta)
  z <- exp(-ut)*(F2^(-1))*ut*(u^(-1))*lnF2^(teta-1);
  z <- ifelse(z<0.9999999999,z,0.9999999999)
  zt <- z/(z-1)
  p <- lnF1^(teta-1)*F1^(-1)*u^(-1)*(1-teta-ut)
  C <- (lnF1^teta)*log(lnF1)+(lnF2^teta)*log(lnF2)
  b <- teta^(-1)*(C/u-log(u)/teta)*(1-ut)+log(lnF2)-C/u
  Ct <- (lnF1^teta)*(log(lnF1))^2+(lnF2^teta)*(log(lnF2))^2
  h <- zt*F2^(-1)
  ph <- dnorm(-eta1)
  d <- lnF2^(-1)*(-F1*lnF1*p-ut)+1
  E <- (teta-1)*F2^(-1)*lnF2^(-2)*(lnF2^teta*u^(-1)*(1-teta-ut+u^(-1)*lnF2^teta*(teta+ut))-1)
  B <- h*f2*( lnF2^(teta-1)*u^(-1)*
	(log(lnF2)*(1-teta-ut)+(teta-1)*teta^(-1)*C*u^(-1)*(teta+ut)-1+ut*log(u)*teta^(-2)) - b*d*(z-1)^(-1) + lnF2^(-1) )
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
  df2kbyf2 <- -psi+log(k*i2)+1-eta2-ye;  df2kbyf2[i1==0]=0;
  d2F2k <-i2*k^(-1)*f2*(2*df2kbyf2+ye-1-1/k) + (1-F2)*(psiprim-psi^2) + 2*psi*Gk*G^(-1)-G2k/G 


  settings <- cbind(F1,F2,f2,lnF1,lnF2,u,ut,z,zt,p,C,Ct,d,b,B,E,h,ye,Gk,G2k,dF2k,df2kbyf2,d2F2k)
  if(any(abs(settings)=='Inf') | any(settings=='NaN') | any(settings=='NA')) stop("Ill-conditioned task.")


# Likelihood:

  l.par <- weights*( i0*log(F1) + i1*(log(1-z)+log(f2)) ) 
  

# Gradient:

  dl.dk.st <- weights*i1*( -dF2k*d*h + df2kbyf2 )*k
  dl.dteta.st <- weights*i1*zt*b*(teta-1)


# (Minus) Hessian:

  d2l.k.st.k.st    <- -weights*i1*(-dF2k*h*d+df2kbyf2+1-k*psiprim-h*k*(d*d2F2k+dF2k^2*(-E-d*F2^(-1)*(1-d/(z-1)))))*k
  d2l.k.st.teta.st <- weights*i1*dF2k*B*f2^(-1)*(teta-1)*k

  d2l.teta.st.teta.st <- -weights*i1*zt*(teta-1)*( b + (teta-1)*(
			teta^(-1)*(1-ut)*(Ct/u-C^2/u^2-2*teta^(-1)*(C/u-log(u)/teta)) - 
			ut*teta^(-2)*(C/u-log(u)/teta)^2 - b^2/(z-1) - (Ct/u-C^2/u^2)
			) )

  H <- rbind( 
              cbind( sum(d2l.k.st.k.st), sum(d2l.k.st.teta.st) ) ,
              cbind( sum(d2l.k.st.teta.st), sum(d2l.teta.st.teta.st) )
            ) 

  
 res <- -sum(l.par)
 G   <- c( -sum( dl.dk.st ) ,  
           -sum( dl.dteta.st )   )

 list(value=res, gradient=G, hessian=H, l=res, dat1=dat1, dat2=dat2) 

    
}
