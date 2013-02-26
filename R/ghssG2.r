ghssG2 <- function(params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, params.f, sp=NULL, qu.mag=NULL){

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
  F1 <- pnorm(-eta1);
  F1 <- ifelse(F1>0.000000001,F1,0.000000001)
  F1 <- ifelse(F1<0.99999999,F1,0.99999999)
  F2 <- pnorm(e2);
  F2 <- ifelse(F2>0.000000001,F2,0.000000001)
  f2 <- dnorm(e2)/sqv;
  f2 <- ifelse(f2>0.000000001,f2,0.000000001)
  lnF1 <- -log(F1)
  lnF2 <- -log(F2)
  u <- lnF1^teta + lnF2^teta
  u <- ifelse(u>0.000000001,u,0.000000001)
  ut <- u^(1/teta)
  z <- exp(-ut)*(F2^(-1))*ut*(u^(-1))*lnF2^(teta-1);	
  z <- ifelse(z<0.99999999,z,0.99999999)
  zt <- z/(z-1)
  p <- lnF1^(teta-1)*F1^(-1)*u^(-1)*(1-teta-ut)
  C <- (lnF1^teta)*log(lnF1)+(lnF2^teta)*log(lnF2)
  b <- teta^(-1)*(C/u-log(u)/teta)*(1-ut)+log(lnF2)-C/u
  Ct <- (lnF1^teta)*(log(lnF1))^2+(lnF2^teta)*(log(lnF2))^2
  h <- zt*F2^(-1)*f2
  d <- lnF2^(-1)*(-F1*lnF1*p-ut)+1
  E <- d*( F2^(-1)*f2*(1-d/(z-1)) + e2/sqv ) + 
	 (teta-1)*(lnF2)^(-2)*F2^(-1)*f2*( lnF2^teta*u^(-1)*(1-teta-ut+u^(-1)*lnF2^teta*(teta+ut))-1 )
  B <- h*( lnF2^(teta-1)*u^(-1)*
	(log(lnF2)*(1-teta-ut)+(teta-1)*teta^(-1)*C*u^(-1)*(teta+ut)-1+ut*log(u)*teta^(-2)) - b*d*(z-1)^(-1) + lnF2^(-1) )


  settings <- cbind(e2,F1,F2,f2,lnF1,lnF2,u,ut,z,zt,p,C,Ct,d,b,B,E,h)
  if(any(abs(settings)=='Inf') | any(settings=='NaN') | any(settings=='NA')) stop("Ill-conditioned task.")
 #{ for(i in 1:length(i1)) if(any(abs(settings[i,])=='Inf') | any(settings[i,]=='NaN') | any(settings[i,]=='NA')) { print(teta); print(c(i1[i],settings[i,]))}
 #stop("Ill-conditioned task.")
 # }


# Likelihood:

  l.par <- weights*( i0*log(F1) + i1*(log(1-z)+log(f2)) ) 
  

# Gradient:

  dl.dsqv.st <- weights*i1*(h*d*e2*sqv+e2^2-1)
  dl.dteta.st <- weights*i1*zt*b*(teta-1)


# Hessian:

  d2l.sqv.st.sqv.st <- -weights*sqv*e2*i1*( h*(E*sqv*e2-d) - 2*e2/sqv )
  d2l.sqv.st.teta.st <- -weights*i1*e2*B*sqv*(teta-1)

  d2l.teta.st.teta.st <- -weights*( i1*zt*(teta-1)*( b+(teta-1)*(
			teta^(-1)*(1-ut)*(Ct/u-C^2/u^2-2*teta^(-1)*(C/u-log(u)/teta)) - 
			ut*teta^(-2)*(C/u-log(u)/teta)^2 - b^2/(z-1) - (Ct/u-C^2/u^2)
			) ) )



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


