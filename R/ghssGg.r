ghssGg <- function(params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, sp=NULL, l.sp1, l.sp2, gp1, gp2, fp, qu.mag=NULL){

  eta1 <- dat1%*%params[1:X1.d2]  
  eta2 <- dat2%*%params[(X1.d2+1):(X1.d2+X2.d2)]   
  k.st <- params[(X1.d2+X2.d2+1)]
  teta.st <- params[(X1.d2+X2.d2+2)]
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
  z <- exp(-ut)*F2^(-1)*ut*u^(-1)*lnF2^(teta-1);
  z <- ifelse(z<0.9999999999,z,0.9999999999)
  zt <- z/(z-1)
  p <- lnF1^(teta-1)*F1^(-1)*u^(-1)*(1-teta-ut)
  C <- lnF1^teta*log(lnF1)+lnF2^teta*log(lnF2)
  b <- teta^(-1)*(C/u-log(u)/teta)*(1-ut)+log(lnF2)-C/u
  Ct <- (lnF1^teta)*(log(lnF1))^2+(lnF2^teta)*(log(lnF2))^2
  h <- zt*F2^(-1)
  ph <- dnorm(-eta1)
  d <- lnF2^(-1)*(-F1*lnF1*p-ut)+1
  E <- (teta-1)*F2^(-1)*lnF2^(-2)*( lnF2^teta*u^(-1)*(1-teta-ut+u^(-1)*lnF2^teta*(teta+ut)) - 1 )
  B <- h*f2*( lnF2^(teta-1)*u^(-1)*
	( log(lnF2)*(1-teta-ut) + (teta-1)*teta^(-1)*C*u^(-1)*(teta+ut) - 1 + ut*log(u)*teta^(-2) ) - b*d*(z-1)^(-1) + lnF2^(-1) )
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


  settings <- cbind(F1,F2,f2,lnF1,lnF2,u,ut,z,zt,p,C,Ct,ph,d,b,B,E,h,ye,Gk,G2k,dF2k,df2kbyf2,d2F2k)
  if(any(abs(settings)=='Inf') | any(settings=='NaN') | any(settings=='NA')) stop("Ill-conditioned task.")


# Likelihood:

  l.par <- weights*(i0*log(F1) + i1*(log(1-z)+log(f2))) 


# Gradient:

  dl.dbe1     <- weights*( -i0*F1^(-1) + i1*zt*p )*ph
  dl.dbe2     <- weights*i1*( i2*f2*d*h - k*(1-ye) )
  dl.dk.st    <- weights*i1*( -dF2k*d*h + df2kbyf2 )*k
  dl.dteta.st <- weights*i1*zt*b*(teta-1)


# (Minus) Hessian:

  d2l.be1.be1     <- -weights*ph*( -i0*(ph/F1-eta1)/F1 +
		i1*zt*p*(
			F1^(-1)*ph*( 1+lnF1^(-1)*(teta-1-lnF1^teta*u^(-1)*(teta+ut*(1-teta-ut)^(-1))) - p*F1*(z-1)^(-1) )
			-eta1
		))
  d2l.be1.be2     <- -weights*i1*i2*f2*h*p*( (teta-1)*(teta+ut)*u^(-1)*lnF2^(teta-1)*(1-teta-ut)^(-1) - d/(z-1) )*ph      
  d2l.be1.k.st    <- weights*i1*dF2k*h*p*( (teta-1)*(teta+ut)*u^(-1)*lnF2^(teta-1)*(1-teta-ut)^(-1) - d/(z-1) )*ph*k     
  d2l.be1.teta.st <- -weights*i1*(teta-1)*zt*p*(
	log(lnF1)-(1-teta-ut)^(-1)*((1-teta)*teta^(-1)*C*u^(-1)*(teta+ut)+1-ut*log(u)*teta^(-2)) - b/(z-1)
	)*ph

  d2l.be2.be2     <- -weights*i1*i2*( f2*h*(i2*f2*E + d*(i2*f2*F2^(-1)*(1-d/(z-1))-k*(1-ye))) - k*exp(-eta2) )
  d2l.be2.k.st    <- -weights*i1*( i2*f2*h*( d*(df2kbyf2-dF2k*F2^(-1)*(1-d/(z-1))) - dF2k*E ) - 1 + ye )*k
  d2l.be2.teta.st <- -weights*i1*i2*B*(teta-1)

  d2l.k.st.k.st    <- -weights*i1*( -dF2k*h*d + df2kbyf2 + 1 - k*psiprim - 
                                      h*k*( d*d2F2k + dF2k^2*(-E-d*F2^(-1)*(1-d/(z-1))) ) )*k
  d2l.k.st.teta.st <- weights*i1*dF2k*B*f2^(-1)*(teta-1)*k

  d2l.teta.st.teta.st <- -weights*i1*zt*(teta-1)*( b + (teta-1)*(
			teta^(-1)*(1-ut)*(Ct/u-C^2/u^2-2*teta^(-1)*(C/u-log(u)/teta)) - 
			ut*teta^(-2)*(C/u-log(u)/teta)^2 - b^2/(z-1) - Ct/u + C^2/u^2   
			) )


  be1.be1     <- crossprod(dat1*c(d2l.be1.be1),dat1)
  be1.be2     <- crossprod(dat1*c(d2l.be1.be2),dat2) 
  be1.k.st    <- t(t(rowSums(t(dat1*c(d2l.be1.k.st)))))
  be1.teta.st <- t(t(rowSums(t(dat1*c(d2l.be1.teta.st)))))

  be2.be2     <- crossprod(dat2*c(d2l.be2.be2),dat2) 
  be2.k.st    <- t(t(rowSums(t(dat2*c(d2l.be2.k.st)))))
  be2.teta.st <- t(t(rowSums(t(dat2*c(d2l.be2.teta.st)))))


  H <- rbind( cbind( be1.be1    , be1.be2    , be1.k.st  ,  be1.teta.st ), 
              cbind( t(be1.be2) , be2.be2    , be2.k.st  ,  be2.teta.st ),
              cbind( t(be1.k.st) , t(be2.k.st) , sum(d2l.k.st.k.st), sum(d2l.k.st.teta.st) ) ,
              cbind( t(be1.teta.st) , t(be2.teta.st) , sum(d2l.k.st.teta.st), sum(d2l.teta.st.teta.st) )
            ) 



  res <- -sum(l.par)
  G   <- c( -colSums( c(dl.dbe1)*dat1 ) ,
            -colSums( c(dl.dbe2)*dat2 )    ,
            -sum( dl.dk.st ) ,  
            -sum( dl.dteta.st )   )


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
       dl.dsqv.st=dl.dk.st,
       dl.dcor.st=dl.dteta.st, 
       d2l.be1.be1=d2l.be1.be1, d2l.be1.be2=d2l.be1.be2, d2l.be2.be2=d2l.be2.be2,
       d2l.be1.sqv.st=d2l.be1.k.st,
       d2l.be1.cor.st=d2l.be1.teta.st,
       d2l.be2.sqv.st=d2l.be2.k.st, 
       d2l.be2.cor.st=d2l.be2.teta.st,
       d2l.sqv.st.sqv.st=d2l.k.st.k.st,
       d2l.sqv.st.cor.st=d2l.k.st.teta.st,    
       d2l.cor.st.cor.st=d2l.teta.st.teta.st )

}


