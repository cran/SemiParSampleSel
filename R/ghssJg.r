ghssJg <- function(params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, sp=NULL, l.sp1, l.sp2, gp1, gp2, fp, qu.mag=NULL){


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
  F1 <- ifelse(F1<0.999999999,F1,0.999999999)

  F2 <- pgamma(i2,shape=k,rate=k*exp(-eta2));		
  #F2 <- ifelse(F2>0.0000000001,F2,0.0000000001)
  F2 <- ifelse(F2<0.999999999,F2,0.999999999)

  f2 <- dgamma(i2,shape=k,rate=k*exp(-eta2));		
  f2[i2==0] <- 0
  f2 <- ifelse(f2>0.0000000001,f2,0.0000000001)

  F1b <- 1-F1;
  F2b <- 1-F2;
  lnF1b <- log(F1b)
  lnF2b <- log(F2b)

  teta1 <- teta-1
  u <- F1b^teta+F2b^teta-(F1b*F2b)^teta
  ut <- u^(1/teta-1)
  z <- (1-F1b^teta)*(F2b^teta1)*ut;  		
  z <- ifelse(z<0.9999999999, z,0.9999999999)
  z <- ifelse(z>0.00000000001,z,0.00000000001)
  zt <- (z-1)^(-1)
  p <- F2b^teta1*ut*(u+teta1)*teta^(-1)*u^(-1)
  C <- (F1b^teta)*lnF1b+(F2b^teta)*lnF2b-((F1b*F2b)^teta)*(lnF1b+lnF2b)
  b <- lnF2b-(teta1/teta)*(C/u)-log(u)/(teta^2)
  B <- z*b-(F2b^teta1)*ut*(F1b^teta)*lnF1b
  Ct <- (F1b^teta)*(lnF1b^2)+(F2b^teta)*(lnF2b^2)-((F1b*F2b)^teta)*((lnF1b+lnF2b)^2)
  h <- z*zt*teta1*F2b^(-1)*F1b^teta*u^(-1)
  ph <- dnorm(-eta1)
  E <- F2b^(-1)*( 1 + teta*(1-u^(-1)*F1b^teta) + h*F2b*z^(-1) ) 
  Bt <- i1*h*(1-teta1*(C/u-lnF1b+B*z^(-1)*zt))
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


  settings <- cbind(i1,i2,F1b,F2b,f2,lnF1b,lnF2b,u,ut,z,zt,p,C,Ct,ph,h,b,B,Bt,E,ye,Gk,G2k,dF2k,df2kbyf2,d2F2k)
  if(any(abs(settings)=='Inf') | any(settings=='NaN') | any(settings=='NA')) stop("Ill-conditioned task.") 


# Likelihood:

  l.par <- weights*( i0*log(F1) + i1*(log(1-z)+log(f2)) )


# Gradient:

  dl.dbe1     <- weights*( -i0*F1^(-1) - i1*teta*zt*F1b^teta1*p )*ph
  dl.dbe2     <- weights*i1*( i2*f2*h - k*(1-ye) )
  dl.dk.st    <- weights*i1*( -dF2k*h + df2kbyf2 )*k
  dl.dteta.st <- weights*i1*zt*B*teta1


# (Minus) Hessian:

  d2l.be1.be1     <- -weights*ph*( -i0*(ph/F1-eta1)/F1 -
		i1*teta*zt*(F1b^(teta-2))*(
			ph*( p*(teta1+teta*p*F1b^teta*zt) - teta1*u^(-1)*(u-F2b^teta)*(2*p+z*(teta*u)^(-1)*(1-F2b^teta)) )
			-F1b*p*eta1
		))
  d2l.be1.be2     <- -weights*teta*i1*i2*f2*h*F1b^teta1*( p*zt*z^(-1) + F2b^teta*F1b^(-teta)*u^(-1) )*ph      
  d2l.be1.k.st    <- weights*teta*i1*h*F1b^teta1*dF2k*( p*zt*z^(-1) + F2b^teta*F1b^(-teta)*u^(-1) )*ph*k     
  d2l.be1.teta.st <- weights*i1*teta1*zt*F1b^teta1*p*(
		1 + teta*(lnF1b-B*zt+b+((1-u)/teta-teta1*C*u^(-1))/(u+teta1))
	)*ph

  d2l.be2.be2     <- weights*i1*i2*( f2*h*(i2*f2*E-k*ye+k) + k*exp(-eta2) )
  d2l.be2.k.st    <- -weights*i1*( i2*f2*h*(dF2k*E+df2kbyf2) - 1 + ye )*k
  d2l.be2.teta.st <- -weights*i2*f2*Bt

  d2l.k.st.k.st    <- -weights*i1*( -dF2k*h + df2kbyf2 - h*k*(d2F2k+dF2k^2*E) + 1 - k*psiprim )*k
  d2l.k.st.teta.st <- weights*dF2k*Bt*k

  d2l.teta.st.teta.st <- -weights*i1*zt*teta1*( B - teta1*(
			zt*B^2-B*b+z*(teta1*teta^(-1)*(Ct/u-C^2/u^2) + 2*teta^(-2)*(C/u-log(u)/teta)) + (z*b-B)*(lnF1b+b)
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


