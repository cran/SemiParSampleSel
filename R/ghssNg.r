ghssNg <- function(params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, sp=NULL, l.sp1, l.sp2, gp1, gp2, fp, qu.mag=NULL){

  eta1 <- dat1%*%params[1:X1.d2]  
  eta2 <- dat2%*%params[(X1.d2+1):(X1.d2+X2.d2)]   
  k.st <- params[(X1.d2+X2.d2+1)]
  teta.st <- params[(X1.d2+X2.d2+2)]
  k <- exp(k.st)
  teta <- tanh(teta.st)

  
# Settings:

  i0 <- 1-dat[, 1]
  i1 <- dat[, 1]
  i2 <- dat[, 2]

  F1 <- pnorm(-eta1);	
  F1 <- ifelse(F1>0.000000001,F1,0.000000001)
  # F1 <- ifelse(F1<0.99999999,F1,0.99999999)

  expeta2 <- exp(-eta2)

  F2 <- pgamma(i2,shape=k,rate=k*expeta2)		
  F2 <- ifelse(F2>0.000000001,F2,0.000000001)
  F2 <- ifelse(F2<0.99999999,F2,0.99999999)

  f2 <- dgamma(i2,shape=k,rate=k*expeta2)		
  f2[i1==0] <- 0
  f2 <- ifelse(f2>0.000000001,f2,0.000000001)

  q <- qnorm(F2)
  tetat <- 1/sqrt(1-teta^2)
  d  <- (eta1+teta*q)*tetat
  b  <- (teta*eta1+q)*tetat
  l1 <- ifelse(eta1<37.5,dnorm(-eta1)/pnorm(-eta1),eta1)
  l2 <- ifelse(d>-37.5,dnorm(d)/pnorm(d),-d)
  h <- teta*tetat*l2*(dnorm(q))^(-1)  
  E <- (teta*tetat*d-q)*(dnorm(q))^(-1)+h
  ye <- i2*expeta2

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
  df2kbyf2 <- -psi+log(k*i2)+1-eta2-ye
  df2kbyf2[i1==0]=0;  
  d2F2k <- i2*k^(-1)*f2*(2*df2kbyf2+ye-1-1/k) + (1-F2)*(psiprim-psi^2) + 2*psi*Gk*G^(-1)-G2k/G


  settings <- cbind(F1,F2,f2,q,d,b,l1,l2,E,h,ye,Gk,G2k,dF2k,df2kbyf2,d2F2k)
  if(any(abs(settings)=='Inf') | any(settings=='NaN') | any(settings=='NA')) 
  { for(i in 1:length(i1)) if(any(abs(settings[i,])=='Inf') | any(settings[i,]=='NaN') | any(settings[i,]=='NA')) print(c(i1[i],settings[i,]))
 stop("Ill-conditioned task.")
  }


# Likelihood:

  l.par <- weights*( i0*log(F1) + i1*log(pnorm(d)*f2) )


# Gradient:

  dl.dbe1     <- weights*( -i0*l1 + i1*l2*tetat )
  dl.dbe2     <- -weights*i1*( i2*f2*h + k*(1-ye) )
  dl.dk.st    <- weights*i1*( dF2k*h + df2kbyf2 )*k
  dl.dteta.st <- weights*i1*l2*b


# (Minus) Hessian:

  d2l.be1.be1     <- weights*( -i0*l1*(eta1-l1) + i1*tetat^2*l2*(d+l2) )
  d2l.be1.be2     <- -weights*i1*i2*f2*tetat*h*(d+l2)
  d2l.be1.k.st    <- weights*i1*dF2k*tetat*h*(d+l2)*k
  d2l.be1.teta.st <- -weights*i1*tetat*l2*( -b*(d+l2) + teta )

  d2l.be2.be2     <- weights*i1*i2*( f2*h*(i2*f2*E+k*ye-k) + k*exp(-eta2) )
  d2l.be2.k.st    <- weights*i1*( -i2*f2*h*(dF2k*E-df2kbyf2) + 1 - ye )*k
  d2l.be2.teta.st <- weights*i1*i2*f2*h*( -b*(d+l2) + teta^(-1) )

  d2l.k.st.k.st    <- -weights*i1*( dF2k*h + df2kbyf2 + h*k*(d2F2k-dF2k^2*E) + 1 - k*psiprim )*k
  d2l.k.st.teta.st <- -weights*i1*dF2k*h*( -b*(d+l2) + teta^(-1) )*k

  d2l.teta.st.teta.st <- -weights*i1*l2*( b*(-b*(d+l2)+teta) + eta1/tetat )


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
            -colSums( c(dl.dbe2)*dat2 ) ,
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

