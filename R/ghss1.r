ghss1 <- function(BivD, margins, params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, sp=NULL, l.sp1, l.sp2, gp1, gp2, exl1, exl2, qu.mag=NULL){


  eta1 <- dat1%*%params[1:X1.d2]  
  eta2 <- dat2%*%params[(X1.d2+1):(X1.d2+X2.d2)]   
  teta.st <- exl2
  eps <- .Machine$double.eps*10^6

  if(BivD %in% c("N","FGM","AMH") ){ teta <- tanh(teta.st); if(teta %in% c(-1,1)) teta <- sign(teta)*0.9999999 }
  if(BivD=="F")                teta <- teta.st + eps
  if(BivD %in% c("C", "rC") )  teta <- exp(teta.st) + eps
  if(BivD %in% c("J", "rJ") )  teta <- exp(teta.st) + 1 + eps 
  if(BivD %in% c("G", "rG") )  teta <- exp(teta.st) + 1 + eps

  i0 <- 1-dat[,1]
  i1 <- dat[,1]
  ind=i1==0

  F1 <- pnorm(-eta1)				
  F1 <- ifelse(F1>0.00000001,F1,0.00000001)
  F1 <- ifelse(F1<0.99999999,F1,0.99999999)
  ph <- dnorm(-eta1)



if(margins[2]=="N"){


    sqv.st <- exl1
    sqv <- exp(sqv.st)


    if(BivD=="N"){ 
 
      a  <- sqrt(1-teta^2)
      e2 <- dat[,2]-eta2  
      A  <- (eta1+(teta/sqv)*e2)/a
      l1 <- ifelse(eta1<37.5,dnorm(-eta1)/pnorm(-eta1),eta1)
      l2 <- ifelse(A>-37.5,dnorm(A)/pnorm(A),-A)
      ec <- exp(2*teta.st) 

      PA    <- -(pnorm(A)*dnorm(A)*A+dnorm(A)^2)/pnorm(A)^2 
      Peta1 <- -(pnorm(-eta1)*dnorm(-eta1)*(-eta1)+dnorm(-eta1)^2)/pnorm(-eta1)^2  


      l.par <- weights*(i0*log(pnorm(-eta1)) + i1*(log(pnorm(A))-1/2*log(2*pi)-log(sqv)-1/2*(e2/sqv)^2) )
 
      dl.1 <- -weights*(i0*l1 - i1*(l2/a)) 
      dl.2 <- weights*i1*( e2/sqv^2 - l2*( teta/(sqv*a) )   )
  
      d2l.11 <- -weights*( i0*Peta1 + i1*(PA/a^2) )  
      d2l.12 <- -weights*( -i1*PA*(teta/(sqv*a^2)) )   
      d2l.22 <- -weights*( i1*( PA*(teta/(sqv*a))^2 - 1/sqv^2 ) ) 

    }
    else{


        e2 <- (dat[,2]-eta2)/sqv
        F2 <- pnorm(e2)
        F2 <- ifelse(F2>0.00000001,F2,0.00000001)
        F2 <- ifelse(F2<0.99999999,F2,0.99999999)
        f2 <- dnorm(e2)/sqv;
        f2 <- ifelse(f2>0.00000001,f2,0.00000001)


        bits <- bitsgHs(cop=BivD,margin=margins[2],i1=i1,F1=F1,F2=F2,f2=f2,eta1=eta1,ph=ph,teta=teta,e2=e2,sqv=sqv,ver=1)

        z <- bits$z; p <- bits$p; h <- bits$h; P <- bits$P; A <- bits$A; E <- bits$E 
        z[ind]=0; p[ind]=0; h[ind]=0; P[ind]=0; A[ind]=0; E[ind]=0;


        l.par <- weights*(i0*log(F1) + i1*(log(1-z)+log(f2)))

        # Gradient:

        dl.1  <- weights*( -i0*F1^(-1) + i1*p )*ph                  # dl.dbe1
        dl.2  <- weights*i1*( h + e2/sqv )                          # dl.dbe2

        # (Minus) Hessian:

        d2l.11 <- -weights*ph*( -i0*(ph/F1-eta1)/F1 + i1*P )        # d2l.be1.be1
        d2l.12 <- -weights*A                                        # d2l.be1.be2

        d2l.22 <- -weights*i1*( h*E - sqv^(-2) )                    # d2l.be2.be2

    }

}
else if(margins[2]=="G"){


    k.st <- exl1
    k <- exp(k.st)
    i2 <- dat[,2]
    F2 <- pgamma(i2,shape=k,rate=k*exp(-eta2));		
    F2 <- ifelse(F2>0.00000001,F2,0.00000001)
    F2 <- ifelse(F2<0.99999999,F2,0.99999999)
    f2 <- dgamma(i2,shape=k,rate=k*exp(-eta2));		
    f2[i2==0] <- 0
    f2 <- ifelse(f2>0.00000001,f2,0.00000001)

    ye <- i2*exp(-eta2)


    bits <- bitsgHs(cop=BivD,margin=margins[2],i1=i1,F1=F1,F2=F2,f2=f2,eta1=eta1,ph=ph,teta=teta,ver=1)

    z <- bits$z; p <- bits$p; h <- bits$h; P <- bits$P; A <- bits$A; E <- bits$E
    z[ind]=0; p[ind]=0; h[ind]=0; P[ind]=0; A[ind]=0; E[ind]=0

    l.par <- weights*(i0*log(F1) + i1*(log(1-z)+log(f2)))


    # Gradient:

    dl.1  <- weights*( -i0*F1^(-1) + i1*p )*ph                   # dl.dbe1
    dl.2  <- weights*i1*( i2*f2*h - k*(1-ye) )                   # dl.dbe2

 
    # (Minus) Hessian:

    d2l.11  <- -weights*ph*( -i0*(ph/F1-eta1)/F1 + i1*P )        # d2l.be1.be1
    d2l.12  <- -weights*i2*f2*A                                  # d2l.be1.be2

    d2l.22  <- -weights*i1*i2*( f2*h*(i2*f2*E+k*ye-k) - k*exp(-eta2) )    # d2l.be2.be2

} 


  H11 <- crossprod(dat1*c(d2l.11),dat1)
  H12 <- crossprod(dat1*c(d2l.12),dat2) 
  H22 <- crossprod(dat2*c(d2l.22),dat2) 


  H <- rbind( cbind( H11    , H12 ), 
              cbind( t(H12) , H22 )   ) 

  res <- -sum(l.par)

  G   <- c( -colSums( c(dl.1)*dat1 ) ,
            -colSums( c(dl.2)*dat2 )   )



      S <- mapply("*", qu.mag$Ss, sp, SIMPLIFY=FALSE)
      S <- do.call(adiag, lapply(S, unlist))

      if(l.sp1!=0 && l.sp2!=0) S.h <- adiag(matrix(0,gp1,gp1),
                                            S[1:(X1.d2-gp1),1:(X1.d2-gp1)],
                                            matrix(0,gp2,gp2),
                                            S[(X1.d2-(gp1-1)):dim(S)[2],(X1.d2-(gp1-1)):dim(S)[2]])

      if(l.sp1==0 && l.sp2!=0) S.h <- adiag(matrix(0,gp1,gp1), matrix(0,gp2,gp2), S)
      if(l.sp1!=0 && l.sp2==0) S.h <- adiag(matrix(0,gp1,gp1), S, matrix(0,gp2,gp2))

  S.h1 <- 0.5*crossprod(params,S.h)%*%params
  S.h2 <- S.h%*%params

         
  S.res <- res
  res <- S.res + S.h1
  G   <- G + S.h2
  H   <- H + S.h  
 

  list(value=res, gradient=G, hessian=H, S.h=S.h, l=S.res, eta1=eta1, eta2=eta2, dat1=dat1, dat2=dat2, 
       dl.dbe1=dl.1, dl.dbe2=dl.2, 
       d2l.be1.be1=d2l.11, d2l.be1.be2=d2l.12, d2l.be2.be2=d2l.22)    

}


