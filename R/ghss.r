ghss <- function(BivD, margins, params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, sp=NULL, l.sp1, l.sp2, gp1, gp2, fp, qu.mag=NULL){

  eta1 <- dat1%*%params[1:X1.d2]  
  eta2 <- dat2%*%params[(X1.d2+1):(X1.d2+X2.d2)]   
  teta.st <- params[(X1.d2+X2.d2+2)]
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


    sqv.st <- params[(X1.d2+X2.d2+1)]
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


      l.par <- weights*(i0*log(pnorm(-eta1)) + i1*(log(pnorm(A))-1/2*log(2*pi)-log(sqv)-1/2*(e2/sqv)^2)) 
 

      dl.1 <- -weights*( i0*l1 - i1*(l2/a) ) 
      dl.2 <- weights*i1*( e2/sqv^2 - l2*( teta/(sqv*a) )   )
      dl.3 <- weights*i1*( -1 + (e2/sqv)^2 - l2*((teta*e2)/(a*sqv)) ) 
      dl.4 <- weights*i1*( l2*( ((2*ec*e2)/((ec+1)*sqv) - (2*teta*e2*ec)/((ec+1)*sqv))/a  - 0.5*( (eta1+(teta*e2)/sqv)*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1) ) )/a^3 )  )   


      d2l.11 <- -weights*( i0*Peta1 + i1*(PA/a^2) )  
      d2l.12 <- -weights*( -i1*PA*(teta/(sqv*a^2)) )   
      d2l.13 <- -weights*( -i1*PA*((teta*e2)/(sqv*a^2)) ) 
      d2l.14 <- -weights*( i1*(       (PA/a)*(  ( (2*ec*e2)/((ec+1)*sqv) - (2*(ec-1)*e2*ec)/((ec+1)^2*sqv) )/a      
               				- 0.5*(eta1+(teta*e2)/sqv)*((-4*teta*ec)/(ec+1)+(4*teta^2*ec)/(ec+1) )/a^3 
         				  )  
   				 -0.5*(l2/a^3)*( (-4*teta*ec)/(ec+1) + (4*teta^2*ec)/(ec+1) )           
                     )    )
      d2l.22 <- -weights*( i1*( PA*(teta/(sqv*a))^2 - 1/sqv^2 ) ) 
      d2l.23 <- -weights*( i1*( PA*(teta^2*e2)/(sqv^2*a^2) + (teta*l2)/(sqv*a) - (2*e2)/sqv^2 ) )
      d2l.24 <- -weights*( i1*(     -PA*(    
    				( (  (2*e2*ec)/((ec+1)*sqv) - (2*teta*e2*ec)/((ec+1)*sqv)  )/a 
      				   -0.5/a^3*( ( eta1+(teta*e2)/sqv )*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1) ))   )*(ec-1)  
     				)/((ec+1)*sqv*a)
   				- (2*l2*ec)/(sqv*a*(ec+1))      
  				 + (2*l2*teta*ec)/(sqv*a*(ec+1))  
  				 + 0.5*l2*teta/(sqv*a^3)*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1)   ) 
                     )    )
      d2l.33 <- -weights*( i1*(  -2*(e2/sqv)^2 + l2*((e2*teta)/(a*sqv)) + PA*((teta*e2)/(sqv*a))^2  )   )
      d2l.34 <- -weights*( i1*(   -PA*(    (( (2*e2*ec)/((ec+1)*sqv) -( 2*teta*e2*ec)/((ec+1)*sqv) )/a 
         				 -0.5/a^3*( ( eta1+(teta*e2)/sqv )*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1) )))*(ec-1)*e2  
      				)/((ec+1)*sqv*a)
   				- (2*l2*ec*e2)/(sqv*a*(ec+1))      
  				 + (2*l2*teta*e2*ec)/(sqv*a*(ec+1))  
  				 + 0.5*l2*teta*e2/(sqv*a^3)*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1)   ) 
                     )    )
      d2l.44 <- -weights*( i1*(   PA*(       (
        				 ( (2*e2*ec)/((ec+1)*sqv) -(2*teta*e2*ec)/((ec+1)*sqv) )/a 
       				 -0.5/a^3*( ( eta1+(teta*e2)/sqv )*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1) )) 
          				   )^2
     				 )
 				+ l2*(
 		 1/a*( 4*ec*e2/((ec+1)*sqv) - 8*ec^2*e2/((ec+1)^2*sqv) + 8*teta*ec^2*e2/((ec+1)^2*sqv)  - 4*teta*ec*e2/((ec+1)*sqv)  ) 
  		-1/a^3*( (2*ec*e2/((ec+1)*sqv) - 2*ec*e2*teta/((ec+1)*sqv) )*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1) )  ) 
  		+3/(4*a^5)*( ( eta1+(teta*e2)/sqv )*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1) )^2 ) 
  		-0.5/a^3*( (eta1+(teta*e2)/sqv )*(  -8*ec^2/(ec+1)^2 + 32*teta*ec^2/(ec+1)^2 - 8*teta*ec/(ec+1) - 24*teta^2*ec^2/(ec+1)^2  + 8*teta^2*ec/(ec+1)    ) 
       		    )
     		 )
                 )    )

    }
    else{


        e2 <- (dat[,2]-eta2)/sqv
        F2 <- pnorm(e2)
        F2 <- ifelse(F2>0.000000001,F2,0.000000001)
        F2 <- ifelse(F2<0.999999999,F2,0.999999999)
        f2 <- dnorm(e2)/sqv;
        f2 <- ifelse(f2>0.000000001,f2,0.000000001)


        bits <- bitsgHs(cop=BivD,margin=margins[2],i1=i1,F1=F1,F2=F2,f2=f2,eta1=eta1,ph=ph,teta=teta,e2=e2,sqv=sqv,ver=0)

        z <- bits$z; p <- bits$p; h <- bits$h; b <- bits$b; P <- bits$P; A <- bits$A; h14 <- bits$h14; E <- bits$E; B <- bits$B; h44 <- bits$h44 

        z[ind]=0; p[ind]=0; h[ind]=0; b[ind]=0; P[ind]=0; A[ind]=0; h14[ind]=0; E[ind]=0; B[ind]=0; h44[ind]=0


        l.par <- weights*(i0*log(F1) + i1*(log(1-z)+log(f2)))

        # Gradient:

        dl.1  <- weights*( -i0*F1^(-1) + i1*p )*ph                  # dl.dbe1
        dl.2  <- weights*i1*( h + e2/sqv )                          # dl.dbe2
        dl.3  <- weights*i1*( h*e2*sqv + e2^2 - 1 )                 # dl.dsqv.st
        dl.4  <- weights*i1*b                                       # dl.dteta.st


        # (Minus) Hessian:

        d2l.11 <- -weights*ph*( -i0*(ph/F1-eta1)/F1 + i1*P )        # d2l.be1.be1
        d2l.12 <- -weights*A                                        # d2l.be1.be2
        d2l.13 <- -weights*A*e2*sqv                                 # d2l.be1.sqv.st
        d2l.14 <- -weights*h14                                      # d2l.be1.teta.st

        d2l.22 <- -weights*i1*( h*E - sqv^(-2) )                    # d2l.be2.be2
        d2l.23 <- -weights*i1*( h*(E*e2*sqv-1) - 2*e2/sqv )         # d2l.be2.sqv.st
        d2l.24 <- -weights*B                                        # d2l.be2.teta.st 

        d2l.33 <- -weights*i1*e2*(h*(E*e2*sqv^2-sqv)-2*e2)          # d2l.sqv.st.sqv.st
        d2l.34 <- -weights*e2*B*sqv                                 # d2l.sqv.st.teta.st

        d2l.44 <- -weights*h44                                      # d2l.teta.st.teta.st

    }

}
else if(margins[2]=="G"){

    k.st <- params[(X1.d2+X2.d2+1)]
    k <- exp(k.st)
    i2 <- dat[,2]
    F2 <- pgamma(i2,shape=k,rate=k*exp(-eta2));		
    F2 <- ifelse(F2>0.000000001,F2,0.000000001)
    F2 <- ifelse(F2<0.999999999,F2,0.999999999)
    f2 <- dgamma(i2,shape=k,rate=k*exp(-eta2));		
    f2[i2==0] <- 0
    f2 <- ifelse(f2>0.000000001,f2,0.000000001)

    ye <- i2*exp(-eta2)

    G <- gamma(k)
    psi <- psigamma(k)
    psiprim <- psigamma(k,deriv=1)
    fun1 <- function(t,k) { t^(k-1)*exp(-t)*log(t) }
    Gkf  <- function(k, y) { sapply(y, function(y) integrate(fun1,y,Inf,k)$value ) }
    Gk   <- Gkf(k,k*ye)
    fun2 <- function(t,k){t^(k-1)*exp(-t)*(log(t))^2}
    G2kf <- function(k, y) { sapply(y, function(y) integrate(fun2,y,Inf,k)$value ) }
    G2k  <- G2kf(k,k*ye)
    dF2k <- i2*f2*k^(-1)+psi*(1-F2)-Gk/G
    df2kbyf2 <- -psi+log(k*i2)+1-eta2-ye;  
    df2kbyf2[i1==0]=0;
    d2F2k <-i2*k^(-1)*f2*(2*df2kbyf2+ye-1-1/k) + (1-F2)*(psiprim-psi^2) + 2*psi*Gk*G^(-1)-G2k/G 



    bits <- bitsgHs(cop=BivD,margin=margins[2],i1=i1,F1=F1,F2=F2,f2=f2,eta1=eta1,ph=ph,teta=teta,ver=0)

    z <- bits$z; p <- bits$p; h <- bits$h; b <- bits$b; P <- bits$P; A <- bits$A; h14 <- bits$h14; E <- bits$E; B <- bits$B; h44 <- bits$h44 

    z[ind]=0; p[ind]=0; h[ind]=0; b[ind]=0; P[ind]=0; A[ind]=0; h14[ind]=0; E[ind]=0; B[ind]=0; h44[ind]=0


    l.par <- weights*(i0*log(F1) + i1*(log(1-z)+log(f2)))


    # Gradient:

    dl.1  <- weights*( -i0*F1^(-1) + i1*p )*ph                   # dl.dbe1
    dl.2  <- weights*i1*( i2*f2*h - k*(1-ye) )                   # dl.dbe2
    dl.3  <- weights*i1*( -dF2k*h + df2kbyf2 )*k                 # dl.dk.st 
    dl.4  <- weights*i1*b                                        # dl.dteta.st

 
    # (Minus) Hessian:

    d2l.11  <- -weights*ph*( -i0*(ph/F1-eta1)/F1 + i1*P )        # d2l.be1.be1
    d2l.12  <- -weights*i2*f2*A                                  # d2l.be1.be2
    d2l.13  <- weights*A*dF2k*k                                  # d2l.be1.k.st
    d2l.14  <- -weights*h14                                      # be1.teta.st

    d2l.22  <- -weights*i1*i2*( f2*h*(i2*f2*E+k*ye-k) - k*exp(-eta2) )    # d2l.be2.be2
    d2l.23  <- -weights*i1*( -i2*f2*h*(dF2k*E-df2kbyf2) - 1 + ye )*k      # d2l.be2.k.st
    d2l.24  <- -weights*i2*B*f2                                           # d2l.be2.teta.st 

    d2l.33   <- -weights*i1*( -dF2k*h + df2kbyf2 - h*k*(d2F2k-dF2k^2*E) + 1 - k*psiprim )*k    # d2l.k.st.k.st
    d2l.34 <- weights*dF2k*B*k                                                                 # d2l.k.st.teta.st

    d2l.44 <- -weights*h44                                       # d2l.teta.st.teta.st

} 




  H11 <- crossprod(dat1*c(d2l.11),dat1)
  H12 <- crossprod(dat1*c(d2l.12),dat2) 
  H13 <- t(t(rowSums(t(dat1*c(d2l.13)))))
  H14 <- t(t(rowSums(t(dat1*c(d2l.14)))))

  H22 <- crossprod(dat2*c(d2l.22),dat2) 
  H23 <- t(t(rowSums(t(dat2*c(d2l.23)))))
  H24 <- t(t(rowSums(t(dat2*c(d2l.24)))))

  H <- rbind( cbind( H11    , H12    , H13  ,  H14 ), 
              cbind( t(H12) , H22    , H23  ,  H24 ),
              cbind( t(H13) , t(H23) , sum(d2l.33), sum(d2l.34) ) ,
              cbind( t(H14) , t(H24) , sum(d2l.34), sum(d2l.44) )
            ) 

  res <- -sum(l.par)

  G   <- c( -colSums( c(dl.1)*dat1 ) ,
            -colSums( c(dl.2)*dat2 )    ,
            -sum( dl.3 ) ,  
            -sum( dl.4 )   )


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
       dl.dbe1=dl.1, dl.dbe2=dl.2, l.par=l.par,
       dl.dsqv.st=dl.3,
       dl.dcor.st=dl.4, 
       d2l.be1.be1=d2l.11, d2l.be1.be2=d2l.12, d2l.be2.be2=d2l.22,
       d2l.be1.sqv.st=d2l.13,
       d2l.be1.cor.st=d2l.14,
       d2l.be2.sqv.st=d2l.23, 
       d2l.be2.cor.st=d2l.24,
       d2l.sqv.st.sqv.st=d2l.33,
       d2l.sqv.st.cor.st=d2l.34,    
       d2l.cor.st.cor.st=d2l.44 )

}






