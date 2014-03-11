ghss2 <- function(BivD, margins, params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, params.f, sp=NULL, qu.mag=NULL){

  eta1 <- dat1%*%params.f[1:X1.d2]  
  eta2 <- dat2%*%params.f[(X1.d2+1):(X1.d2+X2.d2)]   
  teta.st <- params[2]
  eps <- .Machine$double.xmin*10^16

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


    sqv.st <- params[1]
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
 

      dl.3 <- weights*i1*( -1 + (e2/sqv)^2 - l2*((teta*e2)/(a*sqv)) ) 
      dl.4 <- weights*( i1*( l2*( ((2*ec*e2)/((ec+1)*sqv) - (2*teta*e2*ec)/((ec+1)*sqv))/a  - 0.5*( (eta1+(teta*e2)/sqv)*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1) ) )/a^3 )  )  )  


      d2l.33 <- -weights*( i1*(  -2*(e2/sqv)^2 + l2*((e2*teta)/(a*sqv)) + PA*((teta*e2)/(sqv*a))^2  )   )
      d2l.34 <- -weights*( ( i1*(   -PA*(    (( (2*e2*ec)/((ec+1)*sqv) -( 2*teta*e2*ec)/((ec+1)*sqv) )/a 
         				 -0.5/a^3*( ( eta1+(teta*e2)/sqv )*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1) )))*(ec-1)*e2  
      				)/((ec+1)*sqv*a)
   				- (2*l2*ec*e2)/(sqv*a*(ec+1))      
  				 + (2*l2*teta*e2*ec)/(sqv*a*(ec+1))  
  				 + 0.5*l2*teta*e2/(sqv*a^3)*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1)   ) 
                     )    ) )
      d2l.44 <- -weights*( ( i1*(   PA*(       (
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
                 )    ) )

    }
    else{


        e2 <- (dat[,2]-eta2)/sqv
        F2 <- pnorm(e2)
        F2 <- ifelse(F2>0.00000001,F2,0.00000001)
        F2 <- ifelse(F2<0.99999999,F2,0.99999999)
        f2 <- dnorm(e2)/sqv;
        f2 <- ifelse(f2>0.00000001,f2,0.00000001)


        bits <- bitsgHs(cop=BivD,margin=margins[2],i1=i1,F1=F1,F2=F2,f2=f2,eta1=eta1,ph=ph,teta=teta,e2=e2,sqv=sqv,ver=2)

        z <- bits$z; h <- bits$h; b <- bits$b; E <- bits$E; B <- bits$B; h44 <- bits$h44 
        z[ind]=0; h[ind]=0; b[ind]=0; E[ind]=0; B[ind]=0; h44[ind]=0

        l.par <- weights*(i0*log(F1) + i1*(log(1-z)+log(f2)))

        # Gradient:

        dl.3  <- weights*i1*( h*e2*sqv + e2^2 - 1 )                 # dl.dsqv.st
        dl.4  <- weights*i1*b                                       # dl.dteta.st


        # (Minus) Hessian:

        d2l.33 <- -weights*i1*e2*(h*(E*e2*sqv^2-sqv)-2*e2)          # d2l.sqv.st.sqv.st
        d2l.34 <- -weights*e2*B*sqv                                 # d2l.sqv.st.teta.st

        d2l.44 <- -weights*h44                                      # d2l.teta.st.teta.st

    }

}
else if(margins[2]=="G"){

    k.st <- params[1]
    k <- exp(k.st)
    i2 <- dat[,2]
    F2 <- pgamma(i2,shape=k,rate=k*exp(-eta2));		
    F2 <- ifelse(F2>0.00000001,F2,0.00000001)
    F2 <- ifelse(F2<0.99999999,F2,0.99999999)
    f2 <- dgamma(i2,shape=k,rate=k*exp(-eta2));		
    f2[i2==0] <- 0
    f2 <- ifelse(f2>0.00000001,f2,0.00000001)

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


    bits <- bitsgHs(cop=BivD,margin=margins[2],i1=i1,F1=F1,F2=F2,f2=f2,eta1=eta1,ph=ph,teta=teta,ver=2)

    z <- bits$z; h <- bits$h; b <- bits$b; E <- bits$E; B <- bits$B; h44 <- bits$h44 
    z[ind]=0; h[ind]=0; b[ind]=0; E[ind]=0; B[ind]=0; h44[ind]=0

    l.par <- weights*(i0*log(F1) + i1*(log(1-z)+log(f2)))


    # Gradient:

    dl.3  <- weights*i1*( -dF2k*h + df2kbyf2 )*k                 # dl.dk.st 
    dl.4  <- weights*i1*b                                        # dl.dteta.st

 
    # (Minus) Hessian:

    d2l.33   <- -weights*i1*( -dF2k*h + df2kbyf2 - h*k*(d2F2k-dF2k^2*E) + 1 - k*psiprim )*k    # d2l.k.st.k.st
    d2l.34 <- weights*dF2k*B*k                                                                 # d2l.k.st.teta.st

    d2l.44 <- -weights*h44                                       # d2l.teta.st.teta.st

} 


  H <- rbind( 
              cbind( sum(d2l.33), sum(d2l.34) ) ,
              cbind( sum(d2l.34), sum(d2l.44) )
            ) 

  
  res <- -sum(l.par)

  G   <- c( -sum( dl.3 ) ,  
            -sum( dl.4 )   )


  list(value=res, gradient=G, hessian=H, l=res, dat1=dat1, dat2=dat2) 
     
}
