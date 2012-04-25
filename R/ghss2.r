ghss2 <- function(params, dat, X1.d2, X2.d2, gam1, gam2, fp, params.f, dat1, dat2){

  #dat1 <- as.matrix(dat[,3:(X1.d2+2)])
  #dat2 <- as.matrix(dat[,(X1.d2+3):(X1.d2+X2.d2+2)])
  eta1 <- dat1%*%params.f[1:X1.d2]  
  eta2 <- dat2%*%params.f[(X1.d2+1):(X1.d2+X2.d2)]   
  sqv.st <- params[1]
  cor.st <- params[2]
  sqv <- exp(sqv.st)
  cor <- tanh(cor.st)
  
  a  <- sqrt(1-cor^2)
  e2 <- dat[,2]-eta2  
  A  <- (eta1+(cor/sqv)*e2)/a
  l1 <- dnorm(-eta1)/pnorm(-eta1)
  l2 <- dnorm(A)/pnorm(A)
  ec <- exp(2*cor.st) 

  PA    <- -(pnorm(A)*dnorm(A)*A+dnorm(A)^2)/pnorm(A)^2 
  Peta1 <- -(pnorm(-eta1)*dnorm(-eta1)*(-eta1)+dnorm(-eta1)^2)/pnorm(-eta1)^2  

  i0 <- 1-dat[, 1]
  i1 <- dat[, 1] 

  l.par <- i0*log(pnorm(-eta1)) + i1*(log(pnorm(A))-1/2*log(2*pi)-log(sqv)-1/2*(e2/sqv)^2) 
 

  dl.dsqv.st <- i1*( -1 + (e2/sqv)^2 - l2*((cor*e2)/(a*sqv)) ) 
  dl.dcor.st <- i1*( l2*( ((2*ec*e2)/((ec+1)*sqv) - (2*cor*e2*ec)/((ec+1)*sqv))/a  - 0.5*( (eta1+(cor*e2)/sqv)*( -4*cor*ec/(ec+1) + 4*cor^2*ec/(ec+1) ) )/a^3 )  )   



  d2l.sqv.st.sqv.st <- -( i1*(  -2*(e2/sqv)^2 + l2*((e2*cor)/(a*sqv)) + PA*((cor*e2)/(sqv*a))^2  )   )
  d2l.sqv.st.cor.st <- -( i1*(   -PA*(    (( (2*e2*ec)/((ec+1)*sqv) -( 2*cor*e2*ec)/((ec+1)*sqv) )/a 
         				 -0.5/a^3*( ( eta1+(cor*e2)/sqv )*( -4*cor*ec/(ec+1) + 4*cor^2*ec/(ec+1) )))*(ec-1)*e2  
      				)/((ec+1)*sqv*a)
   				- (2*l2*ec*e2)/(sqv*a*(ec+1))      
  				 + (2*l2*cor*e2*ec)/(sqv*a*(ec+1))  
  				 + 0.5*l2*cor*e2/(sqv*a^3)*( -4*cor*ec/(ec+1) + 4*cor^2*ec/(ec+1)   ) 
                     )    )
  d2l.cor.st.cor.st <- -( i1*(   PA*(       (
        				 ( (2*e2*ec)/((ec+1)*sqv) -(2*cor*e2*ec)/((ec+1)*sqv) )/a 
       				 -0.5/a^3*( ( eta1+(cor*e2)/sqv )*( -4*cor*ec/(ec+1) + 4*cor^2*ec/(ec+1) )) 
          				   )^2
     				 )
 				+ l2*(
 		 1/a*( 4*ec*e2/((ec+1)*sqv) - 8*ec^2*e2/((ec+1)^2*sqv) + 8*cor*ec^2*e2/((ec+1)^2*sqv)  - 4*cor*ec*e2/((ec+1)*sqv)  ) 
  		-1/a^3*( (2*ec*e2/((ec+1)*sqv) - 2*ec*e2*cor/((ec+1)*sqv) )*( -4*cor*ec/(ec+1) + 4*cor^2*ec/(ec+1) )  ) 
  		+3/(4*a^5)*( ( eta1+(cor*e2)/sqv )*( -4*cor*ec/(ec+1) + 4*cor^2*ec/(ec+1) )^2 ) 
  		-0.5/a^3*( (eta1+(cor*e2)/sqv )*(  -8*ec^2/(ec+1)^2 + 32*cor*ec^2/(ec+1)^2 - 8*cor*ec/(ec+1) - 24*cor^2*ec^2/(ec+1)^2  + 8*cor^2*ec/(ec+1)    ) 
       		    )
     		 )
                     )    )
                                    



  H <- rbind( 
              cbind( sum(d2l.sqv.st.sqv.st), sum(d2l.sqv.st.cor.st) ) ,
              cbind( sum(d2l.sqv.st.cor.st), sum(d2l.cor.st.cor.st) )
            ) 

  
         res <- -sum(l.par)
         G   <- c( -sum( dl.dsqv.st ) ,  
                   -sum( dl.dcor.st )   )
         H   <- H

         list(value=res, gradient=G, hessian=H, l=res, dat1=dat1, dat2=dat2) 
  
                
     
}


