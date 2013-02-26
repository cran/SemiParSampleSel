ghssN <- function(params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, sp=NULL, l.sp1, l.sp2, gp1, gp2, fp, qu.mag=NULL){

  eta1 <- dat1%*%params[1:X1.d2]  
  eta2 <- dat2%*%params[(X1.d2+1):(X1.d2+X2.d2)]   
  sqv.st <- params[(X1.d2+X2.d2+1)]
  cor.st <- params[(X1.d2+X2.d2+2)]
  sqv <- exp(sqv.st)
  cor <- tanh(cor.st)
  
  a  <- sqrt(1-cor^2)
  e2 <- dat[,2]-eta2  
  A  <- (eta1+(cor/sqv)*e2)/a
  l1 <- ifelse(eta1<37.5,dnorm(-eta1)/pnorm(-eta1),eta1)
  l2 <- ifelse(A>-37.5,dnorm(A)/pnorm(A),-A)
  ec <- exp(2*cor.st) 

  PA    <- -(pnorm(A)*dnorm(A)*A+dnorm(A)^2)/pnorm(A)^2 
  Peta1 <- -(pnorm(-eta1)*dnorm(-eta1)*(-eta1)+dnorm(-eta1)^2)/pnorm(-eta1)^2  

  i0 <- 1-dat[, 1]
  i1 <- dat[, 1] 

  l.par <- weights*(i0*log(pnorm(-eta1)) + i1*(log(pnorm(A))-1/2*log(2*pi)-log(sqv)-1/2*(e2/sqv)^2)) 
 
  dl.dbe1 <- -weights*( i0*l1 - i1*(l2/a) ) 
  dl.dbe2 <- weights*i1*( e2/sqv^2 - l2*( cor/(sqv*a) )   )
  dl.dsqv.st <- weights*i1*( -1 + (e2/sqv)^2 - l2*((cor*e2)/(a*sqv)) ) 
  dl.dcor.st <- weights*i1*( l2*( ((2*ec*e2)/((ec+1)*sqv) - (2*cor*e2*ec)/((ec+1)*sqv))/a  - 0.5*( (eta1+(cor*e2)/sqv)*( -4*cor*ec/(ec+1) + 4*cor^2*ec/(ec+1) ) )/a^3 )  )   


  d2l.be1.be1 <- -weights*( i0*Peta1 + i1*(PA/a^2) )  
  d2l.be1.be2 <- -weights*( -i1*PA*(cor/(sqv*a^2)) )   
  d2l.be1.sqv.st <- -weights*( -i1*PA*((cor*e2)/(sqv*a^2)) ) 
  d2l.be1.cor.st <- -weights*( i1*(       (PA/a)*(  ( (2*ec*e2)/((ec+1)*sqv) - (2*(ec-1)*e2*ec)/((ec+1)^2*sqv) )/a      
               				- 0.5*(eta1+(cor*e2)/sqv)*((-4*cor*ec)/(ec+1)+(4*cor^2*ec)/(ec+1) )/a^3 
         				  )  
   				 -0.5*(l2/a^3)*( (-4*cor*ec)/(ec+1) + (4*cor^2*ec)/(ec+1) )           
                     )    )
  d2l.be2.be2 <- -weights*( i1*( PA*(cor/(sqv*a))^2 - 1/sqv^2 ) ) 
  d2l.be2.sqv.st <- -weights*( i1*( PA*(cor^2*e2)/(sqv^2*a^2) + (cor*l2)/(sqv*a) - (2*e2)/sqv^2 ) )
  d2l.be2.cor.st <- -weights*( i1*(     -PA*(    
    				( (  (2*e2*ec)/((ec+1)*sqv) - (2*cor*e2*ec)/((ec+1)*sqv)  )/a 
      				   -0.5/a^3*( ( eta1+(cor*e2)/sqv )*( -4*cor*ec/(ec+1) + 4*cor^2*ec/(ec+1) ))   )*(ec-1)  
     				)/((ec+1)*sqv*a)
   				- (2*l2*ec)/(sqv*a*(ec+1))      
  				 + (2*l2*cor*ec)/(sqv*a*(ec+1))  
  				 + 0.5*l2*cor/(sqv*a^3)*( -4*cor*ec/(ec+1) + 4*cor^2*ec/(ec+1)   ) 
                     )    )

  d2l.sqv.st.sqv.st <- -weights*( i1*(  -2*(e2/sqv)^2 + l2*((e2*cor)/(a*sqv)) + PA*((cor*e2)/(sqv*a))^2  )   )
  d2l.sqv.st.cor.st <- -weights*( i1*(   -PA*(    (( (2*e2*ec)/((ec+1)*sqv) -( 2*cor*e2*ec)/((ec+1)*sqv) )/a 
         				 -0.5/a^3*( ( eta1+(cor*e2)/sqv )*( -4*cor*ec/(ec+1) + 4*cor^2*ec/(ec+1) )))*(ec-1)*e2  
      				)/((ec+1)*sqv*a)
   				- (2*l2*ec*e2)/(sqv*a*(ec+1))      
  				 + (2*l2*cor*e2*ec)/(sqv*a*(ec+1))  
  				 + 0.5*l2*cor*e2/(sqv*a^3)*( -4*cor*ec/(ec+1) + 4*cor^2*ec/(ec+1)   ) 
                     )    )
  d2l.cor.st.cor.st <- -weights*( i1*(   PA*(       (
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
                                    
  be1.be1 <- crossprod(dat1*c(d2l.be1.be1),dat1)
  be1.be2 <- crossprod(dat1*c(d2l.be1.be2),dat2) 
  be1.sqv.st <- t(t(rowSums(t(dat1*c(d2l.be1.sqv.st))))) 
  be1.cor.st <- t(t(rowSums(t(dat1*c(d2l.be1.cor.st)))))
  
  be2.be2 <- crossprod(dat2*c(d2l.be2.be2),dat2) 
  be2.sqv.st <- t(t(rowSums(t(dat2*c(d2l.be2.sqv.st)))))
  be2.cor.st <- t(t(rowSums(t(dat2*c(d2l.be2.cor.st)))))


  H <- rbind( cbind( be1.be1    , be1.be2    , be1.sqv.st  ,  be1.cor.st ), 
              cbind( t(be1.be2) , be2.be2    , be2.sqv.st  ,  be2.cor.st ),
              cbind( t(be1.sqv.st) , t(be2.sqv.st) , sum(d2l.sqv.st.sqv.st), sum(d2l.sqv.st.cor.st) ) ,
              cbind( t(be1.cor.st) , t(be2.cor.st) , sum(d2l.sqv.st.cor.st), sum(d2l.cor.st.cor.st) )

            ) 

  res <- -sum(l.par)
  G   <- c( -colSums( c(dl.dbe1)*dat1 ) ,
            -colSums( c(dl.dbe2)*dat2 )    ,
            -sum( dl.dsqv.st ) ,  
            -sum( dl.dcor.st )   )



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
       dl.dsqv.st=dl.dsqv.st,
       dl.dcor.st=dl.dcor.st, 
       d2l.be1.be1=d2l.be1.be1, d2l.be1.be2=d2l.be1.be2, d2l.be2.be2=d2l.be2.be2,
       d2l.be1.sqv.st=d2l.be1.sqv.st,
       d2l.be1.cor.st=d2l.be1.cor.st,
       d2l.be2.sqv.st=d2l.be2.sqv.st, 
       d2l.be2.cor.st=d2l.be2.cor.st,
       d2l.sqv.st.sqv.st=d2l.sqv.st.sqv.st,
       d2l.sqv.st.cor.st=d2l.sqv.st.cor.st,    
       d2l.cor.st.cor.st=d2l.cor.st.cor.st )

}

