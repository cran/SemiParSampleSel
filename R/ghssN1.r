ghssN1 <- function(params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, sp=NULL, l.sp1, l.sp2, gp1, gp2, exl1, exl2, qu.mag=NULL){

  eta1 <- dat1%*%params[1:X1.d2]  
  eta2 <- dat2%*%params[(X1.d2+1):(X1.d2+X2.d2)]   
  sqv.st <- exl1
  cor.st <- exl2
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

  l.par <- weights*(i0*log(pnorm(-eta1)) + i1*(log(pnorm(A))-1/2*log(2*pi)-log(sqv)-1/2*(e2/sqv)^2) )
 
  dl.dbe1 <- -weights*(i0*l1 - i1*(l2/a)) 
  dl.dbe2 <- weights*i1*( e2/sqv^2 - l2*( cor/(sqv*a) )   )
  
  d2l.be1.be1 <- -weights*( i0*Peta1 + i1*(PA/a^2) )  
  d2l.be1.be2 <- -weights*( -i1*PA*(cor/(sqv*a^2)) )   
  d2l.be2.be2 <- -weights*( i1*( PA*(cor/(sqv*a))^2 - 1/sqv^2 ) ) 

    be1.be1 <- crossprod(dat1*c(d2l.be1.be1),dat1)
    be1.be2 <- crossprod(dat1*c(d2l.be1.be2),dat2) 
    be2.be2 <- crossprod(dat2*c(d2l.be2.be2),dat2) 
    

  H <- rbind( cbind( be1.be1    , be1.be2 ), 
              cbind( t(be1.be2) , be2.be2 )   ) 

  
  res <- -sum(l.par)
  G   <- c( -colSums( c(dl.dbe1)*dat1 ) ,
            -colSums( c(dl.dbe2)*dat2 )   )


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
       dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, 
       d2l.be1.be1=d2l.be1.be1, d2l.be1.be2=d2l.be1.be2, d2l.be2.be2=d2l.be2.be2)
        

}


