ghss1 <- function(params, dat, X1.d2, X2.d2, S=0, gam1, gam2, fp, exl1, exl2, dat1, dat2){

  # dat1 <- as.matrix(dat[,3:(X1.d2+2)])
  # dat2 <- as.matrix(dat[,(X1.d2+3):(X1.d2+X2.d2+2)])
  eta1 <- dat1%*%params[1:X1.d2]  
  eta2 <- dat2%*%params[(X1.d2+1):(X1.d2+X2.d2)]   
  sqv.st <- exl1
  cor.st <- exl2
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
 
  dl.dbe1 <- -i0*l1 + i1*(l2/a) 
  dl.dbe2 <- i1*( e2/sqv^2 - l2*( cor/(sqv*a) )   )
  
  d2l.be1.be1 <- -( i0*Peta1 + i1*(PA/a^2) )  
  d2l.be1.be2 <- -( -i1*PA*(cor/(sqv*a^2)) )   
  d2l.be2.be2 <- -( i1*( PA*(cor/(sqv*a))^2 - 1/sqv^2 ) ) 
                               
  be1.be1 <- t(dat1*c(d2l.be1.be1))%*%dat1
  be1.be2 <- t(dat1*c(d2l.be1.be2))%*%dat2 
  be2.be2 <- t(dat2*c(d2l.be2.be2))%*%dat2 

  H <- rbind( cbind( be1.be1    , be1.be2 ), 
              cbind( t(be1.be2) , be2.be2 )   ) 

  if( ( length(gam1$smooth)==0 && length(gam2$smooth)==0 ) || fp==TRUE){
         res <- -sum(l.par)
         G   <- c( -colSums( c(dl.dbe1)*dat1 ) ,
                   -colSums( c(dl.dbe2)*dat2 )  
                     )
         H   <- H
         list(value=res, gradient=G, hessian=H, l=res, eta1=eta1, eta2=eta2) 
  }else{
  
         S.h <- adiag(matrix(0,gam1$nsdf,gam1$nsdf),
                      S[1:(X1.d2-gam1$nsdf),1:(X1.d2-gam1$nsdf)],
                      matrix(0,gam2$nsdf,gam2$nsdf),
                      S[(X1.d2-(gam1$nsdf-1)):dim(S)[2],(X1.d2-(gam1$nsdf-1)):dim(S)[2]])  
                              
         S.res <- -sum(l.par)
         res <- S.res + (1/2)*(t(params)%*%S.h%*%params)
         G   <- c( -colSums( c(dl.dbe1)*dat1 ) ,
                   -colSums( c(dl.dbe2)*dat2 )   ) + S.h%*%params

         H   <- H + S.h  
         list(value=res, gradient=G, hessian=H, S.h=S.h, l=S.res, eta1=eta1, eta2=eta2, dat1=dat1, dat2=dat2, 
             dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, 
             d2l.be1.be1=d2l.be1.be1, d2l.be1.be2=d2l.be1.be2, d2l.be2.be2=d2l.be2.be2)
                
   }     

}


