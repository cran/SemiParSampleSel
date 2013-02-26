ghssFGM <- function(params, dat, dat1, dat2, weights=weights, X1.d2, X2.d2, sp=NULL, l.sp1, l.sp2, gp1, gp2, fp, qu.mag=NULL){


  eta1 <- dat1%*%params[1:X1.d2]  
  eta2 <- dat2%*%params[(X1.d2+1):(X1.d2+X2.d2)]   
  sqv.st <- params[(X1.d2+X2.d2+1)]
  teta.st <- params[(X1.d2+X2.d2+2)]
  sqv <- exp(sqv.st)
  teta <- tanh(teta.st)

  
# Settings:

  i0 <- 1-dat[, 1]
  i1 <- dat[, 1]

  e2 <- (dat[,2]-eta2)/sqv
  F1 <- pnorm(-eta1)
  F2 <- pnorm(e2)
  f2 <- dnorm(e2)/sqv
  u <- 1-teta*F1*(1-2*F2)
  ph <- dnorm(-eta1)
  h <- 2*teta*F1*f2*u^(-1)


# Likelihood:

  l.par <- weights*(i0*log(F1) + i1*(log(1-F1)+log(1-teta*F1*(1-2*F2))+log(f2))) 


# Gradient:

  dl.dbe1 <- weights*(-i0*F1^(-1)+i1*((1-F1)^(-1) + teta*(1-2*F2)/u))*ph
  dl.dbe2 <- weights*i1*(-h+e2/sqv)
  dl.dsqv.st <- weights*i1*(-h*e2*sqv+e2^2-1)
  dl.dteta.st <- weights*i1*(1-u^(-1))*(1-teta^2)*teta^(-1)


# Hessian:

  d2l.be1.be1 <- -weights*( ph*(-i0*(ph/F1-eta1)/F1- i1*( 
		ph*((1-F1)^(-2) + (teta*(1-2*F2)/u)^2) + 
		eta1*((1-F1)^(-1) + teta*(1-2*F2)/u )
		) ) )
  d2l.be1.be2 <- -weights*i1*2*teta*u^(-2)*f2*ph
  d2l.be1.sqv.st <- -weights*i1*2*teta*u^(-2)*f2*ph*e2*sqv
  d2l.be1.teta.st <- -weights*i1*(1-2*F2)*u^(-2)*ph*(1-teta^2)

  d2l.be2.be2 <- -weights*i1*(-h*u^(-1)*(e2*u*sqv^(-1) + h*u)-sqv^(-2))
  d2l.be2.sqv.st <- -weights*i1*(-h*u^(-1)*((e2^2-1)*u*sqv^(-1) + h*u*e2)-2*e2/sqv^2)*sqv
  d2l.be2.teta.st <- weights*i1*2*F1*f2*u^(-2)*(1-teta^2)

  d2l.sqv.st.sqv.st <- -weights*i1*(-h*u^(-1)*((e2^2-1)*u*sqv^(-1) + h*u*e2)*sqv^2-2*e2)*e2
  d2l.sqv.st.teta.st <- weights*i1*2*F1*f2*u^(-2)*e2*sqv*(1-teta^2)

  d2l.teta.st.teta.st <- -weights*i1*(u-1)*u^(-1)*((1-teta^2)*teta^(-2)*u^(-1)-(1+teta^2)*teta^(-2))*(1-teta^2)


  be1.be1 <- crossprod(dat1*c(d2l.be1.be1),dat1)
  be1.be2 <- crossprod(dat1*c(d2l.be1.be2),dat2) 
  be1.sqv.st <- t(t(rowSums(t(dat1*c(d2l.be1.sqv.st)))))
  be1.teta.st <- t(t(rowSums(t(dat1*c(d2l.be1.teta.st)))))

  be2.be2 <- crossprod(dat2*c(d2l.be2.be2),dat2) 
  be2.sqv.st <- t(t(rowSums(t(dat2*c(d2l.be2.sqv.st)))))
  be2.teta.st <- t(t(rowSums(t(dat2*c(d2l.be2.teta.st)))))


  H <- rbind( cbind( be1.be1    , be1.be2    , be1.sqv.st  ,  be1.teta.st ), 
              cbind( t(be1.be2) , be2.be2    , be2.sqv.st  ,  be2.teta.st ),
              cbind( t(be1.sqv.st) , t(be2.sqv.st) , sum(d2l.sqv.st.sqv.st), sum(d2l.sqv.st.teta.st) ) ,
              cbind( t(be1.teta.st) , t(be2.teta.st) , sum(d2l.sqv.st.teta.st), sum(d2l.teta.st.teta.st) )
            ) 


  res <- -sum(l.par)
  G   <- c( -colSums( c(dl.dbe1)*dat1 ) ,
            -colSums( c(dl.dbe2)*dat2 )    ,
            -sum( dl.dsqv.st ) ,  
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
       dl.dsqv.st=dl.dsqv.st,
       dl.dcor.st=dl.dteta.st, 
       d2l.be1.be1=d2l.be1.be1, d2l.be1.be2=d2l.be1.be2, d2l.be2.be2=d2l.be2.be2,
       d2l.be1.sqv.st=d2l.be1.sqv.st,
       d2l.be1.cor.st=d2l.be1.teta.st,
       d2l.be2.sqv.st=d2l.be2.sqv.st, 
       d2l.be2.cor.st=d2l.be2.teta.st,
       d2l.sqv.st.sqv.st=d2l.sqv.st.sqv.st,
       d2l.sqv.st.cor.st=d2l.sqv.st.teta.st,    
       d2l.cor.st.cor.st=d2l.teta.st.teta.st )
   
}
