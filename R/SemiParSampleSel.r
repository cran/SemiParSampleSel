SemiParSampleSel <- function(formula.eq1, formula.eq2, data=list(), 
                             iterlimSP=50, pr.tol=1e-6,
                             gamma=1, aut.sp=TRUE, fp=FALSE, start.v=NULL, rinit=1, rmax=100, 
                             fterm=sqrt(.Machine$double.eps), mterm=sqrt(.Machine$double.eps), 
        		     control=list(maxit=50,tol=1e-6,step.half=25,rank.tol=sqrt(.Machine$double.eps)) ){

  conv.sp <- NULL; bs.mgfit <- wor.c <- j.it <- NULL

  #rinit=1;rmax=100;fterm=sqrt(.Machine$double.eps);mterm=sqrt(.Machine$double.eps)
  #gamma=1;control=list(maxit=50,tol=1e-6,step.half=25,rank.tol=.Machine$double.eps^0.5)
  #pr.tol=1e-10; data=list()
  #aut.sp=TRUE; fp=FALSE; start.v=NULL; gcv=FALSE; iterlimFS=1; iterlimSP=25

  #library(mgcv)
  #library(trust)
  #library(mvtnorm)

  #source("ghss.r")
  #source("ghss1.r")
  #source("ghss2.r")
  #source("spS.r")
  #source("S.m.r")
  #source("working.comp.r")

  #f1   <- function(x) -(-1+x-1.6*x^2+sin(5*x))    #sin(x)+1.5*exp(-10*x^2) #(0.08*(x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10))
  #f2   <- function(x) 4*x 
  #f3   <- function(x) (0.08*(x^11*(10*(1-x))^6+10*(10*x)^3*(1-x)^10))

  #eps <- rmvnorm(1000, c(0, 0), matrix(c(1, -0.7, -0.7, 1), 2, 2))
  #xs <- runif(1000)
  #ys <- f1(xs) + eps[, 1] > 0
  #xo <- runif(1000)
  #yoX <- f3(xo) + eps[, 2]
  #yo <- yoX * (ys > 0)

  #data <- as.data.frame(cbind(ys,yo,xs,xo))

  #formula.eq1 <- ys ~ s(xs) 
  #formula.eq2 <- yo ~ s(xo)


  ################################

  gam1   <- gam(formula.eq1, binomial(link="probit"), gamma=gamma, data=data)
  inde <- gam1$y > 0
  environment(formula.eq2) <- environment(NULL)
  gam2  <- gam(formula.eq2, gamma=gamma, data=data, subset=inde)
  environment(gam2$formula) <- environment(gam1$formula)

  X1 <- model.matrix(gam1); X1.d2 <- dim(X1)[2]; X2.d2 <- length(coef(gam2))
  X2 <- matrix(0,length(inde),X2.d2,dimnames = list(c(1:length(inde)),c(names(coef(gam2)))) )
  X2[inde, ] <- model.matrix(gam2)
  y2 <- rep(0,length(inde)); y2[inde] <- gam2$y
  l.sp1 <- length(gam1$smooth); l.sp2 <- length(gam2$smooth)
  dat <- cbind(gam1$y,y2,X1,X2); n <- length(dat[,1])

  p.g1 <- predict(gam1)
  imr <- dnorm(p.g1)/pnorm(p.g1)

  formula.eq2.1 <- update.formula(formula.eq2, ~. + imr)
  environment(formula.eq2.1) <- environment(NULL)

  gam2.1 <- gam(formula.eq2.1, gamma=gamma, data=data, subset=inde) 
  environment(gam2.1$formula) <- environment(gam1$formula)
  
  sigma <- sqrt(mean(resid(gam2.1)^2)+mean(imr[inde]*(imr[inde]+p.g1[inde]))*gam2.1$coef["imr"]^2)[[1]]; names(sigma) <- "sigma"
  co  <- (gam2.1$coef["imr"]/sigma)[[1]]

  rho <- ifelse( abs(co) > 0.99, sign(co)*0.95, co); names(rho) <- "rho"

  sp <- c(gam1$sp,gam2.1$sp)

  if(is.null(start.v)) start.v <- c(coef(gam1),coef(gam2.1)[names(coef(gam2.1))!="imr"],log(sigma),atanh(rho))

  if(l.sp1!=0 && l.sp2!=0 && fp==FALSE){S <- spS(sp,gam1,gam2); qu.mag <- S.m(gam1,gam2)}

  fit  <- trust(ghss, start.v, rinit=rinit, rmax=rmax, dat=dat,
                X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, 
                fterm=fterm, mterm=mterm, iterlim=1e+4)


    if(aut.sp==TRUE){

     l.o <- fit$l; l.n <- 0 

     exl <- (length(start.v)-1):length(start.v)
     coef.v <- c(fit$argument[exl]) 
     coef.p <- c(fit$argument[-exl])

      if(l.sp1!=0 && l.sp2!=0){

       j.it <- 1; conv.sp <- TRUE 

	  while( abs(l.n-l.o)/abs(l.o) > pr.tol){   

                 fit  <- trust(ghss1, coef.p, rinit=rinit, rmax=rmax, dat=dat, exl1=coef.v[1], exl2=coef.v[2], dat1=fit$dat1, dat2=fit$dat2,
                     X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, iterlim=1e+4, 
                     fterm=fterm, mterm=mterm)

                 So <- spS(sp,gam1,gam2); coef.p <- fit$argument

		 wor.c    <- try(working.comp(fit,X1,X2,X1.d2,X2.d2,n))
                 if(class(wor.c)=="try-error") break
             
   
                 bs.mgfit <- try(magic(y=wor.c$rW.Z,X=wor.c$rW.X,sp=sp,S=qu.mag$Ss,
                                    off=qu.mag$off,rank=qu.mag$rank,n.score=2*n,
                                    gcv=FALSE,gamma=gamma,control=control))
                 if(class(bs.mgfit)=="try-error") {conv.sp <- FALSE; break} 

                 sp <- bs.mgfit$sp
                   
             S <- spS(sp,gam1,gam2)
             l.o <- fit$l

             fit  <- try(trust(ghss1, coef.p, rinit=rinit, rmax=rmax, dat=dat, exl1=coef.v[1], exl2=coef.v[2], dat1=fit$dat1, dat2=fit$dat2,
                           X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, iterlim=1e+4, 
                           fterm=fterm, mterm=mterm),silent=TRUE)

             if(class(fit)=="try-error"){ 
               fit  <- trust(ghss1, coef.p, rinit=rinit, rmax=rmax, dat=dat, exl1=coef.v[1], exl2=coef.v[2], dat1=fit$dat1, dat2=fit$dat2,
                           X1.d2=X1.d2, X2.d2=X2.d2, S=So, gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, iterlim=1e+4, 
                           fterm=fterm, mterm=mterm)
               conv.sp <- FALSE
               break
              }

             coef.p <- fit$argument
             l.n <- fit$l

             fit  <- trust(ghss2, coef.v, rinit=rinit, rmax=rmax, dat=dat, params.f=coef.p, dat1=fit$dat1, dat2=fit$dat2, 
                           X1.d2=X1.d2, X2.d2=X2.d2, gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, iterlim=1e+4, 
                           fterm=fterm, mterm=mterm)

             coef.v <- fit$argument


             if(j.it>iterlimSP){
              conv.sp <- FALSE
              break
             }
             j.it <- j.it + 1       
           
        }
      }
      
  up.start.v <- c(coef.p,coef.v) 
  fit  <- trust(ghss, up.start.v, rinit=rinit, rmax=rmax, dat=dat,
                X1.d2=X1.d2, X2.d2=X2.d2, S=S, gam1=gam1, gam2=gam2, fp=fp, blather=TRUE, 
                fterm=fterm, mterm=mterm, iterlim=1e+4)
   }


    He <- fit$hessian
    He.eig <- eigen(He)
    k.e <- sum(as.numeric(He.eig$val<sqrt(.Machine$double.eps)))
    
     if(k.e!=0){
      ind.e <- (length(He.eig$val)-(k.e-1)):length(He.eig$val)
      min.e <- min(He.eig$val[1:(ind.e[1]-1)])
      for(i in 1:k.e) He.eig$val[ind.e[i]] <- min.e/10^i  
      Vb <- He.eig$vec%*%diag(1/He.eig$val)%*%t(He.eig$vec)      
     }else{
      Vb <- He.eig$vec%*%diag(1/He.eig$val)%*%t(He.eig$vec) 
    }
             
  if(l.sp1!=0 && l.sp2!=0 && fp==FALSE){HeSh <- He - fit$S.h; F <- Vb%*%HeSh} else{HeSh <- He; F <- Vb%*%HeSh}      
  t.edf <- sum(diag(F))

L <- list(fit=fit, gam1=gam1, gam2=gam2, gam2.1=gam2.1, sp=sp, iter.sp=j.it, start.v=start.v,
          sigma=exp(fit$argument[length(fit$argument)-1]), rho=tanh(fit$argument[length(fit$argument)]), 
          n=n, n.sel=length(gam2$y), X1=X1, X2=X2, X1.d2=X1.d2, X2.d2=X2.d2, 
          l.sp1=l.sp1, l.sp2=l.sp2, He=He, HeSh=HeSh, Vb=Vb, F=F, 
          t.edf=t.edf, bs.mgfit=bs.mgfit, conv.sp=conv.sp, wor.c=wor.c, eta1=fit$eta1, eta2=fit$eta2, dat=dat)


class(L) <- "SemiParSampleSel"

L

}
