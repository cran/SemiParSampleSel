fit.SemiParSampleSel <- function( dat, X1, X2, qu.mag, sp, X1.d2, X2.d2, gp1, gp2, n, l.sp1, l.sp2, weights, iterlimSP, pr.tol,
                                  BivD, margins, gamma, aut.sp, fp, start.v, rinit, rmax, fterm, mterm, control ) {



  fit  <- trust(ghss, start.v, rinit=rinit, rmax=rmax, BivD=BivD, margins=margins, dat=dat, dat1=X1, dat2=X2, qu.mag=qu.mag, sp=sp, X1.d2=X1.d2, X2.d2=X2.d2,
                l.sp1=l.sp1, l.sp2=l.sp2, gp1=gp1, gp2=gp2, fp=fp, blather=TRUE, fterm=fterm, mterm=mterm, iterlim=1e+4, weights=weights)



  iter.if <- fit$iterations
  conv.sp <- iter.sp <- bs.mgfit <- wor.c <- iter.fi <- NULL


  if(aut.sp==TRUE && fp==FALSE && (l.sp1!=0 || l.sp2!=0)){

      l.o <- fit$l; l.n <- 0 
      exl <- (length(start.v)-1):length(start.v)
      coef.v <- c(fit$argument[exl]) 
      coef.p <- c(fit$argument[-exl])

      conv.sp <- TRUE; iter.sp <- 0  

      myf <- function(x) list( rbind( c(x[1],x[2]), 
                                      c(x[2],x[3])  ) )

      while( abs(l.n-l.o)/abs(l.o) > pr.tol){   

                 fit  <- trust(ghss1, coef.p, rinit=rinit, rmax=rmax, BivD=BivD, margins=margins, dat=dat, exl1=coef.v[1], exl2=coef.v[2], dat1=X1, dat2=X2, qu.mag=qu.mag,
                      sp=sp, l.sp1=l.sp1, l.sp2=l.sp2, X1.d2=X1.d2, X2.d2=X2.d2, gp1=gp1, gp2=gp2, blather=TRUE, iterlim=1e+4, fterm=fterm, 
                      mterm=mterm, weights=weights)
	  
                 spo <- sp; coef.p <- fit$argument

		 wor.c    <- try(working.comp(fit,X1,X2,X1.d2,X2.d2,n,myf))
                 if(class(wor.c)=="try-error") break
             
   
                 bs.mgfit <- try(magic(y=wor.c$rW.Z, X=wor.c$rW.X, sp=sp, S=qu.mag$Ss, off=qu.mag$off,
                                       rank=qu.mag$rank, gcv=FALSE, gamma=gamma, control=control))
                 if(class(bs.mgfit)=="try-error") { conv.sp <- FALSE; break } 

                 sp <- bs.mgfit$sp;  iter.sp <- iter.sp + 1 
                   
                 l.o <- fit$l

                 fit  <- try(trust(ghss1, coef.p, rinit=rinit, rmax=rmax, BivD=BivD, margins=margins, dat=dat, exl1=coef.v[1], exl2=coef.v[2], dat1=X1, dat2=X2,
                           qu.mag=qu.mag, sp=sp, l.sp1=l.sp1, l.sp2=l.sp2, X1.d2=X1.d2, X2.d2=X2.d2, gp1=gp1, gp2=gp2, blather=TRUE, iterlim=1e+4, 
                           fterm=fterm, mterm=mterm, weights=weights), silent=TRUE)

             	 if(class(fit)=="try-error"){ 
               	           fit  <- trust(ghss1, coef.p, rinit=rinit, rmax=rmax, BivD=BivD, margins=margins, dat=dat, exl1=coef.v[1], exl2=coef.v[2], dat1=X1, dat2=X2, qu.mag=qu.mag,
                	                 sp=spo, l.sp1=l.sp1, l.sp2=l.sp2, X1.d2=X1.d2, X2.d2=X2.d2, gp1=gp1, gp2=gp2, blather=TRUE, iterlim=1e+4, 
                        	         fterm=fterm, mterm=mterm, weights=weights)
               	           conv.sp <- FALSE
               	           break
              	 }

                 coef.p <- fit$argument
                 l.n <- fit$l

                 fit  <- trust(ghss2, coef.v, rinit=rinit, rmax=rmax, BivD=BivD, margins=margins, dat=dat, params.f=coef.p, dat1=X1, dat2=X2, qu.mag=qu.mag, sp=sp, 
                               X1.d2=X1.d2, X2.d2=X2.d2, blather=TRUE, iterlim=1e+4, fterm=fterm, mterm=mterm, weights=weights)

                 coef.v <- fit$argument

                 if(iter.sp>iterlimSP){ conv.sp <- FALSE; break}     
           
        }

        up.start.v <- c(coef.p,coef.v) 
  
        fit  <- trust(ghss, up.start.v, rinit=rinit, rmax=rmax, BivD=BivD, margins=margins, dat=dat, dat1=X1, dat2=X2, qu.mag=qu.mag, sp=sp, X1.d2=X1.d2, X2.d2=X2.d2, 
                      l.sp1=l.sp1, l.sp2=l.sp2, gp1=gp1, gp2=gp2, fp=fp, blather=TRUE, fterm=fterm, mterm=mterm, iterlim=1e+4, weights=weights)
        iter.fi <- fit$iterations

  }

  list(fit=fit, bs.mgfit=bs.mgfit, conv.sp=conv.sp, wor.c=wor.c, iter.sp=iter.sp, iter.if=iter.if, iter.fi=iter.fi, sp=sp )

}




