fit.SemiParSampleSel <- function(funcop, start.v, dat, VC, qu.mag, sp, iterlimsp, pr.tolsp, rinit, rmax) {

  fit <- trust(funcop, start.v, dat=dat, VC=VC, qu.mag=qu.mag, sp=sp, rinit=rinit, rmax=rmax, blather=TRUE, iterlim=1e+4, parscale=VC$parscale,)

  iter.if <- fit$iterations

  conv.sp <- iter.sp <- iter.inner <- bs.mgfit <- wor.c <- magpp <- NULL


  if(VC$fp==FALSE && (VC$l.sp1!=0 || VC$l.sp2!=0)){

      myf <- function(x) list( rbind( c(x[1],x[2]), 
                                      c(x[2],x[3])  ) )
                                      
      stoprule.SP <- 1; conv.sp <- TRUE; iter.inner <- iter.sp <- 0                                  
                                      

      while( stoprule.SP > pr.tolsp ){   
      
                 coefo <- fit$argument; spo <- sp 

		 wor.c <- try(working.comp(fit, VC, myf)); if(class(wor.c)=="try-error") break
             
                 bs.mgfit <- try(magic(y=wor.c$rW.Z, X=wor.c$rW.X, sp=sp, S=qu.mag$Ss, off=qu.mag$off,
                                       rank=qu.mag$rank, gcv=FALSE, gamma=VC$gamma))
                 if(class(bs.mgfit)=="try-error") { conv.sp <- FALSE; break } 

                 sp <- bs.mgfit$sp; iter.sp <- iter.sp + 1; names(sp) <- names(spo) 
                   
                 o.ests <- c(fit$argument)  

                 fit  <- try(trust(funcop, o.ests, dat=dat, VC=VC, qu.mag=qu.mag, parscale=VC$parscale, 
                                   sp=sp, rinit=rinit, rmax=rmax, blather=TRUE, iterlim=1e+4), silent=TRUE)

             	 if(class(fit)=="try-error"){ 
               	           fit  <- trust(funcop, coefo, dat=dat, VC=VC, qu.mag=qu.mag, parscale=VC$parscale,
                                   sp=sp, rinit=rinit, rmax=rmax, blather=TRUE, iterlim=1e+4)
               	           conv.sp <- FALSE
               	           break
              	 }

                 iter.inner <- iter.inner + fit$iterations  
                 
                 
             n.ests <- c(fit$argument)

             if(iter.sp > iterlimsp){conv.sp <- FALSE; break }

             stoprule.SP <- max(abs(o.ests-n.ests))   
           
        }
        
  magpp <- magic.post.proc(wor.c$rW.X, bs.mgfit)      

  }

                  list(fit = fit, 
                       iter.if = iter.if, 
                       conv.sp = conv.sp, 
                       iter.sp = iter.sp, 
                       iter.inner = iter.inner, 
                       bs.mgfit = bs.mgfit, 
                       wor.c = wor.c, 
                       sp = sp, magpp = magpp )

}




