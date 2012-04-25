plot.SemiParSampleSel <- function(x, eq, select, rug=TRUE, se=TRUE, se.l=1.95996, seWithMean=FALSE, n=100,
                                 xlab = NULL, ylab=NULL, zlab=NULL, xlim=NULL, ylim = NULL, main=NULL, trans = I, n2 = 40, 
                                 theta = 30, phi = 30, too.far = 0.1, ...){
                                                                                         
  sub.edf <- function(lab, edf){ 
      pos <- regexpr(":", lab)[1]
      if(pos < 0){pos <- nchar(lab) - 1
                  lab <- paste(substr(lab, start = 1, stop = pos),",", round(edf, digits = 2), ")", sep = "")
      }
      else{lab1 <- substr(lab, start = 1, stop = pos - 2)
           lab2 <- substr(lab, start = pos - 1, stop = nchar(lab))
           lab <- paste(lab1, ",", round(edf, digits = 2), lab2,sep = "")
      }
  lab
  }

  if(x$l.sp1==0 && x$l.sp2==0) stop("The model is fully parametric; no smooth components to plot")

  if(eq==1) if(select>x$l.sp1) stop("No more smooth component to plot")
  if(eq==2) if(select>x$l.sp2) stop("No more smooth component to plot")
  
  if(eq==1) if(x$gam1$smooth[[select]]$dim>2) stop("No plotting for smooths of more than two variables")
  if(eq==2) if(x$gam2$smooth[[select]]$dim>2) stop("No plotting for smooths of more than two variables")

  if(eq==1) ind <- x$gam1$smooth[[select]]$first.para:x$gam1$smooth[[select]]$last.para
  if(eq==2) ind <- (x$gam2$smooth[[select]]$first.para:x$gam2$smooth[[select]]$last.para)+x$X1.d2

  est.par.c1 <- x$fit$argument
  est.par.c2 <- est.par.c1[-c(1:x$X1.d2)]
  Vb         <- x$Vb
  d.F        <- diag(x$F)

  if(eq==1){edf <- round(sum(d.F[ind]),2)
         if(x$gam1$smooth[[select]]$dim==1){raw <- x$gam1$model[x$gam1$smooth[[select]]$term] 
                                                 xx  <- seq(min(raw), max(raw), length=n) 
                                                      if(x$gam1$smooth[[select]]$by!= "NA"){by <- rep(1, n)
                                                                                                 d  <- data.frame(x = xx, by = by)
                                                                                                 names(d) <- c(x$gam1$smooth[[select]]$term,x$gam1$smooth[[select]]$by) 
                                                      }
                                                      else{d <- data.frame(x = xx)
                                                           names(d) <- x$gam1$smooth[[select]]$term
                                                      }
         }
         else if(x$gam1$smooth[[select]]$dim==2){xterm <- x$gam1$smooth[[select]]$term[1]
                                                      yterm <- x$gam1$smooth[[select]]$term[2]      
                 						      raw <- data.frame(x = as.numeric(x$gam1$model[xterm][[1]]),
                              			                        y = as.numeric(x$gam1$model[yterm][[1]]) )
               						      n2 <- max(10, n2)
               						      xm <- seq(min(raw$x), max(raw$x), length = n2)
               						  	ym <- seq(min(raw$y), max(raw$y), length = n2)
              						 	xx <- rep(xm, n2)
              						      yy <- rep(ym, rep(n2, n2))
                					       	   if(too.far > 0) exclude <- exclude.too.far(xx, yy, raw$x, raw$y,dist = too.far)
                                                	   else exclude <- rep(FALSE, n2 * n2)
                                                 	   if(x$gam1$smooth[[select]]$by != "NA"){by <- rep(1, n2^2)
                                                 	                                               d  <- data.frame(x = xx, y = yy, by = by)
                                                 	                                               names(d) <- c(xterm, yterm, x$gam1$smooth[[select]]$by)
                                                  	   }
                                                   	   else{d <- data.frame(x = xx, y = yy)
                                                    	        names(d) <- c(xterm, yterm)
                                                         } 
              }
  X <- PredictMat(x$gam1$smooth[[select]], d) 
  f <- X%*%est.par.c1[ind]
  }

  if(eq==2){edf <- round(sum(d.F[ind]),2)
         if(x$gam2$smooth[[select]]$dim==1){raw <- x$gam2$model[x$gam2$smooth[[select]]$term] 
                                                 xx  <- seq(min(raw), max(raw), length=n) 
                                                      if(x$gam2$smooth[[select]]$by!= "NA"){by <- rep(1, n)
                                                                                                 d  <- data.frame(x = xx, by = by)
                                                                                                 names(d) <- c(x$gam2$smooth[[select]]$term,x$gam2$smooth[[select]]$by) 
                                                      }
                                                      else{d <- data.frame(x = xx)
                                                           names(d) <- x$gam2$smooth[[select]]$term
                                                      }
         }
         else if(x$gam2$smooth[[select]]$dim==2){xterm <- x$gam2$smooth[[select]]$term[1]
                                                      yterm <- x$gam2$smooth[[select]]$term[2]      
                 						      raw <- data.frame(x = as.numeric(x$gam2$model[xterm][[1]]),
                              			                        y = as.numeric(x$gam2$model[yterm][[1]]) )
               						      n2 <- max(10, n2)
               						      xm <- seq(min(raw$x), max(raw$x), length = n2)
               						  	ym <- seq(min(raw$y), max(raw$y), length = n2)
              						 	xx <- rep(xm, n2)
              						      yy <- rep(ym, rep(n2, n2))
                					       	   if(too.far > 0) exclude <- exclude.too.far(xx, yy, raw$x, raw$y,dist = too.far)
                                                	   else exclude <- rep(FALSE, n2 * n2)
                                                 	   if(x$gam2$smooth[[select]]$by != "NA"){by <- rep(1, n2^2)
                                                 	                                               d  <- data.frame(x = xx, y = yy, by = by)
                                                 	                                               names(d) <- c(xterm, yterm, x$gam2$smooth[[select]]$by)
                                                  	   }
                                                   	   else{d <- data.frame(x = xx, y = yy)
                                                    	        names(d) <- c(xterm, yterm)
                                                         } 
              }
  X <- PredictMat(x$gam2$smooth[[select]], d)
  f <- X%*%est.par.c2[ind-x$X1.d2]
  }

  if(se){
      if(eq==1)
             if(seWithMean){X1 <- matrix(c(x$gam1$cmX,rep(0,length(est.par.c2))), nrow(X), ncol(Vb), byrow = TRUE)
                            X1[,ind] <- X
                            se.fit <- sqrt(rowSums((X1 %*% Vb) * X1))                           
             }
             else se.fit <- sqrt(rowSums((X %*% Vb[ind,ind]) * X))
         
      else   if(seWithMean){X1 <- matrix(c(rep(0,length(x$gam1$cmX)),x$gam2$cmX,0,0), nrow(X), ncol(Vb), byrow = TRUE)
                            X1[,ind] <- X
                            se.fit <- sqrt(rowSums((X1 %*% Vb) * X1))                           
             }
             else se.fit <- sqrt(rowSums((X %*% Vb[ind,ind]) * X))      
  ub <- (f+se.l*se.fit)
  lb <- (f-se.l*se.fit) 
  if(eq==1 && x$gam1$smooth[[select]]$dim==1) if(is.null(ylim)) ylim <- c(min(lb),max(ub))
  if(eq==2 && x$gam2$smooth[[select]]$dim==1) if(is.null(ylim)) ylim <- c(min(lb),max(ub))

  } 
              
  if(eq==1){
         if(x$gam1$smooth[[select]]$dim==1){
       	   if(is.null(xlab)) x.lab <- x$gam1$smooth[[select]]$term else x.lab <- xlab  
       	   if(is.null(ylab)) y.lab <- sub.edf(x$gam1$smooth[[select]]$label,edf) else y.lab <- ylab  

      	 plot(xx,f,type="l",xlab=x.lab,ylim=ylim,ylab=y.lab, main=main, ...)
      	   if(se){lines( xx , ub, lty=2 )
       	          lines( xx , lb, lty=2 )
       	   }
      	   if(rug) rug(as.numeric(x$gam1$model[x$gam1$smooth[[select]]$term][[1]]))
         }
         else{if(x$gam1$smooth[[select]]$dim==2){     if(is.null(zlab)) zlabel <- sub.edf(x$gam1$smooth[[select]]$label, edf) else zlabel <- zlab
                                                      if(is.null(xlab)) xlabel <- xterm else xlabel <- xlab
                                                      if(is.null(ylab)) ylabel <- yterm else ylabel <- ylab


                pd <- list(fit = f, dim = 2, xm = xm, ym = ym, ylab = ylabel, xlab = xlabel, zlab = zlabel, raw = raw)
                if (is.null(ylim)) pd$ylim <- range(ym) else pd$ylim <- ylim
                if (is.null(xlim)) pd$xlim <- range(xm) else pd$xlim <- xlim
                         persp(pd$xm, pd$ym, matrix(trans(pd$fit), n2, n2), xlab = pd$xlab, 
                               ylab = pd$ylab, zlab = pd$zlab, 
                               ylim = pd$ylim, xlim = pd$xlim,
                               theta = theta, phi = phi, main = main, ...)
                
              }
          }
  }
          
  if(eq==2){
         if(x$gam2$smooth[[select]]$dim==1){
       	   if(is.null(xlab)) x.lab <- x$gam2$smooth[[select]]$term else x.lab <- xlab  
               if(is.null(ylab)) y.lab <- sub.edf(x$gam2$smooth[[select]]$label,edf) else y.lab <- ylab  

      	 plot(xx,f,type="l",xlab=x.lab,ylim=ylim,ylab=y.lab, main=main, ...)
      	   if(se){lines( xx , ub, lty=2 )
       	          lines( xx , lb, lty=2 )
       	   }
      	   if(rug) rug(as.numeric(x$gam2$model[x$gam2$smooth[[select]]$term][[1]]))
         }
         else{if(x$gam2$smooth[[select]]$dim==2){     if(is.null(zlab)) zlabel <- sub.edf(x$gam2$smooth[[select]]$label, edf) else zlabel <- zlab
                                                      if(is.null(xlab)) xlabel <- xterm else xlabel <- xlab
                                                      if(is.null(ylab)) ylabel <- yterm else ylabel <- ylab

                pd <- list(fit = f, dim = 2, xm = xm, ym = ym, ylab = ylabel, xlab = xlabel, zlab = zlabel, raw = raw)
                if (is.null(ylim)) pd$ylim <- range(ym) else pd$ylim <- ylim
                if (is.null(xlim)) pd$xlim <- range(xm) else pd$xlim <- xlim
                         persp(pd$xm, pd$ym, matrix(trans(pd$fit), n2, n2), xlab = pd$xlab, 
                               ylab = pd$ylab, zlab = pd$zlab, 
                               ylim = pd$ylim, xlim = pd$xlim,
                               theta = theta, phi = phi, main = main, ...)
                
              }
          }
  }


}










