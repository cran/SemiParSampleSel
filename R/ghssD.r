ghssD <- function(params, dat, VC, qu.mag=NULL, sp=NULL){

  eta1 <- VC$X1%*%params[1:VC$X1.d2]  
  eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]   
  
  if(VC$margins[2]=="P")                 theta.star <- params[(VC$X1.d2+VC$X2.d2+1)]
  if(VC$margins[2] %in% c("NB","PIG") )  theta.star <- params[(VC$X1.d2+VC$X2.d2+2)]
  if(VC$margins[2] %in% c("D","S") )     theta.star <- params[(VC$X1.d2+VC$X2.d2+3)]
    
  eps <- sqrt(.Machine$double.eps)


  i1 <- dat[,1]
  i0 <- 1-i1
  ind <- i1==0 

  F1 <- pnorm(-eta1)  			
  F1 <- ifelse(F1>0.00000001,F1,0.00000001)
  F1 <- ifelse(F1<0.99999999,F1,0.99999999)
  
  precision<-10^(-8)
  
#-------------------------------------------------------
  dF1 <- -dnorm(-eta1)
  d2F1ddelta1delta1 <- -dF1*eta1
  
  
  fd.prec<-10^(-7) 
  
  
 if(VC$margins[2]=="P") {
    i2 <- dat[,2]
    F2 <- ppois(i2, exp(eta2))
    F2 <- ifelse(F2<(1-precision), F2, 1-precision)
    F2 <- ifelse(F2>precision, F2, precision)    
    f2 <- dpois(i2, exp(eta2))
    f2 <- ifelse(f2>precision, f2, precision)
    F22 <- F2 - f2
    F22 <- ifelse(F22<(1-precision), F22, 1-precision)
    F22 <- ifelse(F22>precision, F22, precision)
    
    # df2, dF2, dF22  
    gamma <- exp(-exp(eta2))*exp(eta2)^(i2)/(dpois(i2, exp(eta2)))
    df2 <- as.vector(1/gamma*(-exp(-exp(eta2))*exp(eta2)^(i2)+exp(-exp(eta2))*i2*exp(eta2)^(i2-1))*exp(eta2))
    
    
    df2.int <- function(y, mu) { 
      gamma <- exp(-mu)*mu^(y)/(dpois(y, mu))
      as.vector(1/gamma*(-exp(-mu)*mu^(y)+exp(-mu)*y*mu^(y-1))*mu) 
    }
    
    
    ly <- length(i2)
    cdf.f <- rep(0, ly)
    nmu <- rep(exp(eta2), length = ly) 
    
    for (i in 1:ly) {
      
      y.y <- i2[i]
      mm <- nmu[i]
      allval <- seq(0, y.y)
      pdfall <- df2.int(allval, mu = mm)
      cdf.f[i] <- sum(pdfall)
      
    }  
    
    dF2 <- cdf.f
    dF22 <- dF2-df2
    
    
    
    # Hessian derivative components
    
    gamma.plus <- exp(-exp(eta2+fd.prec))*exp(eta2+fd.prec)^(i2)/(dpois(i2, exp(eta2+fd.prec)))
    df2.plus <- as.vector(1/gamma.plus*(-exp(-exp(eta2+fd.prec))*exp(eta2+fd.prec)^(i2)+exp(-exp(eta2+fd.prec))*i2*exp(eta2+fd.prec)^(i2-1))*exp(eta2+fd.prec))
    d2f2delta22 <- (df2.plus - df2)/fd.prec
    
    
    ly <- length(i2)
    cdf.f.plus <- rep(0, ly)
    nmu.plus <- rep(exp(eta2 + fd.prec), length = ly) 
    
    for (i in 1:ly) {
      
      y.y <- i2[i]
      mm <- nmu.plus[i]
      allval <- seq(0, y.y)
      pdfall <- df2.int(allval, mu = mm)
      cdf.f.plus[i] <- sum(pdfall)
      
    }  
    
    dF2.plus <- cdf.f.plus
    d2F2ddelta22 <- (dF2.plus - dF2) / fd.prec
    
    
    d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    
    
    
    
    
    
    
  } else { if (VC$margins[2]=="NB") {
    i2 <- dat[,2]
    sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
    sigma <- exp(sigma.star)
    sigma <- ifelse(sigma>10^-8, sigma, 10^-8)
    sigma <- ifelse(sigma<10^7, sigma, 10^7)
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^7, mu, 10^7)
    F2 <- pNBI(i2, sigma = sigma, mu = mu)
    F2 <- ifelse(F2<(1-precision), F2, 1-precision)
    F2 <- ifelse(F2>precision, F2, precision)
    f2 <- dNBI(i2, sigma = sigma,  mu = mu)
    f2 <- ifelse(f2>precision, f2, precision)
    F22 <- F2 - f2
    F22 <- ifelse(F22<(1-precision), F22, 1-precision)
    F22 <- ifelse(F22>precision, F22, precision)
    
    
    
    
    df2 <- f2*(i2*(mu*sigma)^(-1)*sigma*mu-(i2+1/sigma)*(mu*sigma+1)^(-1)*sigma*mu)
    df2 <- as.vector(df2)
    
    df2.sigma <- f2*(digamma(i2+1/sigma)*(-1/sigma^2)+i2*(mu*sigma)^(-1)*mu-(digamma(1/sigma)*(-1/sigma^2)+(-1/sigma^2)*log(mu*sigma+1)+(i2+1/sigma)*(1/(mu*sigma+1))*mu))
    
    
    df2.mu.int<-function(y, mu, sigma, nu) {
      f2 <- dNBI(y, sigma = sigma,  mu = mu)
      f2 <- ifelse(f2>precision, f2, precision)
      df2 <- f2*(y*(mu*sigma)^(-1)*sigma*mu-(y+1/sigma)*(mu*sigma+1)^(-1)*sigma*mu) 
      df2
    }
    
    df2.sigma.int<-function(y, mu, sigma, nu) {
      f2 <- dNBI(y, sigma = sigma,  mu = mu)
      f2 <- ifelse(f2>precision, f2, precision)
      df2 <- f2*(digamma(y+1/sigma)*(-1/sigma^2)+y*(mu*sigma)^(-1)*mu-(digamma(1/sigma)*(-1/sigma^2)+(-1/sigma^2)*log(mu*sigma+1)+(y+1/sigma)*(1/(mu*sigma+1))*mu)) 
      df2
    }
    
    
    ly <- length(i2)
    cdf.f <- rep(0, ly)
    nmu <- rep(mu, length = ly) 
    nsigma <- rep(sigma, length = ly)
    
    for (i in 1:ly) {
      y.y <- i2[i]
      mm <- nmu[i]
      ms <- nsigma[i]
      allval <- seq(0, y.y)
      pdfall <- df2.mu.int(allval, mm, ms)
      cdf.f[i] <- sum(pdfall)
    }  
    
    dF2 <- cdf.f
    dF22 <- dF2 - df2
    
    
    ly <- length(i2)
    cdf.f.sigma <- rep(0, ly)
    nmu <- rep(mu, length = ly) 
    nsigma <- rep(sigma, length = ly)
    
    for (i in 1:ly) {
      y.y <- i2[i]
      mm <- nmu[i]
      ms <- nsigma[i]
      allval <- seq(0, y.y)
      pdfall <- df2.sigma.int(allval, mm, ms)
      cdf.f.sigma[i] <- sum(pdfall)
    } 
    
    dF2.sigma <- cdf.f.sigma
    dF22.sigma <- dF2.sigma - df2.sigma
    
    
    
    # Hessian derivative components
    
    
    mu.plus <- ifelse(exp(eta2+fd.prec)>0.0001, exp(eta2+fd.prec), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    
    f2.plus <- dNBI(i2, sigma = sigma,  mu = mu.plus)
    f2.plus <- ifelse(f2.plus>precision, f2.plus, precision)
    df2.plus <- f2.plus*(i2*(mu.plus*sigma)^(-1)*sigma*mu.plus-(i2+1/sigma)*(mu.plus*sigma+1)^(-1)*sigma*mu.plus)
    df2.plus <- as.vector(df2.plus)
    d2f2delta22 <- (df2.plus - df2) / fd.prec
    
    
    sigma.plus <- sigma+fd.prec
    f2.plus.sigma <- dNBI(i2, sigma = sigma.plus,  mu = mu)
    f2.plus.sigma <- ifelse(f2.plus.sigma>precision, f2.plus.sigma, precision)
    df2.sigma.plus <- f2.plus.sigma*(digamma(i2+1/sigma.plus)*(-1/sigma.plus^2)+i2*(mu*sigma.plus)^(-1)*mu-(digamma(1/sigma.plus)*(-1/sigma.plus^2)+(-1/sigma.plus^2)*log(mu*sigma.plus+1)+(i2+1/sigma.plus)*(1/(mu*sigma.plus+1))*mu))
    d2f2sigma2 <- (df2.sigma.plus - df2.sigma) / fd.prec
    
    
    df2.dsigma.plus <- f2.plus.sigma*(i2*(mu*sigma.plus)^(-1)*sigma.plus*mu-(i2+1/sigma.plus)*(mu*sigma.plus+1)^(-1)*sigma.plus*mu)
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec
    
    
    
    ly <- length(i2)
    cdf.f.plus <- rep(0, ly)
    nmu <- rep(mu.plus, length = ly) 
    nsigma <- rep(sigma, length = ly)
    
    for (i in 1:ly) {
      y.y <- i2[i]
      mm <- nmu[i]
      ms <- nsigma[i]
      allval <- seq(0, y.y)
      pdfall <- df2.mu.int(allval, mm, ms)
      cdf.f.plus[i] <- sum(pdfall)
    }  
    
    dF2.mu.plus <- cdf.f.plus
    d2F2ddelta22 <- (dF2.mu.plus - dF2) / fd.prec
    d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    
    ly <- length(i2)
    cdf.f.sigma.plus <- rep(0, ly)
    nmu <- rep(mu, length = ly) 
    nsigma <- rep(sigma.plus, length = ly)
    
    for (i in 1:ly) {
      y.y <- i2[i]
      mm <- nmu[i]
      ms <- nsigma[i]
      allval <- seq(0, y.y)
      pdfall <- df2.sigma.int(allval, mm, ms)
      cdf.f.sigma.plus[i] <- sum(pdfall)
    } 
    
    dF2.sigma.plus <- cdf.f.sigma.plus
    d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec
    d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
    
    ly <- length(i2)
    cdf.f.dsigma.plus <- rep(0, ly)
    nmu <- rep(mu, length = ly) 
    nsigma <- rep(sigma.plus, length = ly)
    
    for (i in 1:ly) {
      y.y <- i2[i]
      mm <- nmu[i]
      ms <- nsigma[i]
      allval <- seq(0, y.y)
      pdfall <- df2.mu.int(allval, mm, ms)
      cdf.f.dsigma.plus[i] <- sum(pdfall)
    } 
    
    dF2.dsigma.plus <- cdf.f.dsigma.plus
    d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec
    d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma   
    
    
    
    
  } else { if (VC$margins[2]=="D") {
    i2 <- dat[,2]
    sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
    nu.star <- params[(VC$X1.d2+VC$X2.d2+2)]
    nu <- plogis(nu.star)
    nu <- ifelse(nu<precision, precision, nu)
    nu <- ifelse(nu>(1-precision), 1-precision, nu)
    sigma <- exp(sigma.star)
    sigma <- ifelse(sigma>precision, sigma, precision)
    sigma <- ifelse(sigma<10^7, sigma, 10^7)
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^7, mu, 10^7)
    F2 <- pDEL(i2, mu=mu, sigma=sigma, nu=nu)
    F2 <- ifelse(F2<(1-precision), F2, 1-precision)
    F2 <- ifelse(F2>precision, F2, precision)  
    f2 <- dDEL(i2, mu=mu, sigma=sigma, nu=nu)
    f2 <- ifelse(f2>precision, f2, precision)
    F22 <- F2 - f2
    F22 <- ifelse(F22<(1-precision), F22, 1-precision)
    F22 <- ifelse(F22>precision, F22, precision)
    
    
    tofyDELPORT2 <- function (y, mu, sigma, nu, what=1) 
   {
         ly <- max(length(y),length(mu),length(sigma),length(nu)) 
          y <- rep(y, length = ly)    
      sigma <- rep(sigma, length = ly)
         mu <- rep(mu, length = ly)   
         nu <- rep(nu, length = ly) 
         ty <- sumlty <- rep(0,length(y)) 
    for (i in 1:length(y))
    {
         iy <- dum <- 0
         iy <- y[i]+1 
      tofyn <- rep(0,iy)
       tofyn[1] <- mu[i]*nu[i]+mu[i]*(1-nu[i])/(1+mu[i]*sigma[i]*(1-nu[i])) 
        if (iy<2)
        {
          ty[i] <- tofyn[iy] 
      sumlty[i] <- 0
        }
        else
        {
        for (j in 2:iy)
         {
           dum = (1+(1/(mu[i]*sigma[i]*(1-nu[i]))))
             tofyn[j] <- ((j-1)+mu[i]*nu[i]+(1/(sigma[i]*(1-nu[i])))-(mu[i]*nu[i]*(j-1))/ tofyn[j-1])/dum      
         }
              
            ty[i] <- tofyn[iy] 
            sumlty[i] <- sum(log(tofyn[1:iy-1]))
       }     
    }
   result <- sumlty
   result        
   }  
    
    
    
    S_prime_mu <- (exp(tofyDELPORT2(i2, mu+fd.prec, sigma, nu))-exp(tofyDELPORT2(i2, mu, sigma, nu)))/fd.prec
    S_prime_eta2 <- mu*S_prime_mu
    S <- exp(tofyDELPORT2(i2, mu, sigma, nu))
    df2 <- as.vector(f2*(-mu*nu+(-1/sigma)*(1+mu*sigma*(1-nu))^(-1)*mu*sigma*(1-nu)+S_prime_eta2/S))
    
    S_prime_sigma <- (exp(tofyDELPORT2(i2, mu, sigma+fd.prec, nu))-exp(tofyDELPORT2(i2, mu, sigma, nu)))/fd.prec
    d_prime <- ((1/sigma^2)*log(1+mu*sigma*(1-nu))-(1/sigma)*(1/(1+mu*sigma*(1-nu)))*mu*(1-nu))*(1+mu*sigma*(1-nu))^(-1/sigma)
    df2.sigma <- f2*(d_prime*(1+mu*sigma*(1-nu))^(1/sigma)+S_prime_sigma/S)
    
    S_prime_nu <- (exp(tofyDELPORT2(i2, mu, sigma, nu+fd.prec))-exp(tofyDELPORT2(i2, mu, sigma, nu)))/fd.prec
    df2.nu <- f2*(-mu+((-1/sigma)*(1+mu*sigma*(1-nu))^(-1)*(-mu*sigma)+S_prime_nu/S))
    
    
    
    
    df2.mu.int <- function(y, mu, sigma, nu) {
      f2 <- dDEL(y, mu=mu, sigma=sigma, nu=nu)
      f2 <- ifelse(f2>precision, f2, precision)
      
      S_prime_mu <- (exp(tofyDELPORT2(y, mu+fd.prec, sigma, nu))-exp(tofyDELPORT2(y, mu, sigma, nu)))/fd.prec
      S_prime_eta2 <- mu*S_prime_mu
      S <- exp(tofyDELPORT2(y, mu, sigma, nu))
      
      df2 <- f2*(-mu*nu+(-1/sigma)*(1+mu*sigma*(1-nu))^(-1)*mu*sigma*(1-nu)+S_prime_eta2/S)
      df2
      
    }
    
    df2.sigma.int<-function(y, mu, sigma, nu) {
      f2 <- dDEL(y, mu=mu, sigma=sigma, nu=nu)
      f2 <- ifelse(f2>precision, f2, precision)
      
      S_prime_sigma <- (exp(tofyDELPORT2(y, mu, sigma+fd.prec, nu))-exp(tofyDELPORT2(y, mu, sigma, nu)))/fd.prec
      d_prime <- ((1/sigma^2)*log(1+mu*sigma*(1-nu))-(1/sigma)*(1/(1+mu*sigma*(1-nu)))*mu*(1-nu))*(1+mu*sigma*(1-nu))^(-1/sigma)
      S <- exp(tofyDELPORT2(y, mu, sigma, nu))
      
      df2 <- f2*(d_prime*(1+mu*sigma*(1-nu))^(1/sigma)+S_prime_sigma/S)
      df2
      
    }
    
    df2.nu.int <- function(y, mu, sigma, nu) {
      f2 <- dDEL(y, mu=mu, sigma=sigma, nu=nu)
      f2 <- ifelse(f2>precision, f2, precision)
      
      S_prime_nu <- (exp(tofyDELPORT2(y, mu, sigma, nu+fd.prec))-exp(tofyDELPORT2(y, mu, sigma, nu)))/fd.prec
      S <- exp(tofyDELPORT2(y, mu, sigma, nu))
      
      df2 <- f2*(-mu+((-1/sigma)*(1+mu*sigma*(1-nu))^(-1)*(-mu*sigma)+S_prime_nu/S))
      df2
      
    }
    
    
    ly <- length(i2)
    cdf.f <- rep(0, ly)
    nmu <- rep(mu, length = ly) 
    nsigma <- rep(sigma, length = ly)
    nnu <- rep(nu, length = ly)
    
    for (i in 1:ly) {
      y.y <- i2[i]
      mm <- nmu[i]
      ms <- nsigma[i]
      mn<- nnu[i]
      allval <- seq(0, y.y)
      pdfall <- df2.mu.int(allval, mm, ms, mn)
      cdf.f[i] <- sum(pdfall)
    }  
    
    dF2<-cdf.f
    dF22<-dF2-df2
    
    
    ly <- length(i2)
    cdf.f.sigma <- rep(0, ly)
    nmu <- rep(mu, length = ly) 
    nsigma <- rep(sigma, length = ly)
    nnu <- rep(nu, length = ly)
    
    for (i in 1:ly) {
      y.y <- i2[i]
      mm <- nmu[i]
      ms <- nsigma[i]
      mn<- nnu[i]
      allval <- seq(0, y.y)
      pdfall <- df2.sigma.int(allval, mm, ms, mn)
      cdf.f.sigma[i] <- sum(pdfall)
    } 
    
    dF2.sigma <- cdf.f.sigma
    dF22.sigma <- dF2.sigma-df2.sigma
    
    
    ly <- length(i2)
    cdf.f.nu <- rep(0, ly)
    nmu <- rep(mu, length = ly) 
    nsigma <- rep(sigma, length = ly)
    nnu <- rep(nu, length = ly)
    
    for (i in 1:ly) {
      y.y <- i2[i]
      mm <- nmu[i]
      ms <- nsigma[i]
      mn <- nnu[i]
      allval <- seq(0, y.y)
      pdfall <- df2.nu.int(allval, mm, ms, mn)
      cdf.f.nu[i] <- sum(pdfall)
    } 
    
    dF2.nu <- cdf.f.nu
    dF22.nu <- dF2.nu-df2.nu
    
    
    
    
    # Hessian derivative components
    
    fd.prec2 <- 10^-3
    
    
    #delta2delta2
    mu.plus <- ifelse(exp(eta2+fd.prec2)>0.0001, exp(eta2+fd.prec2), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    f2.plus <- dDEL(i2, sigma = sigma,  mu = mu.plus, nu=nu)
    f2.plus <- ifelse(f2.plus>precision, f2.plus, precision)
    S_prime_mu.plus <- (exp(tofyDELPORT2(i2, mu.plus+fd.prec, sigma, nu))-exp(tofyDELPORT2(i2, mu.plus, sigma, nu)))/fd.prec
    S_prime_eta2.plus <- mu.plus*S_prime_mu.plus
    S.plus <- exp(tofyDELPORT2(i2, mu.plus, sigma, nu))
    df2.plus <- as.vector(f2.plus*(-mu.plus*nu+(-1/sigma)*(1+mu.plus*sigma*(1-nu))^(-1)*mu.plus*sigma*(1-nu)+S_prime_eta2.plus/S.plus))
    df2.plus <- as.vector(df2.plus)
    d2f2delta22 <- (df2.plus - df2) / fd.prec2
    
    #sigma^2
    sigma.plus <- sigma + fd.prec2
    f2.plus.sigma <- dDEL(i2, sigma = sigma.plus,  mu = mu, nu = nu)
    f2.plus.sigma <- ifelse(f2.plus.sigma>precision, f2.plus.sigma, precision)
    S_prime_sigma <- (exp(tofyDELPORT2(i2, mu, sigma.plus+fd.prec, nu))-exp(tofyDELPORT2(i2, mu, sigma.plus, nu)))/fd.prec
    d_prime <- ((1/sigma.plus^2)*log(1+mu*sigma.plus*(1-nu))-(1/sigma.plus)*(1/(1+mu*sigma.plus*(1-nu)))*mu*(1-nu))*(1+mu*sigma.plus*(1-nu))^(-1/sigma.plus)
    S.plus.sigma <- exp(tofyDELPORT2(i2, mu, sigma.plus, nu))
    df2.sigma.plus <- f2.plus.sigma*(d_prime*(1+mu*sigma.plus*(1-nu))^(1/sigma.plus)+S_prime_sigma/S.plus.sigma)
    d2f2sigma2 <- (df2.sigma.plus - df2.sigma) / fd.prec2
    
    #delta2sigma
    S_prime_mu.dsigma <- (exp(tofyDELPORT2(i2, mu+fd.prec, sigma.plus, nu))-exp(tofyDELPORT2(i2, mu, sigma.plus, nu)))/fd.prec
    S_prime_eta2.dsigma <- mu*S_prime_mu.dsigma
    df2.dsigma.plus <- as.vector(f2.plus.sigma*(-mu*nu+(-1/sigma.plus)*(1+mu*sigma.plus*(1-nu))^(-1)*mu*sigma.plus*(1-nu)+S_prime_eta2.dsigma/S.plus.sigma))
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec2
    
    
    #nu^2
    nu.plus <- nu + fd.prec2
    f2.plus.nu <- dDEL(i2, sigma = sigma,  mu = mu, nu = nu.plus)
    f2.plus.nu <- ifelse(f2.plus.nu>precision, f2.plus.nu, precision)
    S_prime_nu <- (exp(tofyDELPORT2(i2, mu, sigma, nu.plus+fd.prec))-exp(tofyDELPORT2(i2, mu, sigma, nu.plus)))/fd.prec
    S.plus.nu <- exp(tofyDELPORT2(i2, mu, sigma, nu.plus))
    df2.nu.plus <- f2.plus.nu*(-mu+((-1/sigma)*(1+mu*sigma*(1-nu.plus))^(-1)*(-mu*sigma)+S_prime_nu/S.plus.nu))
    d2f2nu2 <- (df2.nu.plus - df2.nu) / fd.prec2
    
    
    #delta2nu
    S_prime_mu.dnu <- (exp(tofyDELPORT2(i2, mu+fd.prec, sigma, nu.plus))-exp(tofyDELPORT2(i2, mu, sigma, nu.plus)))/fd.prec
    S_prime_eta2.dnu <- mu*S_prime_mu.dnu
    S.dnu <- exp(tofyDELPORT2(i2, mu, sigma, nu.plus))
    df2.dnu.plus <- as.vector(f2.plus.nu*(-mu*nu.plus+(-1/sigma)*(1+mu*sigma*(1-nu.plus))^(-1)*mu*sigma*(1-nu.plus)+S_prime_eta2.dnu/S.dnu))
    df2.dnu.plus <- as.vector(df2.dnu.plus)
    d2f2delta2nu <- (df2.dnu.plus - df2) / fd.prec2
    
    
    #nusigma
    f2.plus.nu.sigma <- dDEL(i2, sigma = sigma.plus,  mu = mu, nu = nu)
    f2.plus.nu.sigma <- ifelse(f2.plus.nu.sigma>precision, f2.plus.nu.sigma, precision)
    S_prime_nu.plus <- (exp(tofyDELPORT2(i2, mu, sigma.plus, nu+fd.prec))-exp(tofyDELPORT2(i2, mu, sigma.plus, nu)))/fd.prec
    df2.nu.sigma.plus <- f2.plus.nu.sigma*(-mu+((-1/sigma.plus)*(1+mu*sigma.plus*(1-nu))^(-1)*(-mu*sigma.plus)+S_prime_nu.plus/S.plus.sigma))
    d2f2nusigma <- (df2.nu.sigma.plus - df2.nu) / fd.prec2
    
    
    
    ly <- length(i2)
    cdf.f.plus <- rep(0, ly)
    nmu <- rep(mu.plus, length = ly) 
    nsigma <- rep(sigma, length = ly)
    nnu <- rep(nu, length = ly)
    
    for (i in 1:ly) {
      y.y <- i2[i]
      mm <- nmu[i]
      ms <- nsigma[i]
      mn<- nnu[i]
      allval <- seq(0, y.y)
      pdfall <- df2.mu.int(allval, mm, ms, mn)
      cdf.f.plus[i] <- sum(pdfall)
    }  
    
    dF2.mu.plus <- cdf.f.plus
    d2F2ddelta22 <- (dF2.mu.plus - dF2) / fd.prec2
    d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    
    ly <- length(i2)
    cdf.f.sigma.plus <- rep(0, ly)
    nmu <- rep(mu, length = ly) 
    nsigma <- rep(sigma.plus, length = ly)
    nnu <- rep(nu, length = ly)
    
    for (i in 1:ly) {
      y.y <- i2[i]
      mm <- nmu[i]
      ms <- nsigma[i]
      mn<- nnu[i]
      allval <- seq(0, y.y)
      pdfall <- df2.sigma.int(allval, mm, ms, mn)
      cdf.f.sigma.plus[i] <- sum(pdfall)
    }  
    
    dF2.sigma.plus <- cdf.f.sigma.plus
    d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec2
    d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
    
    ly <- length(i2)
    cdf.f.dsigma.plus <- rep(0, ly)
    nmu <- rep(mu, length = ly) 
    nsigma <- rep(sigma.plus, length = ly)
    nnu <- rep(nu, length = ly)
    
    for (i in 1:ly) {
      y.y <- i2[i]
      mm <- nmu[i]
      ms <- nsigma[i]
      mn<- nnu[i]
      allval <- seq(0, y.y)
      pdfall <- df2.mu.int(allval, mm, ms, mn)
      cdf.f.dsigma.plus[i] <- sum(pdfall)
    }   
    
    dF2.dsigma.plus <- cdf.f.dsigma.plus
    d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec2
    d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma   
    
    
    ly <- length(i2)
    cdf.f.nu.plus <- rep(0, ly)
    nmu <- rep(mu, length = ly) 
    nsigma <- rep(sigma, length = ly)
    nnu <- rep(nu.plus, length = ly)
    
    for (i in 1:ly) {
      y.y <- i2[i]
      mm <- nmu[i]
      ms <- nsigma[i]
      mn<- nnu[i]
      allval <- seq(0, y.y)
      pdfall <- df2.nu.int(allval, mm, ms, mn)
      cdf.f.nu.plus[i] <- sum(pdfall)
    } 
    
    dF2.nu.plus <- cdf.f.nu.plus
    d2F2dnu2 <- (dF2.nu.plus - dF2.nu) / fd.prec2
    d2F22dnu2 <- d2F2dnu2 - d2f2nu2
    
    
    ly <- length(i2)
    cdf.f.dnu.plus <- rep(0, ly)
    nmu <- rep(mu, length = ly) 
    nsigma <- rep(sigma, length = ly)
    nnu <- rep(nu.plus, length = ly)
    
    for (i in 1:ly) {
      y.y <- i2[i]
      mm <- nmu[i]
      ms <- nsigma[i]
      mn<- nnu[i]
      allval <- seq(0, y.y)
      pdfall <- df2.mu.int(allval, mm, ms, mn)
      cdf.f.dnu.plus[i] <- sum(pdfall)
    }  
    
    dF2.dnu.plus <- cdf.f.dnu.plus
    d2F2ddelta2dnu <- (dF2.dnu.plus - dF2) / fd.prec2
    d2F22ddelta2dnu <- d2F2ddelta2dnu - d2f2delta2nu
    
    
    ly <- length(i2)
    cdf.f.dnu.dsigma.plus <- rep(0, ly)
    nmu <- rep(mu, length = ly) 
    nsigma <- rep(sigma.plus, length = ly)
    nnu <- rep(nu, length = ly)
    
    for (i in 1:ly) {
      y.y <- i2[i]
      mm <- nmu[i]
      ms <- nsigma[i]
      mn<- nnu[i]
      allval <- seq(0, y.y)
      pdfall <- df2.nu.int(allval, mm, ms, mn)
      cdf.f.dnu.dsigma.plus[i] <- sum(pdfall)
    }  
    
    dF2.dnu.dsigma.plus <- cdf.f.dnu.dsigma.plus
    d2F2dnudsigma <- (dF2.dnu.dsigma.plus - dF2.nu) / fd.prec2
    d2F22dnudsigma <- d2F2dnudsigma - d2f2nusigma
    
    
     
  } else { if (VC$margins[2]=="PIG") {
    i2 <- dat[,2]
    sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
    sigma <- exp(sigma.star)
    sigma <- ifelse(sigma>precision, sigma, precision)
    sigma <- ifelse(sigma<10^7, sigma, 10^7)
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^7, mu, 10^7)
    F2 <- pPIG(i2, mu=mu, sigma=sigma)
    F2 <- ifelse(F2<(1-precision), F2, 1-precision)
    F2 <- ifelse(F2>precision, F2, precision)    
    f2 <- dPIG(i2, mu=mu, sigma=sigma)
    f2 <- ifelse(f2>precision, f2, precision)
    F22 <- F2 - f2
    F22 <- ifelse(F22<(1-precision), F22, 1-precision)
    F22 <- ifelse(F22>precision, F22, precision)
    
    
    f2.fd.mu <- dPIG(i2, mu = (mu+fd.prec), sigma = sigma)
    f2.fd.mu <- ifelse(f2.fd.mu>precision, f2.fd.mu, precision)
    df2 <- (f2.fd.mu - f2) / (fd.prec)
    df2 <- as.vector(df2)*as.vector(mu)
    
    f2.fd.sigma <- dPIG(i2, mu=mu, sigma=(sigma+fd.prec))
    f2.fd.sigma <- ifelse(f2.fd.sigma>precision, f2.fd.sigma, precision)
    df2.sigma <- (f2.fd.sigma - f2) / (fd.prec)
    
    F2.fd.mu <- pPIG(i2, mu = (mu+fd.prec), sigma = sigma)
    F2.fd.mu <- ifelse(F2.fd.mu<(1-precision), F2.fd.mu, 1-precision)
    dF2 <- (F2.fd.mu-F2) / (fd.prec)
    dF2 <- as.vector(dF2)*as.vector(mu) 
    dF22 <- dF2 - df2
    
    F2.fd.sigma <- pPIG(i2, mu = mu, sigma = (sigma+fd.prec))
    F2.fd.sigma <- ifelse(F2.fd.sigma<(1-precision), F2.fd.sigma, 1-precision)
    dF2.sigma <- (F2.fd.sigma - F2) / (fd.prec)
    dF22.sigma <- dF2.sigma - df2.sigma
    
    
    
    
    # Hessian derivative components
    
    fd.prec2 <- 10^-3
    
    mu.plus <- ifelse(exp(eta2+fd.prec2)>0.0001, exp(eta2+fd.prec2), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    f2.plus <- dPIG(i2, mu = mu.plus, sigma = sigma)
    f2.plus <- ifelse(f2.plus>precision, f2.plus, precision)
    f2.fd.mu.plus <- dPIG(i2, mu = (mu.plus+fd.prec), sigma = sigma)
    f2.fd.mu.plus <- ifelse(f2.fd.mu.plus>precision, f2.fd.mu.plus, precision)
    df2.plus <- (f2.fd.mu.plus - f2.plus) / (fd.prec)
    df2.plus <- as.vector(df2.plus)*as.vector(mu.plus)
    df2.plus <- as.vector(df2.plus)
    d2f2delta22 <- (df2.plus - df2) / fd.prec2
    
    
    sigma.plus <- sigma + fd.prec2
    f2.plus.sigma <- dPIG(i2, mu = mu, sigma = sigma.plus)
    f2.plus.sigma <- ifelse(f2.plus.sigma>precision, f2.plus.sigma, precision)
    f2.fd.sigma.plus <- dPIG(i2, mu = mu, sigma = (sigma.plus+fd.prec))
    f2.fd.sigma.plus <- ifelse(f2.fd.sigma.plus>precision, f2.fd.sigma.plus, precision)
    df2.sigma.plus <- (f2.fd.sigma.plus - f2.plus.sigma) / (fd.prec)
    d2f2sigma2 <- (df2.sigma.plus - df2.sigma) / fd.prec2
    
    
    f2.dsigma.plus <- dPIG(i2, mu = mu, sigma = sigma.plus)
    f2.dsigma.plus <- ifelse(f2.dsigma.plus>precision, f2.dsigma.plus, precision)
    f2.fd.mu.sigma.plus <- dPIG(i2, mu = (mu+fd.prec), sigma = sigma.plus)
    f2.fd.mu.sigma.plus <- ifelse(f2.fd.mu.sigma.plus>precision, f2.fd.mu.sigma.plus, precision)
    df2.dsigma.plus <- (f2.fd.mu.sigma.plus - f2.dsigma.plus) / (fd.prec)
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)*as.vector(mu)
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec2
    
    
    
    F2.plus <- pPIG(i2, mu = mu.plus, sigma = sigma)
    F2.fd.mu.plus <- pPIG(i2, mu = (mu.plus+fd.prec), sigma = sigma)
    F2.fd.mu.plus <- ifelse(F2.fd.mu.plus<(1-precision), F2.fd.mu.plus, 1-precision)
    dF2.plus <- (F2.fd.mu.plus - F2.plus) / (fd.prec)
    dF2.plus <- as.vector(dF2.plus)*as.vector(mu.plus) 
    d2F2ddelta22 <- (dF2.plus - dF2) / fd.prec2
    d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    
    F2.sigma.plus <- pPIG(i2, mu = mu, sigma = sigma.plus)
    F2.fd.sigma.plus <- pPIG(i2, mu = mu, sigma = (sigma.plus+fd.prec))
    F2.fd.sigma.plus <- ifelse(F2.fd.sigma.plus<(1-precision), F2.fd.sigma.plus, 1-precision)
    dF2.sigma.plus <- (F2.fd.sigma.plus - F2.sigma.plus) / (fd.prec)
    d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) /fd.prec2
    d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
    
    F2.fd.dsigma.plus <- pPIG(i2, mu = (mu+fd.prec), sigma = sigma.plus)
    F2.fd.dsigma.plus <- ifelse(F2.fd.dsigma.plus<(1-precision), F2.fd.dsigma.plus, 1-precision)
    dF2.dsigma.plus <- (F2.fd.dsigma.plus - F2.sigma.plus) / (fd.prec)
    dF2.dsigma.plus <- as.vector(dF2.dsigma.plus)*as.vector(mu) 
    d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec2
    d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma   
    
    
    
    
  } else { if (VC$margins[2]=="S") {
    i2 <- dat[,2]
    sigma.star <- params[(VC$X1.d2+VC$X2.d2+1)]
    nu.star <- params[(VC$X1.d2+VC$X2.d2+2)]
    nu <- nu.star
    sigma <- exp(sigma.star)
    sigma <- ifelse(sigma>precision, sigma, precision)
    sigma <- ifelse(sigma<10^6, sigma, 10^6)
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    mu <- ifelse(mu<10^6, mu, 10^6)    
    sigma <- exp(sigma.star)
    sigma <- ifelse(sigma>precision, sigma, precision)
    mu <- ifelse(exp(eta2)>0.0001, exp(eta2), 0.0001)
    F2 <- pSI(i2, mu=mu, sigma=sigma, nu=nu)
    F2 <- ifelse(F2>precision, F2, precision)
    F2 <- ifelse(F2<(1-precision), F2, 1-precision)    
    f2 <- dSI(i2, mu=mu, sigma=sigma, nu=nu)
    f2 <- ifelse(f2>precision, f2, precision)
    F22 <- F2 - f2
    F22 <- ifelse(F22<(1-precision), F22, 1-precision)
    F22 <- ifelse(F22>precision, F22, precision)
    
    
    
    f2.fd.mu <- dSI(i2, mu = (mu +fd.prec), sigma = sigma, nu = nu)
    f2.fd.mu <- ifelse(f2.fd.mu>precision, f2.fd.mu, precision)
    df2 <- (f2.fd.mu - f2) / (fd.prec)
    df2 <- as.vector(df2)*as.vector(mu)
    
    f2.fd.sigma <- dSI(i2, mu = mu, sigma = (sigma+fd.prec), nu = nu)
    f2.fd.sigma <- ifelse(f2.fd.sigma>precision, f2.fd.sigma, precision)
    df2.sigma <- (f2.fd.sigma - f2) / (fd.prec)
    
    f2.fd.nu <- dSI(i2, mu = mu, sigma = sigma, nu =(nu+fd.prec))
    f2.fd.nu <- ifelse(f2.fd.nu>precision, f2.fd.nu, precision)
    df2.nu <- (f2.fd.nu - f2) / (fd.prec)
    
    F2.fd.mu <- pSI(i2, mu = (mu+fd.prec), sigma = sigma, nu = nu)
    F2.fd.mu <- ifelse(F2.fd.mu<(1-precision), F2.fd.mu, 1-precision)
    dF2 <- (F2.fd.mu - F2) / (fd.prec)
    dF2 <- as.vector(dF2)*as.vector(mu) 
    dF22 <- dF2 - df2
    
    F2.fd.sigma <- pSI(i2, mu = mu, sigma = (sigma+fd.prec), nu = nu)
    F2.fd.sigma <- ifelse(F2.fd.sigma<(1-precision), F2.fd.sigma, 1-precision)
    dF2.sigma <- (F2.fd.sigma - F2) / (fd.prec)
    dF22.sigma <- dF2.sigma - df2.sigma
    
    F2.fd.nu <- pSI(i2, mu = mu, sigma = sigma, nu = (nu+fd.prec))
    F2.fd.nu <- ifelse(F2.fd.nu<(1-precision), F2.fd.nu, 1-precision)
    dF2.nu <- (F2.fd.nu - F2) / (fd.prec)
    dF22.nu <- dF2.nu - df2.nu
    
    
    
    
    # Hessian derivative components
    
    #second derivative fd precision
    fd.prec2 <- 10^-3
    
    
    #pmfs
    
    #delta2^2
    mu.plus <- ifelse(exp(eta2+fd.prec2)>0.0001, exp(eta2+fd.prec2), 0.0001)
    mu.plus <- ifelse(mu.plus<10^7, mu.plus, 10^7)
    f2.plus <- dSI(i2, mu = mu.plus, sigma = sigma, nu = nu)
    f2.plus <- ifelse(f2.plus>precision, f2.plus, precision)
    f2.fd.mu.plus <- dSI(i2, mu = (mu.plus+fd.prec), sigma = sigma, nu = nu)
    f2.fd.mu.plus <- ifelse(f2.fd.mu.plus>precision, f2.fd.mu.plus, precision)
    df2.plus <- (f2.fd.mu.plus - f2.plus) / (fd.prec)
    df2.plus <- as.vector(df2.plus)*as.vector(mu.plus)
    df2.plus <- as.vector(df2.plus)
    d2f2delta22 <- (df2.plus - df2) / fd.prec2
    
    #sigma^2
    sigma.plus <- sigma+fd.prec2
    f2.plus.sigma <- dSI(i2, mu = mu, sigma = sigma.plus, nu = nu)
    f2.plus.sigma <- ifelse(f2.plus.sigma>precision, f2.plus.sigma, precision)
    f2.fd.sigma.plus <- dSI(i2, mu = mu, sigma = (sigma.plus+fd.prec), nu = nu)
    f2.fd.sigma.plus <- ifelse(f2.fd.sigma.plus>precision, f2.fd.sigma.plus, precision)
    df2.sigma.plus <- (f2.fd.sigma.plus - f2.plus.sigma) / (fd.prec)
    d2f2sigma2 <- (df2.sigma.plus - df2.sigma) / fd.prec2
    
    #delta2sigma
    f2.dsigma.plus <- dSI(i2, mu = mu, sigma = sigma.plus, nu = nu)
    f2.dsigma.plus <- ifelse(f2.dsigma.plus>precision, f2.dsigma.plus, precision)
    f2.fd.mu.sigma.plus <- dSI(i2, mu = (mu+fd.prec), sigma = sigma.plus, nu = nu)
    f2.fd.mu.sigma.plus <- ifelse(f2.fd.mu.sigma.plus>precision, f2.fd.mu.sigma.plus, precision)
    df2.dsigma.plus <- (f2.fd.mu.sigma.plus - f2.dsigma.plus) / (fd.prec)
    df2.dsigma.plus <- as.vector(df2.dsigma.plus)*as.vector(mu)
    d2f2delta2sigma <- (df2.dsigma.plus - df2) / fd.prec2
    
    #nu^2
    nu.plus <- nu+fd.prec2
    f2.plus.nu <- dSI(i2, mu = mu, sigma = sigma, nu = nu.plus)
    f2.plus.nu <- ifelse(f2.plus.nu>precision, f2.plus.nu, precision)
    f2.fd.nu.plus <- dSI(i2, mu = mu, sigma = sigma, nu = nu.plus+fd.prec)
    f2.fd.nu.plus <- ifelse(f2.fd.nu.plus>precision, f2.fd.nu.plus, precision)
    df2.nu.plus <- (f2.fd.nu.plus - f2.plus.nu) / (fd.prec)
    d2f2nu2 <- (df2.nu.plus - df2.nu) / fd.prec2
    
    #delta2nu
    f2.dnu.plus <- dSI(i2, mu = mu, sigma = sigma, nu = nu.plus)
    f2.dnu.plus <- ifelse(f2.dnu.plus>precision, f2.dnu.plus, precision)
    f2.fd.mu.nu.plus <- dSI(i2, mu = (mu+fd.prec), sigma = sigma, nu = nu.plus)
    f2.fd.mu.nu.plus <- ifelse(f2.fd.mu.nu.plus>precision, f2.fd.mu.nu.plus, precision)
    df2.dnu.plus <- (f2.fd.mu.nu.plus - f2.dnu.plus) / (fd.prec)
    df2.dnu.plus <- as.vector(df2.dnu.plus)*as.vector(mu)
    d2f2delta2nu <-(df2.dnu.plus - df2) / fd.prec2
    
    #nusigma
    f2.fd.sigma.nu.plus <- dSI(i2, mu = mu, sigma = (sigma+fd.prec), nu = nu.plus)
    f2.fd.sigma.nu.plus <- ifelse(f2.fd.sigma.nu.plus>precision, f2.fd.sigma.nu.plus, precision)
    df2.dnu.sigma.plus <- (f2.fd.sigma.nu.plus - f2.plus.nu) / (fd.prec)
    df2.dnu.sigma.plus <- as.vector(df2.dnu.sigma.plus)
    d2f2sigmanu <- (df2.dnu.sigma.plus - df2.sigma) / fd.prec2
    
    
    
    #cdfs
    
    #delta2^2
    F2.plus <- pSI(i2, mu = mu.plus, sigma = sigma, nu = nu)
    F2.fd.mu.plus <- pSI(i2, mu = (mu.plus+fd.prec), sigma = sigma, nu = nu)
    F2.fd.mu.plus <- ifelse(F2.fd.mu.plus<(1-precision), F2.fd.mu.plus, 1-precision)
    dF2.plus <- (F2.fd.mu.plus - F2.plus) / (fd.prec)
    dF2.plus <- as.vector(dF2.plus)*as.vector(mu.plus) 
    d2F2ddelta22 <-(dF2.plus - dF2) / fd.prec2
    d2F22ddelta22 <- d2F2ddelta22 - d2f2delta22
    
    #sigma^2
    F2.sigma.plus <- pSI(i2, mu = mu, sigma = sigma.plus, nu = nu)
    F2.fd.sigma.plus <- pSI(i2, mu = mu, sigma = (sigma.plus+fd.prec), nu=nu)
    F2.fd.sigma.plus <- ifelse(F2.fd.sigma.plus<(1-precision), F2.fd.sigma.plus, 1-precision)
    dF2.sigma.plus <- (F2.fd.sigma.plus - F2.sigma.plus) / (fd.prec)
    d2F2dsigma2 <- (dF2.sigma.plus - dF2.sigma) / fd.prec2
    d2F22dsigma2 <- d2F2dsigma2 - d2f2sigma2
    
    #delta2sigma
    F2.fd.dsigma.plus <- pSI(i2, mu = (mu+fd.prec), sigma = sigma.plus, nu)
    F2.fd.dsigma.plus <- ifelse(F2.fd.dsigma.plus<(1-precision), F2.fd.dsigma.plus, 1-precision)
    dF2.dsigma.plus <- (F2.fd.dsigma.plus - F2.sigma.plus) / (fd.prec)
    dF2.dsigma.plus <- as.vector(dF2.dsigma.plus)*as.vector(mu) 
    d2F2ddelta2dsigma <- (dF2.dsigma.plus - dF2) / fd.prec2
    d2F22ddelta2dsigma <- d2F2ddelta2dsigma - d2f2delta2sigma
    
    #nu^2
    F2.nu.plus<- pSI(i2, mu = mu, sigma = sigma, nu = nu.plus)
    F2.fd.nu.plus <- pSI(i2, mu = mu, sigma = sigma, nu = nu.plus+fd.prec)
    F2.fd.nu.plus <- ifelse(F2.fd.nu.plus<(1-precision), F2.fd.nu.plus, 1-precision)
    dF2.nu.plus <- (F2.fd.nu.plus - F2.nu.plus) / (fd.prec)
    d2F2dnu2 <- (dF2.nu.plus - dF2.nu) / fd.prec2
    d2F22dnu2 <- d2F2dnu2 - d2f2nu2
    
    #delta2nu
    F2.fd.dnu.plus <- pSI(i2, mu = (mu+fd.prec), sigma = sigma, nu.plus)
    F2.fd.dnu.plus <- ifelse(F2.fd.dnu.plus<(1-precision), F2.fd.dnu.plus, 1-precision)
    dF2.dnu.plus <- (F2.fd.dnu.plus - F2.nu.plus) / (fd.prec)
    dF2.dnu.plus <- as.vector(dF2.dnu.plus)*as.vector(mu) 
    d2F2ddelta2dnu <- (dF2.dnu.plus - dF2) / fd.prec2
    d2F22ddelta2dnu <- d2F2ddelta2dnu - d2f2delta2nu
    
    #sigmanu
    F2.fd.dnu.sigma.plus <- pSI(i2, mu = mu, sigma = (sigma+fd.prec), nu.plus)
    F2.fd.dnu.sigma.plus <- ifelse(F2.fd.dnu.sigma.plus<(1-precision), F2.fd.dnu.sigma.plus, 1-precision)
    dF2.dnu.sigma.plus <- (F2.fd.dnu.sigma.plus - F2.nu.plus) / (fd.prec)
    dF2.dnu.sigma.plus <- as.vector(dF2.dnu.sigma.plus)
    d2F2dsigmadnu <- (dF2.dnu.sigma.plus - dF2.sigma) / fd.prec2
    d2F22dsigmadnu <- d2F2dsigmadnu - d2f2sigmanu
    
    
    
    
    
  }}}}}
  
  
  if(VC$BivD=="FGM") {
    theta <- tanh(theta.star)
    FGM1 <- F1*F2*(1+theta*(1-F1)*(1-F2))
    FGM2 <- F1*F22*(1+theta*(1-F1)*(1-F22))
    
    FGM1 <- ifelse(FGM1>precision, FGM1, precision)
    FGM1 <- ifelse(FGM1<(1-precision), FGM1, 1-precision)
    FGM2 <- ifelse(FGM2>precision, FGM2, precision)
    FGM2 <- ifelse(FGM2<(1-precision), FGM2, 1-precision)
    lx <- f2 - FGM1 + FGM2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- F2*((1+theta*(1-F1)*(1-F2))-F1*theta*(1-F2))
    dcop11 <- F22*((1+theta*(1-F1)*(1-F22))-F1*theta*(1-F22))
    
    dcop2 <- F1*((1+theta*(1-F2)*(1-F1))-F2*theta*(1-F1))
    dcop22 <- F1*((1+theta*(1-F22)*(1-F1))-F22*theta*(1-F1))
    
    dcop.theta1 <- (1-F1)*(1-F2)*F1*F2
    dcop.theta2 <- (1-F1)*(1-F22)*F1*F22
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- F2*(-2*theta*(1-F2))
    d2Cdcop112 <- F22*(-2*theta*(1-F22))
    
    d2Cdcop22 <- F1*(-2*theta*(1-F1))
    d2Cdcop222 <- F1*(-2*theta*(1-F1))
    
    d2Cdcop1dcop2 <- (1+theta*(1-F1)*(1-F2)-F1*theta*(1-F2))+F2*(-theta*(1-F1)+F1*theta)
    d2Cdcop11dcop22 <- (1+theta*(1-F1)*(1-F22)-F1*theta*(1-F22))+F22*(-theta*(1-F1)+F1*theta)
    
    d2Cdcop.theta12 <- F2*((1-F1)*(1-F2)-F1*(1-F2))
    d2Cdcop.theta22 <- F22*((1-F1)*(1-F22)-F1*(1-F22))
    
    d2Cdcop1dcop.theta1 <- F1*((1-F2)*(1-F1)-F2*(1-F1))
    d2Cdcop11dcop.theta2 <- F1*((1-F22)*(1-F1)-F22*(1-F1))
    
    d2Cdcop2dcop.theta1 <- 0
    d2Cdcop22dcop.theta2 <- 0
    
    
  } else { if (VC$BivD=="N") {
    theta <- tanh(theta.star)
    BN1 <- pbinorm(qnorm(F1), qnorm(F2), mean1 = 0, mean2 = 0, var1 = 1, var2 = 1, cov12 = theta) 
    BN2 <- pbinorm(qnorm(F1), qnorm(F22), mean1 = 0, mean2 = 0, var1 = 1, var2 = 1, cov12 = theta)
    
    BN1 <- ifelse(BN1>precision, BN1, precision)
    BN1 <- ifelse(BN1<(1-precision), BN1, 1-precision)
    BN2 <- ifelse(BN2>precision, BN2, precision)
    BN2 <- ifelse(BN2<(1-precision), BN2, 1-precision)
    
    lx <- f2 - BN1 + BN2
    lx <- ifelse(lx>precision, lx, precision)
    
    
    dcop1 <- pnorm((qnorm(F2)-theta*qnorm(F1))/sqrt(1-theta^2))
    dcop11 <- pnorm((qnorm(F22)-theta*qnorm(F1))/sqrt(1-theta^2))
    
    dcop2 <- pnorm((qnorm(F1)-theta*qnorm(F2))/sqrt(1-theta^2))
    dcop22 <- pnorm((qnorm(F1)-theta*qnorm(F22))/sqrt(1-theta^2))
    
    dcop.theta1 <- dbinorm(qnorm(F1), qnorm(F2), mean1 = 0, mean2 = 0, var1 = 1, var2 = 1, cov12 = theta, log = FALSE)
    dcop.theta2 <- dbinorm(qnorm(F1), qnorm(F22), mean1 = 0, mean2 = 0, var1 = 1, var2 = 1, cov12 = theta, log = FALSE)
    
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- dnorm((qnorm(F2)-theta*qnorm(F1))/sqrt(1-theta^2))*(-theta*(dnorm(qnorm(F1)))^(-1)*(1-theta^2)^(-1/2))
    d2Cdcop112 <- dnorm((qnorm(F22)-theta*qnorm(F1))/sqrt(1-theta^2))*(-theta*(dnorm(qnorm(F1)))^(-1)*(1-theta^2)^(-1/2))
    
    d2Cdcop22 <- dnorm((qnorm(F1)-theta*qnorm(F2))/sqrt(1-theta^2))*(-theta*(dnorm(qnorm(F2)))^(-1)*(1-theta^2)^(-1/2))
    d2Cdcop222 <- dnorm((qnorm(F1)-theta*qnorm(F22))/sqrt(1-theta^2))*(-theta*(dnorm(qnorm(F22)))^(-1)*(1-theta^2)^(-1/2))
    
    d2Cdcop1dcop2 <- dnorm((qnorm(F2)-theta*qnorm(F1))/sqrt(1-theta^2))*(dnorm(qnorm(F2)))^(-1)*(1-theta^2)^(-1/2)
    d2Cdcop11dcop22 <- dnorm((qnorm(F22)-theta*qnorm(F1))/sqrt(1-theta^2))*(dnorm(qnorm(F22)))^(-1)*(1-theta^2)^(-1/2)
    
    d2Cdcop.theta12 <- dnorm((qnorm(F2)-theta*qnorm(F1))/sqrt(1-theta^2))*(((-qnorm(F1))*(1-theta^2)^(1/2)+theta*(qnorm(F2)-theta*qnorm(F1))*(1-theta^2)^(-1/2))/(1-theta^2))
    d2Cdcop.theta22 <- dnorm((qnorm(F22)-theta*qnorm(F1))/sqrt(1-theta^2))*(((-qnorm(F1))*(1-theta^2)^(1/2)+theta*(qnorm(F22)-theta*qnorm(F1))*(1-theta^2)^(-1/2))/(1-theta^2))
    
    d2Cdcop1dcop.theta1 <- dnorm((qnorm(F1)-theta*qnorm(F2))/sqrt(1-theta^2))*(((-qnorm(F2))*(1-theta^2)^(1/2)+theta*(qnorm(F1)-theta*qnorm(F2))*(1-theta^2)^(-1/2))/(1-theta^2))
    d2Cdcop11dcop.theta2 <- dnorm((qnorm(F1)-theta*qnorm(F22))/sqrt(1-theta^2))*(((-qnorm(F22))*(1-theta^2)^(1/2)+theta*(qnorm(F1)-theta*qnorm(F22))*(1-theta^2)^(-1/2))/(1-theta^2))
    
    d2Cdcop2dcop.theta1  <- (1-theta^2)^(-1)*theta*dcop.theta1+dcop.theta1*(-(-2*qnorm(F1)*qnorm(F2)*(1-theta^2)+(qnorm(F1)^2-2*theta*qnorm(F1)*qnorm(F2)+qnorm(F2)^2)*(2*theta))/(2*(1-theta^2)^2))
    d2Cdcop22dcop.theta2  <- (1-theta^2)^(-1)*theta*dcop.theta2+dcop.theta2*(-(-2*qnorm(F1)*qnorm(F22)*(1-theta^2)+(qnorm(F1)^2-2*theta*qnorm(F1)*qnorm(F22)+qnorm(F22)^2)*(2*theta))/(2*(1-theta^2)^2))
    
    
    
  } else { if (VC$BivD=="AMH") {
    theta <- tanh(theta.star)
    AMH1 <- (F1*F2)/(1-theta*(1-F1)*(1-F2))
    AMH2 <- (F1*F22)/(1-theta*(1-F1)*(1-F22))
    AMH1 <- ifelse(AMH1>precision, AMH1, precision)
    AMH1 <- ifelse(AMH1<(1-precision), AMH1, 1-precision)
    AMH2 <- ifelse(AMH2>precision, AMH2, precision)
    AMH2 <- ifelse(AMH2<(1-precision), AMH2, 1-precision)
    lx <- f2 - AMH1 + AMH2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- (F2*(1-theta*(1-F1)*(1-F2))-F1*F2*(theta*(1-F2)))/(1-theta*(1-F1)*(1-F2))^2
    dcop11 <- (F22*(1-theta*(1-F1)*(1-F22))-F1*F22*(theta*(1-F22)))/(1-theta*(1-F1)*(1-F22))^2
    
    dcop2 <- (F1*(1-theta*(1-F1)*(1-F2))-F1*F2*(theta*(1-F1)))/(1-theta*(1-F1)*(1-F2))^2
    dcop22 <- (F1*(1-theta*(1-F1)*(1-F22))-F1*F22*(theta*(1-F1)))/(1-theta*(1-F1)*(1-F22))^2
    
    dcop.theta1 <- ((F1*F2)/(1-theta*(1-F1)*(1-F2))^2)*(1-F1)*(1-F2)
    dcop.theta2 <- ((F1*F22)/(1-theta*(1-F1)*(1-F22))^2)*(1-F1)*(1-F22)
    
    
    
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- (-(F2*(1-theta*(1-F1)*(1-F2))-F1*F2*(theta*(1-F2)))*2*(1-theta*(1-F1)*(1-F2))*theta*(1-F2))/(1-theta*(1-F1)*(1-F2))^4
    d2Cdcop112 <- (-(F22*(1-theta*(1-F1)*(1-F22))-F1*F22*(theta*(1-F22)))*2*(1-theta*(1-F1)*(1-F22))*theta*(1-F22))/(1-theta*(1-F1)*(1-F22))^4
    
    d2Cdcop22 <- (-(F1*(1-theta*(1-F1)*(1-F2))-F1*F2*(theta*(1-F1)))*2*(1-theta*(1-F1)*(1-F2))*theta*(1-F1))/(1-theta*(1-F1)*(1-F2))^4
    d2Cdcop222 <- (-(F1*(1-theta*(1-F1)*(1-F22))-F1*F22*(theta*(1-F1)))*2*(1-theta*(1-F1)*(1-F22))*theta*(1-F1))/(1-theta*(1-F1)*(1-F22))^4
    
    d2Cdcop1dcop2 <- (((1-theta*(1-F1)*(1-F2))+F2*theta*(1-F1)-(F1*theta-2*theta*F1*F2))*(1-theta*(1-F1)*(1-F2))^2-(F2*(1-theta*(1-F1)*(1-F2))-F1*F2*(theta*(1-F2)))*(1-theta*(1-F1)*(1-F2))*2*theta*(1-F1))/(1-theta*(1-F1)*(1-F2))^4
    d2Cdcop11dcop22 <- (((1-theta*(1-F1)*(1-F22))+F22*theta*(1-F1)-(F1*theta-2*theta*F1*F22))*(1-theta*(1-F1)*(1-F22))^2-(F22*(1-theta*(1-F1)*(1-F22))-F1*F22*(theta*(1-F22)))*(1-theta*(1-F1)*(1-F22))*2*theta*(1-F1))/(1-theta*(1-F1)*(1-F22))^4
    
    d2Cdcop.theta12 <- ((-F2*(1-F1)*(1-F2)-F1*F2*(1-F2))*(1-theta*(1-F1)*(1-F2))^2-(F2*(1-theta*(1-F1)*(1-F2))-F1*F2*(theta*(1-F2)))*(1-theta*(1-F1)*(1-F2))*2*(-(1-F1)*(1-F2)))/(1-theta*(1-F1)*(1-F2))^4
    d2Cdcop.theta22 <- ((-F22*(1-F1)*(1-F22)-F1*F22*(1-F22))*(1-theta*(1-F1)*(1-F22))^2-(F22*(1-theta*(1-F1)*(1-F22))-F1*F22*(theta*(1-F22)))*(1-theta*(1-F1)*(1-F22))*2*(-(1-F1)*(1-F22)))/(1-theta*(1-F1)*(1-F22))^4
    
    d2Cdcop1dcop.theta1 <- ((-F1*(1-F1)*(1-F2)-F1*F2*(1-F1))*(1-theta*(1-F1)*(1-F2))^2-(F1*(1-theta*(1-F1)*(1-F2))-F1*F2*(theta*(1-F1)))*(1-theta*(1-F1)*(1-F2))*2*(-(1-F1)*(1-F2)))/(1-theta*(1-F1)*(1-F2))^4
    d2Cdcop11dcop.theta2 <- ((-F1*(1-F1)*(1-F22)-F1*F22*(1-F1))*(1-theta*(1-F1)*(1-F22))^2-(F1*(1-theta*(1-F1)*(1-F22))-F1*F22*(theta*(1-F1)))*(1-theta*(1-F1)*(1-F22))*2*(-(1-F1)*(1-F22)))/(1-theta*(1-F1)*(1-F22))^4
    
    d2Cdcop2dcop.theta1 <- F1*F2*(1-F1)*(1-F2)*2*(1-theta*(1-F1)*(1-F2))^(-3)*(1-F1)*(1-F2)
    d2Cdcop22dcop.theta2 <- F1*F22*(1-F1)*(1-F22)*2*(1-theta*(1-F1)*(1-F22))^(-3)*(1-F1)*(1-F22)
    
    
    
    
  } else { if (VC$BivD=="C0") {
    theta <- exp(theta.star) + precision
    Clayton1 <- (F1^(-theta)+F2^(-theta)-1)^(-1/theta)
    Clayton2 <- (F1^(-theta)+F22^(-theta)-1)^(-1/theta)
    Clayton1 <- ifelse(Clayton1>precision, Clayton1, precision)
    Clayton1 <- ifelse(Clayton1<(1-precision), Clayton1, 1-precision)
    Clayton2 <- ifelse(Clayton2>precision, Clayton2, precision)
    Clayton2 <- ifelse(Clayton2<(1-precision), Clayton2, 1-precision)
    lx <- f2 - Clayton1+Clayton2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- F1^(-theta-1)*(F1^(-theta)+F2^(-theta)-1)^(-(1/theta)-1)
    dcop11 <- F1^(-theta-1)*(F1^(-theta)+F22^(-theta)-1)^(-(1/theta)-1)
    
    dcop2 <- F2^(-theta-1)*(F2^(-theta)+F1^(-theta)-1)^(-(1/theta)-1)
    dcop22 <- F22^(-theta-1)*(F22^(-theta)+F1^(-theta)-1)^(-(1/theta)-1)
    
    dcop.theta1 <- (1/theta^2)*log(F1^(-theta)+F2^(-theta)-1)*(F1^(-theta)+F2^(-theta)-1)^(-1/theta)+(1/theta)*(F1^(-theta)*log(F1)+F2^(-theta)*log(F2))*(F1^(-theta)+F2^(-theta)-1)^(-1/theta-1)
    dcop.theta2 <- (1/theta^2)*log(F1^(-theta)+F22^(-theta)-1)*(F1^(-theta)+F22^(-theta)-1)^(-1/theta)+(1/theta)*(F1^(-theta)*log(F1)+F22^(-theta)*log(F22))*(F1^(-theta)+F22^(-theta)-1)^(-1/theta-1)
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- (-theta-1)*F1^(-theta-2)*(F1^(-theta)+F2^(-theta)-1)^(-(1/theta)-1)+F1^(-theta-1)*(F1^(-theta)+F2^(-theta)-1)^(-(1/theta)-2)*F1^(-theta-1)*(-theta)*(-1/theta-1)
    d2Cdcop112 <- (-theta-1)*F1^(-theta-2)*(F1^(-theta)+F22^(-theta)-1)^(-(1/theta)-1)+F1^(-theta-1)*(F1^(-theta)+F22^(-theta)-1)^(-(1/theta)-2)*F1^(-theta-1)*(-theta)*(-1/theta-1)
    
    d2Cdcop22 <- (-theta-1)*F2^(-theta-2)*(F1^(-theta)+F2^(-theta)-1)^(-(1/theta)-1)+F2^(-theta-1)*(F1^(-theta)+F2^(-theta)-1)^(-(1/theta)-2)*F2^(-theta-1)*(-theta)*(-1/theta-1)
    d2Cdcop222 <- (-theta-1)*F22^(-theta-2)*(F1^(-theta)+F22^(-theta)-1)^(-(1/theta)-1)+F22^(-theta-1)*(F1^(-theta)+F22^(-theta)-1)^(-(1/theta)-2)*F22^(-theta-1)*(-theta)*(-1/theta-1)
    
    d2Cdcop1dcop2 <- F1^(-theta-1)*(F1^(-theta)+F2^(-theta)-1)^(-(1/theta)-2)*(-1/theta-1)*(-theta)*F2^(-theta-1)
    d2Cdcop11dcop22 <- F1^(-theta-1)*(F1^(-theta)+F22^(-theta)-1)^(-(1/theta)-2)*(-1/theta-1)*(-theta)*F22^(-theta-1)
    
    d2Cdcop.theta12 <- F1^(-theta-1)*log(F1)*(-1)*(F1^(-theta)+F2^(-theta)-1)^(-(1/theta)-1)+F1^(-theta-1)*(dcop.theta1*(F1^(-theta)+F2^(-theta)-1)^(-1)+(F1^(-theta)+F2^(-theta)-1)^(-(1/theta))*(F1^(-theta)+F2^(-theta)-1)^(-2)*(F1^(-theta)*log(F1)+F2^(-theta)*log(F2)))
    d2Cdcop.theta22 <- F1^(-theta-1)*log(F1)*(-1)*(F1^(-theta)+F22^(-theta)-1)^(-(1/theta)-1)+F1^(-theta-1)*(dcop.theta2*(F1^(-theta)+F22^(-theta)-1)^(-1)+(F1^(-theta)+F22^(-theta)-1)^(-(1/theta))*(F1^(-theta)+F22^(-theta)-1)^(-2)*(F1^(-theta)*log(F1)+F22^(-theta)*log(F22)))
    
    d2Cdcop1dcop.theta1 <- F2^(-theta-1)*log(F2)*(-1)*(F1^(-theta)+F2^(-theta)-1)^(-(1/theta)-1)+F2^(-theta-1)*(dcop.theta1*(F1^(-theta)+F2^(-theta)-1)^(-1)+(F1^(-theta)+F2^(-theta)-1)^(-(1/theta))*(F1^(-theta)+F2^(-theta)-1)^(-2)*(F1^(-theta)*log(F1)+F2^(-theta)*log(F2)))
    d2Cdcop11dcop.theta2 <- F22^(-theta-1)*log(F22)*(-1)*(F1^(-theta)+F22^(-theta)-1)^(-(1/theta)-1)+F22^(-theta-1)*(dcop.theta2*(F1^(-theta)+F22^(-theta)-1)^(-1)+(F1^(-theta)+F22^(-theta)-1)^(-(1/theta))*(F1^(-theta)+F22^(-theta)-1)^(-2)*(F1^(-theta)*log(F1)+F22^(-theta)*log(F22)))
    
    d2Cdcop2dcop.theta1 <- -2/theta^3*log(F1^(-theta)+F2^(-theta)-1)*(F1^(-theta)+F2^(-theta)-1)^(-(1/theta))+1/theta^2*((F1^(-theta)+F2^(-theta)-1)^(-(1/theta)-1)*(-F1^(-theta)*log(F1)-F2^(-theta)*log(F2))+log(F1^(-theta)+F2^(-theta)-1)*dcop.theta1)+(-1/theta^2*(F1^(-theta)*log(F1)+F2^(-theta)*log(F2))+1/theta*(-F1^(-theta)*(log(F1))^2-F2^(-theta)*(log(F2))^2))*(F1^(-theta)+F2^(-theta)-1)^(-(1/theta)-1)+1/theta*(F1^(-theta)*log(F1)+F2^(-theta)*log(F2))*(dcop.theta1*(F1^(-theta)+F2^(-theta)-1)^(-1)+(F1^(-theta)+F2^(-theta)-1)^(-(1/theta)-2)*(F1^(-theta)*log(F1)+F2^(-theta)*log(F2)))
    d2Cdcop22dcop.theta2 <- -2/theta^3*log(F1^(-theta)+F22^(-theta)-1)*(F1^(-theta)+F22^(-theta)-1)^(-(1/theta))+1/theta^2*((F1^(-theta)+F22^(-theta)-1)^(-(1/theta)-1)*(-F1^(-theta)*log(F1)-F22^(-theta)*log(F22))+log(F1^(-theta)+F22^(-theta)-1)*dcop.theta2)+(-1/theta^2*(F1^(-theta)*log(F1)+F22^(-theta)*log(F22))+1/theta*(-F1^(-theta)*(log(F1))^2-F22^(-theta)*(log(F22))^2))*(F1^(-theta)+F22^(-theta)-1)^(-(1/theta)-1)+1/theta*(F1^(-theta)*log(F1)+F22^(-theta)*log(F22))*(dcop.theta2*(F1^(-theta)+F22^(-theta)-1)^(-1)+(F1^(-theta)+F22^(-theta)-1)^(-(1/theta)-2)*(F1^(-theta)*log(F1)+F22^(-theta)*log(F22)))
    
    
  } else { if (VC$BivD=="C90") {
    theta <- exp(theta.star) + precision
    
    Clayton1 <- F2 -((1-F1)^(-theta)+F2^(-theta)-1)^(-1/theta)
    Clayton2 <- F22 -((1-F1)^(-theta)+F22^(-theta)-1)^(-1/theta)
    Clayton1 <- ifelse(Clayton1>precision, Clayton1, precision)
    Clayton1 <- ifelse(Clayton1<(1-precision), Clayton1, 1-precision)
    Clayton2 <- ifelse(Clayton2>precision, Clayton2, precision)
    Clayton2 <- ifelse(Clayton2<(1-precision), Clayton2, 1-precision)
    lx <- f2 - Clayton1 + Clayton2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- (1-F1)^(-theta-1)*((1-F1)^(-theta)+F2^(-theta)-1)^(-(1/theta)-1)
    dcop11 <- (1-F1)^(-theta-1)*((1-F1)^(-theta)+F22^(-theta)-1)^(-(1/theta)-1)
    
    dcop2 <- 1-F2^(-theta-1)*(F2^(-theta)+(1-F1)^(-theta)-1)^(-(1/theta)-1)
    dcop22 <- 1-F22^(-theta-1)*(F22^(-theta)+(1-F1)^(-theta)-1)^(-(1/theta)-1)
    
    dcop.theta1 <- ( (1/theta^2)*log((1-F1)^(-theta)+F2^(-theta)-1)*((1-F1)^(-theta)+F2^(-theta)-1)^(-1/theta)+(1/theta)*((1-F1)^(-theta)*log((1-F1))+F2^(-theta)*log(F2))*((1-F1)^(-theta)+F2^(-theta)-1)^(-1/theta-1)  )
    dcop.theta2 <- ( (1/theta^2)*log((1-F1)^(-theta)+F22^(-theta)-1)*((1-F1)^(-theta)+F22^(-theta)-1)^(-1/theta)+(1/theta)*((1-F1)^(-theta)*log((1-F1))+F22^(-theta)*log(F22))*((1-F1)^(-theta)+F22^(-theta)-1)^(-1/theta-1)  )
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- -( (-theta-1)*(1-F1)^(-theta-2)*((1-F1)^(-theta)+F2^(-theta)-1)^(-(1/theta)-1)+(1-F1)^(-theta-1)*((1-F1)^(-theta)+F2^(-theta)-1)^(-(1/theta)-2)*(1-F1)^(-theta-1)*(-theta)*(-1/theta-1) )
    d2Cdcop112 <- -(  (-theta-1)*(1-F1)^(-theta-2)*((1-F1)^(-theta)+F22^(-theta)-1)^(-(1/theta)-1)+(1-F1)^(-theta-1)*((1-F1)^(-theta)+F22^(-theta)-1)^(-(1/theta)-2)*(1-F1)^(-theta-1)*(-theta)*(-1/theta-1) )
    
    d2Cdcop22 <- -(  (-theta-1)*F2^(-theta-2)*((1-F1)^(-theta)+F2^(-theta)-1)^(-(1/theta)-1)+F2^(-theta-1)*((1-F1)^(-theta)+F2^(-theta)-1)^(-(1/theta)-2)*F2^(-theta-1)*(-theta)*(-1/theta-1) )
    d2Cdcop222 <- -( (-theta-1)*F22^(-theta-2)*((1-F1)^(-theta)+F22^(-theta)-1)^(-(1/theta)-1)+F22^(-theta-1)*((1-F1)^(-theta)+F22^(-theta)-1)^(-(1/theta)-2)*F22^(-theta-1)*(-theta)*(-1/theta-1) )
    
    d2Cdcop1dcop2 <- (1-F1)^(-theta-1)*((1-F1)^(-theta)+F2^(-theta)-1)^(-(1/theta)-2)*(-1/theta-1)*(-theta)*F2^(-theta-1)
    d2Cdcop11dcop22 <- (1-F1)^(-theta-1)*((1-F1)^(-theta)+F22^(-theta)-1)^(-(1/theta)-2)*(-1/theta-1)*(-theta)*F22^(-theta-1)
    
    d2Cdcop.theta12 <- -( (1-F1)^(-theta-1)*log((1-F1))*(-1)*((1-F1)^(-theta)+F2^(-theta)-1)^(-(1/theta)-1)+(1-F1)^(-theta-1)*(dcop.theta1*((1-F1)^(-theta)+F2^(-theta)-1)^(-1)+((1-F1)^(-theta)+F2^(-theta)-1)^(-(1/theta))*((1-F1)^(-theta)+F2^(-theta)-1)^(-2)*((1-F1)^(-theta)*log((1-F1))+F2^(-theta)*log(F2)))   )
    d2Cdcop.theta22 <- -( (1-F1)^(-theta-1)*log((1-F1))*(-1)*((1-F1)^(-theta)+F22^(-theta)-1)^(-(1/theta)-1)+(1-F1)^(-theta-1)*(dcop.theta2*((1-F1)^(-theta)+F22^(-theta)-1)^(-1)+((1-F1)^(-theta)+F22^(-theta)-1)^(-(1/theta))*((1-F1)^(-theta)+F22^(-theta)-1)^(-2)*((1-F1)^(-theta)*log((1-F1))+F22^(-theta)*log(F22)))    )
    
    d2Cdcop1dcop.theta1 <- ( F2^(-theta-1)*log(F2)*(-1)*((1-F1)^(-theta)+F2^(-theta)-1)^(-(1/theta)-1)+F2^(-theta-1)*(dcop.theta1*((1-F1)^(-theta)+F2^(-theta)-1)^(-1)+((1-F1)^(-theta)+F2^(-theta)-1)^(-(1/theta))*((1-F1)^(-theta)+F2^(-theta)-1)^(-2)*((1-F1)^(-theta)*log((1-F1))+F2^(-theta)*log(F2)))  )
    d2Cdcop11dcop.theta2 <- ( F22^(-theta-1)*log(F22)*(-1)*((1-F1)^(-theta)+F22^(-theta)-1)^(-(1/theta)-1)+F22^(-theta-1)*(dcop.theta2*((1-F1)^(-theta)+F22^(-theta)-1)^(-1)+((1-F1)^(-theta)+F22^(-theta)-1)^(-(1/theta))*((1-F1)^(-theta)+F22^(-theta)-1)^(-2)*((1-F1)^(-theta)*log((1-F1))+F22^(-theta)*log(F22)))  )
    
    d2Cdcop2dcop.theta1 <- -(  -2/theta^3*log((1-F1)^(-theta)+F2^(-theta)-1)*((1-F1)^(-theta)+F2^(-theta)-1)^(-(1/theta))+1/theta^2*(((1-F1)^(-theta)+F2^(-theta)-1)^(-(1/theta)-1)*(-(1-F1)^(-theta)*log((1-F1))-F2^(-theta)*log(F2))+log((1-F1)^(-theta)+F2^(-theta)-1)*dcop.theta1)+(-1/theta^2*((1-F1)^(-theta)*log((1-F1))+F2^(-theta)*log(F2))+1/theta*(-(1-F1)^(-theta)*(log((1-F1)))^2-F2^(-theta)*(log(F2))^2))*((1-F1)^(-theta)+F2^(-theta)-1)^(-(1/theta)-1)+1/theta*((1-F1)^(-theta)*log((1-F1))+F2^(-theta)*log(F2))*(dcop.theta1*((1-F1)^(-theta)+F2^(-theta)-1)^(-1)+((1-F1)^(-theta)+F2^(-theta)-1)^(-(1/theta)-2)*((1-F1)^(-theta)*log((1-F1))+F2^(-theta)*log(F2)))   )
    d2Cdcop22dcop.theta2 <- -(  -2/theta^3*log((1-F1)^(-theta)+F22^(-theta)-1)*((1-F1)^(-theta)+F22^(-theta)-1)^(-(1/theta))+1/theta^2*(((1-F1)^(-theta)+F22^(-theta)-1)^(-(1/theta)-1)*(-(1-F1)^(-theta)*log((1-F1))-F22^(-theta)*log(F22))+log((1-F1)^(-theta)+F22^(-theta)-1)*dcop.theta2)+(-1/theta^2*((1-F1)^(-theta)*log((1-F1))+F22^(-theta)*log(F22))+1/theta*(-(1-F1)^(-theta)*(log((1-F1)))^2-F22^(-theta)*(log(F22))^2))*((1-F1)^(-theta)+F22^(-theta)-1)^(-(1/theta)-1)+1/theta*((1-F1)^(-theta)*log((1-F1))+F22^(-theta)*log(F22))*(dcop.theta2*((1-F1)^(-theta)+F22^(-theta)-1)^(-1)+((1-F1)^(-theta)+F22^(-theta)-1)^(-(1/theta)-2)*((1-F1)^(-theta)*log((1-F1))+F22^(-theta)*log(F22)))    )
    
    
    
  } else { if (VC$BivD=="C180") {
    theta <- exp(theta.star) + precision
    Clayton1 <- F1 + F2 - 1 + ((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-1/theta)
    Clayton2 <- F1 + F22 - 1 + ((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-1/theta)
    Clayton1 <- ifelse(Clayton1>precision, Clayton1, precision)
    Clayton1 <- ifelse(Clayton1<(1-precision), Clayton1, 1-precision)
    Clayton2 <-ifelse(Clayton2>precision, Clayton2, precision)
    Clayton2 <- ifelse(Clayton2<(1-precision), Clayton2, 1-precision)
    lx <- f2 - Clayton1 + Clayton2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- 1-(1-F1)^(-theta-1)*((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-1)
    dcop11 <-1-(1-F1)^(-theta-1)*((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-1)
    
    dcop2 <- 1-(1-F2)^(-theta-1)*((1-F2)^(-theta)+(1-F1)^(-theta)-1)^(-(1/theta)-1)
    dcop22 <- 1-(1-F22)^(-theta-1)*((1-F22)^(-theta)+(1-F1)^(-theta)-1)^(-(1/theta)-1)
    
    dcop.theta1 <- (  (1/theta^2)*log((1-F1)^(-theta)+(1-F2)^(-theta)-1)*((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-1/theta)+(1/theta)*((1-F1)^(-theta)*log((1-F1))+(1-F2)^(-theta)*log((1-F2)))*((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-1/theta-1)   )
    dcop.theta2 <- (  (1/theta^2)*log((1-F1)^(-theta)+(1-F22)^(-theta)-1)*((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-1/theta)+(1/theta)*((1-F1)^(-theta)*log((1-F1))+(1-F22)^(-theta)*log((1-F22)))*((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-1/theta-1)    )
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- ( (-theta-1)*(1-F1)^(-theta-2)*((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-1)+(1-F1)^(-theta-1)*((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-2)*(1-F1)^(-theta-1)*(-theta)*(-1/theta-1) )
    d2Cdcop112 <- (  (-theta-1)*(1-F1)^(-theta-2)*((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-1)+(1-F1)^(-theta-1)*((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-2)*(1-F1)^(-theta-1)*(-theta)*(-1/theta-1) )
    
    d2Cdcop22 <- (  (-theta-1)*(1-F2)^(-theta-2)*((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-1)+(1-F2)^(-theta-1)*((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-2)*(1-F2)^(-theta-1)*(-theta)*(-1/theta-1) )
    d2Cdcop222 <- ( (-theta-1)*(1-F22)^(-theta-2)*((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-1)+(1-F22)^(-theta-1)*((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-2)*(1-F22)^(-theta-1)*(-theta)*(-1/theta-1) )
    
    d2Cdcop1dcop2 <-  ( (1-F1)^(-theta-1)*((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-2)*(-1/theta-1)*(-theta)*(1-F2)^(-theta-1)  )
    d2Cdcop11dcop22 <-  ( (1-F1)^(-theta-1)*((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-2)*(-1/theta-1)*(-theta)*(1-F22)^(-theta-1)  )
    
    d2Cdcop.theta12 <- -(  (1-F1)^(-theta-1)*log((1-F1))*(-1)*((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-1)+(1-F1)^(-theta-1)*(dcop.theta1*((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-1)+((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta))*((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-2)*((1-F1)^(-theta)*log((1-F1))+(1-F2)^(-theta)*log((1-F2))))   )
    d2Cdcop.theta22 <-   -(  (1-F1)^(-theta-1)*log((1-F1))*(-1)*((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-1)+(1-F1)^(-theta-1)*(dcop.theta2*((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-1)+((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta))*((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-2)*((1-F1)^(-theta)*log((1-F1))+(1-F22)^(-theta)*log((1-F22))))   )
    
    d2Cdcop1dcop.theta1 <- -( (1-F2)^(-theta-1)*log((1-F2))*(-1)*((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-1)+(1-F2)^(-theta-1)*(dcop.theta1*((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-1)+((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta))*((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-2)*((1-F1)^(-theta)*log((1-F1))+(1-F2)^(-theta)*log((1-F2))))     )
    d2Cdcop11dcop.theta2 <- -( (1-F22)^(-theta-1)*log((1-F22))*(-1)*((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-1)+(1-F22)^(-theta-1)*(dcop.theta2*((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-1)+((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta))*((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-2)*((1-F1)^(-theta)*log((1-F1))+(1-F22)^(-theta)*log((1-F22))))     )
    
    d2Cdcop2dcop.theta1 <- (  -2/theta^3*log((1-F1)^(-theta)+(1-F2)^(-theta)-1)*((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta))+1/theta^2*(((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-1)*(-(1-F1)^(-theta)*log((1-F1))-(1-F2)^(-theta)*log((1-F2)))+log((1-F1)^(-theta)+(1-F2)^(-theta)-1)*dcop.theta1)+(-1/theta^2*((1-F1)^(-theta)*log((1-F1))+(1-F2)^(-theta)*log((1-F2)))+1/theta*(-(1-F1)^(-theta)*(log((1-F1)))^2-(1-F2)^(-theta)*(log((1-F2)))^2))*((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-1)+1/theta*((1-F1)^(-theta)*log((1-F1))+(1-F2)^(-theta)*log((1-F2)))*(dcop.theta1*((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-1)+((1-F1)^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-2)*((1-F1)^(-theta)*log((1-F1))+(1-F2)^(-theta)*log((1-F2))))   )
    d2Cdcop22dcop.theta2 <- (  -2/theta^3*log((1-F1)^(-theta)+(1-F22)^(-theta)-1)*((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta))+1/theta^2*(((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-1)*(-(1-F1)^(-theta)*log((1-F1))-(1-F22)^(-theta)*log((1-F22)))+log((1-F1)^(-theta)+(1-F22)^(-theta)-1)*dcop.theta2)+(-1/theta^2*((1-F1)^(-theta)*log((1-F1))+(1-F22)^(-theta)*log((1-F22)))+1/theta*(-(1-F1)^(-theta)*(log((1-F1)))^2-(1-F22)^(-theta)*(log((1-F22)))^2))*((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-1)+1/theta*((1-F1)^(-theta)*log((1-F1))+(1-F22)^(-theta)*log((1-F22)))*(dcop.theta2*((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-1)+((1-F1)^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-2)*((1-F1)^(-theta)*log((1-F1))+(1-F22)^(-theta)*log((1-F22))))    )
    
    
    
    
  } else { if (VC$BivD=="C270") {
    theta <- exp(theta.star) + precision
    Clayton1 <- F1 - (F1^(-theta)+(1-F2)^(-theta)-1)^(-1/theta)
    Clayton2 <- F1 - (F1^(-theta)+(1-F22)^(-theta)-1)^(-1/theta)
    Clayton1 <- ifelse(Clayton1>precision, Clayton1, precision)
    Clayton1 <- ifelse(Clayton1<(1-precision), Clayton1, 1-precision)
    Clayton2 <- ifelse(Clayton2>precision, Clayton2, precision)
    Clayton2 <- ifelse(Clayton2<(1-precision), Clayton2, 1-precision)
    lx <- f2 - Clayton1 + Clayton2
    lx <- ifelse(lx>precision, lx, precision)
    
    
    dcop1 <- 1-F1^(-theta-1)*(F1^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-1)
    dcop11 <- 1-F1^(-theta-1)*(F1^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-1)
    
    dcop2 <- (1-F2)^(-theta-1)*((1-F2)^(-theta)+F1^(-theta)-1)^(-(1/theta)-1)
    dcop22 <- (1-F22)^(-theta-1)*((1-F22)^(-theta)+F1^(-theta)-1)^(-(1/theta)-1)
    
    dcop.theta1 <-  (  (1/theta^2)*log(F1^(-theta)+(1-F2)^(-theta)-1)*(F1^(-theta)+(1-F2)^(-theta)-1)^(-1/theta)+(1/theta)*(F1^(-theta)*log(F1)+(1-F2)^(-theta)*log((1-F2)))*(F1^(-theta)+(1-F2)^(-theta)-1)^(-1/theta-1)   )    
    dcop.theta2 <-    (   (1/theta^2)*log(F1^(-theta)+(1-F22)^(-theta)-1)*(F1^(-theta)+(1-F22)^(-theta)-1)^(-1/theta)+(1/theta)*(F1^(-theta)*log(F1)+(1-F22)^(-theta)*log((1-F22)))*(F1^(-theta)+(1-F22)^(-theta)-1)^(-1/theta-1)    )    
    
    
    
    # Hessian derivative components
    
    d2Cdcop12 <-  -( (-theta-1)*F1^(-theta-2)*(F1^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-1)+F1^(-theta-1)*(F1^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-2)*F1^(-theta-1)*(-theta)*(-1/theta-1)   )
    d2Cdcop112 <- -(  (-theta-1)*F1^(-theta-2)*(F1^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-1)+F1^(-theta-1)*(F1^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-2)*F1^(-theta-1)*(-theta)*(-1/theta-1)   )
    
    d2Cdcop22 <- -( (-theta-1)*(1-F2)^(-theta-2)*(F1^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-1)+(1-F2)^(-theta-1)*(F1^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-2)*(1-F2)^(-theta-1)*(-theta)*(-1/theta-1)   )
    d2Cdcop222 <- -( (-theta-1)*(1-F22)^(-theta-2)*(F1^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-1)+(1-F22)^(-theta-1)*(F1^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-2)*(1-F22)^(-theta-1)*(-theta)*(-1/theta-1)    )
    
    d2Cdcop1dcop2 <- F1^(-theta-1)*(F1^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-2)*(-1/theta-1)*(-theta)*(1-F2)^(-theta-1)
    d2Cdcop11dcop22 <- F1^(-theta-1)*(F1^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-2)*(-1/theta-1)*(-theta)*(1-F22)^(-theta-1)
    
    d2Cdcop.theta12 <-   (  F1^(-theta-1)*log(F1)*(-1)*(F1^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-1)+F1^(-theta-1)*(dcop.theta1*(F1^(-theta)+(1-F2)^(-theta)-1)^(-1)+(F1^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta))*(F1^(-theta)+(1-F2)^(-theta)-1)^(-2)*(F1^(-theta)*log(F1)+(1-F2)^(-theta)*log((1-F2))))    )       
    d2Cdcop.theta22 <-   (  F1^(-theta-1)*log(F1)*(-1)*(F1^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-1)+F1^(-theta-1)*(dcop.theta2*(F1^(-theta)+(1-F22)^(-theta)-1)^(-1)+(F1^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta))*(F1^(-theta)+(1-F22)^(-theta)-1)^(-2)*(F1^(-theta)*log(F1)+(1-F22)^(-theta)*log((1-F22))))     )        
    
    d2Cdcop1dcop.theta1 <-  -(  (1-F2)^(-theta-1)*log((1-F2))*(-1)*(F1^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-1)+(1-F2)^(-theta-1)*(dcop.theta1*(F1^(-theta)+(1-F2)^(-theta)-1)^(-1)+(F1^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta))*(F1^(-theta)+(1-F2)^(-theta)-1)^(-2)*(F1^(-theta)*log(F1)+(1-F2)^(-theta)*log((1-F2))))       )
    d2Cdcop11dcop.theta2 <-   -(  (1-F22)^(-theta-1)*log((1-F22))*(-1)*(F1^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-1)+(1-F22)^(-theta-1)*(dcop.theta2*(F1^(-theta)+(1-F22)^(-theta)-1)^(-1)+(F1^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta))*(F1^(-theta)+(1-F22)^(-theta)-1)^(-2)*(F1^(-theta)*log(F1)+(1-F22)^(-theta)*log((1-F22))))       )     
    
    d2Cdcop2dcop.theta1 <-    -  (  -2/theta^3*log(F1^(-theta)+(1-F2)^(-theta)-1)*(F1^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta))+1/theta^2*((F1^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-1)*(-F1^(-theta)*log(F1)-(1-F2)^(-theta)*log((1-F2)))+log(F1^(-theta)+(1-F2)^(-theta)-1)*dcop.theta1)+(-1/theta^2*(F1^(-theta)*log(F1)+(1-F2)^(-theta)*log((1-F2)))+1/theta*(-F1^(-theta)*(log(F1))^2-(1-F2)^(-theta)*(log((1-F2)))^2))*(F1^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-1)+1/theta*(F1^(-theta)*log(F1)+(1-F2)^(-theta)*log((1-F2)))*(dcop.theta1*(F1^(-theta)+(1-F2)^(-theta)-1)^(-1)+(F1^(-theta)+(1-F2)^(-theta)-1)^(-(1/theta)-2)*(F1^(-theta)*log(F1)+(1-F2)^(-theta)*log((1-F2))))    )      
    d2Cdcop22dcop.theta2 <-      - ( -2/theta^3*log(F1^(-theta)+(1-F22)^(-theta)-1)*(F1^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta))+1/theta^2*((F1^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-1)*(-F1^(-theta)*log(F1)-(1-F22)^(-theta)*log((1-F22)))+log(F1^(-theta)+(1-F22)^(-theta)-1)*dcop.theta2)+(-1/theta^2*(F1^(-theta)*log(F1)+(1-F22)^(-theta)*log((1-F22)))+1/theta*(-F1^(-theta)*(log(F1))^2-(1-F22)^(-theta)*(log((1-F22)))^2))*(F1^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-1)+1/theta*(F1^(-theta)*log(F1)+(1-F22)^(-theta)*log((1-F22)))*(dcop.theta2*(F1^(-theta)+(1-F22)^(-theta)-1)^(-1)+(F1^(-theta)+(1-F22)^(-theta)-1)^(-(1/theta)-2)*(F1^(-theta)*log(F1)+(1-F22)^(-theta)*log((1-F22))))    )      
    
    
    
    
  } else { if (VC$BivD=="F") {
    theta <- theta.star + precision
    
    Frank1 <- -theta^(-1)*log(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))
    Frank2 <- -theta^(-1)*log(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))
    
    Frank1 <- ifelse(Frank1>precision, Frank1, precision)
    Frank1 <- ifelse(Frank1<(1-precision), Frank1, 1-precision)
    Frank2 <- ifelse(Frank2>precision, Frank2, precision)
    Frank2 <- ifelse(Frank2<(1-precision), Frank2, 1-precision)
    lx <- f2 - Frank1 + Frank2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- -(1/theta)*(1/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1)))*((exp(-theta*F2)-1)/(exp(-theta)-1))*exp(-theta*F1)*(-theta)
    dcop11 <- -(1/theta)*(1/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1)))*((exp(-theta*F22)-1)/(exp(-theta)-1))*exp(-theta*F1)*(-theta) 
    
    dcop2 <- -(1/theta)*(1/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1)))*((exp(-theta*F1)-1)/(exp(-theta)-1))*exp(-theta*F2)*(-theta)
    dcop22 <- -(1/theta)*(1/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1)))*((exp(-theta*F1)-1)/(exp(-theta)-1))*exp(-theta*F22)*(-theta)
    
    dcop.theta1 <- (1/theta^2)*log(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))-(1/theta)*(1/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1)))*((-exp(-theta*F1)*F1*(exp(-theta*F2)-1)-exp(-theta*F2)*F2*(exp(-theta*F1)-1))*(exp(-theta)-1)+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)*exp(-theta))/(exp(-theta)-1)^2
    dcop.theta2 <- (1/theta^2)*log(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))-(1/theta)*(1/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1)))*((-exp(-theta*F1)*F1*(exp(-theta*F22)-1)-exp(-theta*F22)*F22*(exp(-theta*F1)-1))*(exp(-theta)-1)+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)*exp(-theta))/(exp(-theta)-1)^2
    
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- ((exp(-theta*F2)-1)/(exp(-theta)-1)*exp(-theta*F1)*(-theta)*(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))-(exp(-theta*F2)-1)/(exp(-theta)-1)*exp(-theta*F1)*(exp(-theta*F2)-1)/(exp(-theta)-1)*exp(-theta*F1)*(-theta))/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))^2
    d2Cdcop112 <- ((exp(-theta*F22)-1)/(exp(-theta)-1)*exp(-theta*F1)*(-theta)*(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))-(exp(-theta*F22)-1)/(exp(-theta)-1)*exp(-theta*F1)*(exp(-theta*F22)-1)/(exp(-theta)-1)*exp(-theta*F1)*(-theta))/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))^2
    
    d2Cdcop22 <- ((exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F2)*(-theta)*(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))-(exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F2)*(exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F2)*(-theta))/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))^2
    d2Cdcop222 <- ((exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F22)*(-theta)*(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))-(exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F22)*(exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F22)*(-theta))/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))^2
    
    d2Cdcop1dcop2 <- (exp(-theta*F1)/(exp(-theta)-1)*exp(-theta*F2)*(-theta)*(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))-(exp(-theta*F2)-1)/(exp(-theta)-1)*exp(-theta*F1)*(exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F2)*(-theta))/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))^2
    d2Cdcop11dcop22 <- (exp(-theta*F1)/(exp(-theta)-1)*exp(-theta*F22)*(-theta)*(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))-(exp(-theta*F22)-1)/(exp(-theta)-1)*exp(-theta*F1)*(exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F22)*(-theta))/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))^2
    
    d2Cdcop.theta12 <- (((exp(-theta)-1)*((exp(-theta*F2)-1)*exp(-theta*F1)*(-F1)+exp(-theta*F1)*exp(-theta*F2)*(-F2))+(exp(-theta*F2)-1)*exp(-theta*F1)*exp(-theta))/(exp(-theta)-1)^2*(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))-(exp(-theta*F2)-1)/(exp(-theta)-1)*exp(-theta*F1)*((exp(-theta)-1)*((exp(-theta*F2)-1)*exp(-theta*F1)*(-F1)+(exp(-theta*F1)-1)*exp(-theta*F2)*(-F2))+(exp(-theta*F2)-1)*(exp(-theta*F1)-1)*exp(-theta))/(exp(-theta)-1)^2 )/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))^2
    d2Cdcop.theta22 <- (((exp(-theta)-1)*((exp(-theta*F22)-1)*exp(-theta*F1)*(-F1)+exp(-theta*F1)*exp(-theta*F22)*(-F22))+(exp(-theta*F22)-1)*exp(-theta*F1)*exp(-theta))/(exp(-theta)-1)^2*(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))-(exp(-theta*F22)-1)/(exp(-theta)-1)*exp(-theta*F1)*((exp(-theta)-1)*((exp(-theta*F22)-1)*exp(-theta*F1)*(-F1)+(exp(-theta*F1)-1)*exp(-theta*F22)*(-F22))+(exp(-theta*F22)-1)*(exp(-theta*F1)-1)*exp(-theta))/(exp(-theta)-1)^2 )/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))^2
    
    d2Cdcop1dcop.theta1 <- (((exp(-theta)-1)*((exp(-theta*F1)-1)*exp(-theta*F2)*(-F2)+exp(-theta*F2)*exp(-theta*F1)*(-F1))+(exp(-theta*F1)-1)*exp(-theta*F2)*exp(-theta))/(exp(-theta)-1)^2*(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))-(exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F2)*((exp(-theta)-1)*((exp(-theta*F2)-1)*exp(-theta*F1)*(-F1)+(exp(-theta*F1)-1)*exp(-theta*F2)*(-F2))+(exp(-theta*F2)-1)*(exp(-theta*F1)-1)*exp(-theta))/(exp(-theta)-1)^2 )/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))^2
    d2Cdcop11dcop.theta2 <- (((exp(-theta)-1)*((exp(-theta*F1)-1)*exp(-theta*F22)*(-F22)+exp(-theta*F22)*exp(-theta*F1)*(-F1))+(exp(-theta*F1)-1)*exp(-theta*F22)*exp(-theta))/(exp(-theta)-1)^2*(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))-(exp(-theta*F1)-1)/(exp(-theta)-1)*exp(-theta*F22)*((exp(-theta)-1)*((exp(-theta*F22)-1)*exp(-theta*F1)*(-F1)+(exp(-theta*F1)-1)*exp(-theta*F22)*(-F22))+(exp(-theta*F22)-1)*(exp(-theta*F1)-1)*exp(-theta))/(exp(-theta)-1)^2 )/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))^2
    
    d2Cdcop2dcop.theta1 <- -2/theta^3*log(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))+2/theta^2*1/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))*((-exp(-theta*F1)*F1*(exp(-theta*F2)-1)-exp(-theta*F2)*F2*(exp(-theta*F1)-1))*(exp(-theta)-1)+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)*exp(-theta))/(exp(-theta)-1)^2-1/theta*(-1/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))^2*((exp(-theta)-1)*(exp(-theta*F1)*(-F1)*(exp(-theta*F2)-1)+(exp(-theta*F1)-1)*exp(-theta*F2)*(-F2))+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)*exp(-theta) )^2/(exp(-theta)-1)^4+1/(1+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)/(exp(-theta)-1))*( (exp(-theta)-1)^2*((exp(-theta*F1)*F1^2*(exp(-theta*F2)-1)+exp(-theta*F1)*F1*exp(-theta*F2)*F2+exp(-theta*F2)*F2^2*(exp(-theta*F1)-1)+exp(-theta*F2)*F2*exp(-theta*F1)*F1)*(exp(-theta)-1)-(-exp(-theta*F1)*F1*(exp(-theta*F2)-1)-exp(-theta*F2)*F2*(exp(-theta*F1)-1))*exp(-theta)+exp(-theta*F1-theta*F2-theta)*(-F1-F2-1)-exp(-theta*F1-theta)*(-F1-1)-exp(-theta*F2-theta)*(-F2-1)-exp(-theta))+2*(exp(-theta)-1)*exp(-theta)*((-exp(-theta*F1)*F1*(exp(-theta*F2)-1)-exp(-theta*F2)*F2*(exp(-theta*F1)-1))*(exp(-theta)-1)+(exp(-theta*F1)-1)*(exp(-theta*F2)-1)*exp(-theta) ) )/(exp(-theta)-1)^4         )  
    d2Cdcop22dcop.theta2 <- -2/theta^3*log(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))+2/theta^2*1/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))*((-exp(-theta*F1)*F1*(exp(-theta*F22)-1)-exp(-theta*F22)*F22*(exp(-theta*F1)-1))*(exp(-theta)-1)+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)*exp(-theta))/(exp(-theta)-1)^2-1/theta*(-1/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))^2*((exp(-theta)-1)*(exp(-theta*F1)*(-F1)*(exp(-theta*F22)-1)+(exp(-theta*F1)-1)*exp(-theta*F22)*(-F22))+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)*exp(-theta) )^2/(exp(-theta)-1)^4+1/(1+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)/(exp(-theta)-1))*( (exp(-theta)-1)^2*((exp(-theta*F1)*F1^2*(exp(-theta*F22)-1)+exp(-theta*F1)*F1*exp(-theta*F22)*F22+exp(-theta*F22)*F22^2*(exp(-theta*F1)-1)+exp(-theta*F22)*F22*exp(-theta*F1)*F1)*(exp(-theta)-1)-(-exp(-theta*F1)*F1*(exp(-theta*F22)-1)-exp(-theta*F22)*F22*(exp(-theta*F1)-1))*exp(-theta)+exp(-theta*F1-theta*F22-theta)*(-F1-F22-1)-exp(-theta*F1-theta)*(-F1-1)-exp(-theta*F22-theta)*(-F22-1)-exp(-theta))+2*(exp(-theta)-1)*exp(-theta)*((-exp(-theta*F1)*F1*(exp(-theta*F22)-1)-exp(-theta*F22)*F22*(exp(-theta*F1)-1))*(exp(-theta)-1)+(exp(-theta*F1)-1)*(exp(-theta*F22)-1)*exp(-theta) ) )/(exp(-theta)-1)^4    )  
    
    
    
  } else { if (VC$BivD=="G0") {
    theta <- 1 + exp(theta.star)
    Gumbel1 <- exp(-((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta))
    Gumbel2 <- exp(-((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta))
    
    Gumbel1 <- ifelse(Gumbel1>precision, Gumbel1, precision)
    Gumbel1 <- ifelse(Gumbel1<(1-precision), Gumbel1, 1-precision)
    Gumbel2 <- ifelse(Gumbel2>precision, Gumbel2, precision)
    Gumbel2 <- ifelse(Gumbel2<(1-precision), Gumbel2, 1-precision)
    lx <- f2 - Gumbel1 + Gumbel2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- exp(-((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta))*((-1/theta)*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1) 
    dcop11 <- exp(-((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta))*((-1/theta)*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1)  
    
    dcop2 <- exp(-((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta))*((-1/theta)*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(-log(F2))^(theta-1))*(-1/F2)
    dcop22 <- exp(-((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta))*((-1/theta)*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(-log(F22))^(theta-1))*(-1/F22) 
    
    dcop.theta1 <- -exp(-((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta))*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta)*((-1/theta^2)*log((-log(F1))^(theta)+(-log(F2))^(theta)))+(1/theta)*(1/((-log(F1))^(theta)+(-log(F2))^(theta)))*((-log(F1))^(theta)*log(-log(F1))+(-log(F2))^(theta)*log(-log(F2)))*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta)*(-exp(-((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta)))
    dcop.theta2 <- -exp(-((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta))*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta)*((-1/theta^2)*log((-log(F1))^(theta)+(-log(F22))^(theta)))+(1/theta)*(1/((-log(F1))^(theta)+(-log(F22))^(theta)))*((-log(F1))^(theta)*log(-log(F1))+(-log(F22))^(theta)*log(-log(F22)))*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta)*(-exp(-((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta)))
    
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- dcop1*((-1/theta)*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1)+exp(-((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-2)*theta^2*(-log(F1))^(2*theta-2)*(-1/F1)^2+((-1/theta)*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log(F1))^(theta-2)*(-1/F1)^2+theta*(-log(F1))^(theta-1)*(-1/F1)^2) )
    d2Cdcop112 <- dcop11*((-1/theta)*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1)+exp(-((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-2)*theta^2*(-log(F1))^(2*theta-2)*(-1/F1)^2+((-1/theta)*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log(F1))^(theta-2)*(-1/F1)^2+theta*(-log(F1))^(theta-1)*(-1/F1)^2) )
    
    d2Cdcop22 <- dcop2*((-1/theta)*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(-log(F2))^(theta-1))*(-1/F2)+exp(-((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-2)*theta^2*(-log(F2))^(2*theta-2)*(-1/F2)^2+((-1/theta)*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log(F2))^(theta-2)*(-1/F2)^2+theta*(-log(F2))^(theta-1)*(-1/F2)^2) )
    d2Cdcop222 <- dcop22*((-1/theta)*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(-log(F22))^(theta-1))*(-1/F22)+exp(-((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-2)*theta^2*(-log(F22))^(2*theta-2)*(-1/F22)^2+((-1/theta)*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log(F22))^(theta-2)*(-1/F22)^2+theta*(-log(F22))^(theta-1)*(-1/F22)^2) )
    
    d2Cdcop1dcop2 <- dcop2*((-1/theta)*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1)+exp(-((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-2)*theta*(-log(F2))^(theta-1)*(-1/F2))*theta*(-log(F1))^(theta-1)*(-1/F1)
    d2Cdcop11dcop22 <- dcop22*((-1/theta)*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1)+exp(-((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-2)*theta*(-log(F22))^(theta-1)*(-1/F22))*theta*(-log(F1))^(theta-1)*(-1/F1)
    
    d2Cdcop.theta12 <- dcop.theta1*(-1/theta*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1)+(-1/F1)*exp(-((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta))*((1/theta^2*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-1)-1/theta*(dcop.theta1/(-exp(-((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta)))*((-log(F1))^(theta)+(-log(F2))^(theta))^(-1)-((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-2)*((-log(F1))^(theta)*log(-log(F1))+(-log(F2))^(theta)*log(-log(F2)))))*(theta*(-log(F1))^(theta-1))+(-1/theta)*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-1)*((-log(F1))^(theta-1)+theta*(-log(F1))^(theta-1)*(log(-log(F1)))) )
    d2Cdcop.theta22 <- dcop.theta2*(-1/theta*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1)+(-1/F1)*exp(-((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta))*((1/theta^2*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-1)-1/theta*(dcop.theta2/(-exp(-((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta)))*((-log(F1))^(theta)+(-log(F22))^(theta))^(-1)-((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-2)*((-log(F1))^(theta)*log(-log(F1))+(-log(F22))^(theta)*log(-log(F22)))))*(theta*(-log(F1))^(theta-1))+(-1/theta)*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-1)*((-log(F1))^(theta-1)+theta*(-log(F1))^(theta-1)*(log(-log(F1)))) )
    
    d2Cdcop1dcop.theta1 <- dcop.theta1*(-1/theta*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(-log(F2))^(theta-1))*(-1/F2)+(-1/F2)*exp(-((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta))*((1/theta^2*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-1)-1/theta*(dcop.theta1/(-exp(-((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta)))*((-log(F1))^(theta)+(-log(F2))^(theta))^(-1)-((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-2)*((-log(F1))^(theta)*log(-log(F1))+(-log(F2))^(theta)*log(-log(F2)))))*(theta*(-log(F2))^(theta-1))+(-1/theta)*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta-1)*((-log(F2))^(theta-1)+theta*(-log(F2))^(theta-1)*(log(-log(F2)))) )
    d2Cdcop11dcop.theta2 <- dcop.theta2*(-1/theta*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(-log(F22))^(theta-1))*(-1/F22)+(-1/F22)*exp(-((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta))*((1/theta^2*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-1)-1/theta*(dcop.theta2/(-exp(-((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta)))*((-log(F1))^(theta)+(-log(F22))^(theta))^(-1)-((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-2)*((-log(F1))^(theta)*log(-log(F1))+(-log(F22))^(theta)*log(-log(F22)))))*(theta*(-log(F22))^(theta-1))+(-1/theta)*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta-1)*((-log(F22))^(theta-1)+theta*(-log(F22))^(theta-1)*(log(-log(F22)))) )
    
    d2Cdcop2dcop.theta1 <- -dcop.theta1*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta)*(-1/theta^2*log(((-log(F1))^(theta)+(-log(F2))^(theta))))-exp(-((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta))*(dcop.theta1/(-exp(-((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta)))*(-1/theta^2)*log(((-log(F1))^(theta)+(-log(F2))^(theta)))+((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta)*(2/theta^3*log(((-log(F1))^(theta)+(-log(F2))^(theta)))-1/theta^2*((-log(F1))^(theta)+(-log(F2))^(theta))^(-1)*((-log(F1))^(theta)*log(-log(F1))+(-log(F2))^(theta)*log(-log(F2)))))+(-1/theta^2*((-log(F1))^(theta)+(-log(F2))^(theta))^(-1)-1/theta*((-log(F1))^(theta)+(-log(F2))^(theta))^(-2)*((-log(F1))^(theta)*log(-log(F1))+(-log(F2))^(theta)*log(-log(F2))))*((-log(F1))^(theta)*log(-log(F1))+(-log(F2))^(theta)*log(-log(F2)))*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta)*(-Gumbel1)+(-Gumbel1*(((-log(F1))^(theta)*log(-log(F1))^2+(-log(F2))^(theta)*log(-log(F2))^2)*((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta)+dcop.theta1/(-exp(-((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta)))*((-log(F1))^(theta)*log(-log(F1))+(-log(F2))^(theta)*log(-log(F2))))+((-log(F1))^(theta)+(-log(F2))^(theta))^(1/theta)*((-log(F1))^(theta)*log(-log(F1))+(-log(F2))^(theta)*log(-log(F2)))*(-dcop.theta1))*1/theta*((-log(F1))^(theta)+(-log(F2))^(theta))^(-1)
    d2Cdcop22dcop.theta2 <- -dcop.theta2*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta)*(-1/theta^2*log(((-log(F1))^(theta)+(-log(F22))^(theta))))-exp(-((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta))*(dcop.theta2/(-exp(-((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta)))*(-1/theta^2)*log(((-log(F1))^(theta)+(-log(F22))^(theta)))+((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta)*(2/theta^3*log(((-log(F1))^(theta)+(-log(F22))^(theta)))-1/theta^2*((-log(F1))^(theta)+(-log(F22))^(theta))^(-1)*((-log(F1))^(theta)*log(-log(F1))+(-log(F22))^(theta)*log(-log(F22)))))+(-1/theta^2*((-log(F1))^(theta)+(-log(F22))^(theta))^(-1)-1/theta*((-log(F1))^(theta)+(-log(F22))^(theta))^(-2)*((-log(F1))^(theta)*log(-log(F1))+(-log(F22))^(theta)*log(-log(F22))))*((-log(F1))^(theta)*log(-log(F1))+(-log(F22))^(theta)*log(-log(F22)))*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta)*(-Gumbel2)+(-Gumbel2*(((-log(F1))^(theta)*log(-log(F1))^2+(-log(F22))^(theta)*log(-log(F22))^2)*((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta)+dcop.theta2/(-exp(-((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta)))*((-log(F1))^(theta)*log(-log(F1))+(-log(F22))^(theta)*log(-log(F22))))+((-log(F1))^(theta)+(-log(F22))^(theta))^(1/theta)*((-log(F1))^(theta)*log(-log(F1))+(-log(F22))^(theta)*log(-log(F22)))*(-dcop.theta2))*1/theta*((-log(F1))^(theta)+(-log(F22))^(theta))^(-1)
    
    
  } else { if (VC$BivD=="G90") {
    theta<- 1 + exp(theta.star)
    Gumbel1 <- F2 - exp(-((-log(1-F1))^(theta)+(-log(F2))^(theta))^(1/theta))
    Gumbel2 <- F22 - exp(-((-log(1-F1))^(theta)+(-log(F22))^(theta))^(1/theta))
    Gumbel1 <- ifelse(Gumbel1>precision, Gumbel1, precision)
    Gumbel1 <- ifelse(Gumbel1<(1-precision), Gumbel1, 1-precision)
    Gumbel2 <- ifelse(Gumbel2>precision, Gumbel2, precision)
    Gumbel2 <- ifelse(Gumbel2<(1-precision), Gumbel2, 1-precision)
    lx <- f2 - Gumbel1 + Gumbel2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- exp(-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1)) 
    dcop11 <- exp(-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1))  
    
    dcop2 <- 1-exp(-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(-log(F2))^(theta-1))*(-1/F2)
    dcop22 <- 1-exp(-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(-log(F22))^(theta-1))*(-1/F22) 
    
    dcop.theta1 <-  ( -exp(-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta))*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta)*((-1/theta^2)*log((-log((1-F1)))^(theta)+(-log(F2))^(theta)))+(1/theta)*(1/((-log((1-F1)))^(theta)+(-log(F2))^(theta)))*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log(F2))^(theta)*log(-log(F2)))*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta)*(-exp(-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta)))  )
    dcop.theta2 <-  ( -exp(-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta))*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta)*((-1/theta^2)*log((-log((1-F1)))^(theta)+(-log(F22))^(theta)))+(1/theta)*(1/((-log((1-F1)))^(theta)+(-log(F22))^(theta)))*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log(F22))^(theta)*log(-log(F22)))*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta)*(-exp(-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta)))  )
    
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- -( dcop1*((-1/theta)*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1))+exp(-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-2)*theta^2*(-log((1-F1)))^(2*theta-2)*(-1/(1-F1))^2+((-1/theta)*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log((1-F1)))^(theta-2)*(-1/(1-F1))^2+theta*(-log((1-F1)))^(theta-1)*(-1/(1-F1))^2) )      )
    d2Cdcop112 <- -( dcop11*((-1/theta)*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1))+exp(-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-2)*theta^2*(-log((1-F1)))^(2*theta-2)*(-1/(1-F1))^2+((-1/theta)*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log((1-F1)))^(theta-2)*(-1/(1-F1))^2+theta*(-log((1-F1)))^(theta-1)*(-1/(1-F1))^2) )      )
    
    d2Cdcop22 <- -( exp(-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(-log(F2))^(theta-1))*(-1/F2)*((-1/theta)*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(-log(F2))^(theta-1))*(-1/F2)+exp(-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-2)*theta^2*(-log(F2))^(2*theta-2)*(-1/F2)^2+((-1/theta)*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log(F2))^(theta-2)*(-1/F2)^2+theta*(-log(F2))^(theta-1)*(-1/F2)^2) )     )
    d2Cdcop222 <- -(  exp(-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(-log(F22))^(theta-1))*(-1/F22)*((-1/theta)*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(-log(F22))^(theta-1))*(-1/F22)+exp(-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-2)*theta^2*(-log(F22))^(2*theta-2)*(-1/F22)^2+((-1/theta)*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log(F22))^(theta-2)*(-1/F22)^2+theta*(-log(F22))^(theta-1)*(-1/F22)^2) )    )
    
    d2Cdcop1dcop2 <- exp(-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(-log(F2))^(theta-1))*(-1/F2)*((-1/theta)*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1))+exp(-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-2)*theta*(-log(F2))^(theta-1)*(-1/F2))*theta*(-log((1-F1)))^(theta-1)*(-1/(1-F1))
    d2Cdcop11dcop22 <- exp(-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(-log(F22))^(theta-1))*(-1/F22)*((-1/theta)*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1))+exp(-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-2)*theta*(-log(F22))^(theta-1)*(-1/F22))*theta*(-log((1-F1)))^(theta-1)*(-1/(1-F1))
    
    d2Cdcop.theta12 <- -(  dcop.theta1*(-1/theta*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1))+(-1/(1-F1))*exp(-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta))*((1/theta^2*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-1)-1/theta*(dcop.theta1/(-exp(-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta)))*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(-1)-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-2)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log(F2))^(theta)*log(-log(F2)))))*(theta*(-log((1-F1)))^(theta-1))+(-1/theta)*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-1)*((-log((1-F1)))^(theta-1)+theta*(-log((1-F1)))^(theta-1)*(log(-log((1-F1))))) )     )
    d2Cdcop.theta22 <- -(   dcop.theta2*(-1/theta*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1))+(-1/(1-F1))*exp(-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta))*((1/theta^2*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-1)-1/theta*(dcop.theta2/(-exp(-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta)))*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(-1)-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-2)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log(F22))^(theta)*log(-log(F22)))))*(theta*(-log((1-F1)))^(theta-1))+(-1/theta)*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-1)*((-log((1-F1)))^(theta-1)+theta*(-log((1-F1)))^(theta-1)*(log(-log((1-F1))))) )     )
    
    d2Cdcop1dcop.theta1 <- ( dcop.theta1*(-1/theta*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-1))*(theta*(-log(F2))^(theta-1))*(-1/F2)+(-1/F2)*exp(-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta))*((1/theta^2*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-1)-1/theta*(dcop.theta1/(-exp(-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta)))*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(-1)-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-2)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log(F2))^(theta)*log(-log(F2)))))*(theta*(-log(F2))^(theta-1))+(-1/theta)*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta-1)*((-log(F2))^(theta-1)+theta*(-log(F2))^(theta-1)*(log(-log(F2)))) )        )
    d2Cdcop11dcop.theta2 <- ( dcop.theta2*(-1/theta*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-1))*(theta*(-log(F22))^(theta-1))*(-1/F22)+(-1/F22)*exp(-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta))*((1/theta^2*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-1)-1/theta*(dcop.theta2/(-exp(-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta)))*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(-1)-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-2)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log(F22))^(theta)*log(-log(F22)))))*(theta*(-log(F22))^(theta-1))+(-1/theta)*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta-1)*((-log(F22))^(theta-1)+theta*(-log(F22))^(theta-1)*(log(-log(F22)))) )       )
    
    d2Cdcop2dcop.theta1 <- -(  -dcop.theta1*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta)*(-1/theta^2*log(((-log((1-F1)))^(theta)+(-log(F2))^(theta))))-exp(-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta))*(dcop.theta1/(-exp(-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta)))*(-1/theta^2)*log(((-log((1-F1)))^(theta)+(-log(F2))^(theta)))+((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta)*(2/theta^3*log(((-log((1-F1)))^(theta)+(-log(F2))^(theta)))-1/theta^2*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(-1)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log(F2))^(theta)*log(-log(F2)))))+(-1/theta^2*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(-1)-1/theta*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(-2)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log(F2))^(theta)*log(-log(F2))))*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log(F2))^(theta)*log(-log(F2)))*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta)*(-exp(-((-log(1-F1))^(theta)+(-log(F2))^(theta))^(1/theta)))+(-exp(-((-log(1-F1))^(theta)+(-log(F2))^(theta))^(1/theta))*(((-log((1-F1)))^(theta)*log(-log((1-F1)))^2+(-log(F2))^(theta)*log(-log(F2))^2)*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta)+dcop.theta1/(-exp(-((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta)))*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log(F2))^(theta)*log(-log(F2))))+((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(1/theta)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log(F2))^(theta)*log(-log(F2)))*(-dcop.theta1))*1/theta*((-log((1-F1)))^(theta)+(-log(F2))^(theta))^(-1)      )
    d2Cdcop22dcop.theta2 <-  -( -dcop.theta2*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta)*(-1/theta^2*log(((-log((1-F1)))^(theta)+(-log(F22))^(theta))))-exp(-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta))*(dcop.theta2/(-exp(-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta)))*(-1/theta^2)*log(((-log((1-F1)))^(theta)+(-log(F22))^(theta)))+((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta)*(2/theta^3*log(((-log((1-F1)))^(theta)+(-log(F22))^(theta)))-1/theta^2*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(-1)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log(F22))^(theta)*log(-log(F22)))))+(-1/theta^2*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(-1)-1/theta*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(-2)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log(F22))^(theta)*log(-log(F22))))*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log(F22))^(theta)*log(-log(F22)))*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta)*(-exp(-((-log(1-F1))^(theta)+(-log(F22))^(theta))^(1/theta)))+(-exp(-((-log(1-F1))^(theta)+(-log(F22))^(theta))^(1/theta))*(((-log((1-F1)))^(theta)*log(-log((1-F1)))^2+(-log(F22))^(theta)*log(-log(F22))^2)*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta)+dcop.theta2/(-exp(-((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta)))*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log(F22))^(theta)*log(-log(F22))))+((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(1/theta)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log(F22))^(theta)*log(-log(F22)))*(-dcop.theta2))*1/theta*((-log((1-F1)))^(theta)+(-log(F22))^(theta))^(-1)        )
    
    
    
  } else { if (VC$BivD=="G180") {
    theta <- 1 + exp(theta.star)
    Gumbel1 <- F1 + F2 - 1 + exp(-((-log(1-F1))^(theta)+(-log(1-F2))^(theta))^(1/theta))
    Gumbel2 <- F1 + F22 - 1 +  exp(-((-log(1-F1))^(theta)+(-log(1-F22))^(theta))^(1/theta))
    Gumbel1 <- ifelse(Gumbel1>precision, Gumbel1, precision)
    Gumbel1 <- ifelse(Gumbel1<(1-precision), Gumbel1, 1-precision)
    Gumbel2 <- ifelse(Gumbel2>precision, Gumbel2, precision)
    Gumbel2 <- ifelse(Gumbel2<(1-precision), Gumbel2, 1-precision)
    lx <- f2 - Gumbel1 + Gumbel2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- 1-exp(-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1)) 
    dcop11 <- 1-exp(-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1))  
    
    dcop2 <- 1-exp(-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log((1-F2)))^(theta-1))*(-1/(1-F2))
    dcop22 <- 1-exp(-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log((1-F22)))^(theta-1))*(-1/(1-F22)) 
    
    dcop.theta1 <- ( -exp(-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta)*((-1/theta^2)*log((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta)))+(1/theta)*(1/((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta)))*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log((1-F2)))^(theta)*log(-log((1-F2))))*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta)*(-exp(-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta)))  )
    dcop.theta2 <- ( -exp(-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta)*((-1/theta^2)*log((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta)))+(1/theta)*(1/((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta)))*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log((1-F22)))^(theta)*log(-log((1-F22))))*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta)*(-exp(-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta)))  )
    
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- ( exp(-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1))*((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1))+exp(-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-2)*theta^2*(-log((1-F1)))^(2*theta-2)*(-1/(1-F1))^2+((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log((1-F1)))^(theta-2)*(-1/(1-F1))^2+theta*(-log((1-F1)))^(theta-1)*(-1/(1-F1))^2) )      )
    d2Cdcop112 <- ( exp(-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1))*((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1))+exp(-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-2)*theta^2*(-log((1-F1)))^(2*theta-2)*(-1/(1-F1))^2+((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log((1-F1)))^(theta-2)*(-1/(1-F1))^2+theta*(-log((1-F1)))^(theta-1)*(-1/(1-F1))^2) )      )
    
    d2Cdcop22 <- ( exp(-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log((1-F2)))^(theta-1))*(-1/(1-F2))*((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log((1-F2)))^(theta-1))*(-1/(1-F2))+exp(-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-2)*theta^2*(-log((1-F2)))^(2*theta-2)*(-1/(1-F2))^2+((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log((1-F2)))^(theta-2)*(-1/(1-F2))^2+theta*(-log((1-F2)))^(theta-1)*(-1/(1-F2))^2) )     )
    d2Cdcop222 <- (  exp(-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log((1-F22)))^(theta-1))*(-1/(1-F22))*((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log((1-F22)))^(theta-1))*(-1/(1-F22))+exp(-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-2)*theta^2*(-log((1-F22)))^(2*theta-2)*(-1/(1-F22))^2+((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log((1-F22)))^(theta-2)*(-1/(1-F22))^2+theta*(-log((1-F22)))^(theta-1)*(-1/(1-F22))^2) )    )
    
    d2Cdcop1dcop2 <- ( exp(-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log((1-F2)))^(theta-1))*(-1/(1-F2))*((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1))+exp(-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-2)*theta*(-log((1-F2)))^(theta-1)*(-1/(1-F2)))*theta*(-log((1-F1)))^(theta-1)*(-1/(1-F1))     )
    d2Cdcop11dcop22 <- ( exp(-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log((1-F22)))^(theta-1))*(-1/(1-F22))*((-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1))+exp(-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-2)*theta*(-log((1-F22)))^(theta-1)*(-1/(1-F22)))*theta*(-log((1-F1)))^(theta-1)*(-1/(1-F1))      )
    
    d2Cdcop.theta12 <- -( dcop.theta1*(-1/theta*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1))+(-1/(1-F1))*exp(-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*((1/theta^2*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1)-1/theta*(dcop.theta1/(-exp(-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta)))*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(-1)-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-2)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log((1-F2)))^(theta)*log(-log((1-F2))))))*(theta*(-log((1-F1)))^(theta-1))+(-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1)*((-log((1-F1)))^(theta-1)+theta*(-log((1-F1)))^(theta-1)*(log(-log((1-F1))))) )   )
    d2Cdcop.theta22 <-  -( dcop.theta2*(-1/theta*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log((1-F1)))^(theta-1))*(-1/(1-F1))+(-1/(1-F1))*exp(-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*((1/theta^2*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1)-1/theta*(dcop.theta2/(-exp(-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta)))*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(-1)-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-2)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log((1-F22)))^(theta)*log(-log((1-F22))))))*(theta*(-log((1-F1)))^(theta-1))+(-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1)*((-log((1-F1)))^(theta-1)+theta*(-log((1-F1)))^(theta-1)*(log(-log((1-F1))))) )    )
    
    d2Cdcop1dcop.theta1 <- -( dcop.theta1*(-1/theta*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log((1-F2)))^(theta-1))*(-1/(1-F2))+(-1/(1-F2))*exp(-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*((1/theta^2*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1)-1/theta*(dcop.theta1/(-exp(-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta)))*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(-1)-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-2)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log((1-F2)))^(theta)*log(-log((1-F2))))))*(theta*(-log((1-F2)))^(theta-1))+(-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1)*((-log((1-F2)))^(theta-1)+theta*(-log((1-F2)))^(theta-1)*(log(-log((1-F2))))) )    )    
    d2Cdcop11dcop.theta2 <- -( dcop.theta2*(-1/theta*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log((1-F22)))^(theta-1))*(-1/(1-F22))+(-1/(1-F22))*exp(-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*((1/theta^2*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1)-1/theta*(dcop.theta2/(-exp(-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta)))*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(-1)-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-2)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log((1-F22)))^(theta)*log(-log((1-F22))))))*(theta*(-log((1-F22)))^(theta-1))+(-1/theta)*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1)*((-log((1-F22)))^(theta-1)+theta*(-log((1-F22)))^(theta-1)*(log(-log((1-F22))))) )     )   
    
    d2Cdcop2dcop.theta1 <- (  -dcop.theta1*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta)*(-1/theta^2*log(((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))))-exp(-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*(dcop.theta1/(-exp(-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta)))*(-1/theta^2)*log(((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta)))+((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta)*(2/theta^3*log(((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta)))-1/theta^2*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(-1)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log((1-F2)))^(theta)*log(-log((1-F2))))))+(-1/theta^2*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(-1)-1/theta*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(-2)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log((1-F2)))^(theta)*log(-log((1-F2)))))*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log((1-F2)))^(theta)*log(-log((1-F2))))*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta)*(-exp(-((-log(1-F1))^(theta)+(-log(1-F2))^(theta))^(1/theta)))+(-exp(-((-log(1-F1))^(theta)+(-log(1-F2))^(theta))^(1/theta))*(((-log((1-F1)))^(theta)*log(-log((1-F1)))^2+(-log((1-F2)))^(theta)*log(-log((1-F2)))^2)*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta)+dcop.theta1/(-exp(-((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta)))*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log((1-F2)))^(theta)*log(-log((1-F2)))))+((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(1/theta)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log((1-F2)))^(theta)*log(-log((1-F2))))*(-dcop.theta1))*1/theta*((-log((1-F1)))^(theta)+(-log((1-F2)))^(theta))^(-1)      )
    d2Cdcop22dcop.theta2 <-  ( -dcop.theta2*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta)*(-1/theta^2*log(((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))))-exp(-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*(dcop.theta2/(-exp(-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta)))*(-1/theta^2)*log(((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta)))+((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta)*(2/theta^3*log(((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta)))-1/theta^2*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(-1)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log((1-F22)))^(theta)*log(-log((1-F22))))))+(-1/theta^2*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(-1)-1/theta*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(-2)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log((1-F22)))^(theta)*log(-log((1-F22)))))*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log((1-F22)))^(theta)*log(-log((1-F22))))*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta)*(-exp(-((-log(1-F1))^(theta)+(-log(1-F22))^(theta))^(1/theta)))+(-exp(-((-log(1-F1))^(theta)+(-log(1-F22))^(theta))^(1/theta))*(((-log((1-F1)))^(theta)*log(-log((1-F1)))^2+(-log((1-F22)))^(theta)*log(-log((1-F22)))^2)*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta)+dcop.theta2/(-exp(-((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta)))*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log((1-F22)))^(theta)*log(-log((1-F22)))))+((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(1/theta)*((-log((1-F1)))^(theta)*log(-log((1-F1)))+(-log((1-F22)))^(theta)*log(-log((1-F22))))*(-dcop.theta2))*1/theta*((-log((1-F1)))^(theta)+(-log((1-F22)))^(theta))^(-1)        )
    
    
    
    
    
  } else { if (VC$BivD=="G270") {
    theta <- 1 + exp(theta.star)
    Gumbel1 <- F1 - exp(-((-log(F1))^(theta)+(-log(1-F2))^(theta))^(1/theta))
    Gumbel2 <- F1 - exp(-((-log(F1))^(theta)+(-log(1-F22))^(theta))^(1/theta))
    Gumbel1 <- ifelse(Gumbel1>precision, Gumbel1, precision)
    Gumbel1 <- ifelse(Gumbel1<(1-precision), Gumbel1, 1-precision)
    Gumbel2 <- ifelse(Gumbel2>precision, Gumbel2, precision)
    Gumbel2 <- ifelse(Gumbel2<(1-precision), Gumbel2, 1-precision)
    lx <- f2 - Gumbel1 + Gumbel2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- 1-exp(-((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*((-1/theta)*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1) 
    dcop11 <- 1-exp(-((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*((-1/theta)*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1)  
    
    dcop2 <- exp(-((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*((-1/theta)*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log((1-F2)))^(theta-1))*(-1/(1-F2))
    dcop22 <- exp(-((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*((-1/theta)*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log((1-F22)))^(theta-1))*(-1/(1-F22)) 
    
    dcop.theta1 <- ( -exp(-((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta)*((-1/theta^2)*log((-log(F1))^(theta)+(-log((1-F2)))^(theta)))+(1/theta)*(1/((-log(F1))^(theta)+(-log((1-F2)))^(theta)))*((-log(F1))^(theta)*log(-log(F1))+(-log((1-F2)))^(theta)*log(-log((1-F2))))*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta)*(-exp(-((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta)))   )
    dcop.theta2 <- ( -exp(-((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta)*((-1/theta^2)*log((-log(F1))^(theta)+(-log((1-F22)))^(theta)))+(1/theta)*(1/((-log(F1))^(theta)+(-log((1-F22)))^(theta)))*((-log(F1))^(theta)*log(-log(F1))+(-log((1-F22)))^(theta)*log(-log((1-F22))))*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta)*(-exp(-((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta)))  )
    
    
    # Hessian derivative components
    
    d2Cdcop12 <-  -( exp(-((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*((-1/theta)*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1)*((-1/theta)*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1)+exp(-((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-2)*theta^2*(-log(F1))^(2*theta-2)*(-1/F1)^2+((-1/theta)*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log(F1))^(theta-2)*(-1/F1)^2+theta*(-log(F1))^(theta-1)*(-1/F1)^2) )    )
    d2Cdcop112 <-  -(   exp(-((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*((-1/theta)*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1)*((-1/theta)*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1)+exp(-((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-2)*theta^2*(-log(F1))^(2*theta-2)*(-1/F1)^2+((-1/theta)*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log(F1))^(theta-2)*(-1/F1)^2+theta*(-log(F1))^(theta-1)*(-1/F1)^2) )     )
    
    d2Cdcop22 <- -( dcop2*((-1/theta)*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log((1-F2)))^(theta-1))*(-1/(1-F2))+exp(-((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-2)*theta^2*(-log((1-F2)))^(2*theta-2)*(-1/(1-F2))^2+((-1/theta)*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log((1-F2)))^(theta-2)*(-1/(1-F2))^2+theta*(-log((1-F2)))^(theta-1)*(-1/(1-F2))^2) )    )
    d2Cdcop222 <-  -(  dcop22*((-1/theta)*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log((1-F22)))^(theta-1))*(-1/(1-F22))+exp(-((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-2)*theta^2*(-log((1-F22)))^(2*theta-2)*(-1/(1-F22))^2+((-1/theta)*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(theta-1)*(-log((1-F22)))^(theta-2)*(-1/(1-F22))^2+theta*(-log((1-F22)))^(theta-1)*(-1/(1-F22))^2) )    )
    
    d2Cdcop1dcop2 <- dcop2*((-1/theta)*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1)+exp(-((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-2)*theta*(-log((1-F2)))^(theta-1)*(-1/(1-F2)))*theta*(-log(F1))^(theta-1)*(-1/F1)
    d2Cdcop11dcop22 <- dcop22*((-1/theta)*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1)+exp(-((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*( (-1/theta)*(1/theta-1)*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-2)*theta*(-log((1-F22)))^(theta-1)*(-1/(1-F22)))*theta*(-log(F1))^(theta-1)*(-1/F1)
    
    d2Cdcop.theta12 <- ( dcop.theta1*(-1/theta*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1)+(-1/F1)*exp(-((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*((1/theta^2*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1)-1/theta*(dcop.theta1/(-exp(-((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta)))*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(-1)-((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-2)*((-log(F1))^(theta)*log(-log(F1))+(-log((1-F2)))^(theta)*log(-log((1-F2))))))*(theta*(-log(F1))^(theta-1))+(-1/theta)*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1)*((-log(F1))^(theta-1)+theta*(-log(F1))^(theta-1)*(log(-log(F1)))) )     )
    d2Cdcop.theta22 <-  ( dcop.theta2*(-1/theta*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log(F1))^(theta-1))*(-1/F1)+(-1/F1)*exp(-((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*((1/theta^2*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1)-1/theta*(dcop.theta2/(-exp(-((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta)))*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(-1)-((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-2)*((-log(F1))^(theta)*log(-log(F1))+(-log((1-F22)))^(theta)*log(-log((1-F22))))))*(theta*(-log(F1))^(theta-1))+(-1/theta)*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1)*((-log(F1))^(theta-1)+theta*(-log(F1))^(theta-1)*(log(-log(F1)))) )      )
    
    d2Cdcop1dcop.theta1 <- -( dcop.theta1*(-1/theta*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1))*(theta*(-log((1-F2)))^(theta-1))*(-1/(1-F2))+(-1/(1-F2))*exp(-((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*((1/theta^2*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1)-1/theta*(dcop.theta1/(-exp(-((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta)))*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(-1)-((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-2)*((-log(F1))^(theta)*log(-log(F1))+(-log((1-F2)))^(theta)*log(-log((1-F2))))))*(theta*(-log((1-F2)))^(theta-1))+(-1/theta)*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta-1)*((-log((1-F2)))^(theta-1)+theta*(-log((1-F2)))^(theta-1)*(log(-log((1-F2))))) )    )
    d2Cdcop11dcop.theta2 <- -(  dcop.theta2*(-1/theta*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1))*(theta*(-log((1-F22)))^(theta-1))*(-1/(1-F22))+(-1/(1-F22))*exp(-((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*((1/theta^2*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1)-1/theta*(dcop.theta2/(-exp(-((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta)))*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(-1)-((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-2)*((-log(F1))^(theta)*log(-log(F1))+(-log((1-F22)))^(theta)*log(-log((1-F22))))))*(theta*(-log((1-F22)))^(theta-1))+(-1/theta)*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta-1)*((-log((1-F22)))^(theta-1)+theta*(-log((1-F22)))^(theta-1)*(log(-log((1-F22))))) )    )
    
    d2Cdcop2dcop.theta1 <- -( -dcop.theta1*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta)*(-1/theta^2*log(((-log(F1))^(theta)+(-log((1-F2)))^(theta))))-exp(-((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta))*(dcop.theta1/(-exp(-((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta)))*(-1/theta^2)*log(((-log(F1))^(theta)+(-log((1-F2)))^(theta)))+((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta)*(2/theta^3*log(((-log(F1))^(theta)+(-log((1-F2)))^(theta)))-1/theta^2*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(-1)*((-log(F1))^(theta)*log(-log(F1))+(-log((1-F2)))^(theta)*log(-log((1-F2))))))+(-1/theta^2*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(-1)-1/theta*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(-2)*((-log(F1))^(theta)*log(-log(F1))+(-log((1-F2)))^(theta)*log(-log((1-F2)))))*((-log(F1))^(theta)*log(-log(F1))+(-log((1-F2)))^(theta)*log(-log((1-F2))))*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta)*(-exp(-((-log(F1))^(theta)+(-log(1-F2))^(theta))^(1/theta)))+(-exp(-((-log(F1))^(theta)+(-log(1-F2))^(theta))^(1/theta))*(((-log(F1))^(theta)*log(-log(F1))^2+(-log((1-F2)))^(theta)*log(-log((1-F2)))^2)*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta)+dcop.theta1/(-exp(-((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta)))*((-log(F1))^(theta)*log(-log(F1))+(-log((1-F2)))^(theta)*log(-log((1-F2)))))+((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(1/theta)*((-log(F1))^(theta)*log(-log(F1))+(-log((1-F2)))^(theta)*log(-log((1-F2))))*(-dcop.theta1))*1/theta*((-log(F1))^(theta)+(-log((1-F2)))^(theta))^(-1)   )
    d2Cdcop22dcop.theta2 <- -( -dcop.theta2*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta)*(-1/theta^2*log(((-log(F1))^(theta)+(-log((1-F22)))^(theta))))-exp(-((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta))*(dcop.theta2/(-exp(-((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta)))*(-1/theta^2)*log(((-log(F1))^(theta)+(-log((1-F22)))^(theta)))+((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta)*(2/theta^3*log(((-log(F1))^(theta)+(-log((1-F22)))^(theta)))-1/theta^2*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(-1)*((-log(F1))^(theta)*log(-log(F1))+(-log((1-F22)))^(theta)*log(-log((1-F22))))))+(-1/theta^2*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(-1)-1/theta*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(-2)*((-log(F1))^(theta)*log(-log(F1))+(-log((1-F22)))^(theta)*log(-log((1-F22)))))*((-log(F1))^(theta)*log(-log(F1))+(-log((1-F22)))^(theta)*log(-log((1-F22))))*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta)*(-exp(-((-log(F1))^(theta)+(-log(1-F22))^(theta))^(1/theta)))+(-exp(-((-log(F1))^(theta)+(-log(1-F22))^(theta))^(1/theta))*(((-log(F1))^(theta)*log(-log(F1))^2+(-log((1-F22)))^(theta)*log(-log((1-F22)))^2)*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta)+dcop.theta2/(-exp(-((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta)))*((-log(F1))^(theta)*log(-log(F1))+(-log((1-F22)))^(theta)*log(-log((1-F22)))))+((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(1/theta)*((-log(F1))^(theta)*log(-log(F1))+(-log((1-F22)))^(theta)*log(-log((1-F22))))*(-dcop.theta2))*1/theta*((-log(F1))^(theta)+(-log((1-F22)))^(theta))^(-1)    )
    
    
    
  } else { if (VC$BivD=="J0") {
    theta <- 1 + exp(theta.star) + precision
    Joe1 <- 1-((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta)
    Joe2 <- 1-((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta)
    
    Joe1 <- ifelse(Joe1>precision, Joe1, precision)
    Joe1 <- ifelse(Joe1<(1-precision), Joe1, 1-precision)
    Joe2 <- ifelse(Joe2>precision, Joe2, precision)
    Joe2 <- ifelse(Joe2<(1-precision), Joe2, 1-precision)
    lx <- f2 - Joe1 + Joe2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1<- (-1/theta)*((1-F1)^theta+(1-F2)^theta-(1-F1)^theta*(1-F2)^theta)^(1/theta-1)*(-theta*(1-F1)^(theta-1)+(1-F2)^theta*theta*(1-F1)^(theta-1))
    dcop11<- (-1/theta)*((1-F1)^theta+(1-F22)^theta-(1-F1)^theta*(1-F22)^theta)^(1/theta-1)*(-theta*(1-F1)^(theta-1)+(1-F22)^theta*theta*(1-F1)^(theta-1))   
    
    dcop2<- (-1/theta)*((1-F1)^theta+(1-F2)^theta-(1-F1)^theta*(1-F2)^theta)^(1/theta-1)*(-theta*(1-F2)^(theta-1)+(1-F1)^theta*theta*(1-F2)^(theta-1))
    dcop22<- (-1/theta)*((1-F1)^theta+(1-F22)^theta-(1-F1)^theta*(1-F22)^theta)^(1/theta-1)*(-theta*(1-F22)^(theta-1)+(1-F1)^theta*theta*(1-F22)^(theta-1)) 
    
    dcop.theta1<- -(-1/theta^2)*log((1-F1)^theta+(1-F2)^theta-(1-F1)^theta*(1-F2)^theta)*((1-F1)^theta+(1-F2)^theta-(1-F1)^theta*(1-F2)^theta)^(1/theta)-(1/theta)*((1-F1)^theta+(1-F2)^theta-(1-F1)^theta*(1-F2)^theta)^(1/theta-1)*((1-F1)^theta*log(1-F1)+(1-F2)^theta*log(1-F2)-(1-F1)^theta*log(1-F1)*(1-F2)^theta-(1-F2)^theta*log(1-F2)*(1-F1)^theta) 
    dcop.theta2<- -(-1/theta^2)*log((1-F1)^theta+(1-F22)^theta-(1-F1)^theta*(1-F22)^theta)*((1-F1)^theta+(1-F22)^theta-(1-F1)^theta*(1-F22)^theta)^(1/theta)-(1/theta)*((1-F1)^theta+(1-F22)^theta-(1-F1)^theta*(1-F22)^theta)^(1/theta-1)*((1-F1)^theta*log(1-F1)+(1-F22)^theta*log(1-F22)-(1-F1)^theta*log(1-F1)*(1-F22)^theta-(1-F22)^theta*log(1-F22)*(1-F1)^theta)
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- -1/theta*(1/theta-1)*((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta-2)*(-theta*(1-F1)^(theta-1)+(1-F2)^(theta)*theta*(1-F1)^(theta-1))^2-1/theta*((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta-1)*(theta*(theta-1)*(1-F1)^(theta-2)-(1-F2)^(theta)*theta*(theta-1)*(1-F1)^(theta-2))
    d2Cdcop112 <- -1/theta*(1/theta-1)*((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta-2)*(-theta*(1-F1)^(theta-1)+(1-F22)^(theta)*theta*(1-F1)^(theta-1))^2-1/theta*((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta-1)*(theta*(theta-1)*(1-F1)^(theta-2)-(1-F22)^(theta)*theta*(theta-1)*(1-F1)^(theta-2))
    
    d2Cdcop22 <- -1/theta*(1/theta-1)*((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta-2)*(-theta*(1-F2)^(theta-1)+(1-F1)^(theta)*theta*(1-F2)^(theta-1))^2-1/theta*((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta-1)*(theta*(theta-1)*(1-F2)^(theta-2)-(1-F1)^(theta)*theta*(theta-1)*(1-F2)^(theta-2))
    d2Cdcop222 <- -1/theta*(1/theta-1)*((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta-2)*(-theta*(1-F22)^(theta-1)+(1-F1)^(theta)*theta*(1-F22)^(theta-1))^2-1/theta*((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta-1)*(theta*(theta-1)*(1-F22)^(theta-2)-(1-F1)^(theta)*theta*(theta-1)*(1-F22)^(theta-2))
    
    d2Cdcop1dcop2 <- -1/theta*(1/theta-1)*((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta-2)*(-theta*(1-F2)^(theta-1)+(1-F1)^(theta)*theta*(1-F2)^(theta-1))*(-theta*(1-F1)^(theta-1)+(1-F2)^(theta)*theta*(1-F1)^(theta-1))-1/theta*((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta-1)*(-theta^2*(1-F2)^(theta-1)*(1-F1)^(theta-1))
    d2Cdcop11dcop22 <- -1/theta*(1/theta-1)*((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta-2)*(-theta*(1-F22)^(theta-1)+(1-F1)^(theta)*theta*(1-F22)^(theta-1))*(-theta*(1-F1)^(theta-1)+(1-F22)^(theta)*theta*(1-F1)^(theta-1))-1/theta*((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta-1)*(-theta^2*(1-F22)^(theta-1)*(1-F1)^(theta-1))
    
    d2Cdcop.theta12 <- 1/theta^2*((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta-1)*(-theta*(1-F1)^(theta-1)+(1-F2)^(theta)*theta*(1-F1)^(theta-1))-1/theta*((-theta*(1-F1)^(theta-1)+(1-F2)^(theta)*theta*(1-F1)^(theta-1))*(-dcop.theta1*((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(-1)-((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta-2)*((1-F1)^theta*log(1-F1)+(1-F2)^theta*log(1-F2)-(1-F1)^theta*log(1-F1)*(1-F2)^theta-(1-F2)^theta*log(1-F2)*(1-F1)^theta))+((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta-1)*(-(1-F1)^(theta-1)-theta*(1-F1)^(theta-1)*log(1-F1)+theta*(1-F2)^(theta)*log(1-F2)*(1-F1)^(theta-1)+((1-F1)^(theta-1)+theta*(1-F1)^(theta-1)*log(1-F1))*(1-F2)^theta))
    d2Cdcop.theta22 <- 1/theta^2*((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta-1)*(-theta*(1-F1)^(theta-1)+(1-F22)^(theta)*theta*(1-F1)^(theta-1))-1/theta*((-theta*(1-F1)^(theta-1)+(1-F22)^(theta)*theta*(1-F1)^(theta-1))*(-dcop.theta2*((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(-1)-((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta-2)*((1-F1)^theta*log(1-F1)+(1-F22)^theta*log(1-F22)-(1-F1)^theta*log(1-F1)*(1-F22)^theta-(1-F22)^theta*log(1-F22)*(1-F1)^theta))+((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta-1)*(-(1-F1)^(theta-1)-theta*(1-F1)^(theta-1)*log(1-F1)+theta*(1-F22)^(theta)*log(1-F22)*(1-F1)^(theta-1)+((1-F1)^(theta-1)+theta*(1-F1)^(theta-1)*log(1-F1))*(1-F22)^theta))
    
    d2Cdcop1dcop.theta1 <- 1/theta^2*((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta-1)*(-theta*(1-F2)^(theta-1)+(1-F1)^(theta)*theta*(1-F2)^(theta-1))-1/theta*((-theta*(1-F2)^(theta-1)+(1-F1)^(theta)*theta*(1-F2)^(theta-1))*(-dcop.theta1*((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(-1)-((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta-2)*((1-F1)^theta*log(1-F1)+(1-F2)^theta*log(1-F2)-(1-F1)^theta*log(1-F1)*(1-F2)^theta-(1-F2)^theta*log(1-F2)*(1-F1)^theta))+((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta-1)*(-(1-F2)^(theta-1)-theta*(1-F2)^(theta-1)*log(1-F2)+theta*(1-F1)^(theta)*log(1-F1)*(1-F2)^(theta-1)+((1-F2)^(theta-1)+theta*(1-F2)^(theta-1)*log(1-F2))*(1-F1)^theta))
    d2Cdcop11dcop.theta2 <- 1/theta^2*((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta-1)*(-theta*(1-F22)^(theta-1)+(1-F1)^(theta)*theta*(1-F22)^(theta-1))-1/theta*((-theta*(1-F22)^(theta-1)+(1-F1)^(theta)*theta*(1-F22)^(theta-1))*(-dcop.theta2*((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(-1)-((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta-2)*((1-F1)^theta*log(1-F1)+(1-F22)^theta*log(1-F22)-(1-F1)^theta*log(1-F1)*(1-F22)^theta-(1-F22)^theta*log(1-F22)*(1-F1)^theta))+((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta-1)*(-(1-F22)^(theta-1)-theta*(1-F22)^(theta-1)*log(1-F22)+theta*(1-F1)^(theta)*log(1-F1)*(1-F22)^(theta-1)+((1-F22)^(theta-1)+theta*(1-F22)^(theta-1)*log(1-F22))*(1-F1)^theta))
    
    d2Cdcop2dcop.theta1 <- -2/theta^3*log(((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta)))*((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta)+1/theta^2*(1/((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))*((1-F1)^theta*log(1-F1)+(1-F2)^theta*log(1-F2)-(1-F1)^theta*log(1-F1)*(1-F2)^theta-(1-F2)^theta*log(1-F2)*(1-F1)^theta)*((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta)+log(((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta)))*(-dcop.theta1))+1/theta^2*((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta-1)*((1-F1)^theta*log(1-F1)+(1-F2)^theta*log(1-F2)-(1-F1)^theta*log(1-F1)*(1-F2)^theta-(1-F2)^theta*log(1-F2)*(1-F1)^theta)-1/theta*(((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta-1)*((1-F1)^theta*log(1-F1)^2+(1-F2)^theta*log(1-F2)^2-(1-F1)^theta*log(1-F1)^2*(1-F2)^theta-(1-F2)^theta*log(1-F2)^2*(1-F1)^theta-2*(1-F2)^theta*log(1-F2)*(1-F1)^theta*log(1-F1))+((1-F1)^theta*log(1-F1)+(1-F2)^theta*log(1-F2)-(1-F1)^theta*log(1-F1)*(1-F2)^theta-(1-F2)^theta*log(1-F2)*(1-F1)^theta)*(-dcop.theta1*((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(-1)-((1-F1)^(theta)+(1-F2)^(theta)-(1-F1)^(theta)*(1-F2)^(theta))^(1/theta-2)*((1-F1)^theta*log(1-F1)+(1-F2)^theta*log(1-F2)-(1-F1)^theta*log(1-F1)*(1-F2)^theta-(1-F2)^theta*log(1-F2)*(1-F1)^theta)))
    d2Cdcop22dcop.theta2 <- -2/theta^3*log(((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta)))*((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta)+1/theta^2*(1/((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))*((1-F1)^theta*log(1-F1)+(1-F22)^theta*log(1-F22)-(1-F1)^theta*log(1-F1)*(1-F22)^theta-(1-F22)^theta*log(1-F22)*(1-F1)^theta)*((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta)+log(((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta)))*(-dcop.theta2))+1/theta^2*((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta-1)*((1-F1)^theta*log(1-F1)+(1-F22)^theta*log(1-F22)-(1-F1)^theta*log(1-F1)*(1-F22)^theta-(1-F22)^theta*log(1-F22)*(1-F1)^theta)-1/theta*(((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta-1)*((1-F1)^theta*log(1-F1)^2+(1-F22)^theta*log(1-F22)^2-(1-F1)^theta*log(1-F1)^2*(1-F22)^theta-(1-F22)^theta*log(1-F22)^2*(1-F1)^theta-2*(1-F22)^theta*log(1-F22)*(1-F1)^theta*log(1-F1))+((1-F1)^theta*log(1-F1)+(1-F22)^theta*log(1-F22)-(1-F1)^theta*log(1-F1)*(1-F22)^theta-(1-F22)^theta*log(1-F22)*(1-F1)^theta)*(-dcop.theta2*((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(-1)-((1-F1)^(theta)+(1-F22)^(theta)-(1-F1)^(theta)*(1-F22)^(theta))^(1/theta-2)*((1-F1)^theta*log(1-F1)+(1-F22)^theta*log(1-F22)-(1-F1)^theta*log(1-F1)*(1-F22)^theta-(1-F22)^theta*log(1-F22)*(1-F1)^theta)))
    
    
  } else { if (VC$BivD=="J90") {
    theta <- 1 + exp(theta.star) + precision
    Joe1 <- F2 - (1-((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta))
    Joe2 <- F22 - (1-((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta))
    Joe1 <- ifelse(Joe1>precision, Joe1, precision)
    Joe1 <- ifelse(Joe1<(1-precision), Joe1, 1-precision)
    Joe2 <- ifelse(Joe2>precision, Joe2, precision)
    Joe2 <- ifelse(Joe2<(1-precision), Joe2, 1-precision)
    lx <- f2-Joe1+Joe2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1<- (-1/theta)*((1-(1-F1))^theta+(1-F2)^theta-(1-(1-F1))^theta*(1-F2)^theta)^(1/theta-1)*(-theta*(1-(1-F1))^(theta-1)+(1-F2)^theta*theta*(1-(1-F1))^(theta-1))
    dcop11<- (-1/theta)*((1-(1-F1))^theta+(1-F22)^theta-(1-(1-F1))^theta*(1-F22)^theta)^(1/theta-1)*(-theta*(1-(1-F1))^(theta-1)+(1-F22)^theta*theta*(1-(1-F1))^(theta-1))   
    
    dcop2<- 1-(-1/theta)*((1-(1-F1))^theta+(1-F2)^theta-(1-(1-F1))^theta*(1-F2)^theta)^(1/theta-1)*(-theta*(1-F2)^(theta-1)+(1-(1-F1))^theta*theta*(1-F2)^(theta-1))
    dcop22<- 1-(-1/theta)*((1-(1-F1))^theta+(1-F22)^theta-(1-(1-F1))^theta*(1-F22)^theta)^(1/theta-1)*(-theta*(1-F22)^(theta-1)+(1-(1-F1))^theta*theta*(1-F22)^(theta-1)) 
    
    dcop.theta1<- -( (-1/theta^2)*log((1-(1-F1))^theta+(1-F2)^theta-(1-(1-F1))^theta*(1-F2)^theta)*((1-(1-F1))^theta+(1-F2)^theta-(1-(1-F1))^theta*(1-F2)^theta)^(1/theta)+(1/theta)*((1-(1-F1))^theta+(1-F2)^theta-(1-(1-F1))^theta*(1-F2)^theta)^(1/theta-1)*((1-(1-F1))^theta*log(1-(1-F1))+(1-F2)^theta*log(1-F2)-(1-(1-F1))^theta*log(1-(1-F1))*(1-F2)^theta-(1-F2)^theta*log(1-F2)*(1-(1-F1))^theta)    )
    dcop.theta2<- -( (-1/theta^2)*log((1-(1-F1))^theta+(1-F22)^theta-(1-(1-F1))^theta*(1-F22)^theta)*((1-(1-F1))^theta+(1-F22)^theta-(1-(1-F1))^theta*(1-F22)^theta)^(1/theta)+(1/theta)*((1-(1-F1))^theta+(1-F22)^theta-(1-(1-F1))^theta*(1-F22)^theta)^(1/theta-1)*((1-(1-F1))^theta*log(1-(1-F1))+(1-F22)^theta*log(1-F22)-(1-(1-F1))^theta*log(1-(1-F1))*(1-F22)^theta-(1-F22)^theta*log(1-F22)*(1-(1-F1))^theta)    )
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- -(  -1/theta*(1/theta-1)*((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta-2)*(-theta*(1-(1-F1))^(theta-1)+(1-F2)^(theta)*theta*(1-(1-F1))^(theta-1))^2-1/theta*((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta-1)*(theta*(theta-1)*(1-(1-F1))^(theta-2)-(1-F2)^(theta)*theta*(theta-1)*(1-(1-F1))^(theta-2))     )
    d2Cdcop112 <- -( -1/theta*(1/theta-1)*((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta-2)*(-theta*(1-(1-F1))^(theta-1)+(1-F22)^(theta)*theta*(1-(1-F1))^(theta-1))^2-1/theta*((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta-1)*(theta*(theta-1)*(1-(1-F1))^(theta-2)-(1-F22)^(theta)*theta*(theta-1)*(1-(1-F1))^(theta-2))      )
    
    d2Cdcop22 <-  -( -1/theta*(1/theta-1)*((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta-2)*(-theta*(1-F2)^(theta-1)+(1-(1-F1))^(theta)*theta*(1-F2)^(theta-1))^2-1/theta*((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta-1)*(theta*(theta-1)*(1-F2)^(theta-2)-(1-(1-F1))^(theta)*theta*(theta-1)*(1-F2)^(theta-2))     )
    d2Cdcop222 <- -( -1/theta*(1/theta-1)*((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta-2)*(-theta*(1-F22)^(theta-1)+(1-(1-F1))^(theta)*theta*(1-F22)^(theta-1))^2-1/theta*((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta-1)*(theta*(theta-1)*(1-F22)^(theta-2)-(1-(1-F1))^(theta)*theta*(theta-1)*(1-F22)^(theta-2))      )
    
    d2Cdcop1dcop2 <- -1/theta*(1/theta-1)*((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta-2)*(-theta*(1-F2)^(theta-1)+(1-(1-F1))^(theta)*theta*(1-F2)^(theta-1))*(-theta*(1-(1-F1))^(theta-1)+(1-F2)^(theta)*theta*(1-(1-F1))^(theta-1))-1/theta*((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta-1)*(-theta^2*(1-F2)^(theta-1)*(1-(1-F1))^(theta-1))
    d2Cdcop11dcop22 <- -1/theta*(1/theta-1)*((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta-2)*(-theta*(1-F22)^(theta-1)+(1-(1-F1))^(theta)*theta*(1-F22)^(theta-1))*(-theta*(1-(1-F1))^(theta-1)+(1-F22)^(theta)*theta*(1-(1-F1))^(theta-1))-1/theta*((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta-1)*(-theta^2*(1-F22)^(theta-1)*(1-(1-F1))^(theta-1))
    
    d2Cdcop.theta12 <-  -( 1/theta^2*((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta-1)*(-theta*(1-(1-F1))^(theta-1)+(1-F2)^(theta)*theta*(1-(1-F1))^(theta-1))-1/theta*((-theta*(1-(1-F1))^(theta-1)+(1-F2)^(theta)*theta*(1-(1-F1))^(theta-1))*(-dcop.theta1*((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(-1)-((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta-2)*((1-(1-F1))^theta*log(1-(1-F1))+(1-F2)^theta*log(1-F2)-(1-(1-F1))^theta*log(1-(1-F1))*(1-F2)^theta-(1-F2)^theta*log(1-F2)*(1-(1-F1))^theta))+((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta-1)*(-(1-(1-F1))^(theta-1)-theta*(1-(1-F1))^(theta-1)*log(1-(1-F1))+theta*(1-F2)^(theta)*log(1-F2)*(1-(1-F1))^(theta-1)+((1-(1-F1))^(theta-1)+theta*(1-(1-F1))^(theta-1)*log(1-(1-F1)))*(1-F2)^theta))   )
    d2Cdcop.theta22 <-  -(  1/theta^2*((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta-1)*(-theta*(1-(1-F1))^(theta-1)+(1-F22)^(theta)*theta*(1-(1-F1))^(theta-1))-1/theta*((-theta*(1-(1-F1))^(theta-1)+(1-F22)^(theta)*theta*(1-(1-F1))^(theta-1))*(-dcop.theta2*((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(-1)-((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta-2)*((1-(1-F1))^theta*log(1-(1-F1))+(1-F22)^theta*log(1-F22)-(1-(1-F1))^theta*log(1-(1-F1))*(1-F22)^theta-(1-F22)^theta*log(1-F22)*(1-(1-F1))^theta))+((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta-1)*(-(1-(1-F1))^(theta-1)-theta*(1-(1-F1))^(theta-1)*log(1-(1-F1))+theta*(1-F22)^(theta)*log(1-F22)*(1-(1-F1))^(theta-1)+((1-(1-F1))^(theta-1)+theta*(1-(1-F1))^(theta-1)*log(1-(1-F1)))*(1-F22)^theta))   )
    
    d2Cdcop1dcop.theta1 <-  1/theta^2*((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta-1)*(-theta*(1-F2)^(theta-1)+(1-(1-F1))^(theta)*theta*(1-F2)^(theta-1))-1/theta*((-theta*(1-F2)^(theta-1)+(1-(1-F1))^(theta)*theta*(1-F2)^(theta-1))*(-dcop.theta1*((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(-1)-((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta-2)*((1-(1-F1))^theta*log(1-(1-F1))+(1-F2)^theta*log(1-F2)-(1-(1-F1))^theta*log(1-(1-F1))*(1-F2)^theta-(1-F2)^theta*log(1-F2)*(1-(1-F1))^theta))+((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta-1)*(-(1-F2)^(theta-1)-theta*(1-F2)^(theta-1)*log(1-F2)+theta*(1-(1-F1))^(theta)*log(1-(1-F1))*(1-F2)^(theta-1)+((1-F2)^(theta-1)+theta*(1-F2)^(theta-1)*log(1-F2))*(1-(1-F1))^theta))     
    d2Cdcop11dcop.theta2 <-  1/theta^2*((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta-1)*(-theta*(1-F22)^(theta-1)+(1-(1-F1))^(theta)*theta*(1-F22)^(theta-1))-1/theta*((-theta*(1-F22)^(theta-1)+(1-(1-F1))^(theta)*theta*(1-F22)^(theta-1))*(-dcop.theta2*((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(-1)-((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta-2)*((1-(1-F1))^theta*log(1-(1-F1))+(1-F22)^theta*log(1-F22)-(1-(1-F1))^theta*log(1-(1-F1))*(1-F22)^theta-(1-F22)^theta*log(1-F22)*(1-(1-F1))^theta))+((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta-1)*(-(1-F22)^(theta-1)-theta*(1-F22)^(theta-1)*log(1-F22)+theta*(1-(1-F1))^(theta)*log(1-(1-F1))*(1-F22)^(theta-1)+((1-F22)^(theta-1)+theta*(1-F22)^(theta-1)*log(1-F22))*(1-(1-F1))^theta))     
    
    d2Cdcop2dcop.theta1 <- -( -2/theta^3*log(((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta)))*((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta)+1/theta^2*(1/((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))*((1-(1-F1))^theta*log(1-(1-F1))+(1-F2)^theta*log(1-F2)-(1-(1-F1))^theta*log(1-(1-F1))*(1-F2)^theta-(1-F2)^theta*log(1-F2)*(1-(1-F1))^theta)*((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta)+log(((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta)))*(-dcop.theta1))+1/theta^2*((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta-1)*((1-(1-F1))^theta*log(1-(1-F1))+(1-F2)^theta*log(1-F2)-(1-(1-F1))^theta*log(1-(1-F1))*(1-F2)^theta-(1-F2)^theta*log(1-F2)*(1-(1-F1))^theta)-1/theta*(((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta-1)*((1-(1-F1))^theta*log(1-(1-F1))^2+(1-F2)^theta*log(1-F2)^2-(1-(1-F1))^theta*log(1-(1-F1))^2*(1-F2)^theta-(1-F2)^theta*log(1-F2)^2*(1-(1-F1))^theta-2*(1-F2)^theta*log(1-F2)*(1-(1-F1))^theta*log(1-(1-F1)))+((1-(1-F1))^theta*log(1-(1-F1))+(1-F2)^theta*log(1-F2)-(1-(1-F1))^theta*log(1-(1-F1))*(1-F2)^theta-(1-F2)^theta*log(1-F2)*(1-(1-F1))^theta)*(-dcop.theta1*((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(-1)-((1-(1-F1))^(theta)+(1-F2)^(theta)-(1-(1-F1))^(theta)*(1-F2)^(theta))^(1/theta-2)*((1-(1-F1))^theta*log(1-(1-F1))+(1-F2)^theta*log(1-F2)-(1-(1-F1))^theta*log(1-(1-F1))*(1-F2)^theta-(1-F2)^theta*log(1-F2)*(1-(1-F1))^theta)))    )
    d2Cdcop22dcop.theta2 <- -( -2/theta^3*log(((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta)))*((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta)+1/theta^2*(1/((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))*((1-(1-F1))^theta*log(1-(1-F1))+(1-F22)^theta*log(1-F22)-(1-(1-F1))^theta*log(1-(1-F1))*(1-F22)^theta-(1-F22)^theta*log(1-F22)*(1-(1-F1))^theta)*((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta)+log(((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta)))*(-dcop.theta2))+1/theta^2*((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta-1)*((1-(1-F1))^theta*log(1-(1-F1))+(1-F22)^theta*log(1-F22)-(1-(1-F1))^theta*log(1-(1-F1))*(1-F22)^theta-(1-F22)^theta*log(1-F22)*(1-(1-F1))^theta)-1/theta*(((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta-1)*((1-(1-F1))^theta*log(1-(1-F1))^2+(1-F22)^theta*log(1-F22)^2-(1-(1-F1))^theta*log(1-(1-F1))^2*(1-F22)^theta-(1-F22)^theta*log(1-F22)^2*(1-(1-F1))^theta-2*(1-F22)^theta*log(1-F22)*(1-(1-F1))^theta*log(1-(1-F1)))+((1-(1-F1))^theta*log(1-(1-F1))+(1-F22)^theta*log(1-F22)-(1-(1-F1))^theta*log(1-(1-F1))*(1-F22)^theta-(1-F22)^theta*log(1-F22)*(1-(1-F1))^theta)*(-dcop.theta2*((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(-1)-((1-(1-F1))^(theta)+(1-F22)^(theta)-(1-(1-F1))^(theta)*(1-F22)^(theta))^(1/theta-2)*((1-(1-F1))^theta*log(1-(1-F1))+(1-F22)^theta*log(1-F22)-(1-(1-F1))^theta*log(1-(1-F1))*(1-F22)^theta-(1-F22)^theta*log(1-F22)*(1-(1-F1))^theta)))    )
    
    
    
  } else { if (VC$BivD=="J180") {
    theta <- 1 + exp(theta.star) + precision
    Joe1 <- F1 + F2 - 1 + (1-((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta))
    Joe2 <- F1 + F22 - 1 + (1-((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta))
    Joe1 <- ifelse(Joe1>precision, Joe1, precision)
    Joe1 <- ifelse(Joe1<(1-precision), Joe1, 1-precision)
    Joe2 <- ifelse(Joe2>precision, Joe2, precision)
    Joe2 <- ifelse(Joe2<(1-precision), Joe2, 1-precision)
    lx <- f2 - Joe1+ Joe2
    lx <- ifelse(lx>precision, lx, precision)
    
    dcop1 <- 1-(-1/theta)*((1-(1-F1))^theta+(1-(1-F2))^theta-(1-(1-F1))^theta*(1-(1-F2))^theta)^(1/theta-1)*(-theta*(1-(1-F1))^(theta-1)+(1-(1-F2))^theta*theta*(1-(1-F1))^(theta-1))
    dcop11 <- 1-(-1/theta)*((1-(1-F1))^theta+(1-(1-F22))^theta-(1-(1-F1))^theta*(1-(1-F22))^theta)^(1/theta-1)*(-theta*(1-(1-F1))^(theta-1)+(1-(1-F22))^theta*theta*(1-(1-F1))^(theta-1))   
    
    dcop2 <- 1-(-1/theta)*((1-(1-F1))^theta+(1-(1-F2))^theta-(1-(1-F1))^theta*(1-(1-F2))^theta)^(1/theta-1)*(-theta*(1-(1-F2))^(theta-1)+(1-(1-F1))^theta*theta*(1-(1-F2))^(theta-1))
    dcop22 <- 1-(-1/theta)*((1-(1-F1))^theta+(1-(1-F22))^theta-(1-(1-F1))^theta*(1-(1-F22))^theta)^(1/theta-1)*(-theta*(1-(1-F22))^(theta-1)+(1-(1-F1))^theta*theta*(1-(1-F22))^(theta-1)) 
    
    dcop.theta1 <- -( (-1/theta^2)*log((1-(1-F1))^theta+(1-(1-F2))^theta-(1-(1-F1))^theta*(1-(1-F2))^theta)*((1-(1-F1))^theta+(1-(1-F2))^theta-(1-(1-F1))^theta*(1-(1-F2))^theta)^(1/theta)+(1/theta)*((1-(1-F1))^theta+(1-(1-F2))^theta-(1-(1-F1))^theta*(1-(1-F2))^theta)^(1/theta-1)*((1-(1-F1))^theta*log(1-(1-F1))+(1-(1-F2))^theta*log(1-(1-F2))-(1-(1-F1))^theta*log(1-(1-F1))*(1-(1-F2))^theta-(1-(1-F2))^theta*log(1-(1-F2))*(1-(1-F1))^theta)    )
    dcop.theta2 <- -( (-1/theta^2)*log((1-(1-F1))^theta+(1-(1-F22))^theta-(1-(1-F1))^theta*(1-(1-F22))^theta)*((1-(1-F1))^theta+(1-(1-F22))^theta-(1-(1-F1))^theta*(1-(1-F22))^theta)^(1/theta)+(1/theta)*((1-(1-F1))^theta+(1-(1-F22))^theta-(1-(1-F1))^theta*(1-(1-F22))^theta)^(1/theta-1)*((1-(1-F1))^theta*log(1-(1-F1))+(1-(1-F22))^theta*log(1-(1-F22))-(1-(1-F1))^theta*log(1-(1-F1))*(1-(1-F22))^theta-(1-(1-F22))^theta*log(1-(1-F22))*(1-(1-F1))^theta)    )
    
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- (  -1/theta*(1/theta-1)*((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta-2)*(-theta*(1-(1-F1))^(theta-1)+(1-(1-F2))^(theta)*theta*(1-(1-F1))^(theta-1))^2-1/theta*((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*(theta*(theta-1)*(1-(1-F1))^(theta-2)-(1-(1-F2))^(theta)*theta*(theta-1)*(1-(1-F1))^(theta-2))     )
    d2Cdcop112 <- ( -1/theta*(1/theta-1)*((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta-2)*(-theta*(1-(1-F1))^(theta-1)+(1-(1-F22))^(theta)*theta*(1-(1-F1))^(theta-1))^2-1/theta*((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*(theta*(theta-1)*(1-(1-F1))^(theta-2)-(1-(1-F22))^(theta)*theta*(theta-1)*(1-(1-F1))^(theta-2))      )
    
    d2Cdcop22 <-  ( -1/theta*(1/theta-1)*((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta-2)*(-theta*(1-(1-F2))^(theta-1)+(1-(1-F1))^(theta)*theta*(1-(1-F2))^(theta-1))^2-1/theta*((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*(theta*(theta-1)*(1-(1-F2))^(theta-2)-(1-(1-F1))^(theta)*theta*(theta-1)*(1-(1-F2))^(theta-2))     )
    d2Cdcop222 <- ( -1/theta*(1/theta-1)*((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta-2)*(-theta*(1-(1-F22))^(theta-1)+(1-(1-F1))^(theta)*theta*(1-(1-F22))^(theta-1))^2-1/theta*((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*(theta*(theta-1)*(1-(1-F22))^(theta-2)-(1-(1-F1))^(theta)*theta*(theta-1)*(1-(1-F22))^(theta-2))      )
    
    d2Cdcop1dcop2 <- ( -1/theta*(1/theta-1)*((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta-2)*(-theta*(1-(1-F2))^(theta-1)+(1-(1-F1))^(theta)*theta*(1-(1-F2))^(theta-1))*(-theta*(1-(1-F1))^(theta-1)+(1-(1-F2))^(theta)*theta*(1-(1-F1))^(theta-1))-1/theta*((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*(-theta^2*(1-(1-F2))^(theta-1)*(1-(1-F1))^(theta-1))   )
    d2Cdcop11dcop22 <- ( -1/theta*(1/theta-1)*((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta-2)*(-theta*(1-(1-F22))^(theta-1)+(1-(1-F1))^(theta)*theta*(1-(1-F22))^(theta-1))*(-theta*(1-(1-F1))^(theta-1)+(1-(1-F22))^(theta)*theta*(1-(1-F1))^(theta-1))-1/theta*((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*(-theta^2*(1-(1-F22))^(theta-1)*(1-(1-F1))^(theta-1))   )
    
    d2Cdcop.theta12 <- -( 1/theta^2*((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*(-theta*(1-(1-F1))^(theta-1)+(1-(1-F2))^(theta)*theta*(1-(1-F1))^(theta-1))-1/theta*((-theta*(1-(1-F1))^(theta-1)+(1-(1-F2))^(theta)*theta*(1-(1-F1))^(theta-1))*(-dcop.theta1*((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(-1)-((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta-2)*((1-(1-F1))^theta*log(1-(1-F1))+(1-(1-F2))^theta*log(1-(1-F2))-(1-(1-F1))^theta*log(1-(1-F1))*(1-(1-F2))^theta-(1-(1-F2))^theta*log(1-(1-F2))*(1-(1-F1))^theta))+((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*(-(1-(1-F1))^(theta-1)-theta*(1-(1-F1))^(theta-1)*log(1-(1-F1))+theta*(1-(1-F2))^(theta)*log(1-(1-F2))*(1-(1-F1))^(theta-1)+((1-(1-F1))^(theta-1)+theta*(1-(1-F1))^(theta-1)*log(1-(1-F1)))*(1-(1-F2))^theta))   )
    d2Cdcop.theta22 <- -( 1/theta^2*((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*(-theta*(1-(1-F1))^(theta-1)+(1-(1-F22))^(theta)*theta*(1-(1-F1))^(theta-1))-1/theta*((-theta*(1-(1-F1))^(theta-1)+(1-(1-F22))^(theta)*theta*(1-(1-F1))^(theta-1))*(-dcop.theta2*((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(-1)-((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta-2)*((1-(1-F1))^theta*log(1-(1-F1))+(1-(1-F22))^theta*log(1-(1-F22))-(1-(1-F1))^theta*log(1-(1-F1))*(1-(1-F22))^theta-(1-(1-F22))^theta*log(1-(1-F22))*(1-(1-F1))^theta))+((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*(-(1-(1-F1))^(theta-1)-theta*(1-(1-F1))^(theta-1)*log(1-(1-F1))+theta*(1-(1-F22))^(theta)*log(1-(1-F22))*(1-(1-F1))^(theta-1)+((1-(1-F1))^(theta-1)+theta*(1-(1-F1))^(theta-1)*log(1-(1-F1)))*(1-(1-F22))^theta))    )
    
    d2Cdcop1dcop.theta1 <-  -( 1/theta^2*((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*(-theta*(1-(1-F2))^(theta-1)+(1-(1-F1))^(theta)*theta*(1-(1-F2))^(theta-1))-1/theta*((-theta*(1-(1-F2))^(theta-1)+(1-(1-F1))^(theta)*theta*(1-(1-F2))^(theta-1))*(-dcop.theta1*((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(-1)-((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta-2)*((1-(1-F1))^theta*log(1-(1-F1))+(1-(1-F2))^theta*log(1-(1-F2))-(1-(1-F1))^theta*log(1-(1-F1))*(1-(1-F2))^theta-(1-(1-F2))^theta*log(1-(1-F2))*(1-(1-F1))^theta))+((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*(-(1-(1-F2))^(theta-1)-theta*(1-(1-F2))^(theta-1)*log(1-(1-F2))+theta*(1-(1-F1))^(theta)*log(1-(1-F1))*(1-(1-F2))^(theta-1)+((1-(1-F2))^(theta-1)+theta*(1-(1-F2))^(theta-1)*log(1-(1-F2)))*(1-(1-F1))^theta))    )
    d2Cdcop11dcop.theta2 <-   -( 1/theta^2*((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*(-theta*(1-(1-F22))^(theta-1)+(1-(1-F1))^(theta)*theta*(1-(1-F22))^(theta-1))-1/theta*((-theta*(1-(1-F22))^(theta-1)+(1-(1-F1))^(theta)*theta*(1-(1-F22))^(theta-1))*(-dcop.theta2*((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(-1)-((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta-2)*((1-(1-F1))^theta*log(1-(1-F1))+(1-(1-F22))^theta*log(1-(1-F22))-(1-(1-F1))^theta*log(1-(1-F1))*(1-(1-F22))^theta-(1-(1-F22))^theta*log(1-(1-F22))*(1-(1-F1))^theta))+((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*(-(1-(1-F22))^(theta-1)-theta*(1-(1-F22))^(theta-1)*log(1-(1-F22))+theta*(1-(1-F1))^(theta)*log(1-(1-F1))*(1-(1-F22))^(theta-1)+((1-(1-F22))^(theta-1)+theta*(1-(1-F22))^(theta-1)*log(1-(1-F22)))*(1-(1-F1))^theta))      )   
    
    d2Cdcop2dcop.theta1 <- ( -2/theta^3*log(((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta)))*((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta)+1/theta^2*(1/((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))*((1-(1-F1))^theta*log(1-(1-F1))+(1-(1-F2))^theta*log(1-(1-F2))-(1-(1-F1))^theta*log(1-(1-F1))*(1-(1-F2))^theta-(1-(1-F2))^theta*log(1-(1-F2))*(1-(1-F1))^theta)*((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta)+log(((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta)))*(-dcop.theta1))+1/theta^2*((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*((1-(1-F1))^theta*log(1-(1-F1))+(1-(1-F2))^theta*log(1-(1-F2))-(1-(1-F1))^theta*log(1-(1-F1))*(1-(1-F2))^theta-(1-(1-F2))^theta*log(1-(1-F2))*(1-(1-F1))^theta)-1/theta*(((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*((1-(1-F1))^theta*log(1-(1-F1))^2+(1-(1-F2))^theta*log(1-(1-F2))^2-(1-(1-F1))^theta*log(1-(1-F1))^2*(1-(1-F2))^theta-(1-(1-F2))^theta*log(1-(1-F2))^2*(1-(1-F1))^theta-2*(1-(1-F2))^theta*log(1-(1-F2))*(1-(1-F1))^theta*log(1-(1-F1)))+((1-(1-F1))^theta*log(1-(1-F1))+(1-(1-F2))^theta*log(1-(1-F2))-(1-(1-F1))^theta*log(1-(1-F1))*(1-(1-F2))^theta-(1-(1-F2))^theta*log(1-(1-F2))*(1-(1-F1))^theta)*(-dcop.theta1*((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(-1)-((1-(1-F1))^(theta)+(1-(1-F2))^(theta)-(1-(1-F1))^(theta)*(1-(1-F2))^(theta))^(1/theta-2)*((1-(1-F1))^theta*log(1-(1-F1))+(1-(1-F2))^theta*log(1-(1-F2))-(1-(1-F1))^theta*log(1-(1-F1))*(1-(1-F2))^theta-(1-(1-F2))^theta*log(1-(1-F2))*(1-(1-F1))^theta)))    )
    d2Cdcop22dcop.theta2 <- ( -2/theta^3*log(((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta)))*((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta)+1/theta^2*(1/((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))*((1-(1-F1))^theta*log(1-(1-F1))+(1-(1-F22))^theta*log(1-(1-F22))-(1-(1-F1))^theta*log(1-(1-F1))*(1-(1-F22))^theta-(1-(1-F22))^theta*log(1-(1-F22))*(1-(1-F1))^theta)*((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta)+log(((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta)))*(-dcop.theta2))+1/theta^2*((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*((1-(1-F1))^theta*log(1-(1-F1))+(1-(1-F22))^theta*log(1-(1-F22))-(1-(1-F1))^theta*log(1-(1-F1))*(1-(1-F22))^theta-(1-(1-F22))^theta*log(1-(1-F22))*(1-(1-F1))^theta)-1/theta*(((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*((1-(1-F1))^theta*log(1-(1-F1))^2+(1-(1-F22))^theta*log(1-(1-F22))^2-(1-(1-F1))^theta*log(1-(1-F1))^2*(1-(1-F22))^theta-(1-(1-F22))^theta*log(1-(1-F22))^2*(1-(1-F1))^theta-2*(1-(1-F22))^theta*log(1-(1-F22))*(1-(1-F1))^theta*log(1-(1-F1)))+((1-(1-F1))^theta*log(1-(1-F1))+(1-(1-F22))^theta*log(1-(1-F22))-(1-(1-F1))^theta*log(1-(1-F1))*(1-(1-F22))^theta-(1-(1-F22))^theta*log(1-(1-F22))*(1-(1-F1))^theta)*(-dcop.theta2*((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(-1)-((1-(1-F1))^(theta)+(1-(1-F22))^(theta)-(1-(1-F1))^(theta)*(1-(1-F22))^(theta))^(1/theta-2)*((1-(1-F1))^theta*log(1-(1-F1))+(1-(1-F22))^theta*log(1-(1-F22))-(1-(1-F1))^theta*log(1-(1-F1))*(1-(1-F22))^theta-(1-(1-F22))^theta*log(1-(1-F22))*(1-(1-F1))^theta)))    )
    
    
    
    
  } else { if (VC$BivD=="J270") {
    theta <- 1 + exp(theta.star) + precision
    Joe1 <- F1 - (1-((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta))
    Joe2 <- F1 - (1-((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta))
    Joe1 <- ifelse(Joe1>precision, Joe1, precision)
    Joe1 <- ifelse(Joe1<(1-precision), Joe1, 1-precision)
    Joe2 <- ifelse(Joe2>precision, Joe2, precision)
    Joe2 <- ifelse(Joe2<(1-precision), Joe2, 1-precision)
    lx <- f2 - Joe1 + Joe2
    lx <- ifelse(lx>precision, lx, precision)  
    
    dcop1 <- 1-(-1/theta)*((1-F1)^theta+(1-(1-F2))^theta-(1-F1)^theta*(1-(1-F2))^theta)^(1/theta-1)*(-theta*(1-F1)^(theta-1)+(1-(1-F2))^theta*theta*(1-F1)^(theta-1))
    dcop11 <- 1-(-1/theta)*((1-F1)^theta+(1-(1-F22))^theta-(1-F1)^theta*(1-(1-F22))^theta)^(1/theta-1)*(-theta*(1-F1)^(theta-1)+(1-(1-F22))^theta*theta*(1-F1)^(theta-1))   
    
    dcop2 <- (-1/theta)*((1-F1)^theta+(1-(1-F2))^theta-(1-F1)^theta*(1-(1-F2))^theta)^(1/theta-1)*(-theta*(1-(1-F2))^(theta-1)+(1-F1)^theta*theta*(1-(1-F2))^(theta-1))
    dcop22 <- (-1/theta)*((1-F1)^theta+(1-(1-F22))^theta-(1-F1)^theta*(1-(1-F22))^theta)^(1/theta-1)*(-theta*(1-(1-F22))^(theta-1)+(1-F1)^theta*theta*(1-(1-F22))^(theta-1)) 
    
    dcop.theta1 <- -(  (-1/theta^2)*log((1-F1)^theta+(1-(1-F2))^theta-(1-F1)^theta*(1-(1-F2))^theta)*((1-F1)^theta+(1-(1-F2))^theta-(1-F1)^theta*(1-(1-F2))^theta)^(1/theta)+(1/theta)*((1-F1)^theta+(1-(1-F2))^theta-(1-F1)^theta*(1-(1-F2))^theta)^(1/theta-1)*((1-F1)^theta*log(1-F1)+(1-(1-F2))^theta*log(1-(1-F2))-(1-F1)^theta*log(1-F1)*(1-(1-F2))^theta-(1-(1-F2))^theta*log(1-(1-F2))*(1-F1)^theta)   )
    dcop.theta2 <- -(  (-1/theta^2)*log((1-F1)^theta+(1-(1-F22))^theta-(1-F1)^theta*(1-(1-F22))^theta)*((1-F1)^theta+(1-(1-F22))^theta-(1-F1)^theta*(1-(1-F22))^theta)^(1/theta)+(1/theta)*((1-F1)^theta+(1-(1-F22))^theta-(1-F1)^theta*(1-(1-F22))^theta)^(1/theta-1)*((1-F1)^theta*log(1-F1)+(1-(1-F22))^theta*log(1-(1-F22))-(1-F1)^theta*log(1-F1)*(1-(1-F22))^theta-(1-(1-F22))^theta*log(1-(1-F22))*(1-F1)^theta)     )
    
    
    
    # Hessian derivative components
    
    d2Cdcop12 <- -(  -1/theta*(1/theta-1)*((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta-2)*(-theta*(1-F1)^(theta-1)+(1-(1-F2))^(theta)*theta*(1-F1)^(theta-1))^2-1/theta*((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*(theta*(theta-1)*(1-F1)^(theta-2)-(1-(1-F2))^(theta)*theta*(theta-1)*(1-F1)^(theta-2))   )
    d2Cdcop112 <- -( -1/theta*(1/theta-1)*((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta-2)*(-theta*(1-F1)^(theta-1)+(1-(1-F22))^(theta)*theta*(1-F1)^(theta-1))^2-1/theta*((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*(theta*(theta-1)*(1-F1)^(theta-2)-(1-(1-F22))^(theta)*theta*(theta-1)*(1-F1)^(theta-2))    )
    
    d2Cdcop22 <- -(  -1/theta*(1/theta-1)*((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta-2)*(-theta*(1-(1-F2))^(theta-1)+(1-F1)^(theta)*theta*(1-(1-F2))^(theta-1))^2-1/theta*((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*(theta*(theta-1)*(1-(1-F2))^(theta-2)-(1-F1)^(theta)*theta*(theta-1)*(1-(1-F2))^(theta-2))    )
    d2Cdcop222 <- -( -1/theta*(1/theta-1)*((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta-2)*(-theta*(1-(1-F22))^(theta-1)+(1-F1)^(theta)*theta*(1-(1-F22))^(theta-1))^2-1/theta*((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*(theta*(theta-1)*(1-(1-F22))^(theta-2)-(1-F1)^(theta)*theta*(theta-1)*(1-(1-F22))^(theta-2))    )
    
    d2Cdcop1dcop2 <- -1/theta*(1/theta-1)*((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta-2)*(-theta*(1-(1-F2))^(theta-1)+(1-F1)^(theta)*theta*(1-(1-F2))^(theta-1))*(-theta*(1-F1)^(theta-1)+(1-(1-F2))^(theta)*theta*(1-F1)^(theta-1))-1/theta*((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*(-theta^2*(1-(1-F2))^(theta-1)*(1-F1)^(theta-1))
    d2Cdcop11dcop22 <- -1/theta*(1/theta-1)*((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta-2)*(-theta*(1-(1-F22))^(theta-1)+(1-F1)^(theta)*theta*(1-(1-F22))^(theta-1))*(-theta*(1-F1)^(theta-1)+(1-(1-F22))^(theta)*theta*(1-F1)^(theta-1))-1/theta*((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*(-theta^2*(1-(1-F22))^(theta-1)*(1-F1)^(theta-1))
    
    d2Cdcop.theta12 <- (  1/theta^2*((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*(-theta*(1-F1)^(theta-1)+(1-(1-F2))^(theta)*theta*(1-F1)^(theta-1))-1/theta*((-theta*(1-F1)^(theta-1)+(1-(1-F2))^(theta)*theta*(1-F1)^(theta-1))*(-dcop.theta1*((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(-1)-((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta-2)*((1-F1)^theta*log(1-F1)+(1-(1-F2))^theta*log(1-(1-F2))-(1-F1)^theta*log(1-F1)*(1-(1-F2))^theta-(1-(1-F2))^theta*log(1-(1-F2))*(1-F1)^theta))+((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*(-(1-F1)^(theta-1)-theta*(1-F1)^(theta-1)*log(1-F1)+theta*(1-(1-F2))^(theta)*log(1-(1-F2))*(1-F1)^(theta-1)+((1-F1)^(theta-1)+theta*(1-F1)^(theta-1)*log(1-F1))*(1-(1-F2))^theta))    )
    d2Cdcop.theta22 <-  (  1/theta^2*((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*(-theta*(1-F1)^(theta-1)+(1-(1-F22))^(theta)*theta*(1-F1)^(theta-1))-1/theta*((-theta*(1-F1)^(theta-1)+(1-(1-F22))^(theta)*theta*(1-F1)^(theta-1))*(-dcop.theta2*((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(-1)-((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta-2)*((1-F1)^theta*log(1-F1)+(1-(1-F22))^theta*log(1-(1-F22))-(1-F1)^theta*log(1-F1)*(1-(1-F22))^theta-(1-(1-F22))^theta*log(1-(1-F22))*(1-F1)^theta))+((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*(-(1-F1)^(theta-1)-theta*(1-F1)^(theta-1)*log(1-F1)+theta*(1-(1-F22))^(theta)*log(1-(1-F22))*(1-F1)^(theta-1)+((1-F1)^(theta-1)+theta*(1-F1)^(theta-1)*log(1-F1))*(1-(1-F22))^theta))       )
    
    d2Cdcop1dcop.theta1 <-  -( 1/theta^2*((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*(-theta*(1-(1-F2))^(theta-1)+(1-F1)^(theta)*theta*(1-(1-F2))^(theta-1))-1/theta*((-theta*(1-(1-F2))^(theta-1)+(1-F1)^(theta)*theta*(1-(1-F2))^(theta-1))*(-dcop.theta1*((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(-1)-((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta-2)*((1-F1)^theta*log(1-F1)+(1-(1-F2))^theta*log(1-(1-F2))-(1-F1)^theta*log(1-F1)*(1-(1-F2))^theta-(1-(1-F2))^theta*log(1-(1-F2))*(1-F1)^theta))+((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*(-(1-(1-F2))^(theta-1)-theta*(1-(1-F2))^(theta-1)*log(1-(1-F2))+theta*(1-F1)^(theta)*log(1-F1)*(1-(1-F2))^(theta-1)+((1-(1-F2))^(theta-1)+theta*(1-(1-F2))^(theta-1)*log(1-(1-F2)))*(1-F1)^theta))   )
    d2Cdcop11dcop.theta2 <-  -(   1/theta^2*((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*(-theta*(1-(1-F22))^(theta-1)+(1-F1)^(theta)*theta*(1-(1-F22))^(theta-1))-1/theta*((-theta*(1-(1-F22))^(theta-1)+(1-F1)^(theta)*theta*(1-(1-F22))^(theta-1))*(-dcop.theta2*((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(-1)-((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta-2)*((1-F1)^theta*log(1-F1)+(1-(1-F22))^theta*log(1-(1-F22))-(1-F1)^theta*log(1-F1)*(1-(1-F22))^theta-(1-(1-F22))^theta*log(1-(1-F22))*(1-F1)^theta))+((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*(-(1-(1-F22))^(theta-1)-theta*(1-(1-F22))^(theta-1)*log(1-(1-F22))+theta*(1-F1)^(theta)*log(1-F1)*(1-(1-F22))^(theta-1)+((1-(1-F22))^(theta-1)+theta*(1-(1-F22))^(theta-1)*log(1-(1-F22)))*(1-F1)^theta))  )
    
    d2Cdcop2dcop.theta1 <- -(   -2/theta^3*log(((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta)))*((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta)+1/theta^2*(1/((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))*((1-F1)^theta*log(1-F1)+(1-(1-F2))^theta*log(1-(1-F2))-(1-F1)^theta*log(1-F1)*(1-(1-F2))^theta-(1-(1-F2))^theta*log(1-(1-F2))*(1-F1)^theta)*((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta)+log(((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta)))*(-dcop.theta1))+1/theta^2*((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*((1-F1)^theta*log(1-F1)+(1-(1-F2))^theta*log(1-(1-F2))-(1-F1)^theta*log(1-F1)*(1-(1-F2))^theta-(1-(1-F2))^theta*log(1-(1-F2))*(1-F1)^theta)-1/theta*(((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta-1)*((1-F1)^theta*log(1-F1)^2+(1-(1-F2))^theta*log(1-(1-F2))^2-(1-F1)^theta*log(1-F1)^2*(1-(1-F2))^theta-(1-(1-F2))^theta*log(1-(1-F2))^2*(1-F1)^theta-2*(1-(1-F2))^theta*log(1-(1-F2))*(1-F1)^theta*log(1-F1))+((1-F1)^theta*log(1-F1)+(1-(1-F2))^theta*log(1-(1-F2))-(1-F1)^theta*log(1-F1)*(1-(1-F2))^theta-(1-(1-F2))^theta*log(1-(1-F2))*(1-F1)^theta)*(-dcop.theta1*((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(-1)-((1-F1)^(theta)+(1-(1-F2))^(theta)-(1-F1)^(theta)*(1-(1-F2))^(theta))^(1/theta-2)*((1-F1)^theta*log(1-F1)+(1-(1-F2))^theta*log(1-(1-F2))-(1-F1)^theta*log(1-F1)*(1-(1-F2))^theta-(1-(1-F2))^theta*log(1-(1-F2))*(1-F1)^theta)))   )
    d2Cdcop22dcop.theta2 <-  -(   -2/theta^3*log(((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta)))*((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta)+1/theta^2*(1/((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))*((1-F1)^theta*log(1-F1)+(1-(1-F22))^theta*log(1-(1-F22))-(1-F1)^theta*log(1-F1)*(1-(1-F22))^theta-(1-(1-F22))^theta*log(1-(1-F22))*(1-F1)^theta)*((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta)+log(((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta)))*(-dcop.theta2))+1/theta^2*((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*((1-F1)^theta*log(1-F1)+(1-(1-F22))^theta*log(1-(1-F22))-(1-F1)^theta*log(1-F1)*(1-(1-F22))^theta-(1-(1-F22))^theta*log(1-(1-F22))*(1-F1)^theta)-1/theta*(((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta-1)*((1-F1)^theta*log(1-F1)^2+(1-(1-F22))^theta*log(1-(1-F22))^2-(1-F1)^theta*log(1-F1)^2*(1-(1-F22))^theta-(1-(1-F22))^theta*log(1-(1-F22))^2*(1-F1)^theta-2*(1-(1-F22))^theta*log(1-(1-F22))*(1-F1)^theta*log(1-F1))+((1-F1)^theta*log(1-F1)+(1-(1-F22))^theta*log(1-(1-F22))-(1-F1)^theta*log(1-F1)*(1-(1-F22))^theta-(1-(1-F22))^theta*log(1-(1-F22))*(1-F1)^theta)*(-dcop.theta2*((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(-1)-((1-F1)^(theta)+(1-(1-F22))^(theta)-(1-F1)^(theta)*(1-(1-F22))^(theta))^(1/theta-2)*((1-F1)^theta*log(1-F1)+(1-(1-F22))^theta*log(1-(1-F22))-(1-F1)^theta*log(1-F1)*(1-(1-F22))^theta-(1-(1-F22))^theta*log(1-(1-F22))*(1-F1)^theta)))    )
    
    
    
  }}}}}}}}}}}}}}}}
  
  
  
  if (VC$BivD=="FGM") {
    theta.append<-1-theta^2
    theta.append.der<- -2*theta*(1-theta^2)
  } else if (VC$BivD=="N") {
    theta.append<-1-theta^2
    theta.append.der<- -2*theta*(1-theta^2)
  } else if (VC$BivD=="AMH") {
    theta.append<-1-theta^2
    theta.append.der<- -2*theta*(1-theta^2)
  } else if (VC$BivD=="C0") {
    theta.append<-theta-precision
    theta.append.der<- theta-precision
  } else if (VC$BivD=="C90") {
    theta.append<- -(theta-precision)
    theta.append.der<- -(theta-precision)
  } else if (VC$BivD=="C180") {
    theta.append<-theta-precision
    theta.append.der<- theta-precision
  } else if (VC$BivD=="C270") {
    theta.append<- -(theta-precision)
    theta.append.der<- -(theta-precision)
  } else if (VC$BivD=="F") {
    theta.append<-1
    theta.append.der<- 0
  } else if (VC$BivD=="G0") {
    theta.append<-theta-1
    theta.append.der<- theta-1
  } else if (VC$BivD=="G90") {
    theta.append<- -(theta-1)
    theta.append.der<- -(theta-1)
  } else if (VC$BivD=="G180") {
    theta.append<-theta-1
    theta.append.der<- theta-1
  } else if (VC$BivD=="G270") {
    theta.append<- -(theta-1)
    theta.append.der<- -(theta-1)
  } else if (VC$BivD=="J0") {
    theta.append<- theta-1-precision
    theta.append.der<- theta-1-precision
  } else if (VC$BivD=="J90") {
    theta.append<- -(theta-1-precision)
    theta.append.der<- -(theta-1-precision)
  } else if (VC$BivD=="J180") {
    theta.append<- theta-1-precision
    theta.append.der<- theta-1-precision
  } else if (VC$BivD=="J270") {
    theta.append<- -(theta-1-precision)
    theta.append.der<- -(theta-1-precision)
  }
  


  
#-------------------------------------------------------

    # Likelihood

    l.par <- VC$weights*(i0*log(F1) + i1*log(lx))


    # Gradient components:

    dl.1  <- VC$weights*(  i0*1/as.vector(F1)*dF1 + i1*as.vector(1/lx)*as.vector(-dcop1+dcop11)*dF1 )                  # dl.dbe1
    dl.2  <- VC$weights*( i1*as.vector(1/lx)*(df2-as.vector(dcop2)*dF2+as.vector(dcop22)*dF22) )                 # dl.dbe2
    dl.5  <- VC$weights*( i1*(1/lx)*(-dcop.theta1+dcop.theta2)*theta.append )                          # dl.dteta.st
     
    
    if (VC$margins[2]=="P") {
     dl.3  <- NULL                  # dl.sigma.st 
     dl.4  <- NULL                    # dl.nu.st
  } else if (VC$margins[2]=="NB") {
      dl.3  <- VC$weights*( i1*as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma)*sigma )         # dl.sigma.st 
      dl.4  <- NULL                   # dl.nu.st
  } else if (VC$margins[2]=="D") {
      dl.3  <- VC$weights*( i1*as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma)*sigma )                    # dl.sigma.st 
      dl.4  <- VC$weights*( i1*as.vector(1/lx)*(df2.nu-as.vector(dcop2)*dF2.nu+as.vector(dcop22)*dF22.nu)*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2) )      # dl.nu.st
  } else if (VC$margins[2]=="PIG") {
     dl.3  <- VC$weights*( (i1*as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma)*sigma) )      # dl.sigma.st 
     dl.4  <- NULL                   # dl.nu.st
  } else if (VC$margins[2]=="S") {
     dl.3  <- VC$weights*( i1*as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma)*sigma )       # dl.sigma.st 
     dl.4  <- VC$weights*( i1*as.vector(1/lx)*(df2.nu-as.vector(dcop2)*dF2.nu+as.vector(dcop22)*dF22.nu)  )                     # dl.nu.st  
  }
    
    
    # (Minus) Hessian components:

  
  
    d2l.11  <- -VC$weights*( i0*as.vector(-1/(F1^2)*(dF1)^2) + i0*as.vector(1/F1*d2F1ddelta1delta1) + i1*as.vector(-1/(lx^2)*(-dcop1+dcop11)^2*(dF1)^2) + i1*as.vector((1/lx)*((-d2Cdcop12+d2Cdcop112)*dF1^2 + (-dcop1+dcop11)*d2F1ddelta1delta1))  )        # d2l.be1.be1
    d2l.12  <- -VC$weights*i1*( as.vector((-1/(lx^2)*(-dcop1*dF1+dcop11*dF1)*(df2-dcop2*dF2+dcop22*dF22))+1/lx*(-d2Cdcop1dcop2*dF2*dF1+d2Cdcop11dcop22*dF22*dF1))  )                                  # d2l.be1.be2
    d2l.15  <- -VC$weights*i1*( as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(-dcop1*dF1+dcop11*dF1)+1/lx*(-d2Cdcop.theta12*dF1+d2Cdcop.theta22*dF1))*(theta.append)  )                                     # be1.teta.st
  
    d2l.22  <- -VC$weights*i1*( as.vector(-1/(lx^2)*(df2-dcop2*dF2+dcop22*dF22)^2) + as.vector(1/lx*(d2f2delta22-(d2Cdcop22*dF2^2+dcop2*d2F2ddelta22)+(d2Cdcop222*dF22^2+dcop22*d2F22ddelta22)))  )  # d2l.be2.be2
    d2l.25  <- -VC$weights*i1*( as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(df2-dcop2*dF2+dcop22*dF22)+1/lx*(-d2Cdcop1dcop.theta1*dF2+d2Cdcop11dcop.theta2*dF22))*(theta.append)  )                                           # d2l.be2.teta.st

    d2l.55  <- -VC$weights*i1*((-1/(lx^2)*(-dcop.theta1+dcop.theta2)^2+1/lx*(-d2Cdcop2dcop.theta1+d2Cdcop22dcop.theta2))*theta.append^2 + (1/lx)*(-dcop.theta1+dcop.theta2)*theta.append.der  )                                      # d2l.teta.st.teta.st

  
  
  
  if (VC$margins[2]=="P") {
    
    
    d2l.13  <-  NULL                                 # d2l.be1.sigma.st
    d2l.14  <-  NULL                                     # d2l.be1.nu.st
    
    d2l.23  <- NULL                                           # d2l.be2.k.st
    d2l.24  <- NULL                                           # d2l.be2.nu.st
    
    d2l.33  <- NULL    # d2l.sigma.st.k.st
    d2l.34  <- NULL                                                                  # d2l.sigma.st.nu.st
    d2l.35  <- NULL                                                                 # d2l.sigma.st.teta.st
    
    d2l.44  <- NULL                                      # d2l.nu.st.nu.st
    d2l.45  <- NULL                                       # d2l.nu.st.teta.st
      
      
    
  } else if (VC$margins[2]=="NB") {
    
    
    d2l.13  <- -VC$weights*i1*(as.vector((-1/(lx^2)*(-dcop1*dF1+dcop11*dF1)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma))+1/lx*(-d2Cdcop1dcop2*dF2.sigma*dF1+d2Cdcop11dcop22*dF22.sigma*dF1))*sigma)                                   # d2l.be1.sigma.st
    d2l.14  <- NULL                                      # d2l.be1.nu.st
    
    d2l.23  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2-dcop2*dF2+dcop22*dF22)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)) + as.vector(1/lx*(d2f2delta2sigma-(d2Cdcop22*dF2*dF2.sigma+dcop2*d2F2ddelta2dsigma)+(d2Cdcop222*dF22*dF22.sigma+dcop22*d2F22ddelta2dsigma))))*sigma)                                            # d2l.be2.k.st
    d2l.24  <-  NULL                                          # d2l.be2.nu.st
    
    d2l.33  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)^2)+as.vector(1/lx*(d2f2sigma2-(d2Cdcop22*dF2.sigma^2+dcop2*d2F2dsigma2)+(d2Cdcop222*dF22.sigma^2+dcop22*d2F22dsigma2))))*sigma^2+(as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma))*sigma )    # d2l.sigma.st.k.st
    d2l.34  <-  NULL                                                                 # d2l.sigma.st.nu.st
    d2l.35  <- -VC$weights*i1*(as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)+1/lx*(-d2Cdcop1dcop.theta1*dF2.sigma+d2Cdcop11dcop.theta2*dF22.sigma))*(theta.append)*sigma)                                                                # d2l.sigma.st.teta.st
    
    d2l.44  <-  NULL                                       # d2l.nu.st.nu.st
    d2l.45  <-  NULL                                     # d2l.nu.st.teta.st
    
    
  } else if (VC$margins[2]=="D") {
    
    d2l.13  <- -VC$weights*i1*(as.vector((-1/(lx^2)*(-dcop1*dF1+dcop11*dF1)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma))+1/lx*(-d2Cdcop1dcop2*dF2.sigma*dF1+d2Cdcop11dcop22*dF22.sigma*dF1))*sigma)                                   # d2l.be1.sigma.st
    d2l.14  <- -VC$weights*i1*(as.vector((-1/(lx^2)*(-dcop1*dF1+dcop11*dF1)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu))+1/lx*(-d2Cdcop1dcop2*dF2.nu*dF1+d2Cdcop11dcop22*dF22.nu*dF1))*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2))                                      # d2l.be1.nu.st
    
    d2l.23  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2-dcop2*dF2+dcop22*dF22)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)) + as.vector(1/lx*(d2f2delta2sigma-(d2Cdcop22*dF2*dF2.sigma+dcop2*d2F2ddelta2dsigma)+(d2Cdcop222*dF22*dF22.sigma+dcop22*d2F22ddelta2dsigma))))*sigma)                                            # d2l.be2.k.st
    d2l.24  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2-dcop2*dF2+dcop22*dF22)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu))+as.vector(1/lx*(d2f2delta2nu-(d2Cdcop22*dF2*dF2.nu+dcop2*d2F2ddelta2dnu)+(d2Cdcop222*dF22*dF22.nu+dcop22*d2F22ddelta2dnu))))*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2) )                                           # d2l.be2.nu.st
    
    d2l.33  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)^2)+as.vector(1/lx*(d2f2sigma2-(d2Cdcop22*dF2.sigma^2+dcop2*d2F2dsigma2)+(d2Cdcop222*dF22.sigma^2+dcop22*d2F22dsigma2))))*sigma^2+(as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma))*sigma )    # d2l.sigma.st.k.st
    d2l.34  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu))+as.vector(1/lx*(d2f2nusigma-(d2Cdcop22*dF2.sigma*dF2.nu+dcop2*d2F2dnudsigma)+(d2Cdcop222*dF22.sigma*dF22.nu+dcop22*d2F22dnudsigma))))*sigma*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2) )                                                                 # d2l.sigma.st.nu.st
    d2l.35  <- -VC$weights*i1*(as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)+1/lx*(-d2Cdcop1dcop.theta1*dF2.sigma+d2Cdcop11dcop.theta2*dF22.sigma))*(theta.append)*sigma)                                                                # d2l.sigma.st.teta.st
    
    d2l.44  <-  -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu)^2)+as.vector(1/lx*(d2f2nu2-(d2Cdcop22*dF2.nu^2+dcop2*d2F2dnu2)+(d2Cdcop222*dF22.nu^2+dcop22*d2F22dnu2))))*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2)^2+(as.vector(1/lx)*(df2.nu-as.vector(dcop2)*dF2.nu+as.vector(dcop22)*dF22.nu))*((nu/(1-nu))*(1+(nu/(1-nu)))^2-2*(1+(nu/(1-nu)))*(nu/(1-nu))^2)/(1+nu/(1-nu))^4  )                                       # d2l.nu.st.nu.st
    d2l.45  <-  -VC$weights*i1*(as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu)+1/lx*(-d2Cdcop1dcop.theta1*dF2.nu+d2Cdcop11dcop.theta2*dF22.nu))*theta.append*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2)  )                                     # d2l.nu.st.teta.st
      
    
  } else if (VC$margins[2]=="PIG") {
    
    d2l.13  <- -VC$weights*i1*(as.vector((-1/(lx^2)*(-dcop1*dF1+dcop11*dF1)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma))+1/lx*(-d2Cdcop1dcop2*dF2.sigma*dF1+d2Cdcop11dcop22*dF22.sigma*dF1))*sigma)                                   # d2l.be1.sigma.st
    d2l.14  <- NULL                                      # d2l.be1.nu.st
    
    d2l.23  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2-dcop2*dF2+dcop22*dF22)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)) + as.vector(1/lx*(d2f2delta2sigma-(d2Cdcop22*dF2*dF2.sigma+dcop2*d2F2ddelta2dsigma)+(d2Cdcop222*dF22*dF22.sigma+dcop22*d2F22ddelta2dsigma))))*sigma)                                            # d2l.be2.k.st
    d2l.24  <-  NULL                                          # d2l.be2.nu.st
    
    d2l.33  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)^2)+as.vector(1/lx*(d2f2sigma2-(d2Cdcop22*dF2.sigma^2+dcop2*d2F2dsigma2)+(d2Cdcop222*dF22.sigma^2+dcop22*d2F22dsigma2))))*sigma^2+(as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma))*sigma )    # d2l.sigma.st.k.st
    d2l.34  <-  NULL                                                                 # d2l.sigma.st.nu.st
    d2l.35  <- -VC$weights*i1*(as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)+1/lx*(-d2Cdcop1dcop.theta1*dF2.sigma+d2Cdcop11dcop.theta2*dF22.sigma))*(theta.append)*sigma)                                                                # d2l.sigma.st.teta.st
    
    d2l.44  <-  NULL                                       # d2l.nu.st.nu.st
    d2l.45  <-  NULL                                     # d2l.nu.st.teta.st
                                             
    
  } else if (VC$margins[2]=="S") {
    
    d2l.13  <- -VC$weights*i1*(as.vector((-1/(lx^2)*(-dcop1*dF1+dcop11*dF1)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma))+1/lx*(-d2Cdcop1dcop2*dF2.sigma*dF1+d2Cdcop11dcop22*dF22.sigma*dF1))*sigma)                                   # d2l.be1.sigma.st
    d2l.14  <- -VC$weights*i1*( as.vector((-1/(lx^2)*(-dcop1*dF1+dcop11*dF1)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu))+1/lx*(-d2Cdcop1dcop2*dF2.nu*dF1+d2Cdcop11dcop22*dF22.nu*dF1)) )                                      # d2l.be1.nu.st
    
    d2l.23  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2-dcop2*dF2+dcop22*dF22)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)) + as.vector(1/lx*(d2f2delta2sigma-(d2Cdcop22*dF2*dF2.sigma+dcop2*d2F2ddelta2dsigma)+(d2Cdcop222*dF22*dF22.sigma+dcop22*d2F22ddelta2dsigma))))*sigma)                                            # d2l.be2.k.st
    d2l.24  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2-dcop2*dF2+dcop22*dF22)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu))+as.vector(1/lx*(d2f2delta2nu-(d2Cdcop22*dF2*dF2.nu+dcop2*d2F2ddelta2dnu)+(d2Cdcop222*dF22*dF22.nu+dcop22*d2F22ddelta2dnu))))  )                                           # d2l.be2.nu.st
    
    d2l.33  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)^2)+as.vector(1/lx*(d2f2sigma2-(d2Cdcop22*dF2.sigma^2+dcop2*d2F2dsigma2)+(d2Cdcop222*dF22.sigma^2+dcop22*d2F22dsigma2))))*sigma^2+(as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma))*sigma )    # d2l.sigma.st.k.st
    d2l.34  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu))+as.vector(1/lx*(d2f2sigmanu-(d2Cdcop22*dF2.sigma*dF2.nu+dcop2*d2F2dsigmadnu)+(d2Cdcop222*dF22.sigma*dF22.nu+dcop22*d2F22dsigmadnu))))*sigma )                                                                 # d2l.sigma.st.nu.st
    d2l.35  <- -VC$weights*i1*(as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)+1/lx*(-d2Cdcop1dcop.theta1*dF2.sigma+d2Cdcop11dcop.theta2*dF22.sigma))*(theta.append)*sigma)                                                                # d2l.sigma.st.teta.st
    
    d2l.44  <-  -VC$weights*i1*(as.vector(-1/(lx^2)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu)^2)+as.vector(1/lx*(d2f2nu2-(d2Cdcop22*dF2.nu^2+dcop2*d2F2dnu2)+(d2Cdcop222*dF22.nu^2+dcop22*d2F22dnu2)))  )                                       # d2l.nu.st.nu.st
    d2l.45  <-  -VC$weights*i1*(as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu)+1/lx*(-d2Cdcop1dcop.theta1*dF2.nu+d2Cdcop11dcop.theta2*dF22.nu))*theta.append )                                     # d2l.nu.st.teta.st
      
    
    
  }


#---------------------------------------------------------------------

  
  res <- -sum(l.par)
  
  
  
  # Creating gradient and Hessian
  
  if (VC$margins[2]=="P") {
     
    H11 <- crossprod(VC$X1*c(d2l.11),VC$X1)
    H12 <- crossprod(VC$X1*c(d2l.12),VC$X2) 
    H15 <- t(t(rowSums(t(VC$X1*c(d2l.15)))))

    H22 <- crossprod(VC$X2*c(d2l.22),VC$X2) 
    H25 <- t(t(rowSums(t(VC$X2*c(d2l.25)))))

    H <- rbind( cbind( H11    , H12    ,  H15  ), 
              cbind( t(H12) , H22    ,  H25  ),
              cbind( t(H15) , t(H25) ,  sum(d2l.55) )
            ) 


    G   <- c( -colSums( c(dl.1)*VC$X1 ) ,
            -colSums( c(dl.2)*VC$X2 )    ,   
            -sum( dl.5 )                 )
    
    
  } 
  
  
  if (VC$margins[2] %in% c("NB","PIG") ) {
      
    H11 <- crossprod(VC$X1*c(d2l.11),VC$X1)
    H12 <- crossprod(VC$X1*c(d2l.12),VC$X2) 
    H13 <- t(t(rowSums(t(VC$X1*c(d2l.13)))))
    H15 <- t(t(rowSums(t(VC$X1*c(d2l.15)))))

    H22 <- crossprod(VC$X2*c(d2l.22),VC$X2) 
    H23 <- t(t(rowSums(t(VC$X2*c(d2l.23)))))
    H25 <- t(t(rowSums(t(VC$X2*c(d2l.25)))))

    H <- rbind( cbind( H11    , H12    , H13  , H15  ), 
              cbind( t(H12) , H22    , H23  ,  H25  ),
              cbind( t(H13) , t(H23) , sum(d2l.33),  sum(d2l.35) ) ,
              cbind( t(H15) , t(H25) , sum(d2l.35),  sum(d2l.55) )
            ) 


    G   <- c( -colSums( c(dl.1)*VC$X1 ) ,
            -colSums( c(dl.2)*VC$X2 )    ,
            -sum( dl.3 ) ,    
            -sum( dl.5 )                 )
    
  } 
  
  if (VC$margins[2] %in% c("D","S") ) {
      
    H11 <- crossprod(VC$X1*c(d2l.11),VC$X1)
    H12 <- crossprod(VC$X1*c(d2l.12),VC$X2) 
    H13 <- t(t(rowSums(t(VC$X1*c(d2l.13)))))
    H14 <- t(t(rowSums(t(VC$X1*c(d2l.14)))))
    H15 <- t(t(rowSums(t(VC$X1*c(d2l.15)))))

    H22 <- crossprod(VC$X2*c(d2l.22),VC$X2) 
    H23 <- t(t(rowSums(t(VC$X2*c(d2l.23)))))
    H24 <- t(t(rowSums(t(VC$X2*c(d2l.24)))))
    H25 <- t(t(rowSums(t(VC$X2*c(d2l.25)))))

    H <- rbind( cbind( H11    , H12    , H13  ,  H14 , H15  ), 
              cbind( t(H12) , H22    , H23  ,  H24 , H25  ),
              cbind( t(H13) , t(H23) , sum(d2l.33), sum(d2l.34), sum(d2l.35) ) ,
              cbind( t(H14) , t(H24) , sum(d2l.34), sum(d2l.44), sum(d2l.45) ) ,
              cbind( t(H15) , t(H25) , sum(d2l.35), sum(d2l.45), sum(d2l.55) )
            ) 


    G   <- c( -colSums( c(dl.1)*VC$X1 ) ,
            -colSums( c(dl.2)*VC$X2 )    ,
            -sum( dl.3 ) ,  
            -sum( dl.4 ) ,   
            -sum( dl.5 )                 )
  
  }
  

  
  
  
  

  
  if( ( VC$l.sp1==0 && VC$l.sp2==0 ) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC)
     
     
    S.res <- res
    res <- S.res + ps$S.h1
    G   <- G + ps$S.h2
    H   <- H + ps$S.h  

  list(value=res, gradient=G, hessian=H, S.h=ps$S.h, S.h2 = ps$S.h2, l=S.res, eta1=eta1, eta2=eta2,
       dl.dbe1=dl.1, dl.dbe2=dl.2, l.par=l.par,
       dl.dsqv.st=dl.3,
       dl.dsh.st=dl.4,
       dl.dcor.st=dl.5, 
       d2l.be1.be1=d2l.11, d2l.be1.be2=d2l.12, d2l.be2.be2=d2l.22,
       d2l.be1.sqv.st=d2l.13,
       d2l.be1.sh.st=d2l.14,
       d2l.be1.cor.st=d2l.15,
       d2l.be2.sqv.st=d2l.23, 
       d2l.be2.sh.st=d2l.24,
       d2l.be2.cor.st=d2l.25,
       d2l.sqv.st.sqv.st=d2l.33,
       d2l.sqv.st.sh.st=d2l.34,
       d2l.sqv.st.cor.st=d2l.35,
       d2l.sh.st.sh.st=d2l.44, 
       d2l.sh.st.cor.st=d2l.45,
       d2l.cor.st.cor.st=d2l.55  )

}





