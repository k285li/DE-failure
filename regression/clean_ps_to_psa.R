library(dplyr)
library(msm)
library(survival)
library(xtable)
library(lubridate)
library(parallel)
library(nleqslv)
library(tibble)
library(MASS)

cores.num <- 1

logLR1.f <- function(par, indata1){
  
  lam12c <- exp(par[1])
  b12c1 <- par[2]
  b12c2 <- par[3]
  b12c3 <- par[4]
  rho1 <- exp(par[5])
  
  L.retro1  <- mclapply(sort(unique(as.character(indata1$id))),
                        function(whichsubj, locdata) {
                          locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
                          a1r <- locdata$age_first_visit
                          a1 <- locdata$AGE_PS
                          A1E <- locdata$A1E
                          
                          X1 <- locdata$AGE_PS-40
                          # X1 <- locdata$AGE_early
                          X2 <- locdata$gender
                          X3 <- locdata$HLAB27
                          
                          lam12c_X <- lam12c*exp(b12c1*X1+b12c2*X2+b12c3*X3)
                          
                          Li <- (rho1 + lam12c_X)*exp(-(rho1 + lam12c_X)*(a1r-a1))/(1-exp(-(rho1 + lam12c_X)*(A1E-a1)))
                          return(Li) 
                        }, 
                        locdata=indata1, mc.cores=cores.num)
  
  
  
  log.retro1 <- log(unlist(L.retro1))
  
  logL <- sum(log.retro1, na.rm = TRUE)
  return( (-1)*logL )
  
}




logLR2.f <- function(par, indata2){
  lam12c <- exp(par[1])
  b12c1 <- par[2]
  b12c2 <- par[3]
  b12c3 <- par[4]
  rho1 <- exp(par[5])
  rho2 <- exp(par[6])
  
  L.retro2  <- mclapply(sort(unique(as.character(indata2$id))),
                        function(whichsubj, locdata) {
                          locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
                          a1 <- locdata$AGE_PS
                          a2 <- locdata$AGE_PSA
                          a2r <- locdata$age
                          A2E <- locdata$A2E
                          
                          X1 <- locdata$AGE_PS-40
                          # X1 <- locdata$AGE_early
                          X2 <- locdata$gender
                          X3 <- locdata$HLAB27
                          
                          lam12c_X <- lam12c*exp(b12c1*X1+b12c2*X2+b12c3*X3)
                          
                          Li <- ( ((rho1 + lam12c_X)-rho2)*exp(-(rho1 + lam12c_X)*(a2-a1))*exp(-rho2*(a2r-a2))  )/( 1/rho2*(1-exp(-rho2*(A2E-a1)) ) - 1/(rho1 + lam12c_X)*(1-exp(-(rho1 + lam12c_X)*(A2E-a1)) ) )
                          return(Li) 
                        }, 
                        locdata=indata2,mc.cores=cores.num)
  
  log.retro2 <- log(unlist(L.retro2))
  
  logL <- sum(log.retro2, na.rm = TRUE)
  return( (-1)*logL )
  
}


logLR.f <- function(par, indata1, indata2){
  lam12c <- exp(par[1])
  b12c1 <- par[2]
  b12c2 <- par[3]
  b12c3 <- par[4]
  rho1 <- exp(par[5])
  rho2 <- exp(par[6])

  L.retro1  <- mclapply(sort(unique(as.character(indata1$id))),
                        function(whichsubj, locdata) {
                          locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
                          a1r <- locdata$age_first_visit
                          a1 <- locdata$AGE_PS
                          A1E <- locdata$A1E
                          
                          X1 <- locdata$AGE_PS-40
                          # X1 <- locdata$AGE_early
                          X2 <- locdata$gender
                          X3 <- locdata$HLAB27
                          
                          lam12c_X <- lam12c*exp(b12c1*X1+b12c2*X2+b12c3*X3)
                          
                          Li <- (rho1 + lam12c_X)*exp(-(rho1 + lam12c_X)*(a1r-a1))/(1-exp(-(rho1 + lam12c_X)*(A1E-a1)))
                          return(Li) 
                        }, 
                        locdata=indata1, mc.cores=cores.num)
  
  log.retro1 <- log(unlist(L.retro1))
  

  L.retro2  <- mclapply(sort(unique(as.character(indata2$id))),
                        function(whichsubj, locdata) {
                          locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
                          a1 <- locdata$AGE_PS
                          a2 <- locdata$AGE_PSA
                          a2r <- locdata$age
                          A2E <- locdata$A2E
                          
                          X1 <- locdata$AGE_PS-40
                          #X1 <- locdata$AGE_early
                          X2 <- locdata$gender
                          X3 <- locdata$HLAB27
                          
                          lam12c_X <- lam12c*exp(b12c1*X1+b12c2*X2+b12c3*X3)
                          
                          Li <- ( ((rho1 + lam12c_X)-rho2)*exp(-(rho1 + lam12c_X)*(a2-a1))*exp(-rho2*(a2r-a2))  )/( 1/rho2*(1-exp(-rho2*(A2E-a1)) ) - 1/(rho1 + lam12c_X)*(1-exp(-(rho1 + lam12c_X)*(A2E-a1)) ) )
                          return(Li) 
                        }, 
                        locdata=indata2,mc.cores=cores.num)
  
  log.retro2 <- log(unlist(L.retro2))
  
  logL <- sum(log.retro1, log.retro2, na.rm = TRUE)
  return( (-1)*logL )
  
}



logLRLP.f <- function(par, indata1, indata2){
  lam12c <- exp(par[1])
  b12c1 <- par[2]
  b12c2 <- par[3]
  b12c3 <- par[4]
  rho1 <- exp(par[5])
  rho2 <- exp(par[6])
  
  L.retro1  <- mclapply(sort(unique(as.character(indata1$id))),
                        function(whichsubj, locdata) {
                          locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
                          a1r <- locdata$age_first_visit
                          a1 <- locdata$AGE_PS
                          A1E <- locdata$A1E
                          
                          X1 <- locdata$AGE_PS-40
                          # X1 <- locdata$AGE_early
                          X2 <- locdata$gender
                          X3 <- locdata$HLAB27
                          
                          lam12c_X <- lam12c*exp(b12c1*X1+b12c2*X2+b12c3*X3)
                          
                          Li <- (rho1 + lam12c_X)*exp(-(rho1 + lam12c_X)*(a1r-a1))/(1-exp(-(rho1 + lam12c_X)*(A1E-a1)))
                          return(Li) 
                        }, 
                        locdata=indata1, mc.cores=cores.num)
  
  log.retro1 <- log(unlist(L.retro1))
  
  L.retro2 <- mclapply(sort(unique(as.character(indata2$id))),
                        function(whichsubj, locdata) {
                          locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
                          a1 <- locdata$AGE_PS
                          a2 <- locdata$AGE_PSA
                          a2r <- locdata$age
                          A2E <- locdata$A2E
                          
                          X1 <- locdata$AGE_PS-40
                          # X1 <- locdata$AGE_early
                          X2 <- locdata$gender
                          X3 <- locdata$HLAB27
                          
                          lam12c_X <- lam12c*exp(b12c1*X1+b12c2*X2+b12c3*X3)
                          
                          Li <- ( ((rho1 + lam12c_X)-rho2)*exp(-(rho1 + lam12c_X)*(a2-a1))*exp(-rho2*(a2r-a2))  )/( 1/rho2*(1-exp(-rho2*(A2E-a1)) ) - 1/(rho1 + lam12c_X)*(1-exp(-(rho1 + lam12c_X)*(A2E-a1)) ) )
                          return(Li) 
                        }, 
                        locdata=indata2, mc.cores=cores.num)
  
  log.retro2 <- log(unlist(L.retro2))
  
  
  L.prosp1  <- mclapply(sort(unique(as.character(indata1$id))),
                        function(whichsubj, locdata) {
                          locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
                          a1r <- locdata$age_first_visit
                          v <- locdata$time
                          delta2 <- locdata$status
                          
                          X1 <- locdata$AGE_PS-40
                          # X1 <- locdata$AGE_early
                          X2 <- locdata$gender
                          X3 <- locdata$HLAB27
                          
                          lam12_X <- lam12c*exp(b12c1*X1+b12c2*X2+b12c3*X3)
                          
                          Li <- (lam12_X^delta2)*exp(-lam12_X*(v-a1r))
                          return(Li)
                        },
                        locdata=indata1, mc.cores=cores.num)
  
  log.prosp1 <- log(unlist(L.prosp1))
  
  logL <- sum(log.retro1, log.retro2, log.prosp1, na.rm = TRUE)
  return( (-1)*logL )
  
}


logLP.f <- function(par, indata1){
  lam12 <- exp(par[1])
  
  b12.1 <- par[2]
  b12.2 <- par[3]
  b12.3 <- par[4]
  
  L.prosp1  <- mclapply(sort(unique(as.character(indata1$id))),
                        function(whichsubj, locdata) {
                          locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
                          a1r <- locdata$age_first_visit
                          v <- locdata$time
                          delta2 <- locdata$status
                          
                          X1 <- locdata$AGE_PS-40
                          # X1 <- locdata$AGE_early
                          X2 <- locdata$gender
                          X3 <- locdata$HLAB27
                          
                          lam12_X <- lam12*exp(b12.1*X1+b12.2*X2+b12.3*X3)
                          
                          Li <- (lam12_X^delta2)*exp(-lam12_X*(v-a1r))
                          return(Li)
                        },
                        locdata=indata1, mc.cores=cores.num)
  
  log.prosp1 <- log(unlist(L.prosp1))
  
  logL <- sum(log.prosp1, na.rm = TRUE)
  return( (-1)*logL )
  
}



logLR2.plot.f <- function(par, indata2, rho1){
  lam12c <- exp(par[1])
  b12c1 <- par[2]
  b12c2 <- par[3]
  b12c3 <- par[4]
  rho1 <- rho1
  rho2 <- exp(par[5])
  
  L.retro2  <- mclapply(sort(unique(as.character(indata2$id))),
                        function(whichsubj, locdata) {
                          locdata <- locdata[as.character(locdata$id) == as.character(whichsubj),]
                          a1 <- locdata$AGE_PS
                          a2 <- locdata$AGE_PSA
                          a2r <- locdata$age
                          A2E <- locdata$A2E
                          
                          X1 <- locdata$AGE_PS-40
                          # X1 <- locdata$AGE_early
                          X2 <- locdata$gender
                          X3 <- locdata$HLAB27
                          
                          lam12c_X <- lam12c*exp(b12c1*X1+b12c2*X2+b12c3*X3)
                          
                          Li <- ( ((rho1 + lam12c_X)-rho2)*exp(-(rho1 + lam12c_X)*(a2-a1))*exp(-rho2*(a2r-a2))  )/( 1/rho2*(1-exp(-rho2*(A2E-a1)) ) - 1/(rho1 + lam12c_X)*(1-exp(-(rho1 + lam12c_X)*(A2E-a1)) ) )
                          return(Li) 
                        }, 
                        locdata=indata2,mc.cores=cores.num)
  
  log.retro2 <- log(unlist(L.retro2))
  
  logL <- sum(log.retro2, na.rm = TRUE)
  return( (-1)*logL )
  
}



