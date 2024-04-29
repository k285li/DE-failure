# =====================================================
# Likelihood Functions
#
# Last Update: August 2, 2022
# =====================================================


expit.f <- function(x) { return( exp(x)/(1 + exp(x)) ) }
logit.f <- function(x) { return( log( x/(1 - x) ) ) }


logL.obs.f <- function(p, indata) {
  lam <- exp(p[1])
  rho <- exp(p[2])

  theta <- lam + rho
  pii   <- lam/theta

  prosp <- indata[indata$prosp == 1,]

  term1 <- log(rho) + ((-1)*rho*prosp$time1)
  term2 <- log(lam) + ((-1)*lam*prosp$time2)
  termRstar <- log(rho) - log(theta) + log(1 - exp((-1)*theta*prosp$Rstar))

  logL <- (-1)*sum(term1 + term2 - termRstar)
  logL <- ifelse(is.finite(logL), logL, NaN)
  
#  print(data.frame(lam, rho, logL))
  return(logL)
}


logL.obs.full.f <- function(p, indata) {
  lam0 <- exp(p[1])
  rho  <- exp(p[2])
  lam  <- exp(p[3])

  theta <- lam0 + rho
  pii   <- lam0/theta

  prosp <- indata[indata$prosp == 1,]

  term1 <- log(rho) + ((-1)*rho*prosp$time1)
  term2 <- log(lam) + ((-1)*lam*prosp$time2)
  termRstar <- log(rho) - log(theta) + log(1 - exp((-1)*theta*prosp$Rstar))

  logL <- (-1)*sum(term1 + term2 - termRstar)
  logL <- ifelse(is.finite(logL), logL, NaN)
  
  print(data.frame(lam0, rho, lam, logL))
  return(logL)
}


logL.comp.f <- function(p, indata, inrej) {
  lam <- exp(p[1])
  rho <- exp(p[2])

  theta <- lam + rho
  pii   <- lam/theta

  prosp <- indata[indata$prosp == 1,]

  term1 <- rho*exp((-1)*rho*prosp$time1)
  term2 <- lam*exp((-1)*lam*prosp$time2)

  termRstar <- 1 - ((1 - pii)*(1 - exp((-1)*theta*inrej)))

  logL <- (-1)*( sum(log(term1)) + sum(log(term2)) + sum(log(termRstar)) )
  logL <- ifelse(is.finite(logL), logL, NaN)
  return(logL)
}


logL.comp.full.f <- function(p, indata, inrej) {
  lam0 <- exp(p[1])
  rho  <- exp(p[2])
  lam  <- exp(p[3])

  prosp <- indata[indata$prosp == 1,]

  term1 <- log(rho) + ((-1)*rho*prosp$time1)
  term2 <- (-1)*lam0*prosp$time1
  term3 <- log(lam) + ((-1)*lam*(prosp$time2 - prosp$time1))

  theta <- lam0 + rho
  pii   <- lam0/theta

  termRstar <- log( (1 - pii)*(1 - exp((-1)*theta*prosp$Rstar)) )

  logL1 <- sum(term1 + term2 + term3) - sum(termRstar)


  termS <- (1 - pii)*(1 - exp((-1)*theta*prosp$Rstar))
  termR <- 1 - ( (1 - pii)*(1 - exp((-1)*theta*inrej)) )

  logL2 <- sum(log(termS)) + sum(log(termR))

  logL <- (-1)*(logL1 + logL2)
  logL <- ifelse(is.finite(logL), logL, NaN)
  return(logL)
}


logL.retro.f <- function(p, indata, inrej) {
  lam <- exp(p[1])
  rho <- exp(p[2])

  theta <- lam + rho
  pii   <- lam/theta

  prosp <- indata[indata$prosp == 1,]

  term <- rho*exp((-1)*theta*prosp$time1)
  termRstar <- 1 - ( (1 - pii)*(1 - exp((-1)*theta*inrej)) )

  logL <- (-1)*( sum(log(term)) + sum(log(termRstar)) )
  logL <- ifelse(is.finite(logL), logL, NaN)
  return(logL)
}


logL.bin.f <- function(p, indata, inrej) {
  lam <- exp(p[1])
  rho <- exp(p[2])

  theta <- lam + rho
  pii   <- lam/theta

  prosp <- indata[indata$prosp == 1,]

  termS <- (1 - pii)*(1 - exp((-1)*theta*prosp$Rstar))
  termR <- 1 - ( (1 - pii)*(1 - exp((-1)*theta*inrej)) )

  logL <- (-1)*( sum(log(termS)) + sum(log(termR)) )
  logL <- ifelse(is.finite(logL), logL, NaN)
  return(logL)
}


#####################################################
## As Checks
#####################################################
#
# logL.obs.theta.and.rho.f <- function(p, indata) {
#  theta <- exp(p[1])
#  rho   <- exp(p[2])
#
#  lam   <- theta - rho
#  pii   <- lam/theta
#
#  prosp <- indata[indata$prosp == 1,]
# 
#  term1 <- log(rho) + ((-1)*rho*prosp$time1)
#  term2 <- log(lam) + ((-1)*lam*prosp$time2)
#  termRstar <- log(rho) - log(theta) + log(1 - exp((-1)*theta*prosp$Rstar))
#
#  logL <- (-1)*sum(term1 + term2 - termRstar)
#  logL <- ifelse(is.finite(logL), logL, NaN)
#  
#  print(data.frame(lam, rho, logL))
#  return(logL)
#}
#
#
#logL.obs.retro.theta.f <- function(p, indata) {
#  theta <- exp(p)
#
#  prosp <- indata[indata$prosp == 1,]
#
#  term1 <- log(theta) + ((-1)*theta*prosp$time1)
#  term2 <- log(1 - exp((-1)*theta*prosp$Rstar))
#  logL <- (-1)*sum(term1 - term2)
#  logL <- ifelse(is.finite(logL), logL, NaN)
#  
#  print(data.frame(lam, rho, logL))
#  return(logL)
#}
#
#
#logL.obs.prosp.lam.f <- function(p, indata) {
#  lam <- exp(p)
#
#  prosp <- indata[indata$prosp == 1,]
#
#  logL <- (-1)*sum( log(lam) + ((-1)*lam*(prosp$time2 - prosp$time1)) )
#  logL <- ifelse(is.finite(logL), logL, NaN)
#  
#  print(data.frame(lam, rho, logL))
#  return(logL)
#}
#
#####################################################

