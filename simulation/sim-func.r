# =====================================================
# Some shared functions
#
# Last Update: July 13, 2022
# =====================================================


library(parallel)
library(rootSolve)
library(msm)
library(MASS)


qmat.f <- function(lams) {
  qmat <- matrix(0, nrow=2, ncol=2)
  qmat[1,2]  <- lams
  diag(qmat) <- (-1)*apply(qmat, 1, sum)
  return(qmat)
}


qmat.expanded.f <- function(lams, rhos, lamsp) {
  qmat <- matrix(0, nrow=4, ncol=4)
  qmat[1,c(2,3)] <- c(lams, rhos)
  qmat[3,4]      <- lamsp
  diag(qmat)     <- (-1)*apply(qmat, 1, sum)
  return(qmat)
}


getstats.f <- function(pi1, pi1o) {
  A <- 1

  func.lams.f <- function(p, A, const) {
    lam <- exp(p)

    qq <- qmat.f(lams=lam)
    PP <- MatrixExp(qq, t=A)
    eq <- PP[1,ncol(PP)] - const
    return(eq)
  }

  detcpoint.f <- function(p, lam, const) {
    tt <- p

    qq <- qmat.f(lams=lam)
    PP <- MatrixExp(qq, t=tt)
    eq <- PP[1,ncol(PP)] - const
    return(eq)
  }

  func.rho.f <- function(p, lam, tau50, tau95, const) {
    rho <- exp(p)

    term <- rho/(lam + rho)
    term95 <- tau95 + (exp((-1)*(lam + rho)*tau95)/(lam + rho))
    term50 <- tau50 + (exp((-1)*(lam + rho)*tau50)/(lam + rho))
    eq <- term*(term95 - term50) - const
    return(eq)
  }


  est <- uniroot(f=func.lams.f, interval=c(-10, 10),
                 A=A, const=pi1)$root
  lam <- exp(est)

  tau50 <- uniroot(detcpoint.f, interval=c(0,1000),
                   lam=lam, const=0.5)$root

  tau95 <- uniroot(detcpoint.f, interval=c(0,1000),
                   lam=lam, const=0.95)$root

  est <- uniroot(func.rho.f, interval=c(-10, 10),
                 lam=lam, tau50=tau50, tau95=tau95, const=pi1o)$root
  rho <- exp(est)

  outstats <- data.frame(lam, rho, A, tau50, tau95)
  print(outstats)
}


generatedata.manual.f <- function(npop, nsamples, lam, lamp, rho, A, tau50, tau95) { 
  theta <- lam + rho
  pii   <- lam/theta

  Ristar <- runif(npop, min=tau50, max=tau95)
  deltai <- rbinom(npop, size=1, prob=((1 - pii)*(1 - exp((-1)*theta*Ristar))))

  matR <- data.frame(pid=c(1:npop), Rstar=Ristar, delta=deltai)

  recruit <- matR[matR$delta == 1,]

  rlen <- length(recruit$pid)
  ui <- runif(rlen, min=0, max=1)
  Ri <- (-1)*(1/theta)*log( 1 - ui*(1 - exp((-1)*theta*recruit$Rstar)) )
  Ti <- Ri + rexp(rlen, rate=lam)

  timeR  <- cbind(rep(0, rlen), Ri)
  stateR <- cbind(rep(1, rlen), rep(3, rlen)) 
  nR     <- rep(2, rlen)
  retro  <- cbind(1:rlen, recruit$pid, 0, timeR, stateR, nR, recruit$Rstar, Ti)

  timeP  <- cbind(Ri, Ti)
  stateP <- cbind(rep(3, rlen), rep(4, rlen)) 
  nP     <- rep(2, rlen)
  prosp  <- cbind(1:rlen, recruit$pid, 1, timeP, stateP, nP, recruit$Rstar, Ti)
  
  events <- data.frame(rbind(retro, prosp))
  dimnames(events)[[2]] <- c("id","pid", "prosp", paste("time",1:2,sep=""),
                             paste("state",1:2,sep=""), "n", "Rstar","tau")
  events <- events[order(events$id, events$prosp),]
  events <- events[events$id <= nsamples,]
 
  cutoff <- max(events$pid)

  rej <- matR[matR$pid <= cutoff,]
  rej <- rej[rej$delta == 0,]
  nrej <- length(rej$pid)
  rej$time2  <- rexp(nrej, rate=lam)
  rej$state2 <- rep(2, nrej)

  tmpdata <- rbind(events[events$prosp == 0, c("time2","state2")],
                   rej[,c("time2","state2")])
  tmpdata <- tmpdata[tmpdata$time2 <= A,]
  nsubj <- length(unique(c(events$pid, rej$pid))) 

  cat("pi1 = ", nrow(tmpdata)/nsubj, sep="", "\n")

  cat("pi1o = ", nsamples/nsubj, sep="", "\n")

  outdata <- NULL
  outdata$recruited <- events
  outdata$rejected  <- rej$Rstar
  return(outdata)
}

