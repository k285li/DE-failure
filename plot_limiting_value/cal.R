# ==============================================
# Limiting value calculation and plot illustration
#
# Last Update: April 29, 2024
# ==============================================

library(msm)
library(nleqslv)
library(dplyr)
library(survival)
library(latex2exp)

cal.f <- function(ebeta, p1, egamma, pi.F.const, pi.S.const){
  
  fn1 <- function(x) {
    Q <- rbind(c(-(x+x*(1-p1)/p1), x, x*(1-p1)/p1),
               c(0, -(x*(1-p1)/p1*ebeta), x*(1-p1)/p1*ebeta),
               c(0, 0, 0)) 
    
    
    pi.F <- MatrixExp(Q, t=1)[1, 3]
    
    eq1 <- pi.F - pi.F.const
    
    return(eq1)
    
  }
  
  out <- nleqslv(c(0), fn1, control = list(allowSingular=TRUE))
  lam01 <- out$x
  print(out$message)
  
  lam02 <- lam01*(1-p1)/p1
  lam12 <- lam02*ebeta
  
  Q3 <- rbind(c(-(lam01+lam01*(1-p1)/p1), lam01, lam01*(1-p1)/p1),
              c(0, -(lam01*(1-p1)/p1*ebeta), lam01*(1-p1)/p1*ebeta),
              c(0, 0, 0)) 
  
  dat <- NULL
  for (i in 1:1000) {
    sim.i <- sim.msm(Q3, maxtime=10)
    dat.i <- data.frame(id=i, states=sim.i$states, time=sim.i$times)
    dat <- rbind(dat, dat.i)
  }
  dat.fail <- dat[dat$states == 3 ,] 
  tau <- quantile(dat.fail$time, 0.9)          
  
  
  fn2 <- function(x) {
    Q <- rbind(c(-(lam01+lam01*(1-p1)/p1+x), lam01, lam01*(1-p1)/p1, x, 0, 0),
               c(0, -(lam01*(1-p1)/p1*ebeta+x*egamma), lam01*(1-p1)/p1*ebeta, 0, x*egamma, 0),
               c(0, 0, 0, 0, 0, 0),
               c(0, 0, 0, -(lam01+lam01*(1-p1)/p1), lam01, lam01*(1-p1)/p1),
               c(0, 0, 0, 0, -lam01*(1-p1)/p1*ebeta, lam01*(1-p1)/p1*ebeta),
               c(0, 0, 0, 0, 0, 0)) 
    
    
    P <- MatrixExp(Q, t=tau)
    pi.S <- P[1, 4]+P[1, 5]+P[1, 6]
    
    eq2 <- pi.S - pi.S.const
    
    return(eq2)
    
  }
  
  out <- nleqslv(c(0), fn2, control = list(allowSingular=TRUE))
  rho0 <- out$x
  print(out$message)
  
  rho1 <- rho0*egamma
  
  Q6 <- rbind(c(-(lam01+lam01*(1-p1)/p1+rho0), lam01, lam01*(1-p1)/p1, rho0, 0, 0),
              c(0, -(lam01*(1-p1)/p1*ebeta+rho0*egamma), lam01*(1-p1)/p1*ebeta, 0, rho0*egamma, 0),
              c(0, 0, 0, 0, 0, 0),
              c(0, 0, 0, -(lam01+lam01*(1-p1)/p1), lam01, lam01*(1-p1)/p1),
              c(0, 0, 0, 0, -lam01*(1-p1)/p1*ebeta, lam01*(1-p1)/p1*ebeta),
              c(0, 0, 0, 0, 0, 0)) 
  
  outstats <- NULL
  outstats[[1]] <- data.frame(lam01, lam02, lam12, rho0, rho1, tau)
  outstats[[2]] <- Q6
  outstats[[3]] <- Q3
  return(outstats)
}



lam.t <- function(t) {
  P11 <- exp(-(lam01+lam02)*t)
  P12 <- lam01/(lam01+lam02-lam12)*( -exp(-(lam01+lam02)*t)+exp(-lam12*t) )
  return( (P11*lam02+P12*lam12)/(P11+P12) )
  #return( (P11*lam02+P12*lam12) )
}


lam0.t <- function(t) {
  P11 <- exp(-(lam01+lam02+rho0)*t)
  P12 <- lam01/(lam01+lam02+rho0-lam12-rho1)*( -exp(-(lam01+lam02+rho0)*t)+exp(-(lam12+rho1)*t) )
  return( (P11*lam02+P12*lam12)/(P11+P12) )
}


lamP.t <- function(t) {
  P14 <- -exp(-(lam01+lam02+rho0)*t)+exp(-(lam01+lam02)*t)
  a <- lam01+lam02-lam12
  P15 <- ( a*exp(-(lam01+lam02+rho0)*t)+(-a-rho0+rho1)*exp(-(lam01+lam02)*t)-a*exp(-(lam12+rho1)*t)+exp(-lam12*t)*(a+rho0-rho1) )*lam01/(a*(a+rho0-rho1))
  return( (P14*lam02+P15*lam12)/(P14+P15) )
}


lamM.t <- function(t){
  P11 <- exp(-(lam01+lam02+rho0)*t)
  P12 <- lam01/(lam01+lam02+rho0-lam12-rho1)*( -exp(-(lam01+lam02+rho0)*t)+exp(-(lam12+rho1)*t) )
  P14 <- -exp(-(lam01+lam02+rho0)*t)+exp(-(lam01+lam02)*t)
  a <- lam01+lam02-lam12
  P15 <- ( a*exp(-(lam01+lam02+rho0)*t)+(-a-rho0+rho1)*exp(-(lam01+lam02)*t)-a*exp(-(lam12+rho1)*t)+exp(-lam12*t)*(a+rho0-rho1) )*lam01/(a*(a+rho0-rho1))
  #return( (P11*lam02+P12*lam12)+ (P14*lam02+P15*lam12) )
  return( ((P11*lam02+P12*lam12)+ (P14*lam02+P15*lam12))/(P11+P12+P14+P15) )
}


lamM2.t <- expression(
  log( exp(-(lam01+lam02+rho0)*t) + lam01/(lam01+lam02+rho0-lam12-rho1)*( -exp(-(lam01+lam02+rho0)*t)+exp(-(lam12+rho1)*t) ) - exp(-(lam01+lam02+rho0)*t)+exp(-(lam01+lam02)*t) +  ( (lam01+lam02-lam12)*exp(-(lam01+lam02+rho0)*t)+(-(lam01+lam02-lam12)-rho0+rho1)*exp(-(lam01+lam02)*t)-(lam01+lam02-lam12)*exp(-(lam12+rho1)*t)+exp(-lam12*t)*(lam01+lam02-lam12+rho0-rho1) )*lam01/((lam01+lam02-lam12)*((lam01+lam02-lam12)+rho0-rho1)) )
) 




# cumulative hazard -----------------------------------------------------------------------

egamma.opt <- c(1, 4, 8)
piS.opt <- c(0.05, 0.25, 0.5)
t.opt <- seq(0.01, 1.5, 0.01)

p1 <- 0.5

out <- NULL
out0 <- NULL
outP <- NULL
outM <- NULL


set.seed(2024)
for (i in 1:length(piS.opt)) {
  pi.S <- piS.opt[i]
  for (j in 1:length(egamma.opt)) {
    egamma <- egamma.opt[j]
    instats <- cal.f(ebeta=4, p1=p1, egamma=egamma, pi.F.const=0.9, pi.S.const=pi.S)
    instats.par <- instats[[1]]
    
    lam01 <- instats.par$lam01; lam02 <- instats.par$lam02; lam12 <- instats.par$lam12; rho0 <- instats.par$rho0; rho1 <- instats.par$rho1
    tau <- instats.par$tau
    
    instats.Q6 <- instats[[2]]
    
    print(round(instats.par, 2))
    
    outt <- c(pi.S, egamma)
    out00 <- c(pi.S, egamma)
    outPP <- c(pi.S, egamma)
    outMM <- c(pi.S, egamma)
    
    for (r in 1:length(t.opt)) {
      # print(r)
      t <- t.opt[r]
      cumlam.t <- integrate(lam.t, lower = 0, upper = t)$value
      cumlam0.t <- integrate(lam0.t, lower = 0, upper = t)$value
      cumlamP.t <- integrate(lamP.t, lower = 0, upper = t)$value
      cumlamM.t <- integrate(lamM.t, lower = 0, upper = t)$value
      
      outt <- c(outt, cumlam.t)
      out00 <- c(out00, cumlam0.t)
      outPP <- c(outPP, cumlamP.t)
      outMM <- c(outMM, cumlamM.t)
    }
    
    out <- rbind(out, outt)
    out0 <- rbind(out0, out00)
    outP <- rbind(outP, outPP)
    outM <- rbind(outM, outMM)
    
    # # simulation
    # k <- k+1
    # dat <- NULL
    # for (i in 1:100000) {
    #   sim.i <- sim.msm(instats.Q6, maxtime=100)
    #   dat.i <- data.frame(id=i, states=sim.i$states, time=sim.i$times)
    #   dat <- rbind(dat, dat.i)
    # }
    # dat.enrol <- dat %>% group_by(id) %>% filter(max(states) >=4 ) %>% filter(states != 1)
    # 
    # dat.sur <- NULL
    # for (g in unique(dat.enrol$id)) {
    #   #print(g)
    #   dat.enrol.i <- dat.enrol %>% filter(id == g)
    #   dat.sur.i <- data.frame(id=g,
    #                           start=ifelse(min(dat.enrol.i$states)==4, dat.enrol.i$time[dat.enrol.i$states==4], dat.enrol.i$time[dat.enrol.i$states==5]),
    #                           stop=max(dat.enrol.i$time),
    #                           status=ifelse(max(dat.enrol.i$states)==6, 1, 0))
    #   dat.sur <- rbind(dat.sur, dat.sur.i)
    # }
    # 
    # fit <- survfit(Surv(start, stop, status) ~ 1, stype = 2, ctype = 1, id=id, data = dat.sur)
    # sim.t[k, 1: length(fit$time)] <- fit$time
    # sim.cumhaz[k, 1: length(fit$cumhaz)] <- fit$cumhaz
    
    # length(fit$time)
    # t <- 0.8
    # P <- MatrixExp(instats.Q6, t)
    # test <- dat %>% filter(time <= t) %>% group_by(id) %>% mutate(state.cur=max(states)) %>% slice(1)
    # nrow(test %>% filter(state.cur == 1))/10000
    # nrow(test %>% filter(state.cur == 2))/10000
    # nrow(test %>% filter(state.cur == 3))/10000
    # nrow(test %>% filter(state.cur == 4))/10000
    # nrow(test %>% filter(state.cur == 5))/10000
    # nrow(test %>% filter(state.cur == 6))/10000
  }
}

survival.prob <- exp(-out[, 3:ncol(out)])
survival.prob0 <- exp(-out0[, 3:ncol(out0)])
survival.probP <- exp(-outP[, 3:ncol(outP)])
survival.probM <- exp(-outM[, 3:ncol(outM)])



pdf("plot_limit_p1_0.5.pdf", 10, 6)
par(mfrow=c(1,2))

plot(1, type="n", xlab="TIME", xlim=c(0, 1.5), ylim=c(0, 5.5), ylab="Limiting Values of Cumulative Hazards", main = expression(pi[S]~"=0.05"))
lines(t.opt, out[1, 3:ncol(out0)], lty=1, lwd=2)

lines(t.opt, out0[1, 3:ncol(out0)], lty=1, col="blue")
lines(t.opt, out0[2, 3:ncol(out0)], lty=2, col="blue")
lines(t.opt, out0[3, 3:ncol(out0)], lty=3, col="blue")

lines(t.opt, outP[1, 3:ncol(outP)], lty=1, col="red")
lines(t.opt, outP[2, 3:ncol(outP)], lty=2, col="red")
lines(t.opt, outP[3, 3:ncol(outP)], lty=3, col="red")
legend("bottomright", bty='n',
       legend=c(expression(exp(gamma)~"=1"), expression(exp(gamma)~"=4"), expression(exp(gamma)~"=8")), lty=c(1:3), cex=1, lwd = 1.5)

plot(1, type="n", xlab="TIME", xlim=c(0, 1.5), ylim=c(0, 5.5), ylab="Limiting Values of Cumulative Hazards", main = expression(pi[S]~"=0.50"))
lines(t.opt, out[7, 3:ncol(out0)], lty=1, lwd=2)

lines(t.opt, out0[7, 3:ncol(out0)], lty=1, col="blue")
lines(t.opt, out0[8, 3:ncol(out0)], lty=2, col="blue")
lines(t.opt, out0[9, 3:ncol(out0)], lty=3, col="blue")

lines(t.opt, outP[7, 3:ncol(outP)], lty=1, col="red")
lines(t.opt, outP[8, 3:ncol(outP)], lty=2, col="red")
lines(t.opt, outP[9, 3:ncol(outP)], lty=3, col="red")
legend("bottomright", bty='n',
       legend=c(expression(exp(gamma)~"=1"), expression(exp(gamma)~"=4"), expression(exp(gamma)~"=8")), lty=c(1:3), cex=1, lwd = 1.5)

dev.off()



pdf("plot_survival_p1_0.5.pdf", 10, 6)
par(mfrow=c(1,2))
plot(1, type="n", xlab="TIME", xlim=c(0, 1.5), ylim=c(0, 1), ylab="Limiting Values of Survivor Functions", main = expression(pi[S]~"=0.05"))
lines(t.opt, survival.prob[1, ], lty=1)
lines(t.opt, survival.probP[2, ], lty=2)
lines(t.opt, survival.probP[3, ], lty=3)
legend("topright", bty='n',
       legend=c(expression(F(t)~"="~F~"*"~(t)~","~exp(gamma)~"=1"), expression(F~"*"~(t)~","~exp(gamma)~"=4"), expression(F~"*"~(t)~","~exp(gamma)~"=8") ), lty=c(1:3), cex=1, lwd = 1.5)

plot(1, type="n", xlab="TIME", xlim=c(0, 1.5), ylim=c(0, 1), ylab="Limiting Values of Survivor Functions", main = expression(pi[S]~"=0.50"))
lines(t.opt, survival.prob[7, ], lty=1)
lines(t.opt, survival.probP[8, ], lty=2)
lines(t.opt, survival.probP[9, ], lty=3)
legend("topright", bty='n',
       legend=c(expression(F(t)~"="~F~"*"~(t)~","~exp(gamma)~"=1"), expression(F~"*"~(t)~","~exp(gamma)~"=4"), expression(F~"*"~(t)~","~exp(gamma)~"=8") ), lty=c(1:3), cex=1, lwd = 1.5)
dev.off()
