# ==============================================
# Scripts for real-data application on likelihood 
# minimization and regression parameters estimation
#
# Last Update: April 29, 2024
# ==============================================



library(HelpersMG)
library(dplyr)
source("clean_ps_to_psa.R")


# minimize likelihood

fit.LR1 <- optim(par=c(-2, 0, 0, 0, -4), fn=logLR1.f, method = "BFGS", hessian = TRUE,
               indata1=indata1)
round(fit.LR1$par, 2)
round(SEfromHessian(fit.LR1$hessian), 5)

fit.LR2 <- optim(par=c(-2, 0, 0, 0, -4, -2), fn=logLR2.f, method = "BFGS", hessian = TRUE,
               indata2=indata2)
round(fit.LR2$par, 2)
round(SEfromHessian(fit.LR2$hessian), 3)

fit.LR <- optim(par=c(-2, 0, 0, 0, -4, -2), fn = logLR.f,  method = "Nelder-Mead", hessian = TRUE, 
                control = list(maxit=10000),
              indata1=indata1, indata2=indata2)
round(fit.LR$par, 2)
round(SEfromHessian(fit.LR$hessian), 3)

fit.LRLP <- optim(par=c(-2, 0, 0, 0, -4, -2), fn=logLRLP.f,  method = "Nelder-Mead", hessian = TRUE,
                  control = list(maxit=10000),
               indata1=indata1, indata2=indata2)
round(fit.LRLP$par, 2)
round(SEfromHessian(fit.LRLP$hessian), 3)

fit.LP <- optim(par=c(-2, 0, 0, 0), fn=logLP.f,  method = "BFGS", hessian = TRUE,
                  indata1=indata1)
round(fit.LP$par, 2)
round(SEfromHessian(fit.LP$hessian), 3)
