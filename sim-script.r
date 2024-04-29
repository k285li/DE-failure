# ==============================================
# Script to perform simulations
#
# Last Update: August 2, 2022
# ==============================================


source('sim-func.r')
source('sim-methods.r')


script.f <- function(pi1, pi1o) {
  file.out <- paste("sim-pi1", pi1, "-pi1o", pi1o, "-nocov.dat", sep="")
  if ( file.exists(as.character(file.out)) ) { unlink(as.character(file.out)) }

  stats <- getstats.f(pi1=pi1, pi1o=pi1o)

  for (nsim in 1:200) {
    cat("nsim = ", nsim, "\n")

    set.seed(999 + nsim)
    simdata <- generatedata.manual.f(npop=80000, nsamples=5000,
      lam=stats$lam, lamp=stats$lam, rho=stats$rho,
      A=stats$A, tau50=stats$tau50, tau95=stats$tau95)

    fit1 <- nlm(f=logL.obs.f,
                p=round(c(log(stats$lam), log(stats$rho)), 1),
                hessian=TRUE, indata=simdata$recruited)
    fit1$cov <- ginv(fit1$hessian)
    print(fit1)

    fit2 <- nlm(f=logL.comp.f,
                p=round(c(log(stats$lam), log(stats$rho)), 1),
                hessian=TRUE, indata=simdata$recruited, inrej=simdata$rejected)
    fit2$cov <- ginv(fit2$hessian)
    print(fit2)

    fit3 <- nlm(f=logL.retro.f,
                p=round(c(log(stats$lam), log(stats$rho)), 1),
                hessian=TRUE, indata=simdata$recruited, inrej=simdata$rejected)
    fit3$cov <- ginv(fit3$hessian)
    print(fit3)

    fit4 <- nlm(f=logL.bin.f,
                p=round(c(log(stats$lam), log(stats$rho)), 1),
                hessian=TRUE, indata=simdata$recruited, inrej=simdata$rejected)
    fit4$cov <- ginv(fit4$hessian)
    print(fit4)

    cat(nsim, fit1$code,  max(fit1$gradient),  fit1$estimate,  diag(fit1$cov),
              fit2$code,  max(fit2$gradient),  fit2$estimate,  diag(fit2$cov),
              fit3$code,  max(fit3$gradient),  fit3$estimate,  diag(fit3$cov),
              fit4$code,  max(fit4$gradient),  fit4$estimate,  diag(fit4$cov),
        "\n", sep=" ", append=TRUE, file=file.out)
  }

  return()
}


script.full.f <- function(pi1, pi1o) {
  file.out <- paste("sim-pi1", pi1, "-pi1o", pi1o, "-nocov-full.dat", sep="")
  if ( file.exists(as.character(file.out)) ) { unlink(as.character(file.out)) }

  stats <- getstats.f(pi1=pi1, pi1o=pi1o)

  for (nsim in 1:200) {
    cat("nsim = ", nsim, "\n")

    set.seed(999 + nsim)
    simdata <- generatedata.manual.f(npop=80000, nsamples=5000,
      lam=stats$lam, lamp=stats$lam, rho=stats$rho,
      A=stats$A, tau50=stats$tau50, tau95=stats$tau95)

    fit1 <- nlm(f=logL.comp.full.f,
                p=round(c(log(stats$lam), log(stats$rho), log(stats$lam)), 1),
                hessian=TRUE, indata=simdata$recruited, inrej=simdata$rejected)
    fit1$cov <- ginv(fit1$hessian)
    print(fit1)


    cat(nsim, fit1$code, max(fit1$gradient), fit1$estimate, diag(fit1$cov),
        "\n", sep=" ", append=TRUE, file=file.out)
  }

  return()
}

#script.f(pi1=0.9, pi1o=0.05)
#script.f(pi1=0.9, pi1o=0.10)
#script.f(pi1=0.9, pi1o=0.25)

script.full.f(pi1=0.9, pi1o=0.10)
#script.full.f(pi1=0.9, pi1o=0.25)
