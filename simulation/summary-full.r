# ==============================================
# Summarize estimates
#
# Last Update: August 2, 2022
# ==============================================


source('../shared.func.r')
source('sim-func.r')


summary.f <- function(pi1, file.out) {
  stats1 <- getstats.f(pi1=pi1, pi1o=0.10)
  stats2 <- getstats.f(pi1=pi1, pi1o=0.25)

  col.tx <- c("nsim","fit1.code","fit1.maxgrad",
                     "fit1.loglam","fit1.logrho",
                     "fit1.var.loglam","fit1.var.logrho",
                     "fit2.code","fit2.maxgrad",
                     "fit2.loglam","fit2.logrho",
                     "fit2.var.loglam","fit2.var.logrho",
                     "fit2d.code","fit2d.maxgrad",
                     "fit2d.loglam","fit2d.logrho","fit2d.loglamp",
                     "fit2d.var.loglam","fit2d.var.logrho","fit2d.var.loglamp",
                     "fit3.code","fit3.maxgrad",
                     "fit3.loglam","fit3.logrho",
                     "fit3.var.loglam","fit3.var.logrho",
                     "fit4.code","fit4.maxgrad",
                     "fit4.loglam","fit4.logrho",
                     "fit4.var.loglam","fit4.var.logrho")

  sim1 <- read.table(paste("sim-pi1", pi1, "-pi1o", 0.10, "-nocov.dat", sep=""), 
                     header=FALSE)
  dimnames(sim1)[[2]] <- as.character(col.tx)

  sim2 <- read.table(paste("sim-pi1", pi1, "-pi1o", 0.25, "-nocov.dat", sep=""), 
                     header=FALSE)
  dimnames(sim2)[[2]] <- as.character(col.tx)
  print(summary(sim2$fit2d.loglam))
  print(summary(sim2$fit2d.loglamp))


  L.tx <- c("$L(\\lambda = \\lambda^{\\circ}, \\rho)$",
            "$L^{AUG}(\\lambda = \\lambda^{\\circ}, \\rho)$",
            "$L^{AUG}_{CD}(\\lambda, \\rho)$",
            "$L_{B}(\\lambda, \\rho)$")

  for (k in 1:4) {
    est11 <- getbias.full.f(x=exp(sim1[, paste("fit", k, ".loglam", sep="")]), 
               var.x=( (exp(sim1[, paste("fit", k, ".loglam", sep="")])^2)*sim1[, paste("fit", k, ".var.loglam", sep="")] ), 
               true.x=stats1$lam)
    est12 <- getbias.full.f(x=exp(sim1[, paste("fit", k, ".logrho", sep="")]), 
               var.x=( (exp(sim1[, paste("fit", k, ".logrho", sep="")])^2)*sim1[, paste("fit", k, ".var.logrho", sep="")] ), 
               true.x=stats1$rho)

    est21 <- getbias.full.f(x=exp(sim2[, paste("fit", k, ".loglam", sep="")]), 
               var.x=( (exp(sim2[, paste("fit", k, ".loglam", sep="")])^2)*sim2[, paste("fit", k, ".var.loglam", sep="")] ), 
               true.x=stats2$lam)
    est22 <- getbias.full.f(x=exp(sim2[, paste("fit", k, ".logrho", sep="")]), 
               var.x=( (exp(sim2[, paste("fit", k, ".logrho", sep="")])^2)*sim2[, paste("fit", k, ".var.logrho", sep="")] ), 
               true.x=stats2$rho)
    cat(as.character(L.tx[k]),
        " && ", as.character(est11$bias), " & ", as.character(est11$ese),
        " && ", as.character(est12$bias), " & ", as.character(est12$ese),
        " && ", as.character(est21$bias), " & ", as.character(est21$ese),
        " && ", as.character(est22$bias), " & ", as.character(est22$ese),
        "\\\\", "\n", sep="", append=TRUE, file=file.out)
  }

  est11 <- getbias.full.f(x=exp(sim1[, "fit2d.loglam"]), 
             var.x=( (exp(sim1[, "fit2d.loglam"])^2)*sim1[, "fit2d.var.loglam"] ), 
             true.x=stats1$lam)
  est12 <- getbias.full.f(x=exp(sim1[, "fit2d.logrho"]), 
             var.x=( (exp(sim1[, "fit2d.logrho"])^2)*sim1[, "fit2d.var.logrho"] ), 
             true.x=stats1$rho)

  est21 <- getbias.full.f(x=exp(sim2[, "fit2d.loglam"]), 
             var.x=( (exp(sim2[, "fit2d.loglam"])^2)*sim2[, "fit2d.var.loglam"] ), 
             true.x=stats2$lam)
  est22 <- getbias.full.f(x=exp(sim2[, "fit2d.logrho"]), 
             var.x=( (exp(sim2[, "fit2d.logrho"])^2)*sim2[, "fit2d.var.logrho"] ), 
             true.x=stats2$rho)
  cat("$L^{AUG}(\\lambda, \\rho, \\lambda^{\\circ})$",
      " && ", as.character(est11$bias), " & ", as.character(est11$ese),
      " && ", as.character(est12$bias), " & ", as.character(est12$ese),
      " && ", as.character(est21$bias), " & ", as.character(est21$ese),
      " && ", as.character(est22$bias), " & ", as.character(est22$ese),
      "\\\\", "\n", sep="", append=TRUE, file=file.out)

  est11 <- getbias.full.f(x=exp(sim1[, "fit2d.loglamp"]), 
             var.x=( (exp(sim1[, "fit2d.loglamp"])^2)*sim1[, "fit2d.var.loglamp"] ), 
             true.x=stats1$lam)

  est21 <- getbias.full.f(x=exp(sim2[, "fit2d.loglamp"]), 
             var.x=( (exp(sim2[, "fit2d.loglamp"])^2)*sim2[, "fit2d.var.loglamp"] ), 
             true.x=stats2$lam)

  cat(" && ", as.character(est11$bias), " & ", as.character(est11$ese),
      " && ", " & ", 
      " && ", as.character(est21$bias), " & ", as.character(est21$ese),
      "\\\\", "\n", sep="", append=TRUE, file=file.out)

  cat("\n", append=TRUE, file=file.out)

  return()
}


file.out <- "summary.out"
if ( file.exists(as.character(file.out)) ) { unlink(as.character(file.out)) }

summary.f(pi1=0.9, file.out=file.out)
