# ============================ CODE METADATA ================================= #
# AUTHOR: Fernando Colchero
# DATE CREATED: 2024-05-29
# DATE MODIFIED: 
# DESCRIPTION: Functions for competing risks MCMC:
# COMMENTS:
# ============================ START CODE ==================================== #
# ========================= #
# ==== USER FUNCTIONS: ====
# ========================= #
# Main competing risks function:
BayesCR <- function(object, ...) UseMethod("BayesCR")

BayesCR.default <- function(object, niter = 11000, burnin = 1001, 
                                 thinning = 10, ncpus = 4, nsim = 4, 
                                 silent = TRUE) {
  # Start time:
  StartTime <- Sys.time()
  
  # Number of records:
  ndat <- nrow(object)
  
  # Vector of unique causes:
  unCause <- sort(unique(object$causes))
  nCause <- length(unCause)
  
  # Theta parameter names:
  thetaNames <- c("b0", "b1")
  nTheta <- length(thetaNames)
  
  # Full parameter names:
  fullThetaNames <- paste(rep(thetaNames, nCause), 
                          rep(unCause, each = nTheta), sep = "_")
  nFullTheta <- length(fullThetaNames)
  
  # Data object:
  dataObj <- list(x = object$ageLast, xt = object$ageFirst)
  
  # Design matrix for cause of death (theta):
  Xc <- model.matrix(~ causes - 1, data = object)
  colnames(Xc) <- unCause
  dataObj$Xc <- Xc
  
  # Design matrix for contraception and cause of death (gamma):
  Xtemp <- model.matrix(~ contra:causes - 1, data = object)
  Xcon <- Xtemp[, grep("contrayes", colnames(Xtemp))]
  colnames(Xcon) <- gsub("contrayes[[:punct:]]{1}causes", "", colnames(Xcon))
  Xcon <- Xcon[, unCause]
  dataObj$Xcon <- Xcon
  
  # Table of causes vs contra:
  causeCont <- table(object$contra, object$causes)
  
  # Find which causes are present in contracepted:
  idConCause <- which(causeCont["yes", ] > 0)
  dataObj$idcc <- idConCause
  
  # Indicator for contracepted individuals:
  indcont <- rep(0, ndat)
  indcont[which(object$contra == "yes")] <- 1
  dataObj$indcont <- indcont
  
  # Indicator for missing causes:
  IndCause <- matrix(1, ndat, nCause)
  IndCause[which(indcont == 1), -idConCause] <- 0
  dataObj$IC <- IndCause
  
  # -------------------- #
  # ---- MCMC PREP: ---- #
  # -------------------- #
  # Theta object:
  thetaDf <- data.frame(Start = rep(c(-5, 0.05), nCause),
                        Low = rep(c(-Inf, 0), nCause),
                        Mean = rep(c(-5, 0.1), nCause),
                        SD = rep(c(5, 1), nCause),
                        Jitter = rep(c(1, 0.2), nCause))
  rownames(thetaDf) <- fullThetaNames
  
  # Gamma object:
  gammaDf <- data.frame(Start = rep(0, nCause),
                        Mean = rep(0, nCause), SD = rep(2, nCause),
                        Jitter = rep(0.5, nCause))
  rownames(gammaDf) <- unCause
  
  # Parameter object:
  paramObj <- list(theta = thetaDf, gamma = gammaDf,
                   nTheta = nTheta, nCause = nCause,
                   thetaName = thetaNames, causeName = unCause,
                   nFull = nFullTheta, fullName = fullThetaNames)
  
  # Vector of indices for parameters:
  keep <- seq(burnin, niter, thinning)
  
  # Jump SD computation:
  if (!silent) {
    cat("\nRunning sequence to find jump SDs... ")
  }
  outJump <- .RunMCMC(sim = 1, dataObj = dataObj, paramObj = paramObj, 
                     .CalcLikePost = .CalcLikePost, .CalcMHratio = .CalcMHratio,
                     jumpUpdate = TRUE, .dtnorm = .dtnorm, .rtnorm = .rtnorm, 
                     .CalcSx = .CalcSx, .Calcmux = .Calcmux,
                     .Calcmux.matrix = .Calcmux.matrix, 
                     .CalcUx = .CalcUx, .PrepJumpObj = .PrepJumpObj, 
                     .CreateUpdObj = .CreateUpdObj, 
                     .UpdateJumps = .UpdateJumps, thetaJump = NA, gammaJump = NA,
                     niter = 10000)
  
  if (!silent) {
    cat(" Done\n")
  }
  
  # Initialize snowfall:
  if (!silent) {
    cat("Multiple MCMC simulations started...\n\n") 
  }
  
  sfInit(parallel = TRUE, cpus = nsim)
  
  # run parallel MCMC:
  outParal <- sfClusterApplyLB(1:nsim, fun = .RunMCMC, dataObj = dataObj, 
                               paramObj = paramObj, .CalcLikePost = .CalcLikePost,
                               .CalcMHratio = .CalcMHratio,
                               jumpUpdate = FALSE, .dtnorm = .dtnorm, 
                               .rtnorm = .rtnorm, .CalcSx = .CalcSx, 
                               .Calcmux = .Calcmux, 
                               .Calcmux.matrix = .Calcmux.matrix, 
                               .CalcUx = .CalcUx,
                               .PrepJumpObj = NA, .CreateUpdObj = NA, 
                               .UpdateJumps = NA, 
                               thetaJump = outJump$thetaJump, 
                               gammaJump = outJump$gammaJump, niter = niter)
  
  # stop snowfall:
  sfStop()
  
  if (!silent) {
    cat("MCMC simulations finished.\n")
  }
  
  # Extract results from parallel runs:
  outTheta <- outParal[[1]]$theta[keep, ]
  outGamma <- outParal[[1]]$gamma[keep, ]
  for (isim in 2:nsim) {
    outTheta <- rbind(outTheta, outParal[[isim]]$theta[keep, ])
    outGamma <- rbind(outGamma, outParal[[1]]$gamma[keep, ])
  }
  
  # Calculate convergence statistics:
  RhatTheta <- .CalcPSRF(object = outParal, keep = keep, nsim = nsim)
  RhatGamma <- .CalcPSRF(object = outParal, keep = keep, nsim = nsim,  
                        partype = "gamma")
  
  # Serial autocorrelation:
  SerCorTheta <- apply(outTheta, 2, function(xth) {
    cor(xth[-1], xth[-length(xth)])
  })
  SerCorGamma <- rep(NA, paramObj$nCause)
  SerCorGamma[dataObj$idcc] <- apply(outGamma[, dataObj$idcc], 2, 
                                     function(xth) {
                                       cor(xth[-1], xth[-length(xth)])
                                     })
  
  # Calculate coefficients:
  coefs <- list(theta = cbind(Mean = apply(outTheta, 2, mean),
                              SE = apply(outTheta, 2, sd),
                              Lower = apply(outTheta, 2, quantile, 0.025),
                              Upper = apply(outTheta, 2, quantile, 0.975),
                              SerCor = SerCorTheta,
                              Rhat = RhatTheta[, "Rhat"]),
                gamma = cbind(Mean = apply(outGamma, 2, mean),
                              SE = apply(outGamma, 2, sd),
                              Lower = apply(outGamma, 2, quantile, 0.025),
                              Upper = apply(outGamma, 2, quantile, 0.975),
                              SerCor = SerCorGamma,
                              Rhat = RhatGamma[, "Rhat"]))
  # Zero coverage for Gamma:
  zerCov <- rep(NA, paramObj$nCause)
  zerCov[dataObj$idcc] <- sapply(dataObj$idcc, function(icc) {
    cof <- coefs$gamma[icc, ]
    zc <- pnorm(q = 0, mean = cof["Mean"], sd = cof["SE"])
    if (zc > 0.5) {
      pz <- 2 * (1 - zc)
    } else {
      pz <- 2 * zc
    }
    return(pz)
  })
  coefs$gamma <- cbind(coefs$gamma, zeroCov = zerCov)
  
  # Kullback-Leibler discrepancies:
  CrossKL <- .CalcCrossKL(coefs = coefs, paramObj = paramObj, 
                         dataObj = dataObj, .CalcKL = .CalcKL)
  
  # Calculate cumulative incidence:  
  dx <- 0.1
  xv <- seq(0, ceiling(max(dataObj$x)), dx)
  nxv <- length(xv)
  Xcv <- matrix(1, nxv, paramObj$nCause)
  
  FxArr <- array(0, dim = c(nxv, paramObj$nCause, 4), 
                 dimnames = list(NULL, paramObj$causeName, 
                                 c("Mean", "SE", "Lower", "Upper")))
  
  Fx <- list(noContr = FxArr, contr = FxArr)
  for (ic in 1:2) {
    if (ic == 1) {
      gamFact <- 0
      idcc <- 1:paramObj$nCause
    } else {
      gamFact <- 1
      idcc <- dataObj$idcc
    }
    nidc <- length(idcc)
    SxMat <- sapply(1:nrow(outTheta), function(iith) {
      thetaMat <- matrix(outTheta[iith, ], paramObj$nCause, paramObj$nTheta,
                         byrow = TRUE, dimnames = list(NULL, paramObj$thetaName))
      Uxmat <- sapply(idcc, function(ith) {
        .CalcUx(theta = thetaMat[ith, ], x = xv, 
               outGamma[iith, ith] * gamFact)
      })
      cumUx <- c(Uxmat %*% rep(1, nidc))
      Sxc <- exp(-cumUx)
      return(Sxc)
    })
    
    for (icc in idcc) {
      FxMat <- sapply(1:nrow(outTheta), function(iith) {
        thetaMat <- matrix(outTheta[iith, ], paramObj$nCause, 
                           paramObj$nTheta, byrow = TRUE, 
                           dimnames = list(NULL, paramObj$thetaName))
        
        muxi <- .Calcmux(theta = thetaMat[icc, ],
                        gamma = outGamma[iith, icc] * gamFact, x = xv)
        Fxi <- cumsum(muxi * SxMat[, iith] * dx)
      })
      Fx[[ic]][, icc, "Mean"] <- apply(FxMat, 1, mean)
      Fx[[ic]][, icc, "SE"] <- apply(FxMat, 1, sd)
      Fx[[ic]][, icc, "Lower"] <- apply(FxMat, 1, quantile, 0.025)
      Fx[[ic]][, icc, "Upper"] <- apply(FxMat, 1, quantile, 0.975)
    }
  }
  
  # MCMC settings:
  settings <- c(niter = niter, burnin = burnin, thinning = thinning, 
                nsim = nsim)
  
  # End time:
  EndTime <- Sys.time()
  CompTime <- round(as.numeric(EndTime - StartTime, 
                               units = units(EndTime - StartTime)), 2)
  if (!silent) {
    cat(sprintf("Total computing time: %.2f %s.\n\n", CompTime, 
                units(EndTime - StartTime)))
  }
  
  # Output object:
  outCompRisks <- list(coefficients = coefs, theta = outTheta, 
                       gamma = outGamma, KL = CrossKL, x = xv,
                       Fx = Fx, data = dataObj, params = paramObj, 
                       settings = settings, runs = outParal,
                       CompTime = CompTime)
  class(outCompRisks) <- "BayesCR"
  return(outCompRisks)
}

# Plot outputs:
plot.BayesCR <- function(x, type = "traces", param = "theta", ...) {
  causeCols <- brewer.pal(9, "Set1")[-6][1:x$params$nCause]
  names(causeCols) <- x$params$causeName
  # Traces:
  if (type == "traces") {
    nsim <- x$settings["nsim"]
    niter <- x$settings["niter"]
    colTrace <- brewer.pal(nsim, "Set1")
    if (param == "theta") {
      namesPars <- x$params$fullName
      nPars <- x$params$nFull
      idPars <- 1:nPars
    } else {
      namesPars <- x$params$causeName
      nPars <- x$params$nCause
      idPars <- x$data$idcc
    }
    runs <- x$runs
    ylim <- apply(runs[[1]][[param]], 2, range)
    for (isim in 2:nsim) {
      ylim <- apply(rbind(ylim, runs[[isim]][[param]]), 2, range)
    }
    
    par(mfrow = c(ceiling(nPars / 2), 2), mar = c(4, 4, 1, 1))
    for (ip in idPars) {
      
      plot(c(1, niter), ylim[, ip], col = NA, main = namesPars[ip],
           xlab = "Iteration", ylab = "Parameter")
      for (isim in 1:nsim) {
        lines(1:niter, runs[[isim]][[param]][, ip], col = colTrace[isim])
      }
    }
    
  }
  
  # Posterior densities:
  if (type == "densities") {
    
    # Posterior Densities for theta:
    thetaDens <- list()
    limTheta <- list(x = cbind(b0 = c(NA, NA), b1 = c(NA, NA)), 
                     y = cbind(b0 = c(0, 0), b1 = c(0, 0)))
    
    for (ith in 1:x$params$nFull) {
      thn <- x$params$fullName[ith]
      thVec <- x$theta[, ith]
      idin <- which(thVec >= x$coefficients$theta[ith, "Lower"] & 
                      thVec <= x$coefficients$theta[ith, "Upper"])
      xx <- thVec[idin]
      thetaDens[[thn]] <- density(xx)
      iGenTh <- substr(thn, 1, 2)
      limTheta$y[2, iGenTh] <- max(c(limTheta$y[2, iGenTh],
                                     thetaDens[[thn]]$y))
      limTheta$x[, iGenTh] <- range(c(xx, limTheta$x[, iGenTh]), 
                                    na.rm = TRUE)
      
    }
    
    # Posterior Densities for gamma:
    gammaDens <- list()
    limGamma <- list(x = c(NA, NA), y = c(0, 0))
    
    
    for (iga in x$data$idcc) {
      gan <- x$params$causeName[iga]
      gaVec <- x$gamma[, iga]
      idin <- which(gaVec >= x$coefficients$gamma[iga, "Lower"] & 
                      gaVec <= x$coefficients$gamma[iga, "Upper"])
      xx <- gaVec[idin]
      gammaDens[[gan]] <- density(xx)
      limGamma$y[2] <- max(c(limGamma$y[2], gammaDens[[gan]]$y))
      limGamma$x <- 
        range(c(xx, limGamma$x), na.rm = TRUE)
      
    }
    
    par(mfrow = c(x$params$nTheta + 1, 1), mar = c(4, 4, 1, 1))
    for (jth in 1:x$params$nTheta) {
      plot(limTheta$x[, jth], limTheta$y[, jth], col = NA, 
           xlab = "", ylab = "")
      for (iith in grep(x$params$thetaName[jth], x$params$fullName)) {
        dens <- thetaDens[[iith]]
        cn <- substr(x$params$fullName[iith], 4, nchar(x$params$fullName[iith]))
        polygon(c(dens$x, rev(dens$x)), c(dens$y, rep(0, length(dens$y))),
                col = adjustcolor(causeCols[cn], alpha.f = 0.25), border = NA)
        lines(dens$x, dens$y, col = causeCols[cn], lwd = 2)
      }
      legend("topleft", legend = x$params$thetaName[jth], cex = 1.5)
    }
    
    plot(limGamma$x, limGamma$y, col = NA, xlab = "Parameter", ylab = "")
    abline(v = 0)
    for (jga in x$data$idcc) {
      cn <- x$params$causeName[jga]
      dens <- gammaDens[[cn]]
      polygon(c(dens$x, rev(dens$x)), c(dens$y, rep(0, length(dens$y))),
              col = adjustcolor(causeCols[cn], alpha.f = 0.25), border = NA)
      lines(dens$x, dens$y, col = causeCols[cn], lwd = 2)
    }
    legend("topleft", "gamma", cex = 1.5)
    legend('topright', legend = x$params$causeName, 
           col = causeCols[x$params$causeName], 
           pt.bg = adjustcolor(causeCols[x$params$causeName], alpha.f = 0.25),
           pch = 21, pt.cex = 2)
  }
  
  # Boxplots of posterior densities:
  if (type == "densbp") {
    par(mfrow = c(x$params$nTheta + 1, 1), mar = c(4, 4, 1, 1))
    for (jth in 1:x$params$nTheta) {
      idth <- grep(x$params$thetaName[jth], x$params$fullName)
      boxplot(x$theta[, idth], col = adjustcolor(causeCols, alpha.f = 0.25), 
              border = causeCols, names = x$params$causeName, notch = TRUE,
              range = 1)
      legend("topleft", x$params$thetaName[jth], cex = 1.5)
    }
    
    boxplot(x$gamma, col = adjustcolor(causeCols, alpha.f = 0.25), 
            border = causeCols, names = x$params$causeName, notch = TRUE, 
            range = 1)
    legend("topleft", "gamma", cex = 1.5)
    abline(h = 0)
  }
  
  if (type == "cumincidence" | type == "cuminc") {
    colcont <- c(noContr = causeCols[1], contr = 'grey40')
    par(mfrow = c(ceiling(x$params$nCause / 2), 2))
    for (icc in 1:x$params$nCause) {
      plot(range(x$x), c(0, 0.35), col = NA, main = x$params$causeName[icc],
           xlab = "Age", ylab = "Cumulative incidence")
      for (ic in 2:1) {
        col <- causeCols[x$params$causeName[icc]]
        if (ic == 2) col <- 'grey40'
        FxArr <- x$Fx[[ic]]
        polygon(c(x$x, rev(x$x)), c(FxArr[, icc, "Lower"], 
                                    rev(FxArr[, icc, "Upper"])),
                col = adjustcolor(colcont[ic], alpha.f = 0.25),
                border = NA)
        lines(x$x, FxArr[, icc, 'Mean'], col = colcont[ic])
      }  
    }
  }
}

# Results functions:
summary.BayesCR <- function(object, ...) {
  args <- list(...)
  if ("digits" %in% names(args)) {
    digits <- args$digits
  } else {
    digits <- 3
  }
  
  # Coefficients:
  cat("------------------------\n")
  cat("Kullback-Leibler discr.:\n")
  cat("------------------------\n")
  for (ipp in 1:length(object$KL)) {
    pname <- names(object$KL)[ipp]
    cat(sprintf("%s:\n", pname))
    print(object$KL[[pname]], digits = digits)
    
  }
  
  # Coefficients:
  cat("\n-------------\n")
  cat("Coefficients:\n")
  cat("-------------\n")
  cat("theta:\n")
  print(object$coefficients$theta, digits = digits)
  cat("\ngamma:\n")
  print(object$coefficients$gamma, digits = digits)
}

print.BayesCR <- function(x, ...) {
  args <- list(...)
  if ("digits" %in% names(args)) {
    digits <- args$digits
  } else {
    digits <- 3
  }
  
  # Coefficients:
  cat("\n-------------\n")
  cat("Coefficients:\n")
  cat("-------------\n")
  cat("theta:\n")
  print(x$coefficients$theta, digits = digits)
  cat("\ngamma:\n")
  print(x$coefficients$gamma, digits = digits)
}

# ============================= #
# ==== INTERNAL FUNCTIONS: ====
# ============================= #
# -------------------------------- #
# ---- Probability functions: ----
# -------------------------------- #
# Truncated normal:
.rtnorm <- function(n, mean, sd, lower = -Inf, upper = Inf) {
  Flow <- pnorm(lower, mean, sd)
  Fup <- pnorm(upper, mean, sd)
  ru <- runif(n, Flow, Fup)
  rx <- qnorm(ru, mean, sd)
  return(rx)
}

.dtnorm <- function(x, mean, sd, lower = -Inf, upper = Inf, log = FALSE) {
  Flow <- pnorm(lower, mean, sd)
  Fup <- pnorm(upper, mean, sd)
  densx <- dnorm(x, mean, sd) / (Fup - Flow)
  if (log) densx <- log(densx)
  return(densx)
}

.ptnorm <- function(q, mean, sd, lower = -Inf, upper = Inf, log = FALSE) {
  p <- (pnorm(q, mean, sd) - pnorm(lower, mean, sd)) / 
    (pnorm(upper, mean, sd) - pnorm(lower, mean, sd))
  if (log) {
    p <- log(p)
  }
  return(p)
}

.qtnorm <- function (p, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  p2 <- (p) * (pnorm(upper, mean, sd) - pnorm(lower, mean, sd)) + 
    pnorm(lower, mean, sd)
  q <- qnorm(p2, mean, sd)
  return(q)
}

# ---------------------------- #
# ---- Survival analysis: ----
# ---------------------------- #
# Cumulative hazards:
.CalcUx <- function(theta, gamma, x) {
  Ui <- exp(theta["b0"]) / theta["b1"] * 
    (exp(theta["b1"] * x) - 1) * exp(gamma)
  return(Ui)
}

# Mortality function:
.Calcmux <- function(theta, ...) UseMethod(".Calcmux")
.Calcmux.matrix <- function(theta, gamma, x) {
  exp(theta[, "b0"] + theta[, "b1"] * x) * exp(gamma)
}
.Calcmux.numeric <- function(theta, gamma, x) {
  exp(theta["b0"] + theta["b1"] * x) * exp(gamma)
}

# Overall Survival:
.CalcSx <- function(thetaMat, gamma, x, dataObj, .CalcUx) {
  nrth <- nrow(thetaMat)
  Uxmat <- sapply(1:nrth, function(ith) {
    .CalcUx(theta = thetaMat[ith, ], x = x, 
           gamma[ith] * dataObj$indcont)
  }) * dataObj$IC
  cumUx <- c(Uxmat %*% rep(1, nrth))
  Sxc <- exp(-cumUx)
  return(Sxc)
}
  
# ------------------------------ #
# ---- Inference functions: ----
# ------------------------------ #
# Likelihood and posterior:
.CalcLikePost <- function(theta, gamma, dataObj, paramObj, .dtnorm,
                         .CalcSx, .Calcmux, .Calcmux.matrix, .CalcUx) {
  # Matrix of parameters per cause:
  thetaMat <- matrix(theta, paramObj$nCause, paramObj$nTheta, 
                     byrow = TRUE, dimnames = list(NULL, paramObj$thetaName))
  thetaCauseMat <- dataObj$Xc %*% thetaMat
  gammaContrVec <- dataObj$Xcon %*% gamma
  Sx <- .CalcSx(thetaMat = thetaMat, gamma = gamma, x = dataObj$x, 
               dataObj = dataObj, .CalcUx = .CalcUx)
  Sxt <- .CalcSx(thetaMat = thetaMat, gamma = gamma, x = dataObj$xt,
                dataObj = dataObj, .CalcUx = .CalcUx)
  mux <- .Calcmux(theta = thetaCauseMat, gamma = gammaContrVec, x = dataObj$x)
  fx <- mux * Sx
  like <- sum(log(fx / Sxt))
  post <- like + sum(.dtnorm(x = theta, mean = paramObj$theta$Mean, 
                            sd = paramObj$theta$SD, 
                            lower = paramObj$theta$Low, log = TRUE)) +
    sum(dnorm(x = gamma, mean = paramObj$gamma$Mean, sd = paramObj$gamma$SD,
              log = TRUE))
  return(list(Like = like, Post = post))
}

# ---------------------------------- #
# ---- MCMC internal functions: ----
# ---------------------------------- #
# Metropolis-Hastings ratio:
.CalcMHratio <- function(thetaNow, thetaNew, paramObj, .dtnorm, pp) {
  mhRatio <- .dtnorm(x = thetaNow[pp], mean = thetaNew[pp], 
                    sd = paramObj$theta$SD[pp], 
                    lower = paramObj$theta$Low[pp], log = TRUE) - 
    .dtnorm(x = thetaNew[pp], mean = thetaNow[pp], 
           sd = paramObj$theta$SD[pp], lower = paramObj$theta$Low[pp], 
           log = TRUE)
  return(mhRatio)
}

# function to create jump object:
.PrepJumpObj <- function(jump, paramObj, partype = "theta") {
  if (partype == "theta") {
    jobj <- list(jump = jump, jumpsMat = matrix(NA, 0, paramObj$nFull), 
                 update = TRUE, updateRate = rep(0, paramObj$nFull))
  } else {
    jobj <- list(jump = jump, jumpsMat = matrix(NA, 0, paramObj$nCause), 
                 update = TRUE, updateRate = rep(0, paramObj$nCause))
  }
  return(jobj)
}

# put this in the MCMC function, before the loop starts:
.CreateUpdObj <- function(len = 50, targ = 0.25, niter, paramObj,
                         partype = "theta") {
  if (partype == "theta") {
    p <- paramObj$nFull
  } else {
    p <- paramObj$nCause
  }
  updObj <- list()
  updObj$len <- len
  updObj$targ <- targ
  updObj$niter <- niter
  updObj$int <- seq(updObj$len, updObj$niter, updObj$len) 
  updObj$updVec <- matrix(0, updObj$niter, p)
  return(updObj)
}

# Function to update jump variances:
.UpdateJumps <- function(jumpObj, updObj, step) {
  jumpObj$updateRate <- 
    apply(matrix(updObj$updVec[step + c(-(updObj$len - 1):0), ], 
                 ncol = length(jumpObj$jump)), 2, function(upd) sum(upd)) / 
    updObj$len
  jumpObj$updateRate[jumpObj$updateRate == 0] <- 1e-2
  jumpObj$jump <- jumpObj$jump * jumpObj$updateRate / updObj$targ
  jumpObj$jumpsMat <- rbind(jumpObj$jumpsMat, jumpObj$jump)
  return(jumpObj)
}

# ----------------------------- #
# ---- Main MCMC function: ----
# ----------------------------- #
# MCMC function:
.RunMCMC <- function(sim = 1, dataObj, paramObj, .CalcLikePost, .CalcMHratio,
                    jumpUpdate = FALSE, .dtnorm, .rtnorm, .CalcSx, .CalcUx,
                    .Calcmux, .Calcmux.matrix, .PrepJumpObj, .CreateUpdObj,
                    .UpdateJumps, thetaJump, gammaJump, niter) {
  # Theta matrix:
  thetaNow <- paramObj$theta$Start
  
  # Gamma vector:
  gammaNow <- paramObj$gamma$Start
  
  # Likelihood and posterior:
  likePostNow <- .CalcLikePost(theta = thetaNow, gamma = gammaNow, 
                              dataObj = dataObj, 
                              paramObj = paramObj, .dtnorm = .dtnorm,
                              .CalcSx = .CalcSx, .Calcmux = .Calcmux, 
                              .Calcmux.matrix = .Calcmux.matrix, 
                              .CalcUx = .CalcUx)
  
  if (sim > 1) {
    JIT <- TRUE
    while(JIT) {
      thetaNow <- .rtnorm(n = paramObj$nFull, mean = paramObj$theta$Start,
                         sd = paramObj$theta$Jitter, 
                         lower = paramObj$theta$Low)
      gammaNow <- paramObj$gamma$Start
      gammaNow[dataObj$idcc] <-
        rnorm(n = length(dataObj$idcc),
              mean = paramObj$gamma$Start[dataObj$idcc],
              sd = paramObj$gamma$Jitter[dataObj$idcc])
      
      # Likelihood and posterior:
      likePostNow <- .CalcLikePost(theta = thetaNow, gamma = gammaNow, 
                                  dataObj = dataObj, 
                                  paramObj = paramObj, .dtnorm = .dtnorm,
                                  .CalcSx = .CalcSx, .Calcmux = .Calcmux, 
                                  .Calcmux.matrix = .Calcmux.matrix, 
                                  .CalcUx = .CalcUx)
      
      if (!is.na(likePostNow$Post)) {
        if (likePostNow$Post != -Inf) {
          JIT <- FALSE
        }
      } 
    }
  } 
  
  # Setup jump update (if needed):
  if (jumpUpdate) {
    updObjTheta <- .CreateUpdObj(len = 50, targ = 0.25, niter = niter, 
                                paramObj = paramObj)
    updObjGamma <- .CreateUpdObj(len = 50, targ = 0.25, niter = niter, 
                                paramObj = paramObj, partype = "gamma")
    
    # Theta jump SD:
    thetaJump <- rep(c(0.1, 0.01), paramObj$nCause)
    names(thetaJump) <- paramObj$fullName
    
    # Gamma jump SD:
    gammaJump <- rep(0.1, paramObj$nCause)
    
    # Prepare jump object:
    jumpObjTheta <- .PrepJumpObj(jump = thetaJump, paramObj = paramObj)
    thetaOut <- NA
    jumpObjGamma <- .PrepJumpObj(jump = gammaJump, paramObj = paramObj, 
                                partype = "gamma")
    thetaOut <- NA
    gammaOut <- NA
  } else {
    # MCMC output matrices:
    thetaOut <- matrix(0, niter, paramObj$nFull, 
                       dimnames = list(NULL, paramObj$fullName))
    gammaOut <- matrix(0, niter, paramObj$nCause, 
                       dimnames = list(NULL, paramObj$causeName))
    thetaOut[1, ] <- thetaNow
    gammaOut[1, ] <- gammaNow
    jumpObjTheta <- NA
    jumpObjGamma <- NA
  }
  
  
  # Run MCMC
  for (iter in 2:niter) {
    for (pp in 1:paramObj$nFull) {
      thetaNew <- thetaNow
      thetaNew[pp] <- .rtnorm(n = 1, mean = thetaNow[pp], sd = thetaJump[pp],
                             lower = paramObj$theta$Low[pp])
      likePostNew <- .CalcLikePost(theta = thetaNew, gamma = gammaNow,
                                  dataObj = dataObj, 
                                  paramObj = paramObj, .dtnorm = .dtnorm,
                                  .CalcSx = .CalcSx, .Calcmux = .Calcmux, 
                                  .Calcmux.matrix = .Calcmux.matrix, 
                                  .CalcUx = .CalcUx)
      mhRatio <- .CalcMHratio(thetaNow = thetaNow, thetaNew = thetaNew,
                             paramObj = paramObj, .dtnorm = .dtnorm, pp = pp)
      postRatio <- exp(likePostNew$Post - likePostNow$Post + mhRatio)
      if (!is.na(postRatio)) {
        if (postRatio >= runif(1)) {
          thetaNow <- thetaNew
          likePostNow <- likePostNew
          if (jumpUpdate) {
            updObjTheta$updVec[iter, pp] <- 1
          }
        }
      } 
    }
    
    # Gamma update:
    for (pp in dataObj$idcc) {
      gammaNew <- gammaNow
      gammaNew[pp] <- rnorm(n = 1, mean = gammaNow[pp], sd = gammaJump[pp])
      likePostNew <- .CalcLikePost(theta = thetaNow, gamma = gammaNew,
                                  dataObj = dataObj,
                                  paramObj = paramObj, .dtnorm = .dtnorm,
                                  .CalcSx = .CalcSx, .Calcmux = .Calcmux,
                                  .Calcmux.matrix = .Calcmux.matrix,
                                  .CalcUx = .CalcUx)
      postRatio <- exp(likePostNew$Post - likePostNow$Post)
      if (!is.na(postRatio)) {
        if (postRatio >= runif(1)) {
          gammaNow <- gammaNew
          likePostNow <- likePostNew
          if (jumpUpdate) {
            updObjGamma$updVec[iter, pp] <- 1
          }
        }
      }
    }
    
    # Refresh jumps:
    if (jumpUpdate) {
      if (iter %in% updObjTheta$int) {
        jumpObjTheta <- .UpdateJumps(jumpObj = jumpObjTheta, 
                                    updObj = updObjTheta, step = iter)
        thetaJump <- jumpObjTheta$jump
        jumpObjGamma <- .UpdateJumps(jumpObj = jumpObjGamma, 
                                    updObj = updObjGamma, step = iter)
        gammaJump <- jumpObjGamma$jump
      } else if (iter == max(updObjTheta$int)) {
        ljupd <- length(updObjTheta$int)
        idbjmean <- floor(ljupd/2):ljupd
        thetaJump <- apply(jumpObjTheta$jumpsMat[idbjmean, ], 2, mean)
        gammaJump <- apply(jumpObjGamma$jumpsMat[idbjmean, ], 2, mean)
      }
      
    }
    
    # Fill up parameter matrix:
    if (!jumpUpdate) {
      thetaOut[iter, ] <- thetaNow
      gammaOut[iter, ] <- gammaNow
    }
  }
  # Final output object:
  outObj <- list(theta = thetaOut, gamma = gammaOut, thetaJump = thetaJump,
                 gammaJump = gammaJump)
  
  return(outObj)
}

# -------------------------------- #
# ---- MCMC output functions: ----
# -------------------------------- #
# Function to calculate convergence statistics 
#      based on Gelman et al. (2014).
.CalcPSRF <- function(object, keep, nsim, partype = "theta") {
  nthin <- length(keep)
  Means <- t(sapply(1:nsim, function(i) {
    apply(object[[i]][[partype]][keep, ], 2, mean)
  }))
  Vars <- t(sapply(1:nsim, function(i) {
    apply(object[[i]][[partype]][keep, ], 2, var)
  }))
  meanall <- apply(Means, 2, mean)
  B <- nthin / (nsim - 1) * apply(t((t(Means) - meanall)^2), 2, sum)
  W <- 1 / nsim * apply(Vars, 2, sum)
  Varpl <- (nthin - 1) / nthin * W + 1 / nthin * B
  Rhat <- sqrt(Varpl / W)
  Rhat[which(Varpl == 0 | Rhat < 1)] <- 1
  conv <- cbind(B, W, Varpl, Rhat)
  rownames(conv) <- colnames(Means)
  return(conv)
}

# Kullback-Leibler discrepancies:
.CalcKL <- function(m1, sd1, m2, sd2, low = -Inf) {
  qf1 <- .qtnorm(c(0.0001, 0.9999), mean = m1, sd = sd1, lower = low)
  qf2 <- .qtnorm(c(0.0001, 0.9999), mean = m2, sd = sd2, lower = low)
  xp <- seq(min(c(qf1[1], qf2[1])), max(c(qf1[2], qf2[2])), length = 1000)
  dxp <- xp[2] - xp[1]
  df1 <- .dtnorm(xp, mean = m1, sd = sd1, lower = low)
  df2 <- .dtnorm(xp, mean = m2, sd = sd2, lower = low)
  idn0 <- which(df1 > 0 & df2 > 0)
  kl12 <- sum(df1 * log(df1 / df2) * dxp, na.rm = TRUE)
  kl21 <- sum(df2 * log(df2 / df1) * dxp, na.rm = TRUE)
  qKl12 <- (1 + (1 - exp(-2 * kl12)^(1 / 2))) / 2
  qKl21 <- (1 + (1 - exp(-2 * kl21)^(1 / 2))) / 2
  mqKl <- (qKl21 + qKl12) / 2
  outKL <- c(kl12 = kl12, kl21 = kl21, qkl12 = qKl12, 
             qkl21 = qKl21, mqKl = mqKl)
  return(outKL)
}

# Cross referenced KL:
.CalcCrossKL <- function(coefs, paramObj, dataObj, .CalcKL) {
  KLs <- list()
  for (ip in 1:(paramObj$nTheta + 1)) {
    if (ip <= paramObj$nTheta) {
      pname <- paramObj$thetaName[ip]
      idp <- grep(pname, paramObj$fullName)
      mvec <- coefs$theta[idp, "Mean"]
      sdvec <- coefs$theta[idp, "SE"]
      low <- paramObj$theta$Low[ip]
      idpall <- 1:paramObj$nCause
    } else {
      pname <- "gamma"
      mvec <- coefs$gamma[, "Mean"]
      sdvec <- coefs$gamma[, "SE"]
      low <- -Inf
      idpall <- dataObj$idcc
    }
    np <- length(idpall)
    incCauseName <- paramObj$causeName[idpall]
    CrossKL <- matrix(NA, np - 1, np - 1,
                      dimnames = list(incCauseName[-np],
                                      incCauseName[-1]))
    for (ii in 1:(np - 1)) {
      for (jj in 2:np) {
        if (jj > ii) {
          KL <- .CalcKL(m1 = mvec[idpall[ii]], sd1 = sdvec[idpall[ii]], 
                       m2 = mvec[idpall[jj]], sd2 = sdvec[idpall[jj]], 
                       low = low)
          CrossKL[ii, jj - 1] <- KL["mqKl"]
        }
      }
    }
    KLs[[pname]] <- CrossKL
  }
  return(KLs)  
}


# ============================== CODE END ==================================== #

