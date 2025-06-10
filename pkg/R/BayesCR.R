# ============================ CODE METADATA ================================= #
# AUTHOR: Fernando Colchero
# DATE CREATED: 2025-06-05
# DATE MODIFIED: 
# DESCRIPTION: Functions for competing risks MCMC:
# COMMENTS:
# ============================ START CODE ==================================== #
# ========================= #
# ==== USER FUNCTIONS: ====
# ========================= #
# Main function:
BayesCR <- function(object, ...) UseMethod("BayesCR")
BayesCR.default <- function(object, niter = 11000, burnin = 1001, 
                            thinning = 20, nsim = 4, ncpus = 4, 
                            UPDJUMP = TRUE, jumpSD = NULL) {
  # start timer:
  StartFull <- Sys.time()
  
  # Algorithm object:
  algObj <- .CreateAlgObj(niter, burnin, thinning, nsim, UPDJUMP, jumpSD)
  
  # data object:
  dataObj <- .CreateDataObj(object = object)
  
  # Parameter object:
  parObj <- .CreateParObj(dataObj = dataObj)
  
  # jump sd:
  jumpSD <- parObj$jump
  
  # Initial theta object:
  thetaObj <- .CreateThetaObj(parObj = parObj, dataObj = dataObj)
  
  # run Jump MCMC:
  cat("\nRunning sequence to find jump SDs... ")
  StartJump <- Sys.time()
  
  jumpMCMC <- .RunMCMC(sim = 1, dataObj, parObj, thetaObj, algObj, .CalcMu, 
                       .CalcMu.matrix, .CalcMu.numeric, .CalcU, .CalcU.numeric,
                       .CalcU.matrix, .CalcS, .JitterTheta, .SampleTheta, 
                       .CalcDemoFuns, .CalcLikePost, .CalcMetHastRatio,
                       .UpdateJumps, .rtnorm, .dtnorm, jumpSD = NULL, 
                       UPDJUMP = TRUE)
  jumpSD <- jumpMCMC$jumpSD
  parObj$jump <- jumpSD
  EndJump <- Sys.time()
  cat("Done\n")
  compTime <- round(as.numeric(EndJump-StartJump, 
                               units = units(EndJump - StartJump)), 2)
  cat(sprintf("Total jump SDs computing time: %.2f %s.\n\n", compTime, 
              units(EndJump - StartJump)))
  
  # Main MCMC:
  cat("Multiple simulations started...\n\n") 
  
  # StartFull parallel computing:
  sfInit(parallel = TRUE, cpus = ncpus)
  
  # Load functions to CPUS:
  sfLibrary(BayesCompRisk)
  # sfSource("devel/code/sourceToCPUS.R")
  
  # Load variables to CPUS:
  # sfExport(list = cpuVars)
  
  # Run MCMC in parallel:
  outMCMC <- sfClusterApplyLB(1:nsim, .RunMCMC, dataObj, parObj, thetaObj, 
                              algObj, .CalcMu, .CalcMu.matrix, .CalcMu.numeric, 
                              .CalcU, .CalcU.numeric, .CalcU.matrix, .CalcS, 
                              .JitterTheta, .SampleTheta, 
                              .CalcDemoFuns, .CalcLikePost, .CalcMetHastRatio,
                              .UpdateJumps, .rtnorm, .dtnorm, jumpSD = jumpSD,
                              UPDJUMP = FALSE)
  
  # Stop cluster:
  sfStop()
  
  cat("Simulations finished.\n")
  
  # Index of kept parameters:
  idKeep <- seq(algObj$burnin, algObj$niter, algObj$thinning)
  
  # Matrices of converged parameters and likelihood-posterior to be stored:
  thetaMat <- outMCMC[[1]]$theta[idKeep, ]
  likePostMat <- outMCMC[[1]]$likePost[idKeep, ]
  for (isim in 2:nsim) {
    thetaMat <- rbind(thetaMat, outMCMC[[isim]]$theta[idKeep, ])
    likePostMat <- rbind(likePostMat, outMCMC[[isim]]$likePost[idKeep, ])
  }
  
  # Calculate convergence statistics:
  Conv <- .CalcPSRF(object = outMCMC, keep = idKeep, nsim = nsim)
  
  # Logical for convergence in all parameters:
  if (all(Conv[, "Rhat"] < 1.05)) {
    convergence <- TRUE
  } else {
    convergence <- FALSE
  }
  
  # Effective sample size:
  Neff <- .CalcNeff(object = outMCMC, keep = idKeep, nsim = nsim, Rhat = Conv)
  
  # coefficients:
  coeffs <- cbind(Mean = apply(thetaMat,  2, mean), 
                  SD = apply(thetaMat, 2, sd),
                  Lower = apply(thetaMat, 2, quantile, 0.025, names = FALSE),
                  Upper = apply(thetaMat, 2, quantile, 0.975, names = FALSE),
                  Neff = Neff, Rhat = Conv[, "Rhat"])
  
  
  # Calculate DIC:
  DIC <- .CalcDIC(likelihood = likePostMat[, "Likelihood"], k = parObj$pc)
  
  # Calculate demography quantiles:
  demoQuant <- .CalcDemoQuants(thetaMat = thetaMat, parObj = parObj, 
                               dataObj = dataObj)
  
  # Product limit estimators:
  PLE <- list()
  for (ic in 1:dataObj$nCause) {
    depType <- rep("C", dataObj$n)
    depType[which(dataObj$iCause == dataObj$causes[ic])] <- "D"
    plei <- .CalcPLE(ageLast = dataObj$x, departType = depType)
    PLE[[dataObj$causes[ic]]] <- plei
  }
  depType <- rep("D", dataObj$n)
  PLE$All <- .CalcPLE(ageLast = dataObj$x, departType = depType)
  
  # end timer:
  EndFull <- Sys.time()
  compTime <- .CalcTimeDiff(Start = StartFull, End = EndFull)
  cat(sprintf("Total MCMC computing time: %s\n\n", compTime))
  
  # Settings:
  settings <- list(niter = niter, burnin = burnin, thinning = thinning, 
                   nsim = nsim, ncpus = ncpus, compTime = compTime)
  
  # store output:
  fullOut <- list(coefficients = coeffs, x = demoQuant$x, 
                  mort = demoQuant$mort, surv = demoQuant$surv, 
                  cuts = demoQuant$cuts, theta = thetaMat, 
                  likePost = likePostMat, DIC = DIC, PLE = PLE,
                  runs = outMCMC, data = dataObj, settings = settings, 
                  keep = idKeep, params = parObj, convergence = convergence)
  
  class(fullOut) <- "BayesCR"
  return(fullOut)
}

# Plotting:
plot.BayesCR <- function(x, type = "traces", noCIs = FALSE, logMort = FALSE,
                         ...) {
  # User par settings:
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  # Cause name and number:
  causes <- x$data$causes
  nCause <- x$data$nCause
  
  # Settings for traces and posterior densities:
  if (type %in% c("traces", "densities", "densComp")) {
    namesFull <- x$params$namesFull
    namesSimp <- x$params$names
    np <- x$params$p
    npc <- x$params$pc
  }
  if (type == "traces") {
    nsim <- x$settings$nsim
    if (nsim <= 12) {
      Palette <- c('#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C',
                   '#FB9A99', '#E31A1C', '#FDBF6F', '#FF7F00',
                   '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928')
    } else {
      Palette <- rainbow(nsim)
    }
    xkeep <- seq(1, x$settings$niter, x$settings$thinning)
    xlim <- c(1, x$settings$niter)
    xrun <- x$runs
    par(mfrow = c(np, nCause), mar = c(4, 4, 1, 1))
    for (ip in 1:npc) {
      thName <- namesFull[ip]
      ylim <- range(sapply(1:nsim, function(is) {
        range(xrun[[is]]$theta[, thName])
      }))
      plot(xlim, ylim, col = NA, xlab = "", ylab = "", 
           main = thName)
      for (is in 1:nsim) {
        lines(xkeep, xrun[[is]]$theta[xkeep, thName], 
              col = Palette[is])
      }
    }
  } else if (type == "densities") {
    par(mfrow = c(np, nCause), mar = c(4, 4, 1, 1))
    
    for (ip in 1:npc) {
      thName <- namesFull[ip]
      xthe <- x$theta[, thName]
      sthe <- x$coefficients[thName, 1:4]
      low <- x$coefficients[thName, "Lower"]
      dens <- density(xthe)
      if (low == 0) {
        id0 <- which(dens$x >= 0)
        dens$x <- dens$x[id0]
        dens$y <- dens$y[id0]
      }
      id95 <- which(dens$x >= sthe["Lower"] & dens$x <= sthe["Upper"])
      ylim <- c(0, max(dens$y))
      xlim <- quantile(xthe, probs = c(0.001, 0.999))
      plot(xlim, ylim, col = NA, xlab = "", ylab = "", 
           main = thName)
      xx95 <- dens$x[id95]
      yy95 <- dens$y[id95]
      polygon(x = c(xx95, rev(xx95)), 
              y = c(yy95, rep(0, length(yy95))), col = "orange",
              border = NA)
      lines(dens$x, dens$y, col = "orange")
    }
  } else if (type == "densComp") {
    if (nCause <= 8) {
      Palette <- c('#1B9E77', '#D95F02', '#7570B3', '#E7298A',
                   '#66A61E', '#E6AB02', '#A6761D', 
                   '#666666')[1:nCause]
    } else {
      Palette <- rainbow(nCause)
    }
    
    par(mfrow = c(ceiling(np / 2), 2), mar = c(4, 4, 1, 1))
    for (ip in 1:np) {
      thName <- namesSimp[ip]
      idth <- which(namesFull %in% sprintf("%s.%s", thName, causes))
      denl <- lapply(1:nCause, function(ic) {
        iname <- namesFull[idth[ic]]
        xthe <- x$theta[, iname]
        xdens <- density(xthe)
        q99 <- quantile(xthe, probs = c(0.001, 0.999))
        sthe <- x$coefficients[iname, 1:4]
        id95 <- which(xdens$x >= sthe["Lower"] & 
                        xdens$x <= sthe["Upper"])
        return(list(dens = xdens, q99 = q99, id95 = id95, xthe = xthe))
      })
      names(denl) <- causes
      xlim <- range(sapply(1:nCause, function(ic) {
        denl[[ic]]$q99
      }))
      ylim <- c(0, max(sapply(1:nCause, function(ic) {
        max(denl[[ic]]$dens$y)
      })))
      plot(xlim, ylim, col = NA, xlab = "", ylab = "", 
           main = thName)
      for (ic in 1:nCause) {
        xic <- denl[[ic]]$dens$x
        yic <- denl[[ic]]$dens$y
        id95 <- denl[[ic]]$id95
        lines(xic, yic, col = Palette[ic])
        polygon(c(xic[id95], rev(xic[id95])), 
                y = c(yic[id95], rep(0, length(id95))),
                col = adjustcolor(Palette[ic], alpha.f = 0.25), 
                border = NA)
      }
    }
    legend("topright", causes, col = Palette, pch = 15)
  } else if (type == "demorates") {
    par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
    if (nCause <= 7) {
      Palette <- c('#1B9E77', '#D95F02', '#7570B3', '#E7298A',
                   '#66A61E', '#E6AB02', '#A6761D')[1:nCause]
    } else {
      Palette <- rainbow(nCause)
    }
    icut <- x$cuts
    xx <- x$x[icut]
    xlim <- range(xx)
    for (idem in c("surv", "mort")) {
      if (idem == "surv") {
        yy <- x$surv
        ylim <- c(0, 1)
        ylab <- "Survival"
      } else {
        yy <- x$mort
        if (logMort) {
          lowy <- min(sapply(1:length(yy), function(ic) {
            min(log(yy[[ic]]["Mean", icut]), na.rm = TRUE)}))
          ylab <- "log-mortality"
          for (ic in 1:length(yy)) {
            yy[[ic]] <- log(yy[[ic]])
          }
        } else {
          lowy <- 0
          ylab <- "Mortality"
        }
        if (noCIs) {
          ylimSource <- "Mean"
        } else {
          ylimSource <- "Upper"
        }
        ylim <- c(lowy, max(sapply(1:length(yy), function(ic) {
          max(yy[[ic]][ylimSource, icut])
        })))
        
      }
      plot(xlim, ylim, col = NA, xlab = "Age", ylab = ylab)
      for (ic in 1:(nCause + 1)) {
        yyi <- yy[[ic]][, icut]
        if (ic <= nCause) {
          cols <- Palette[ic]
          lwd <- 2
        } else {
          cols <- "grey40"
          lwd <- 6
        }
        if (!noCIs) {
          polygon(x = c(xx, rev(xx)), 
                  y = c(yyi["Lower", ], rev(yyi["Upper", ])),
                  col = adjustcolor(cols, alpha.f = 0.25), 
                  border = NA)
        }
        lines(xx, yyi["Mean", ], col = cols, lwd = lwd)
      }
      if (idem == "surv") {
        legend("bottomleft", c(causes, "All"), 
               col = c(Palette, 'grey40'), pch = NA, 
               lwd = c(rep(2, nCause), 6))
        
      }
    }
  } else if (type == "cumInc") {
    PLE <- x$PLE
    xx <- x$x
    par(mfrow = c(ceiling(nCause / 2), 2))
    for (ic in 1:nCause) {
      plot(PLE[[ic]]$Ages, 1 - PLE[[ic]]$ple, type = 's', main = causes[ic],
           ylim = c(0, 1), xlab = "Age", ylab = "Survival")
      yy <- x$surv[[ic]]
      
      polygon(x = c(xx, rev(xx)), 
              y = 1 - c(yy["Lower", ], rev(yy["Upper", ])),
              col = adjustcolor("orange", alpha.f = 0.25), 
              border = NA)
      lines(xx, 1 - yy["Mean", ], col = "orange", lwd = 2)
    }
  } else if (type == "gof") {
    PLE <- x$PLE
    xx <- x$x
    par(mfrow = c(ceiling((nCause + 1) / 2), 2))
    for (ic in 1:nCause) {
      plot(PLE[[ic]]$Ages, PLE[[ic]]$ple, type = 's', main = causes[ic],
           ylim = c(0, 1), xlab = "Age", ylab = "Survival")
      yy <- x$surv[[ic]]
      
      polygon(x = c(xx, rev(xx)), 
              y = c(yy["Lower", ], rev(yy["Upper", ])),
              col = adjustcolor("orange", alpha.f = 0.25), 
              border = NA)
      lines(xx, yy["Mean", ], col = "orange", lwd = 2)
    }
    plot(PLE$All$Ages, PLE$All$ple, type = 's', main = "All",
         ylim = c(0, 1), xlab = "Age", ylab = "Survival")
    yy <- x$surv$All
    
    polygon(x = c(xx, rev(xx)), 
            y = c(yy["Lower", ], rev(yy["Upper", ])),
            col = adjustcolor("orange", alpha.f = 0.25), 
            border = NA)
    lines(xx, yy["Mean", ], col = "orange", lwd = 2)
  }
}

# printing:
print.BayesCR <- function(x, ...) {
  args <- list(...)
  if (is.element('digits', names(args))){
    digits <- args$digits
  } else {
    digits <- 4
  }
  
  # Coefficients:
  cat("\nCoefficients:\n")
  print.default(x$coefficients, digits, ...)
  
  # Convergence:
  cat("\nConvergence:\n")
  if (x$convergence) {
    convMessage <- "Appropriate convergence reached for all parameters.\n"
  } else {
    convMessage <- "Convergence not reached for some parameters.\n"
  }
  cat(convMessage)
  
  # DIC:
  cat("\nModel fit:\n")
  if (is.na(x$convergence)){
    cat("DIC not provided due to lack of convergence.")
  } else {
    cat(sprintf("DIC = %s\n", round(x$DIC["DIC"], 2)))
  }
}

# Summary:
summary.BayesCR <- function(object, ...) {
  args <- list(...)
  if (is.element('digits', names(args))){
    digits <- args$digits
  } else {
    digits <- 4
  }
  
  # Model settings:
  cat("\nMCMC settings:\n")
  cat(paste("Iterations :", object$settings$niter, "\n"))
  cat(paste("Burn-in    :", object$settings$burnin, "\n"))
  cat(paste("Thinning   :", object$settings$thinning, "\n"))
  cat(paste("Simulations:", object$settings$nsim, "\n"))
  cat(paste("Comp. time :", object$settings$compTime, "\n"))
  
  # Data description:
  cat("\nData:\n")
  cat(paste("Number of records :", object$data$n, "\n"))
  cat(paste("Causes of death   :", paste(object$data$causes, collapse = ", "), 
            "\n"))
  cat(paste("Number of causes  :", object$data$nCause, "\n"))
  
  # Coefficients:
  cat("\nCoefficients:\n")
  print.default(object$coefficients, digits, ...)
  
  # Convergence:
  cat("\nConvergence:\n")
  if (object$convergence) {
    convMessage <- "Appropriate convergence reached for all parameters.\n"
  } else {
    convMessage <- "Convergence not reached for some parameters.\n"
  }
  cat(convMessage)
  
  # DIC:
  cat("\nModel fit:\n")
  if (is.na(object$convergence)){
    cat("DIC not provided due to lack of convergence.")
  } else {
    cat(sprintf("DIC = %s\n", round(object$DIC["DIC"], 2)))
  }
  
  # Output if assigned:
  sumout <- list(coefficients = object$coefficients, DIC = object$DIC,
                 convergence = object$convergence, settings = object$settings)
  return(invisible(sumout))
}


# ============================= #
# ==== INTERNAL FUNCTIONS: ====
# ============================= #
# ------------------------------------- #
# ---- Internal object management: ----
# ------------------------------------- #
# Create algorithm object:
.CreateAlgObj <- function(niter, burnin, thinning, nsim, UPDJUMP, jumpSD) {
  return(list(niter = niter, burnin = burnin, thinning = thinning, 
              nsim = nsim, UPDJUMP = UPDJUMP, jumpSD = jumpSD))
}

# Create data object:
.CreateDataObj <- function(object) {
  # Number of observations:
  n <- nrow(object)
  
  # Unique causes:
  uCauses <- sort(unique(object$cause))
  
  # Number of causes:
  nCause <- length(uCauses)
  
  # Create design matrix per cause:
  causeMat <- model.matrix(~cause - 1, data = object)
  
  # long data.frame of causes:
  iCauseLong <- data.frame(cause = rep(uCauses, n))
  
  # Design matrix for long data frame:
  causeLong <- model.matrix(~cause - 1, data = iCauseLong)
  
  # Long ages (repeated by causes):
  xlong <- rep(object$age, each = nCause)
  
  # One vector per cause:
  oneCause <- rep(1, nCause)
  
  # Create data object:
  dataObj <- list(x = object$age, iCause = object$cause, xlong = xlong, n = n, 
                  causes = uCauses, nCause = nCause, causeMat = causeMat, 
                  causeLong = causeLong, oneCause = oneCause)
  return(dataObj)
}

# Create Parameter object:
.CreateParObj <- function(dataObj) {
  thetaNames <- c("a0", "a1", "c", "b0", "b1", "b2")
  p <- length(thetaNames)
  theta <- matrix(c(-1, 2, 0.001, -5, 0.3, 2), nrow = dataObj$nCause,
                  ncol = p, byrow = TRUE)
  jump <- matrix(rep(0.1, 6), nrow = dataObj$nCause,
                 ncol = p, byrow = TRUE)
  priorMean <- matrix(c(-2, 0.01, 0, -3, 0.01, 1e-10), nrow = dataObj$nCause,
                      ncol = p, byrow = TRUE)
  priorSd <- matrix(c(1, 5, 1, 1, 1, 1), nrow = dataObj$nCause,
                    ncol = p, byrow = TRUE)
  low <- matrix(c(-Inf, 0, 0, -Inf, 0, 0), nrow = dataObj$nCause,
                ncol = p, byrow = TRUE)
  jitter <- matrix(c(0.5, 0.2, 0.2, 0.5, 0.2, 0.5), nrow = dataObj$nCause,
                   ncol = p, byrow = TRUE)
  # Create matrix of theta Indices:
  ind <- matrix(1:(dataObj$nCause * p), nrow = dataObj$nCause, 
                ncol = p)
  
  dimnames(theta) <- dimnames(jump) <- dimnames(priorMean) <- 
    dimnames(priorSd) <- dimnames(low) <- dimnames(jitter) <- 
    dimnames(ind) <- list(dataObj$causes, thetaNames)
  
  # Full names:
  namesFull <- paste(rep(thetaNames, each = dataObj$nCause),
                     rep(dataObj$causes, p), sep = ".")
  
  # output:
  parObj <- list(theta = theta, priorMean = priorMean, priorSd = priorSd,
                 jump = jump, low = low, jitter = jitter, ind = ind, p = p, 
                 pc = p * dataObj$nCause, names = thetaNames,
                 namesFull = namesFull)
  return(parObj)
}

# Create analysis theta object:
.CreateThetaObj <- function(parObj, dataObj) {
  # Create matrix of theta parameters per cause:
  thetaMat <- parObj$theta
  
  # Matrix of theta vectors per individual:
  iTheta <- dataObj$causeMat %*% thetaMat
  
  # long matrix of thetaCause:
  iThetaLong <- dataObj$causeLong %*% thetaMat
  
  # Theta indices on n x p matrix:
  iInd <- dataObj$causeMat %*% parObj$ind
  
  # Theta indices on (n x nc) x p matrix:
  iIndLong <- dataObj$causeLong %*% parObj$ind
  
  # Prior:
  thPrior <- .dtnorm(x = c(thetaMat), mean = c(parObj$priorMean),
                     sd = c(parObj$priorSd), lower = c(parObj$low),
                     log = TRUE)
  thetaPrior <- matrix(thPrior, nrow = dataObj$nCause, ncol = parObj$p,
                       dimnames = dimnames(parObj$theta))
  
  # Output:
  thetaObj <- list(theta = thetaMat, iTheta = iTheta, 
                   iThetaLong = iThetaLong, iInd = iInd, iIndLong = iIndLong,
                   prior = thetaPrior)
  
  return(thetaObj)
}

# Sample theta object (Note: add jump object):
.JitterTheta <- function(thetaObj, parObj, dataObj) {
  contJitt <- TRUE
  while(contJitt) {
    thetaJitt <- thetaObj
    # Jitter all parameters:
    thJitt <- .rtnorm(n = parObj$pc, mean = c(thetaObj$theta), 
                      sd = c(parObj$jitter), lower = c(parObj$low))
    thetaJitt$theta <- matrix(thJitt, nrow = dataObj$nCause, ncol = parObj$p,
                              dimnames = dimnames(parObj$theta))
    
    # Find location of theta ip in short table:
    thetaJitt$iTheta <- dataObj$causeMat %*% thetaJitt$theta
    
    # Find location of theta[ip] in long table:
    thetaJitt$iThetaLong <- dataObj$causeLong %*% thetaJitt$theta
    
    # Update prior:
    thPrior <- .dtnorm(x = (thetaJitt$theta), mean = c(parObj$priorMean),
                       sd = c(parObj$priorSd), lower = c(parObj$low),
                       log = TRUE)
    thetaJitt$prior <- matrix(thPrior, nrow = dataObj$nCause, ncol = parObj$p,
                              dimnames = dimnames(parObj$theta))
    
    # Calculate demographic functions:
    demoJitt <- .CalcDemoFuns(thetaObj = thetaJitt, dataObj = dataObj)
    
    # Likelihood and posterior:
    likePostJitt <- .CalcLikePost(demoObj = demoJitt, thetaObj = thetaJitt, 
                                  parObj = parObj)
    
    if (!is.na(likePostJitt$post) & likePostJitt$post > -Inf &
        likePostJitt$post < Inf) {
      contJitt <- FALSE
    }
  }
  return(thetaJitt)
}

# -------------------------------- #
# ---- Demographic functions: ----
# -------------------------------- #
# Mortality:
.CalcMu <- function(theta, ...) UseMethod(".CalcMu")
.CalcMu.numeric <- function(theta, x) {
  exp(theta["a0"] - theta["a1"] * x) + theta["c"] + 
    exp(theta["b0"] + theta["b1"] * x) / 
    (1 + theta["b2"] * exp(theta["b0"]) / 
       theta["b1"] * (exp(theta["b1"] * x) - 1))
}
.CalcMu.matrix <- function(theta, x) {
  exp(theta[, "a0"] - theta[, "a1"] * x) + theta[, "c"] + 
    exp(theta[, "b0"] + theta[, "b1"] * x) / 
    (1 + theta[, "b2"] * exp(theta[, "b0"]) / 
       theta[, "b1"] * (exp(theta[, "b1"] * x) - 1))
}

# Cummulative hazard:
.CalcU <- function(theta, ...) UseMethod(".CalcU")
.CalcU.matrix <- function(theta, x) {
  exp(theta[, "a0"]) / theta[, "a1"] * (1 - exp(-theta[, "a1"] * x)) +
    theta[, "c"] * x + log(1 + theta[, "b2"] * 
                             exp(theta[, "b0"]) / theta[, "b1"] * 
                             (exp(theta[, "b1"] * x) - 1)) *
    (1 / theta[, "b2"])
}
.CalcU.numeric <- function(theta, x) {
  exp(theta["a0"]) / theta["a1"] * (1 - exp(-theta["a1"] * x)) +
    theta["c"] * x + log(1 + theta["b2"] * 
                           exp(theta["b0"]) / theta["b1"] * 
                           (exp(theta["b1"] * x) - 1)) *
    (1 / theta["b2"])
}

# Survival:
.CalcS <- function(theta, x) {
  exp(-.CalcU(theta, x))
}

# Product limit estimator:
.CalcPLE <- function(ageLast, ageFirst = NULL, departType) {
  # Set age first to 0 if NULL:
  if (is.null(ageFirst)) {
    ageFirst <- rep(0, length(ageLast))
  } 
  
  # Eliminate potential rounding errors:
  ageFirst <- round(ageFirst, 8)
  ageLast <- round(ageLast, 8)
  
  # Find records with same first and last age:
  idsame <- which(ageLast == ageFirst)
  
  # Increase last age by one day:
  if (length(idsame) > 0) {
    ageLast[idsame] <- ageLast[idsame] + 1/365.25
  }
  
  # Sort ages:
  idsort <- sort.int(ageLast, index.return = TRUE)$ix
  
  # Create new age vector:
  agev <- unique(ageLast[idsort])
  
  # Number of unique ages:
  nage <- length(agev)
  
  # Cx and delta x vector:
  Cx <- rep(0, nage)
  delx <- rep(0, nage)
  
  # Fill up Cx and delta:
  for (ii in 1:nage) {
    idNx <- which(ageFirst <= agev[ii] & ageLast >= agev[ii])
    Cx[ii] <- length(idNx) / nage
    idd <- which(ageLast == agev[ii] & departType == "D")
    delx[ii] <- length(idd)
  }
  
  # Calculate product limit estimator:
  ple <- cumprod((1 - 1 / (nage * Cx))^delx)
  
  # Add age 0:
  if (agev[1] > 0) {
    agev <- c(0, agev)
    ple <- c(1, ple)
  }
  
  # Finish at the last age at death:
  idd <- which(agev <= max(ageLast[which(departType == "D")]))
  
  # Create data frame:
  pleTab <- data.frame(Ages = agev[idd], ple = ple[idd])
  
  return(pleTab)
}

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

# ------------------------------------------ #
# ---- Likelihood, posterior, sampling: ----
# ------------------------------------------ #
# Sample theta object (Note: add jump object):
.SampleTheta <- function(thetaObj, parObj, dataObj, jumpSD, ip) {
  thetaObj$theta[ip] <- .rtnorm(n = 1, mean = thetaObj$theta[ip], 
                                sd = jumpSD[ip], lower = parObj$low[ip])
  
  # Find location of theta ip in short table:
  idIp <- which(thetaObj$iInd == ip)
  thetaObj$iTheta[idIp] <- thetaObj$theta[ip]
  
  # Find location of theta[ip] in long table:
  idIpl <- which(thetaObj$iIndLong == ip)
  thetaObj$iThetaLong[idIpl] <- thetaObj$theta[ip]
  
  # Update prior:
  thetaObj$prior[ip] <- .dtnorm(x = thetaObj$theta[ip], 
                                mean = parObj$priorMean[ip],
                                sd = parObj$priorSd[ip],
                                lower = parObj$low[ip], log = TRUE)
  
  return(thetaObj)
}

# Calculate demography functions:
.CalcDemoFuns <- function(thetaObj, dataObj) {
  # Calculate individual hazard:
  muix <- .CalcMu(theta = thetaObj$iTheta, x = dataObj$x)
  
  # Calculate cumulative hazards for all causes:
  Uix <- .CalcU(theta = thetaObj$iThetaLong, x = dataObj$xlong)
  
  # Matrix of individual hazards:
  UiMat <- t(matrix(Uix, nrow = dataObj$nCause, ncol = dataObj$n))
  
  # Additive cumulative hazard:
  Ux <- c(UiMat %*% dataObj$oneCause)
  
  # Survival:
  Sx <- exp(- Ux)
  
  # pdf of ages at death:
  fx <- Sx * muix
  
  # output:
  demoList <- list(muix = muix, Uix = Uix, Ux = Ux, Sx = Sx, fx = fx)
  return(demoList)
}

# Likelihood and posterior:
.CalcLikePost <- function(demoObj, thetaObj, parObj) {
  iLike <- -demoObj$Ux + log(demoObj$muix)
  like <- sum(iLike)
  post <- like + sum(thetaObj$prior)
  return(list(iLike = iLike, like = like, post = post))
}

# Metropolis-Hastings ratio:
.CalcMetHastRatio <- function(likePostNow, likePostNew, thetaNow, thetaNew, 
                              parObj, jumpSD, ip) {
  Hratio <- .dtnorm(x = thetaNow$theta[ip], mean = thetaNew$theta[ip],
                    sd = jumpSD[ip], lower = parObj$low[ip], log = TRUE) -
    .dtnorm(x = thetaNew$theta[ip], mean = thetaNow$theta[ip],
            sd = jumpSD[ip], lower = parObj$low[ip], log = TRUE)
  MHratio <- exp(likePostNew$post - likePostNow$post + Hratio)
  return(MHratio)
}

# Function to update jumps:
.UpdateJumps <- function(jumpSD, updMat, iter, iterUpd = 100, updTarg = 0.25,
                         dataObj, parObj) {
  updRate <- apply(updMat[iter - ((iterUpd - 1):0), ], 2, sum) / iterUpd  
  updRate[which(updRate == 0)] <- 1e-2
  jumpSD <- jumpSD * 
    matrix(updRate, dataObj$nCause, parObj$p) / updTarg
  return(jumpSD)
}

# --------------- #
# ---- MCMC: ----
# --------------- #
.RunMCMC <- function(sim, dataObj, parObj, thetaObj, algObj, .CalcMu, 
                     .CalcMu.matrix, .CalcMu.numeric, .CalcU, .CalcU.numeric,
                     .CalcU.matrix, .CalcS, .JitterTheta, .SampleTheta, 
                     .CalcDemoFuns, .CalcLikePost, .CalcMetHastRatio,
                     .UpdateJumps, .rtnorm, .dtnorm, jumpSD = NULL, 
                     UPDJUMP = TRUE) {
  
  # Start theta object for analysis:
  thetaNow <- .JitterTheta(thetaObj = thetaObj, parObj = parObj,
                           dataObj = dataObj)
  
  # Calculate demographic functions:
  demoNow <- .CalcDemoFuns(thetaObj = thetaNow, dataObj = dataObj)
  
  # Likelihood and posterior:
  likePostNow <- .CalcLikePost(demoObj = demoNow, thetaObj = thetaNow, 
                               parObj = parObj)
  
  # Start update jumps or regular analysis:
  if (UPDJUMP) {
    niter <- 10000
    jumpSD <- parObj$jump
    jumpOut <- matrix(jumpSD, nrow = 1, ncol = parObj$pc,
                      dimnames = list(NULL, parObj$namesFull))
    updJumpIter <- seq(50, niter, 50)
    updateMat <- matrix(0, nrow = niter, ncol = parObj$pc,
                        dimnames = list(NULL, parObj$namesFull))
    parOut <- NA
    likePostOut <- NA
  } else {
    niter <- algObj$niter
    parOut <- matrix(NA, nrow = niter, ncol = parObj$pc,
                     dimnames = list(NULL, parObj$namesFull))
    likePostOut <- matrix(NA, nrow = niter, ncol = 2,
                          dimnames = list(NULL, c("Likelihood", "Posterior")))
  }
  
  # Start MCMC clock:
  StartMCMC <- Sys.time()
  for (iter in 1:niter) {
    # ------------------------ #
    # Sample theta parameters:
    # ------------------------ #
    for (ip in 1:parObj$pc) {
      # Draw random theta:
      thetaNew <- .SampleTheta(thetaObj = thetaNow, parObj = parObj, 
                               dataObj = dataObj, jumpSD = jumpSD, ip = ip)
      
      # Calculate demographic functions:
      demoNew <- .CalcDemoFuns(thetaObj = thetaNew, dataObj = dataObj)
      
      # Likelihood and posterior:
      likePostNew <- .CalcLikePost(demoObj = demoNew, thetaObj = thetaNew, 
                                   parObj = parObj)
      
      # Calculate Metropolis-Hastings ratio:
      acceptRatio <- .CalcMetHastRatio(likePostNow, likePostNew, thetaNow, 
                                       thetaNew, parObj, jumpSD, ip)
      # Accept or reject:
      if (acceptRatio > runif(1)) {
        thetaNow <- thetaNew
        demoNow <- demoNew
        likePostNow <- likePostNew
        if (UPDJUMP) updateMat[iter, ip] <- 1
      }
    }
    
    # ---------------- #
    # Update jump SDs:
    # ---------------- #
    if (UPDJUMP) {
      if (iter %in% updJumpIter) {
        jumpSD <- .UpdateJumps(jumpSD, updateMat, iter, iterUpd = 50, 
                               updTarg = 0.25, dataObj, parObj)
        jumpOut <- rbind(jumpOut, c(jumpSD))
      } 
      if (iter == niter) {
        nJumps <- nrow(jumpOut)
        idJumps <- floor(nJumps / 2):nJumps
        jumpSD <- matrix(apply(jumpOut[idJumps, ], 2, mean),
                         nrow = dataObj$nCause, ncol = parObj$p,
                         dimnames = list(dataObj$causes, parObj$names))
      }
    } else {
      parOut[iter, ] <- c(thetaNow$theta)
      likePostOut[iter, ] <- c(likePostNow$like, likePostNow$post)
    }
  }
  
  # Output:
  outList <- list(theta = parOut, likePost = likePostOut, jumpSD = jumpSD)
  return(outList)
}

# ------------------------------------------- #
# ---- MCMC output management functions: ----
# ------------------------------------------- #
# based on Gelman et al. (2014).
.CalcPSRF <- function(object, keep, nsim) {
  nthin <- length(keep)
  Means <- t(sapply(1:nsim, function(i) {
    apply(object[[i]]$theta[keep, ], 2, mean)
  }))
  Vars <- t(sapply(1:nsim, function(i) {
    apply(object[[i]]$theta, 2, var)
  }))
  meanall <- apply(Means, 2, mean)
  B <- nthin / (nsim - 1) * apply(t((t(Means) - meanall)^2), 2, sum)
  W <- 1 / nsim * apply(Vars, 2, sum)
  Varpl <- (nthin - 1) / nthin * W + 1 / nthin * B
  Rhat <- sqrt(Varpl / W)
  Rhat[which(Varpl == 0)] <- 1
  Rhat[which(Rhat < 1)] <- 1
  conv <- cbind(B, W, Varpl, Rhat)
  rownames(conv) <- colnames(Means)
  return(conv)
}

# Calculate effective sample size:
.CalcNeff <- function(object, keep, nsim, Rhat) {
  nthin <- length(keep)
  Varpl <- Rhat[, "Varpl"]
  Tlags <- min(c(200, nthin - 1))
  nTheta <- ncol(object[[1]]$theta)
  rhoHat <- matrix(NA, nrow = Tlags, ncol = nTheta)
  for (tt in 1:Tlags) {
    Vt <- 1 / (nsim * (nthin - tt)) * 
      apply(sapply(1:nsim, function(im) {
        bMat <- object[[im]]$theta[keep, ]
        apply(bMat, 2, function(bi) {
          sum((bi[1:(nthin - tt + 1)] - bi[tt:nthin])^2)
        })
      }), 1, sum)
    rhoHat[tt, ] <- 1 - Vt / Varpl
  }
  
  # Find T:
  Tthresh <- apply(rhoHat, 2, function(rhoi) {
    which(rhoi[-1] + rhoi[-Tlags] < 0)[1]
  })
  
  idna <- which(is.na(Tthresh))
  if (length(idna) > 0) Tthresh[idna] <- Tlags
  
  neff <- floor(sapply(1:nTheta, function(ip) {
    nsim * nthin / (1 + 2 * sum(rhoHat[1:Tthresh[ip], ip]))
  }))
  names(neff) <- rownames(Rhat)
  return(neff)
}

# Calculate DIC:
.CalcDIC <- function(likelihood, k) {
  L <- length(likelihood)
  Dm <- -2 * likelihood
  Dave <- mean(Dm)
  pD <- 1/2 * 1/(L-1) * sum((Dm - Dave)^2)
  DIC <- Dave + pD
  modSel <- c(Dave, pD, k, DIC)
  names(modSel) <- c("D.ave", "pD", "k", "DIC")
  return(modSel)
}

# Calculate demographic quantiles:
.CalcDemoQuants <- function(thetaMat, parObj, dataObj) {
  np <- nrow(thetaMat)
  xMax <- ceiling(max(dataObj$x)) * 4
  dx <- 0.1
  xv <- seq(0, xMax, dx)
  mortList <- list()
  survList <- list()
  mortQuants <- list()
  survQuants <- list()
  for (ic in 1:dataObj$nCause) {
    icause <- dataObj$causes[ic]
    causeCols <- sprintf("%s.%s", parObj$names, icause)
    thetaCause <- thetaMat[, causeCols]
    colnames(thetaCause) <- parObj$names
    mortic <- apply(thetaCause, 1, function(theic) {
      .CalcMu(theic, xv)
    })
    mortList[[icause]] <- mortic
    mortQuants[[icause]] <- rbind(Mean = apply(mortic, 1, mean, na.rm = TRUE),
                                  Lower = apply(mortic, 1, quantile, 0.025,
                                                na.rm = TRUE, names = FALSE),
                                  Upper = apply(mortic, 1, quantile, 0.975,
                                                na.rm = TRUE, names = FALSE))
    if (ic == 1) {
      mortAll <- mortic
    } else {
      mortAll <- mortAll + mortic
    }
    
    survic <- apply(thetaCause, 1, function(theic) {
      exp(-.CalcU(theic, xv))
    })
    survList[[icause]] <- survic
    survQuants[[icause]] <- rbind(Mean = apply(survic, 1, mean, na.rm = TRUE),
                                  Lower = apply(survic, 1, quantile, 0.025,
                                                na.rm = TRUE, names = FALSE),
                                  Upper = apply(survic, 1, quantile, 0.975,
                                                na.rm = TRUE, names = FALSE))
    
    if (ic == 1) {
      survAll <- survic
    } else {
      survAll <- survAll * survic
    }
  }
  mortQuants$All <- rbind(Mean = apply(mortAll, 1, mean, na.rm = TRUE),
                               Lower = apply(mortAll, 1, quantile, 0.025,
                                             na.rm = TRUE, names = FALSE),
                               Upper = apply(mortAll, 1, quantile, 0.975,
                                             na.rm = TRUE, names = FALSE))
  survQuants$All <- rbind(Mean = apply(survAll, 1, mean, na.rm = TRUE),
                                Lower = apply(survAll, 1, quantile, 0.025,
                                              na.rm = TRUE, names = FALSE),
                                Upper = apply(survAll, 1, quantile, 0.975,
                                              na.rm = TRUE, names = FALSE))
  
  # Calculate cuts based on main survival:
  cuts <- which(survQuants$All[1, ] >= 0.01)
  return(list(mort = mortQuants, surv = survQuants, x = xv, cuts = cuts))
}

# ------------------------------------ #
# ---- Other ancillary functions: ----
# ------------------------------------ #
# Calculate time differences:
.CalcTimeDiff <- function(Start, End) {
  timeDiff <- End - Start
  diffUnits <- units(timeDiff)
  if (diffUnits == "secs" & as.numeric(timeDiff) < 60) {
    tdiffUn <- round(as.numeric(End - Start, 
                                units = 'secs'), 2)
    compTime <- sprintf("%s secs", tdiffUn)
  } else if (diffUnits == "secs" & as.numeric(timeDiff) >= 60 | 
             diffUnits == "mins" & as.numeric(timeDiff) < 60) {
    tdiffUn <- round(as.numeric(End - Start, 
                                units = 'mins'), 2)
    compTime <- sprintf("%s mins", tdiffUn)
  } else if (diffUnits == "mins" & as.numeric(timeDiff) >= 60 | 
             diffUnits == "hours") {
    tdiffUn <- round(as.numeric(End - Start, 
                                units = 'hours'), 2)
    compTime <- sprintf("%s hours", tdiffUn)
  } 
  return(compTime)
}
