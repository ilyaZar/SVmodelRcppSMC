sampleGibbsParams <- function(y, x, param, a, b) {
  TT <- (length(y) + 1)
  aSigma <- a + TT/2
  bSigma <- b + 0.5*sum((x[2:TT] - param*x[1:(TT - 1)])^2)
  sigSqX <- 1/rgamma(1, aSigma, bSigma)
  aBeta <- a + TT/2
  bBeta <- b + 0.5*sum(exp(-x[2:TT])*(y^2))
  betSqX = 1/rgamma(1, aBeta, bBeta)

  sigSqX <- sqrt(sigSqX)
  betSqX <- sqrt(betSqX)
  return(c(sigSqX, betSqX))
}
bpfPG <- function(N, yt, phiX, sigmaX, betaY, xR, asSampling = FALSE) {
  # Bookkeeping:
  TT      <- length(yt) + 1
  yt     <- c(NA, yt)
  # The above implies that yt[1] equals yt[2] in the for-loop of
  # the BPF and thus the iteration is goes to TT + 1 instead of TT.
  x <- matrix(0, N, TT)
  a <- matrix(0, N, TT)
  w <- matrix(0, N, TT)
  
  wLog <- numeric(N)
  # Initialization t = 0:
  x[, 1]  <- rnorm(N, mean = 0, sd = (sigmaX/sqrt(1 - phiX^2)))
  x[N, 1] <- xR[1]
  a[, 1]  <- 1:N
  w[, 1]  <- 1/N
  # Iteration t = 1 to t = TT. NOTE: for-loop runs from t=2 to TT+1
  for (t in 2:TT) {
    a[, t] <- sample.int(N, size = N, replace = TRUE, prob = w[, t - 1])
    # SV BPF-proposal:
    meanP <- phiX*x[a[, t], t - 1]
    x[, t] <- rnorm(N, mean = meanP, sd = sigmaX)
    # Conditioning step:
    x[N, t] <- xR[t]
    # ancestor sampling
    # a[N, t] <- N
    # browser()
    a[N, t] <- computeAS(numParticles = N,
                         currentXR = xR[t],
                         currentSigSqX = sigmaX^2,
                         currentCondMean = phiX * x[, t- 1],
                         currentWLog = log(w[, t - 1]),
                         asSampling = asSampling)
    # SV BPF weighting:
    wLog        <- dnorm(yt[t],
                         mean = 0,
                         sd = betaY*exp(0.5*x[, t]),
                         log = TRUE)
    wLogmax <- max(wLog)
    wTilde  <- exp(wLog - wLogmax)
    w[, t]   <- wTilde/sum(wTilde)
  }
  index <- a[, TT]
  for (t in (TT - 1):1) {
    x[, t] <- x[index, t]
    index  <- a[index, t]
  }
  asIndex <- sample.int(N, size = 1, prob = w[, TT])
  xParticleOut <- x[asIndex, ]
  return(xParticleOut)
}
computeAS <- function(numParticles,
                      currentXR,
                      currentSigSqX,
                      currentCondMean,
                      currentWLog,
                      asSampling) {
  if (asSampling) {
    m       <- -1/(2*currentSigSqX)*(currentXR - currentCondMean)^2
    wASLog  <- currentWLog + m
    wLogmax <- max(wASLog)
    wTilde  <- exp(wASLog - wLogmax)
    wAS     <- wTilde/sum(wTilde)
    indOut  <- sample.int(n = numParticles,
                          size = 1,
                          replace = TRUE,
                          prob = wAS)
  } else {
    indOut <- numParticles
  }
  return(indOut)
}
#' Runs a particle Gibbs procedure on the toy SV model
#'
#' @param y measurements i.e. observations of latent state process
#' @param numParticles number of particles
#' @param startingValues vector of three giving the starting values for the 
#'   parameters in the following order: phiX, sigmaX, and betaY
#' @param numIter number of PG iterations
#' @param xRinit initial particle to condition on (reference trajectory)
#' @param asSampling logical; if TRUE, then ancestor sampling is performed
#' @param numProgressOutputs  int giving the number of progress outputs i.e.
#'   if set to 10, then progress output occurs for every additional 10% of 
#'   completion
#'
#' @return Rcpp::List containing the results: parameter samples (sigmaX, betaY)
#' @export
svModelPGr <- function(y,
                       numParticles, 
                       startingValues,
                       numIter, 
                       xRinit,
                       asSampling = FALSE,
                       numProgressOutputs = 10) {
  progressIntervallNum  <- round(numIter / numProgressOutputs)

  TT <- length(y)
  X <- matrix(0, nrow = numIter, ncol = TT + 1)
  a <- 0.01
  b <- 0.01
  sigSqX <- numeric(numIter)
  betSqX <- numeric(numIter)
  ############################################################################
  ### Initialization:
  ### Iteration m=1
  X[1, ] <- bpfPG(N = numParticles, yt = y, 
                  phiX   = startingValues[1],
                  sigmaX = startingValues[2],
                  betaY  = startingValues[3],
                  xR = xRinit,
                  asSampling)
  samples <- sampleGibbsParams(y = y, x = X[1, ],
                               param = startingValues[1],
                               a = a, b = b)
  sigSqX[1] <- samples[1]
  betSqX[1] <- samples[2]
  X[1, ] <- bpfPG(N = numParticles, yt = y,
                  phiX   =  startingValues[1],
                  sigmaX = sigSqX[1],
                  betaY  = betSqX[1],
                  xR = X[1, ],
                  asSampling)
  # particleCond <- xt
  for (m in 2:numIter) {
    samples <- sampleGibbsParams(y = y,
                                 x = X[m - 1, ], 
                                 param = startingValues[1],
                                 a = a, b = b)
    sigSqX[m] <- samples[1]
    betSqX[m] <- samples[2]
    X[m, ] <- bpfPG(N = numParticles, yt = y,
                    phiX   =  startingValues[1],
                    sigmaX = sigSqX[m], 
                    betaY  = betSqX[m],
                    xR = X[m - 1, ],
                    asSampling)
    if (m %% progressIntervallNum == 0) {
      cat(sprintf("########################\n"))
      cat(sprintf(" Iteration: %d of : %d completed.\n \n",
                  m, numIter))
      cat(sprintf(" Current posterior mean:                  %.4f %.4f  \n",
                  mean(sigSqX[0:m]), mean(betSqX[0:m])))
    }
  }
  return(list(samplesSigmaX = sigSqX,
              samplesBetaY  = betSqX))
}
