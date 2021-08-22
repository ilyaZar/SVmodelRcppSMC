#' Boostrap Particle Filter for the toy SV model
#' 
#' Computes, besides other things, an estimate of the log-likelihood.
#'
#' @param NN number of particles
#' @param yt measurements/data of the unobservered latent state process
#' @param phiX phi parameter (see model equations for details)
#' @param sigmaX sigma parameter (see model equations for details)
#' @param betaY beta parameter (see model equations for details)
#'
#' @return a list of one element (was previously returniong more) giving the 
#'   log-likelihood estimate
#' @export
bpfPMMH <- function(NN, yt, phiX, sigmaX, betaY) {
  # Bookkeeping:
  TT      <- length(yt) + 1
  yt     <- c(NA, yt)
  # The above implies that yt[1] equals yt[2] in the for-loop of
  # the BPF and thus the iteration is goes to TT + 1 instead of TT.
  pStore <- matrix(0, NN, TT)
  wStore <- matrix(0, NN, TT)
  aStore <- matrix(0, NN, TT)
  xFltr <- numeric(TT)
  loglik <- 0
  # Initialization t = 0:
  pStore[, 1] <- rnorm(NN, mean = 0, sd = (sigmaX / sqrt(1 - phiX^2)))
  wStore[, 1] <- 1 / NN
  aStore[, 1] <- 1:NN
  xFltr      <- pStore[, 1]
  # Iteration t = 1 to t = TT. NOTE: for-loop runs from t=2 to TT+1
  for (t in 2:TT) {
    aStore[, t] <- sample.int(NN,
                               size = NN,
                               replace = TRUE,
                               prob = wStore[, t - 1])
    # SV BPF-proposal:
    particleMean <- phiX * pStore[aStore[, t], t - 1]
    pStore[, t]  <- rnorm(NN, mean = particleMean, sd = sigmaX)
    # SV BPF weighting:
    logWeights        <- dnorm(yt[t],
                               mean = 0,
                               sd = betaY * exp(0.5 * pStore[, t]),
                               log = TRUE)
    logWeightMax <- max(logWeights)
    logWeightsNormalized <- exp(logWeights - logWeightMax)
    wStore[, t] <- logWeightsNormalized / sum(logWeightsNormalized)
    loglik      <- loglik + logWeightMax + log(sum(logWeightsNormalized))
  }
  loglik <- loglik - (T - 1) * log(NN)
  # xFltr <- apply(pStore * wStore, 2, sum)
  # return(list(xFltr = xFltr, loglik = loglik,
  #             pStore = pStore, wStore = wStore))
  return(list(loglik = loglik))
}
#' Runs a Particle Marginal Metropolis Hastings on the toy SV model
#'
#' @param numParticles number of particles
#' @param rwVCMprop variance-covariance matrix of the RW-MH proposal
#' @param y measurements
#' @param initVals vector of starting values for phi, sigmaX, and betaY
#' @param numIter number of PMMH iterations
#' @param numProgressOutputs  int giving the number of progress outputs i.e.
#'   if set to 10, then progress output occurs for every additional 10% of 
#'   completion
#' @return a list of parameter samples of for sigma-y and beta-y
#' @export
svModelPMMHr <- function(numParticles,
                         rwVCMprop,
                         y,
                         initVals,
                         numIter,
                         numProgressOutputs = 10) {
  progressIntervallNum  <- round(numIter / numProgressOutputs)
  phiXinit  <- initVals[1]
  sigXinit  <- initVals[2]
  betaYinit <- initVals[3]
  # Bookkeeping
  ll      <- 0
  llp     <- 0
  samplesProp <- rep(0, times = 2)
  samples <- matrix(0, nrow = numIter, ncol = 2)
  accepts <- numeric(numIter)
  # Initilaize by PF pre-run:
  samples[1, 1] <- sigXinit
  samples[1, 2] <- betaYinit
  
  outBPF <- bpfPMMH(NN = numParticles,
                    yt = y,
                    phiX = phiXinit,
                    sigmaX = sqrt(samples[1, 1]),
                    betaY  = sqrt(samples[1, 2]))
  ll <- outBPF$loglik
  # Iteration of PMH:
  for (m in 2:numIter) {
    # RW proposal
    samplesProp <- MASS::mvrnorm(1, mu = samples[m - 1, ], Sigma = rwVCMprop)
    if ((samplesProp[1] > 0) && (samplesProp[2] > 0)) {
      # Compute acceptance probability
      outBPF <- bpfPMMH(NN = numParticles,
                        yt = y, phiX = phiXinit,
                        sigmaX = sqrt(samplesProp[1]),
                        betaY  = sqrt(samplesProp[2]))
      llp <- outBPF$loglik
      
      priorRatioSigma <- (-dgamma(samplesProp[1],
                                  shape = 0.01,
                                  rate = 0.01,
                                  log = TRUE)
                           + dgamma(samples[m - 1, 1],
                                    shape = 0.01,
                                    rate = 0.01,
                                    log = TRUE))
      priorRatioBeta <- (-dgamma(samplesProp[2],
                                 shape = 0.01,
                                 rate = 0.01,
                                 log = TRUE)
                           + dgamma(samples[m - 1, 2],
                                    shape = 0.01,
                                    rate = 0.01,
                                    log = TRUE))
      llRatio <- llp - ll
      alpha <- exp(priorRatioSigma + priorRatioBeta + llRatio)
      # Determine next sample
      if (runif(1) <= alpha) {
        samples[m, ] <- samplesProp
        accepts[m]   <- 1
        ll <- llp
      } else {
        samples[m, ] <- samples[m - 1, ]
        accepts[m]   <- 0
      }
    } else {
      samples[m, ] <- samples[m - 1, ]
      accepts[m]   <- 0
    }
    if (m %% progressIntervallNum == 0) {
      cat(sprintf("########################\n"))
      cat(sprintf(" Iteration: %d of : %d completed.\n \n",
                  m, numIter))
      cat(sprintf(" Current posterior mean:                  %.4f %.4f  \n",
                  mean(samples[0:m, 1]), mean(samples[0:m, 2])))
      cat(sprintf(" Current acceptance rate:                 %.4f \n",
                  mean(accepts[0:m])))
      cat(sprintf("########################\n"))
    }
  }
  return(list(samplesSigmaX = samples[, 1],
              samplesBetaY  = samples[, 2]))
}
