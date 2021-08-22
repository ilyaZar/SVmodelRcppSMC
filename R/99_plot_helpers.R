#' Helper function for plotting PMMH and PGAS estimation results
#'
#' @param outputPMCMC output object as provided by the pmmh estimatio functions
#' @param burnin burnin period
#' @param trueVals defaults to NULL e.g. when a real dataset is used; for simulated data
#'   can be set to true parameter values
#'
#' @return invisible return; function called for plot side effects
#' @export
plotPMCMCoutput <- function(outputPMCMC,
                            burnin,
                            trueVals = NULL) {
  if (!is.null(trueVals)) {
    outTypeSimulation = TRUE
    sigmaXtrue <- trueVals[1]
    betaYtrue  <- trueVals[2]
  } else {
    outTypeSimulation = FALSE
  }
  MM <- length(outputPMCMC[[1]])
  samplesSigmaX <- outputPMCMC$samplesSigmaX
  samplesBetaY  <- outputPMCMC$samplesBetaY
  postMeanSigX <- mean(samplesSigmaX[burnin:MM])
  postMeanBetZ <- mean(samplesBetaY[burnin:MM])
  
  graphics::par(mfrow = c(2, 4), oma = c(0, 0, 2, 0))
  
  expressionsSigma <- list(expression(Histogramm~of~sigma),
                            expression(sigma),
                            expression(sigma~trace))
  
  samplesList <- list(samplesSigmaX, samplesBetaY)
  plotParameterPostResults(samples = samplesList,
                           postMean = postMeanSigX,
                           trueVals = sigmaXtrue,
                           expressionList = expressionsSigma,
                           burnin = burnin,
                           MM = MM,
                           outTypeSimulation,
                           fromAboveAfterBurnin = TRUE)
  expressionsBeta <- list(expression(Histogramm~of~beta),
                           expression(beta),
                           expression(beta~trace))
  samplesList <- list(samplesBetaY, samplesSigmaX)
  plotParameterPostResults(samples = samplesList,
                           postMean = postMeanBetZ,
                           trueVals = betaYtrue,
                           expressionList = expressionsBeta,
                           burnin = burnin,
                           MM = MM,
                           outTypeSimulation,
                           fromAboveAfterBurnin = FALSE)
  titleText <- paste0("Posterior Estimatation Results: ",
                      "(red: post. mean; green: true values if provided)")
  graphics::title(main = titleText,
                  outer = TRUE)
  graphics::par(mfrow = c(1, 1))
}
plotParameterPostResults <- function(samples,
                                     postMean,
                                     trueVals,
                                     expressionList,
                                     burnin,
                                     MM,
                                     outTypeSimulation,
                                     fromAboveAfterBurnin) {
  expressionMain   <- expressionList[[1]]
  expressionPar    <- expressionList[[2]]
  expressionTrace  <- expressionList[[3]]
  samplesMain      <- samples[[1]]
  samplesRemaining <- samples[[2]]
  graphics::hist(samplesMain[burnin:MM], 
                 main = expressionMain,
                 xlab = expressionPar)
  graphics::abline(v = postMean,
                   col = "red",
                   lty = 2)
  if (outTypeSimulation) graphics::abline(v = trueVals,
                                          col = "green",
                                          lty = 1)
  
  graphics::plot(samplesMain,
                 main = "Full Trace Plot",
                 type = "l",
                 ylab = expressionTrace,
                 xlab = "Number of MCMC iteration")
  graphics::abline(h = postMean,
                   col = "red",
                   lty = 2)
  if (outTypeSimulation) graphics::abline(h = trueVals,
                                          col = "green",
                                          lty = 1)
  
  graphics::plot(samplesMain[burnin:MM],
                 main = "Trace Plot After Burnin",
                 type = "l",
                 ylab = expressionTrace,
                 xlab = "Number of MCMC iteration")
  graphics::abline(h = postMean,
                   col = "red",
                   lty = 2)
  if (outTypeSimulation) graphics::abline(h = trueVals,
                                          col = "green",
                                          lty = 1)
  if (fromAboveAfterBurnin) {
    samplesMain      <- samplesMain[burnin:MM]
    samplesRemaining <- samplesRemaining[burnin:MM]
    titleFromAbove   <- "Sample pairs after burnin"
  } else {
    titleFromAbove <- "All sample pairs"
  }
  graphics::plot(samplesMain,
                 samplesRemaining,
                 main = titleFromAbove,
                 xlab = expression(sigma), 
                 ylab = expression(beta),
                 col = "darkblue")
}
#' Generates a plot of boxplots of repeated log-likelihood BPF estimates 
#' 
#' Loglikelihood should be obtained via a previous BPF run, that estimates
#' the SV model likelihood. These are passed, along with a sequence for the 
#' autoregressive parameter phi, the number of simulations (for each individual)
#' boxplot.
#'
#' @param outLogLike a matrix of dimension \code{length(parSeqPhi)} (rows)
#'   times \code{numSimul}
#' @param parSeqPhi parameter sequence for phi parameter for which a 
#'   log-likelihood estimate is computed
#' @param numSimul number of BPF runs per element of parSeqPhi
#'
#' @return invisible return; pure side effect of a boxplot generation
#' @export
plotBoxplotLogLikeEstimates <- function(outLogLike,
                                        parSeqPhi,
                                        numSimul) {
  parCol   <- as.character(rep(parSeqPhi, times = numSimul))
  llikeCol <- as.vector(outLogLike)
  dataBoxplots <- data.frame(parCol, llikeCol)
  
  names(dataBoxplots) <- c("parameter", "loglikeValue")
  xExpression <- expression(paste("Parameter values for ", phi))
  
  p <- ggplot2::ggplot(dataBoxplots,
                       ggplot2::aes(x = .data$parameter,
                                    y = .data$loglikeValue)) +
    ggplot2::geom_boxplot(varwidth = TRUE,
                          outlier.colour = "red",
                          outlier.shape = 1) +
    ggplot2::scale_x_discrete(name =  xExpression) +
    ggplot2::scale_y_continuous(name = "Log-likelihood values") +
    ggplot2::geom_jitter(width = 0.2)
  plot(p)
  return(invisible(NULL))
}
