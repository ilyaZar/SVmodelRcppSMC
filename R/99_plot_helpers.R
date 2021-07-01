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
                           expression_list = expressionsSigma,
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
                           expression_list = expressionsBeta,
                           burnin = burnin,
                           MM = MM,
                           outTypeSimulation,
                              fromAboveAfterBurnin = FALSE)
  title_text <- paste0("Posterior Estimatation Results: ",
                       "(red: post. mean; green: true values if provided)")
  graphics::title(main = title_text,
                  outer = TRUE)
  graphics::par(mfrow = c(1, 1))
}
plotParameterPostResults <- function(samples,
                                     postMean,
                                     trueVals,
                                     expression_list,
                                     burnin,
                                     MM,
                                     outTypeSimulation,
                                     fromAboveAfterBurnin) {
  expressionMain    <- expression_list[[1]]
  expressionPar     <- expression_list[[2]]
  expressionTrace   <- expression_list[[3]]
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
#' @param out_log_like a matrix of dimension \code{length(par_seq_phi)} (rows)
#'   times \code{num_simul}
#' @param par_seq_phi parameter sequence for phi parameter for which a 
#'   log-likelihood estimate is computed
#' @param num_simul number of BPF runs per element of par_seq_phi
#'
#' @return invisible return; pure side effect of a boxplot generation
#' @export
plot_boxplot_loglike_estimates <- function(out_log_like,
                                           par_seq_phi,
                                           num_simul) {
  par_col   <- as.character(rep(par_seq_phi, times = num_simul))
  llike_col <- as.vector(out_log_like)
  data_boxplots <- data.frame(par_col, llike_col)
  
  names(data_boxplots) <- c("parameter", "loglike_value")
  x_expression <- expression(paste("Parameter values for ", phi))
  
  p <- ggplot2::ggplot(data_boxplots,
                       ggplot2::aes(x = .data$parameter,
                                    y = .data$loglike_value)) +
    ggplot2::geom_boxplot(varwidth = TRUE,
                          outlier.colour = "red",
                          outlier.shape = 1) +
    ggplot2::scale_x_discrete(name =  x_expression) +
    ggplot2::scale_y_continuous(name = "Log-likelihood values") +
    ggplot2::geom_jitter(width = 0.2)
  plot(p)
  return(invisible(NULL))
}