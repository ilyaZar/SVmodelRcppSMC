#' Helper function for plotting PMMH and PGAS estimation results
#'
#' @param output_pmcmc output object as provided by the pmmh estimatio functions
#' @param burnin burnin period
#' @param true_vals defaults to NULL e.g. when a real dataset is used; for simulated data
#'   can be set to true parameter values
#'
#' @return invisible return; function called for plot side effects
#' @export
plot_pmcmc_output <- function(output_pmcmc,
                              burnin,
                              true_vals = NULL) {
  if (!is.null(true_vals)) {
    out_type_simulation = TRUE
    sigma_x_true <- true_vals[1]
    beta_y_true  <- true_vals[2]
  } else {
    out_type_simulation = FALSE
  }
  MM <- length(output_pmcmc[[1]])
  samples_sigma_x <- output_pmcmc$samples_sigma_x
  samples_beta_y  <- output_pmcmc$samples_beta_y
  post_mean_sig_x <- mean(samples_sigma_x[burnin:MM])
  post_mean_bet_y <- mean(samples_beta_y[burnin:MM])
  
  graphics::par(mfrow = c(2, 4), oma = c(0, 0, 2, 0))
  
  expressions_sigma <- list(expression(Histogramm~of~sigma),
                            expression(sigma),
                            expression(sigma~trace))
  samples_list <- list(samples_sigma_x, samples_beta_y)
  plot_parameter_post_results(samples = samples_list,
                              post_mean = post_mean_sig_x,
                              true_vals = sigma_x_true,
                              expression_list = expressions_sigma,
                              burnin = burnin,
                              MM = MM,
                              out_type_simulation,
                              from_above_after_burnin = TRUE)
  expressions_beta <- list(expression(Histogramm~of~beta),
                           expression(beta),
                           expression(beta~trace))
  samples_list <- list(samples_beta_y, samples_sigma_x)
  plot_parameter_post_results(samples = samples_list,
                              post_mean = post_mean_bet_y,
                              true_vals = beta_y_true,
                              expression_list = expressions_beta,
                              burnin = burnin,
                              MM = MM,
                              out_type_simulation,
                              from_above_after_burnin = FALSE)
  title_text <- paste0("Posterior Estimatation Results: ",
                       "(red: post. mean; green: true values if provided)")
  graphics::title(main = title_text,
                  outer = TRUE)
  graphics::par(mfrow = c(1, 1))
}
plot_parameter_post_results <- function(samples,
                                        post_mean,
                                        true_vals,
                                        expression_list,
                                        burnin,
                                        MM,
                                        out_type_simulation,
                                        from_above_after_burnin) {
  expression_main  <- expression_list[[1]]
  expression_par   <- expression_list[[2]]
  expression_trace <- expression_list[[3]]
  samples_main      <- samples[[1]]
  samples_remaining <- samples[[2]]
  graphics::hist(samples_main[burnin:MM], 
                 main = expression_main,
                 xlab = expression_par)
  graphics::abline(v = post_mean,
                   col = "red",
                   lty = 2)
  if (out_type_simulation) graphics::abline(v = true_vals,
                                            col = "green",
                                            lty = 1)
  
  graphics::plot(samples_main,
                 main = "Full Trace Plot",
                 type = "l",
                 ylab = expression_trace,
                 xlab = "Number of MCMC iteration")
  graphics::abline(h = post_mean,
                   col = "red",
                   lty = 2)
  if (out_type_simulation) graphics::abline(h = true_vals,
                                            col = "green",
                                            lty = 1)
  
  graphics::plot(samples_main[burnin:MM],
                 main = "Trace Plot After Burnin",
                 type = "l",
                 ylab = expression_trace,
                 xlab = "Number of MCMC iteration")
  graphics::abline(h = post_mean,
                   col = "red",
                   lty = 2)
  if (out_type_simulation) graphics::abline(h = true_vals,
                                            col = "green",
                                            lty = 1)
  if (from_above_after_burnin) {
    samples_main      <- samples_main[burnin:MM]
    samples_remaining <- samples_remaining[burnin:MM]
    title_from_above  <- "Sample pairs after burnin"
  } else {
    title_from_above <- "All sample pairs"
  }
  graphics::plot(samples_main,
                 samples_remaining,
                 main = title_from_above,
                 xlab = expression(sigma), 
                 ylab = expression(beta),
                 col = "darkblue")
}