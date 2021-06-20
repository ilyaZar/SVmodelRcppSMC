#' Runs PMMH procedure for the toy SV model
#'
#' @param data measurements or y
#' @param starting_vals starting values as a vector in the follwoing order
#'   phi_x, sigma_x, and beta_y
#' @param particles number of particles
#' @param resample_freq int giving the number of progress outputs i.e.
#'   if set to 10, then progress output occurs for every additional 10% of 
#'   completion
#'
#' @return returns PMMH output
#' 
#' @export
sv_model_al_tracking <- function(data,
                                 starting_vals,
                                 particles = 10,
                                 resample_freq = 10) {
  if (missing(data)) {
    msg_err <- paste0("'data' argument is missing: ",
                      "either real data or simulate e.g. from ",
                      "'SVmodelPMMH::generate_data_simul_sv()'")
    stop(msg_err)
  }
  
  res <- SVmodelRcppSMC::sv_model_al_tracking_impl(as.matrix(data),
                                                   starting_vals,
                                                   particles,
                                                   resample_freq) 
  plot_al <- t(res[["ancestor_lines"]])
  matplot(plot_al, type = c("b"),pch = 1, col = 1:particles) #plot
  legend("topleft", legend = 1:particles, col=1:particles, pch=1)
  title(main = "Particle ancestral lines",
        xlab = "Iteration",
        ylab = "Particle values (Genealogy)")
  return(res)
}
