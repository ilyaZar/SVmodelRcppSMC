#' Runs PMMH procedure for the toy SV model
#'
#' @param data measurements or y
#' @param starting_vals starting values as a vector in the follwoing order
#'   phi_x, sigma_x, and beta_y
#' @param particles number of particles
#' @param resampleFreq frequency at which resampling is performed; if negative,
#'   then resampling is never performed; if between [0,1), then resampling is
#'   performed when  the ESS falls below that proportion of the number of
#'   particles and when it is greater than or equal to 1, resampling is carried
#'   out when the ESS falls below that value (note: if this parameter is larger
#'   than the total number of particles, then resampling will always be
#'   performed!)
#'
#' @return returns PMMH output
#' 
#' @export
sv_model_al_tracking <- function(data,
                                 starting_vals,
                                 particles = 10,
                                 resampleFreq = 10) {
  if (missing(data)) {
    msg_err <- paste0("'data' argument is missing: ",
                      "either real data or simulate e.g. from ",
                      "'SVmodelPMMH::generate_data_simul_sv()'")
    stop(msg_err)
  }
  
  res <- SVmodelRcppSMC::sv_model_al_tracking_impl(as.matrix(data),
                                                   starting_vals,
                                                   particles,
                                                   resampleFreq) 
  plot_al <- t(res[["ancestor_lines"]])
  matplot(plot_al, type = c("b"),pch = 1, col = 1:particles) #plot
  legend("topleft", legend = 1:particles, col=1:particles, pch=1)
  title(main = "Particle ancestral lines (full genealogy)",
        xlab = "Iteration",
        ylab = "Particle values")
  return(invisible(res))
}
