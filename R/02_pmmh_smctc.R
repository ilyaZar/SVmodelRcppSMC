#' Runs PMMH procedure for the toy SV model
#'
#' @param data measurements or y
#' @param starting_vals starting values as a vector in the follwoing order
#'   phi_x, sigma_x, and beta_y
#' @param rw_mh_sd a vector of two elements giving the standard deviations for the RW-MH step
#' @param particles number of particles
#' @param iterations number of PMCMC iterations
#' @param burnin burnin period (number of samples to exclude)
#' @param num_progress_outputs int giving the number of progress outputs i.e.
#'   if set to 10, then progress output occurs for every additional 10% of 
#'   completion
#'
#' @return returns PMMH output
#' 
#' @export
sv_model_pmmh <- function(data,
                          starting_vals,
                          rw_mh_sd,
                          particles = 2000,
                          iterations= 5000,
                          burnin = round(iterations / 2, digits = 0),
                          num_progress_outputs = 10) {
    if (missing(data)) {
         msg_err <- paste0("'data' argument is missing: ",
                           "either real data or simulate e.g. from ",
                           "'SVmodelPMMH::generate_data_simul_sv()'")
         stop(msg_err)
    }
    if (burnin >= iterations) {
        stop("Burn-in must be less than iterations!")
    }

    res <- SVmodelExamples::sv_model_pmmh_cpp(as.matrix(data),
                                              particles,
                                              iterations,
                                              starting_vals,
                                              rw_mh_sd,
                                              num_progress_outputs) 
    return(res)
}
