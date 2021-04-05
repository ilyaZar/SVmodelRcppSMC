#' Helper function simulating from the toy SV model
#'
#' @param TT number of time periods
#' @param phi_x phi parameter of state process
#' @param sigma_x sigma parameter of measurements y
#' @param beta_y beta parameter of measurements
#' @param init_state_x0 state initialization
#'
#' @return a named list of two: simulated state trajectory (states_xt) and
#'  corresponding measurements (measurements_yt) under the model
#' @export
generate_data_simul_sv <- function(TT, phi_x, sigma_x, beta_y, init_state_x0) {
  x0 <- init_state_x0
  xt <- numeric(TT)
  yt <- numeric(TT)
  xt[1] <- phi_x * x0 + rnorm(1, 0, sd = sigma_x)
  for (t in 2:TT) {
    xt[t] <- phi_x * xt[t - 1] + rnorm(1, 0, sd = sigma_x)
  }
  yt <- rnorm(TT, mean = 0, sd = beta_y * exp(0.5 * xt))
  return(list(states_xt = xt, measurements_yt = yt))
}