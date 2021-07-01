#' Helper function simulating from the toy SV model
#'
#' @param TT number of time periods
#' @param phiX phi parameter of state process
#' @param sigmaX sigma parameter of measurements y
#' @param betaY beta parameter of measurements
#' @param initStateX0 state initialization
#'
#' @return a named list of two: simulated state trajectory (statesXt) and
#'  corresponding measurements (measurementsYt) under the model
#' @export
generateDataSimulSV <- function(TT, phiX, sigmaX, betaY, initStateX0) {
  x0 <- initStateX0
  xt <- numeric(TT)
  yt <- numeric(TT)
  xt[1] <- phiX * x0 + rnorm(1, 0, sd = sigmaX)
  for (t in 2:TT) {
    xt[t] <- phiX * xt[t - 1] + rnorm(1, 0, sd = sigmaX)
  }
  yt <- rnorm(TT, mean = 0, sd = betaY * exp(0.5 * xt))
  return(list(statesXt = xt, measurementsYt = yt))
}