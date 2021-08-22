#' Runs PMMH procedure for the toy SV model
#'
#' @param data measurements or y
#' @param initVals starting values as a vector in the follwoing order
#'   phiX, sigmaX, and betaY
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
#' @export
svModelALTracking <- function(data,
                              initVals,
                              particles = 10,
                              resampleFreq = 10) {
  if (missing(data)) {
    msgErr <- paste0("'data' argument is missing: ",
                      "either real data or simulate e.g. from ",
                      "'SVmodelPMMH::generateDataSimulSv()'")
    stop(msgErr)
  }
  
  browser()
  res <- SVmodelRcppSMC::svModelALtrackingImp(as.matrix(data),
                                              initVals,
                                              particles,
                                              resampleFreq) 
  plotAL <- t(res[["ancestorLines"]])
  matplot(plotAL, type = c("b"),pch = 1, col = 1:particles) #plot
  legend("topleft", legend = 1:particles, col=1:particles, pch=1)
  title(main = "Particle ancestral lines (full genealogy)",
        xlab = "Iteration",
        ylab = "Particle values")
  return(invisible(res))
}