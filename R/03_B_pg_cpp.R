#' Runs PG procedure for the toy SV model
#'
#' @param data measurements or y
#' @param startingVals starting values as a vector in the following order
#'   phiX, sigmaX, and betaY
#' @param startingTrajectory initial reference trajectory condition for the
#'   conditional SMC part of the PG algorith
#' @param particles number of particles
#' @param iterations number of PMCMC iterations
#' @param numProgressOutputs int giving the number of progress outputs i.e.
#'   if set to 10, then progress output occurs for every additional 10% of 
#'   completion
#'
#' @return Rcpp::List containing the results: parameter samples (sigmaX, betaY)
#' 
#' @export
svModelPG <- function(data,
                      startingVals,
                      startingTrajectory,
                      particles = 2000,
                      iterations= 5000,
                      numProgressOutputs = 10) {
    if (missing(data)) {
         msgErr <- paste0("'data' argument is missing: ",
                          "either real data or simulate e.g. from ",
                          "'SVmodelPG::generateDataSimulSv()'")
         stop(msgErr)
    }

    res <- SVmodelRcppSMC::svModelPGimpl(as.matrix(data),
                                         particles,
                                         iterations,
                                         startingVals,
                                         startingTrajectory,
                                         numProgressOutputs)
    return(res)
}
