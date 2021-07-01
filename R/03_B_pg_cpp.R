#' Runs PG procedure for the toy SV model
#'
#' @param data measurements or y
#' @param startingVals starting values as a vector in the follwoing order
#'   phi_x, sigma_x, and beta_y
#' @param startingTrajectory
#' @param particles number of particles
#' @param iterations number of PMCMC iterations
#' @param burnin burnin period (number of samples to exclude)
#' @param numProgressOutputs int giving the number of progress outputs i.e.
#'   if set to 10, then progress output occurs for every additional 10% of 
#'   completion
#'
#' @return returns PG output
#' 
#' @export
svmodelPG <- function(data,
                      startingVals,
                      startingTrajectory,
                      particles = 2000,
                      iterations= 5000,
                      burnin = round(iterations / 2, digits = 0),
                      numProgressOutputs = 10) {
    if (missing(data)) {
         msgErr <- paste0("'data' argument is missing: ",
                          "either real data or simulate e.g. from ",
                          "'SVmodelPG::generateDataSimulSv()'")
         stop(msgErr)
    }
    if (burnin >= iterations) {
        stop("Burn-in must be less than iterations!")
    }

    res <- SVmodelRcppSMC::svmodelPGimpl(as.matrix(data),
                                         particles,
                                         iterations,
                                         startingVals,
                                         startingTrajectory,
                                         numProgressOutputs)
    return(res)
}
