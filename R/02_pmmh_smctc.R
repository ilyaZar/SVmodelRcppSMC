#' Runs PMMH procedure for the toy SV model
#'
#' @param data measurements or y
#' @param initVals starting values as a vector in the follwoing order
#'   phiX, sigmaX, and betaY
#' @param rwMHsd a vector of two elements giving the standard deviations for the RW-MH step
#' @param particles number of particles
#' @param iterations number of PMCMC iterations
#' @param burnin burnin period (number of samples to exclude)
#' @param numProgressOutputs int giving the number of progress outputs i.e.
#'   if set to 10, then progress output occurs for every additional 10% of 
#'   completion
#'
#' @return returns PMMH output
#' 
#' @export
svModelPMMH <- function(data,
                        initVals,
                        rwMHsd,
                        particles = 2000,
                        iterations= 5000,
                        burnin = round(iterations / 2, digits = 0),
                        numProgressOutputs = 10) {
    if (missing(data)) {
         msgErr <- paste0("'data' argument is missing: ",
                           "either real data or simulate e.g. from ",
                           "'SVmodelPMMH::generateDataSimulSV()'")
         stop(msgErr)
    }
    if (burnin >= iterations) {
        stop("Burn-in must be less than iterations!")
    }

    res <- SVmodelRcppSMC::svModelPMMHimpl(as.matrix(data),
                                           particles,
                                           iterations,
                                           initVals,
                                           rwMHsd,
                                           numProgressOutputs) 
    return(res)
}
