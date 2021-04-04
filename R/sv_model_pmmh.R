#' Runs PMMH procedure on the simulated data example
#'
#' @param data measurements or y
#' @param particles number of particles
#' @param iterations number of PMCMC iterations
#' @param burnin burnin period (number of samples to exclude)
#' @param plot logical; if TRUE, a plot of the reslults will be returned
#'
#' @return invisibly returns PMMH output and, if plot = TRUE, a plot of the
#'  output
#' @export
sv_model_pmmh <- function(data,
                          particles= 5000,
                          iterations= 10000,
                          burnin = 0,
                          plot = FALSE) {

    if (missing(data)) {
         msg_w <- "'data' argument is missing: using data simulated from model."
         warning(msg_w)
         data <- simNonlin(len = 500,
                           var_init = 5,
                           var_evol = 10,
                           var_obs = 1,
                           cosSeqOffset = 0)$data
    }

    if (burnin >= iterations) {
        stop("Burn-in must be less than iterations!")
    }

    res <- SVmodelExamples::sv_model_pmmh_cpp(as.matrix(data),
                                              particles,
                                              iterations)

    res_plot <- res[burnin + 1:iterations, ]

    if (plot) {
        graphics::par(mfrow = c(2, 3), oma = c(0, 0, 2, 0))
        with(res_plot, graphics::hist(samples_sigv,
                       xlab = expression(sigma_v), ylab = "density",
                       main = NA))
        graphics::abline(v = sqrt(10), lty = 2, col = "darkblue", lwd = 2)
        with(res_plot, graphics::plot(samples_sigw, samples_sigv,
                       xlab = expression(sigma_w), ylab = expression(sigma_v),
                       main = NA, col = "darkblue"))
        with(res_plot, graphics::plot(samples_sigv, type = "l",
                       ylab = expression(sigma_v),
                       main = NA, xlim = c(0, iterations - burnin)))
        with(res_plot, graphics::plot(samples_sigv, samples_sigw,
                       xlab = expression(sigma_v), ylab = expression(sigma_w),
                       main = NA, col = "darkblue"))
        with(res_plot, graphics::hist(samples_sigw,
                       xlab = expression(sigma_w), ylab = "density",
                       main = NA))
        graphics::abline(v = 1, lty = 2, col = "darkblue", lwd = 2)
        with(res_plot, graphics::plot(samples_sigw, type = "l",
                       ylab = expression(sigma_w),
                       main = NA, xlim = c(0, iterations - burnin)))
        graphics::title("Posterior Estimates", outer = TRUE)
        graphics::par(mfrow = c(1, 1))
    }

    invisible(res)
}
