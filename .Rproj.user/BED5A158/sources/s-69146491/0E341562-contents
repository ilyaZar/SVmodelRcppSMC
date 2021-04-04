#' Ricker example 
#'
#' We perform the nested SMC algorithm of Duan and Fulop (2015) on the Ricker
#' example from Fasiolo et al. (2016).
#'
#' @details
#' Fasiolo et al. (2016) consider completing models for the evolution of a
#' population over time.
#'
#' Under one of the models, the population \eqn{N_t} is modelled at time \code{t}
#' under the discrete time stochastic process
#' \deqn{N_{t-1} = r N_{t} exp[-N_{t} + z_{t+1}]}
#' where \eqn{z_t ~ N(0,\sigma^2)}. The true population is observed with error
#' as \eqn{y_t ~ Pois(\phi N_t)}. We wish to obtain approximate samples from the
#' posterior distribution for parameters
#' \eqn{\theta = [log r, log \phi,log \sigma]} and also to estimate
#' the model evidence so that we can perform model choice (other
#' models will be added to the package later).
#' 
#' The prior, denoted \eqn{p(\theta)}, for this model is:
#' \deqn{log r ~ U(2,5)}
#' \deqn{log \phi ~ U(1.61,3)}
#' \deqn{log \sigma ~ U(-3,-0.22).}
#' We can get an unbiased estimate of the likelihood, \eqn{f(y|\theta)} using a
#' bootstrap particle filter.
#'
#' We use the method of Duan and Fulop (2015) to obtain approximate samples
#' from the posterior distribution. We start by drawing an initial population
#' of \eqn{\theta} from the prior and we traverse the particles through a
#' series of distributions with targets
#' \eqn{\pi_t \propto f(y|theta)^{\gamma_t} p(\theta)}, where \eqn{t=0,\ldots,T}
#' and \eqn{0=\gamma_0 \le \ldots \le \gamma_T = 1}. This is a sequential
#' Monte Carlo (SMC) method, so we are performing two levels of SMC.
#' 
#' @references
#'  J.-C. Duan and A. Fulop (2015). Density-tempered Marginalized
#' Sequential Monte Carlo Samplers. Journal of Business & Economic Statistics,
#' 33(2):192-202.
#'
#' M. Fasiolo, N. Pya and S. Wood (2016). A Comparison of Inferential Methods for
#' Highly Nonlinear State Space Models in Ecology and Epidemiology. Statistical
#' Science, 31(1):96-118.
"_PACKAGE"
#> [1] "_PACKAGE"