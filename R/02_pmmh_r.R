#' Boostrap Particle Filter for the toy SV model
#' 
#' Computes, besides other things, an estimate of the log-likelihood.
#'
#' @param NN number of particles
#' @param yt measurements/data of the unobservered latent state process
#' @param phi_x phi parameter (see model equations for details)
#' @param sigma_x sigma parameter (see model equations for details)
#' @param beta_y beta parameter (see model equations for details)
#'
#' @return a list of one element (was previously returniong more) giving the 
#'   log-likelihood estimate
#' @export
bpf_pmmh <- function(NN, yt, phi_x, sigma_x, beta_y) {
  # Bookkeeping:
  TT      <- length(yt) + 1
  yt     <- c(NA, yt)
  # The above implies that yt[1] equals yt[2] in the for-loop of
  # the BPF and thus the iteration is goes to TT + 1 instead of TT.
  p_store <- matrix(0, NN, TT)
  w_store <- matrix(0, NN, TT)
  a_store <- matrix(0, NN, TT)
  x_fltr <- numeric(TT)
  loglik <- 0
  # Initialization t = 0:
  p_store[, 1] <- rnorm(NN, mean = 0, sd = (sigma_x / sqrt(1 - phi_x^2)))
  w_store[, 1] <- 1 / NN
  a_store[, 1] <- 1:NN
  x_fltr      <- p_store[, 1]
  # Iteration t = 1 to t = TT. NOTE: for-loop runs from t=2 to TT+1
  for (t in 2:TT) {
    a_store[, t] <- sample.int(NN,
                               size = NN,
                               replace = TRUE,
                               prob = w_store[, t - 1])
    # SV BPF-proposal:
    particle_mean <- phi_x * p_store[a_store[, t], t - 1]
    p_store[, t] <- rnorm(NN, mean = particle_mean, sd = sigma_x)
    # SV BPF weighting:
    log_weights        <- dnorm(yt[t],
                                mean = 0,
                                sd = beta_y * exp(0.5 * p_store[, t]),
                                log = TRUE)
    log_weight_max <- max(log_weights)
    log_weights_normalized <- exp(log_weights - log_weight_max)
    w_store[, t] <- log_weights_normalized / sum(log_weights_normalized)
    loglik       <- loglik + log_weight_max + log(sum(log_weights_normalized))
  }
  loglik <- loglik - (T - 1) * log(NN)
  # x_fltr <- apply(p_store * w_store, 2, sum)
  # return(list(x_fltr = x_fltr, loglik = loglik,
  #             p_store = p_store, w_store = w_store))
  return(list(loglik = loglik))
}
#' Runs a Particle Marginal Metropolis Hastings on the toy SV model
#'
#' @param num_particles number of particles
#' @param rw_vcm_prop variance-covariance matrix of the RW-MH proposal
#' @param y measurements
#' @param starting_vals vector of starting values for phi, sigma_x, and beta_y
#' @param num_iter number of PMMH iterations
#'
#' @return a list of parameter samples of for sigma-y and beta-y
#' @export
pmmh_sv <- function(num_particles,
                    rw_vcm_prop,
                    y,
                    starting_vals,
                    num_iter,
                    num_progress_outputs = 10) {
  progress_intervall_num  <- round(num_iter / num_progress_outputs)
  phi_x_init  <- starting_vals[1]
  sig_x_init  <- starting_vals[2]
  beta_y_init <- starting_vals[3]
  # Bookkeeping
  ll      <- 0
  llp     <- 0
  samples_prop <- rep(0, times = 2)
  samples <- matrix(0, nrow = num_iter, ncol = 2)
  accepts <- numeric(num_iter)
  # Initilaize by PF pre-run:
  samples[1, 1] <- sig_x_init
  samples[1, 2] <- beta_y_init
  
  out_bpf         <- bpf_pmmh(NN = num_particles,
                              yt = y,
                              phi_x = phi_x_init,
                              sigma_x = sqrt(samples[1, 1]),
                              beta_y  = sqrt(samples[1, 2]))
  ll <- out_bpf$loglik
  # Iteration of PMH:
  for (m in 2:num_iter) {
    # RW proposal
    samples_prop <- MASS::mvrnorm(1, mu = samples[m - 1, ], Sigma = rw_vcm_prop)
    if ((samples_prop[1] > 0) && (samples_prop[2] > 0)) {
      # Compute acceptance probability
      out_bpf <- bpf_pmmh(NN = num_particles,
                          yt = y, phi_x = phi_x_init,
                          sigma_x = sqrt(samples_prop[1]),
                          beta_y  = sqrt(samples_prop[2]))
      llp <- out_bpf$loglik
      
      prio_ratio_sigma <- (-dgamma(samples_prop[1],
                                   shape = 0.01,
                                   rate = 0.01,
                                   log = TRUE)
                           + dgamma(samples[m - 1, 1],
                                    shape = 0.01,
                                    rate = 0.01,
                                    log = TRUE))
      prior_ratio_beta <- (-dgamma(samples_prop[2],
                                   shape = 0.01,
                                   rate = 0.01,
                                   log = TRUE)
                           + dgamma(samples[m - 1, 2],
                                    shape = 0.01,
                                    rate = 0.01,
                                    log = TRUE))
      ll_ratio <- llp - ll
      alpha <- exp(prio_ratio_sigma + prior_ratio_beta + ll_ratio)
      # Determine next sample
      if (runif(1) <= alpha) {
        samples[m, ] <- samples_prop
        accepts[m]   <- 1
        ll <- llp
      } else {
        samples[m, ] <- samples[m - 1, ]
        accepts[m]   <- 0
      }
    } else {
      samples[m, ] <- samples[m - 1, ]
      accepts[m]   <- 0
    }
    if (m %% progress_intervall_num == 0) {
      cat(sprintf("########################\n"))
      cat(sprintf(" Iteration: %d of : %d completed.\n \n",
                  m, num_iter))
      cat(sprintf(" Current posterior mean:                  %.4f %.4f  \n",
                  mean(samples[0:m, 1]), mean(samples[0:m, 2])))
      cat(sprintf(" Current acceptance rate:                 %.4f \n",
                  mean(accepts[0:m])))
      cat(sprintf("########################\n"))
    }
  }
  return(list(samples_sigma_x = samples[, 1],
              samples_beta_y  = samples[, 2]))
}