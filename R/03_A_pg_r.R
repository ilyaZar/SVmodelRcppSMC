sample_gibbs_params <- function(y, x, param, a, b) {
  TT <- (length(y) + 1)
  a_sigma <- a + TT/2
  b_sigma <- b + 0.5*sum((x[2:TT] - param*x[1:(TT - 1)])^2)
  sig_sq_x <- 1/rgamma(1, a_sigma, b_sigma)
  a_beta <- a + TT/2
  b_beta <- b + 0.5*sum(exp(-x[2:TT])*(y^2))
  bet_sq_x = 1/rgamma(1, a_beta, b_beta)

  sig_sq_x <- sqrt(sig_sq_x)
  bet_sq_x <- sqrt(bet_sq_x)
  return(c(sig_sq_x, bet_sq_x))
}
bpf_pg <- function(N, yt, phi_x, sigma_x, beta_y, x_r, as_sampling) {
  # Bookkeeping:
  TT      <- length(yt) + 1
  yt     <- c(NA, yt)
  # The above implies that yt[1] equals yt[2] in the for-loop of
  # the BPF and thus the iteration is goes to TT + 1 instead of TT.
  x <- matrix(0, N, TT)
  a <- matrix(0, N, TT)
  w <- matrix(0, N, TT)
  
  w_log <- numeric(N)
  # Initialization t = 0:
  x[, 1]  <- rnorm(N, mean = 0, sd = (sigma_x/sqrt(1 - phi_x^2)))
  x[N, 1] <- x_r[1]
  a[, 1]  <- 1:N
  w[, 1]  <- 1/N
  # Iteration t = 1 to t = TT. NOTE: for-loop runs from t=2 to TT+1
  for (t in 2:TT) {
    a[, t] <- sample.int(N, size = N, replace = TRUE, prob = w[, t - 1])
    # SV BPF-proposal:
    meanP <- phi_x*x[a[, t], t - 1]
    x[, t] <- rnorm(N, mean = meanP, sd = sigma_x)
    # Conditioning step:
    x[N, t] <- x_r[t]
    # ancestor sampling
    # a[N, t] <- N
    # browser()
    a[N, t] <- compute_as(num_particles = N,
                          current_x_r = x_r[t],
                          current_sig_sq_x = sigma_x^2,
                          current_cond_mean = phi_x * x[, t- 1],
                          current_w_log = log(w[, t - 1]),
                          as_sampling = as_sampling)
    # SV BPF weighting:
    w_log        <- dnorm(yt[t],
                          mean = 0,
                          sd = beta_y*exp(0.5*x[, t]),
                          log = TRUE)
    w_logmax <- max(w_log)
    w_tilde  <- exp(w_log - w_logmax)
    w[, t]   <- w_tilde/sum(w_tilde)
  }
  index <- a[, TT]
  for (t in (TT - 1):1) {
    x[, t] <- x[index, t]
    index  <- a[index, t]
  }
  as_index <- sample.int(N, size = 1, prob = w[, TT])
  x_particle_out <- x[as_index, ]
  return(x_particle_out)
}
compute_as <- function(num_particles,
                       current_x_r,
                       current_sig_sq_x,
                       current_cond_mean,
                       current_w_log,
                       as_sampling) {
  if (as_sampling) {
    m        <- -1/(2*current_sig_sq_x)*(current_x_r - current_cond_mean)^2
    w_as_log <- current_w_log + m
    w_logmax <- max(w_as_log)
    w_tilde  <- exp(w_as_log - w_logmax)
    w_as     <- w_tilde/sum(w_tilde)
    ind_out  <- sample.int(n = num_particles,
                           size = 1,
                           replace = TRUE,
                           prob = w_as)
  } else {
    ind_out <- num_particles
  }
  return(ind_out)
}
#' Runs a particle Gibbs procedure on the toy SV model
#'
#' @param y measurements i.e. observations of latent state process
#' @param num_particles number of particles
#' @param starting_values vector of three giving the starting values for the 
#'   parameters in the following order: phi_x, sigma_x, and beta_y
#' @param num_iter number of PG iterations
#' @param x_r_init initial particle to condition on (reference trajectory)
#'
#' @return 
#' @export
pg_sv_r <- function(y,
                    num_particles, 
                    starting_values,
                    num_iter, 
                    x_r_init,
                    as_sampling,
                    num_progress_outputs = 10) {
  progress_intervall_num  <- round(num_iter / num_progress_outputs)

  TT <- length(y)
  X <- matrix(0, nrow = num_iter, ncol = TT + 1)
  a <- 0.01
  b <- 0.01
  sig_sq_x <- numeric(num_iter)
  bet_sq_x <- numeric(num_iter)
  ############################################################################
  ### Initialization:
  ### Iteration m=1
  X[1, ] <- bpf_pg(N = num_particles, yt = y, 
                   phi_x   = starting_values[1],
                   sigma_x = starting_values[2],
                   beta_y  = starting_values[3],
                   x_r = x_r_init,
                   as_sampling)
  samples <- sample_gibbs_params(y = y, x = X[1, ],
                                 param = starting_values[1],
                                 a = a, b = b)
  sig_sq_x[1] <- samples[1]
  bet_sq_x[1] <- samples[2]
  X[1, ] <- bpf_pg(N = num_particles, yt = y,
                   phi_x   =  starting_values[1],
                   sigma_x = sig_sq_x[1],
                   beta_y  = bet_sq_x[1],
                   x_r = X[1, ],
                   as_sampling)
  # particleCond <- xt
  for (m in 2:num_iter) {
    samples <- sample_gibbs_params(y = y,
                                   x = X[m - 1, ], 
                                   param = starting_values[1],
                                   a = a, b = b)
    sig_sq_x[m] <- samples[1]
    bet_sq_x[m] <- samples[2]
    X[m, ] <- bpf_pg(N = num_particles, yt = y,
                     phi_x   =  starting_values[1],
                     sigma_x = sig_sq_x[m], 
                     beta_y  = bet_sq_x[m],
                     x_r = X[m - 1, ],
                     as_sampling)
    if (m %% progress_intervall_num == 0) {
      cat(sprintf("########################\n"))
      cat(sprintf(" Iteration: %d of : %d completed.\n \n",
                  m, num_iter))
      cat(sprintf(" Current posterior mean:                  %.4f %.4f  \n",
                  mean(sig_sq_x[0:m]), mean(bet_sq_x[0:m])))
      cat(sprintf("########################\n"))
    }
  }
  return(list(samplesSigmaX = sig_sq_x,
              samplesBetaY  = bet_sq_x))
}