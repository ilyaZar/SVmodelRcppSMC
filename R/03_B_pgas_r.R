#' Title
#'
#' @param M 
#' @param y 
#' @param Z 
#' @param par_prior 
#' @param N 
#' @param par_inits 
#'
#' @return
#' @export
#'
#' @examples
pgas <- function(M, y, Z, par_prior, N, par_inits) {
  T <- length(y)
  
  sig_sq_x <- numeric(M)
  bet_sq_y <- numeric(M)
  phi_x    <- numeric(M)
  bet_x    <- matrix(0, nrow = length(par_inits[[4]]), ncol = M)
  
  regs       <- matrix(0, nrow = T - 1, ncol = ncol(Z) + 1) # ncol(Z) + 1)
  Z          <- as.matrix(Z)
  regs[, -1] <- Z[2:T, ]
  
  w <- numeric(N)
  X <- matrix(0, nrow = M, ncol = T)
  #  Initialize the parameters
  sig_sq_x[1] <- par_inits[[1]]
  bet_sq_y[1] <- par_inits[[2]]
  phi_x[1]    <- par_inits[[3]]
  bet_x[, 1]  <- par_inits[[4]]
  #  Initialize priors
  prior_a     <- par_prior[1]
  prior_b     <- par_prior[2]
  # Initialize the states by running a PF
  cpfOut <- cpf_as(y = y, N = N, Z = Z,
                   sig_sq_x = sig_sq_x[1],
                   bet_sq_y = bet_sq_y[1],
                   phi_x = phi_x[1],
                   bet_x = bet_x[, 1, drop = F],
                   x_r = X[1, ])
  w      <- cpfOut[[2]][, T]
  b      <- sample.int(n = N, size = 1, replace = TRUE, prob = w)
  X[1, ] <- cpfOut[[1]][b, ]
  # Run MCMC loop
  for (m in 2:M) {
    how_long(m, M)
    # Run GIBBS PART
    err_sig_sq_x <- X[m - 1, 2:T] - f(xtt = X[m - 1, 1:(T - 1)],
                                      z = Z[2:T, , drop = F],
                                      t = 1:(T - 1), 
                                      phi_x = phi_x[m - 1],
                                      bet_x = bet_x[, m - 1])
    sig_sq_x[m]  <- 1/rgamma(n = 1, prior_a + (T - 1)/2,
                             prior_b + crossprod(err_sig_sq_x)/2)
    # err_bet_sq_y <- y^2 * exp(-X[m - 1,])
    # bet_sq_y[m]  <- 1/rgamma(n = 1, prior_a + T/2,
    #                          prior_b + sum(err_bet_sq_y)/2)
    # phi_x[m]     <- phi_x[m - 1]
    # bet_x[, m]   <- bet_x[, m - 1]
    # sig_sq_x[m]  <- true_sigSQ_x
    bet_sq_y[m]  <- true_betSQ_y
    regs[, 1]    <- X[m - 1, 1:(T - 1)]
    x_lhs        <- X[m - 1, 2:T]
    Omega_xa     <- solve(crossprod(regs, regs)/sig_sq_x[m] + 1)
    mu_xa        <- Omega_xa %*% (crossprod(regs, x_lhs)/sig_sq_x[m])
    beta_xa      <- rmvnorm(n = 1, mean = mu_xa, sigma = Omega_xa)
    phi_x[m]     <- beta_xa[1]  # true_phi_x
    bet_x[, m]   <- beta_xa[-1] # true_bet_x
    while (near(abs(phi_x[m]), 1, tol = 0.01) | abs(phi_x[m]) > 1) {
      beta_xa      <- rmvnorm(n = 1, mean = mu_xa, sigma = Omega_xa)
      phi_x[m]     <- beta_xa[1]  # true_phi_x
      bet_x[, m]   <- beta_xa[-1] # true_bet_x
    }
    # 
    # R <- chol(sigma, pivot = TRUE)
    # R[, order(attr(R, "pivot"))]
    # retval <- matrix(rnorm(n * ncol(sigma)),nrow = n,byrow =!pre0.9_9994)%*%R
    # retval <- sweep(retval, 2, mean, "+")
    # 
    # Run CPF-AS PART
    cpfOut <- cpf_as(y = y, N = N, Z = Z,
                     sig_sq_x = sig_sq_x[m],
                     bet_sq_y = bet_sq_y[m],
                     phi_x = phi_x[m],
                     bet_x = bet_x[, m, drop = F],
                     x_r = X[m - 1,])
    w      <- cpfOut[[2]][, T]
    # draw b
    b <- sample.int(n = N, size = 1, replace = TRUE, prob = w)
    X[m, ] <- cpfOut[[1]][b, ]   
  }
  return(list(sigma_sq_x = sig_sq_x, beta_sq_y = bet_sq_y, 
              phi_x = phi_x, bet_x = bet_x, xtraj = X))
}
cpf_as <- function(y, N, Z, sig_sq_x, bet_sq_y, phi_x, bet_x, x_r) {
  # DATA CONTAINER
  T <- length(y)                      # time period from t = 0 to t = T
  x <- matrix(0, nrow = N, ncol = T)  # Particles
  a <- matrix(0, nrow = N, ncol = T)  # Ancestor indices
  w <- matrix(0, nrow = N, ncol = T)  # Weights
  
  # I. INITIALIZATION (t = 0)
  # Deterministic initial condition
  # x[, 1]  <- 8*cos(1.2) # 0 # -0.22  
  # corresponding weights
  # w[, 1]  <- 1/N
  # Sampling initial condition from prior 
  x[, 1]  <- rnorm(n = N, mean = Z[1, , drop = F] %*% bet_x/(1 - phi_x),
                   sd = sqrt(sig_sq_x/(1 - phi_x^2)))
  # corresponding weights
  w[, 1]  <- 1/N
  # II. FIRST PERIOD APPROXIMATION (t = 1)
  # resampling
  a[, 1]     <- sample.int(n = N, replace = TRUE, prob = w[, 1])
  # propagation
  eval_f     <- f(xtt = x[, 1], z = Z[1, , drop = F], t = 1,
                  phi_x = phi_x, bet_x = bet_x)
  x[, 1]     <- eval_f[a[, 1]] + sqrt(sig_sq_x)*rnorm(N)
  # corresponding weights
  w_log   <- -1/(2*bet_sq_y*exp(x[, 1]))*(y[1])^2 - 0.5 * x[, 1]
  w_max   <- max(w_log)
  w_tilde <- exp(w_log - w_max)
  w[, 1]  <- w_tilde/sum(w_tilde)
  # conditioning
  x[N, 1] <- x_r[1]  
  
  # II. FOR t = 2,..,T
  for (t in 2:T) {
    # resampling
    a[, t]     <- sample.int(n = N, replace = TRUE, prob = w[, t - 1])
    # propagation
    eval_f     <- f(xtt = x[, t - 1], z = Z[t, , drop = F], t = t - 1,
                    phi_x = phi_x, bet_x = bet_x)
    x[, t]     <- eval_f[a[, t]] + sqrt(sig_sq_x)*rnorm(N)
    # conditioning
    x[N, t]    <- x_r[t]
    # ancestor sampling
    m          <- exp(-1/(2*sig_sq_x)*(x_r[t] - eval_f)^2)
    w_as       <- w[, t - 1]*m
    w_as       <- w_as/sum(w_as)
    a[N, t]    <- sample.int(n = N, size = 1, replace = TRUE, prob = w_as)
    # weighting
    w_log      <- -1/(2*bet_sq_y*exp(x[, t]))*(y[t])^2 - 0.5 * x[, t]
    w_max      <- max(w_log)           
    w_tilde    <- exp(w_log - w_max)
    w[, t]     <- w_tilde/sum(w_tilde)
  }
  # trajectories
  ind <- a[, T]
  for (t in (T - 1):1) {
    x[, t] <- x[ind, t]
    ind    <- a[ind, t]
  }
  return(list(x, w))
}
how_long <- function(current, total, len = 5) {
  print_iter <- total/len
  if ((current %% print_iter) == 0) {
    cat(sprintf("###############################################################
              \n"))
    cat(sprintf(" Iteration: %d of : %d completed.\n",
                current, total))
    cat(sprintf(" %.2f%% completed.\n",
                current*100/total))
  }
}
f <- function(xtt, z, t, phi_x, bet_x) {
  # xt <- phi_x*xtt
  xt <- phi_x*xtt + z %*% bet_x
  # xt <- phi_x*xtt + 8*cos(1.2*t)
  # xt <- phi_x*xtt + 25*xtt/(1 + xtt^2) 
  # xt <- phi_x*xtt + 25*xtt/(1 + xtt^2) + 8*cos(1.2*t)
  return(xt)
}
# g <- function(xt) {
#   yt <- xt^2/20
#   return(yt)
# }
