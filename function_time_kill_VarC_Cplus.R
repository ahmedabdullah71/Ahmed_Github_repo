# ============================================================
# Piecewise survival function with fixed work c0
# c = c0 + exponential variable work
# ============================================================
d_time_kill_VarC2_0 <- function(tau, lambda, mu, c0) {
  if (any(lambda < 0, na.rm = TRUE) ||
      any(mu < 0, na.rm = TRUE) ||
      any(c0 < 0, na.rm = TRUE)) {
    stop("lambda, mu, and c0 must be > 0.")
  }
  
  if (any(tau < 0, na.rm = TRUE)) {
    stop("tau must be >= 0.")
  }
  
  # avoid division by zero at tau = 0
  tau <- pmax(tau, 1e-12)
  
  ifelse(
    tau < c0,
    1 - (mu * tau / (lambda + mu * tau)) *
      exp(-lambda * (c0 / tau - 1)),
    
    (lambda * exp(mu * c0 - mu * tau)) /
      (lambda + mu * tau)
  )
}


# ============================================================
# Interval survival prediction
# prediction = S0(tau) - S0(tau + delta_t)
# ============================================================
d_time_kill_VarC2 <- function(tau, lambda, mu, delta_t, c0) {
  
  out <- d_time_kill_VarC2_0(tau, lambda, mu, c0) -
    d_time_kill_VarC2_0(tau + delta_t, lambda, mu, c0)
  
  out
}


uniformly_divide_0To1 <- function(NN) {
  seq(1, 1 / NN, length = NN) - 1 / (2 * NN)
}


# ============================================================
# Negative log-likelihood function
# ============================================================
ll_function_VarC2 <- function(delta_t, lambda, mu, sd, c0, t, SF,
                              min_delta_t = 0,
                              epsilon = 1e-12,
                              min_c0 = 0) {
  
  c0 <- min_c0 + exp(c0)
  delta_t <- min_delta_t + exp(delta_t)
  lambda <- exp(lambda)
  mu <- exp(mu)
  
  SF_pred <- d_time_kill_VarC2(
    tau = t,
    delta_t = delta_t,
    lambda = lambda,
    mu = mu,
    c0 = c0
  )
  
  SF_pred <- pmax(pmin(SF_pred, 1 - epsilon), epsilon)
  
  y <- -sum(dnorm(
    log(SF),
    mean = log(SF_pred + epsilon),
    sd = sd,
    log = TRUE
  ))
  
  y
}


# ============================================================
# GA fitting function
# ============================================================
fit_VarC2_GA <- function(t, SF,
                         fixed = list(),
                         lower = list(
                           delta_t = log(1e-6),
                           lambda = log(1e-6),
                           mu = log(1e-6),
                           sd = 1e-6,
                           c0 = log(1e-6)
                         ),
                         upper = list(
                           delta_t = log(1e6),
                           lambda = log(1e6),
                           mu = log(1e2),
                           sd = 10,
                           c0 = log(1e6)
                         ),
                         start = list(
                           delta_t = log(50),
                           lambda = log(10),
                           mu = log(0.01),
                           sd = 1.4,
                           c0 = log(500)
                         ),
                         min_delta_t = 0,
                         epsilon = 1e-12,
                         min_c0 = 0,
                         popSize = 100,
                         maxiter = 200,
                         run = 50,
                         seed = 1,
                         parallel = FALSE) {
  
  if (!requireNamespace("GA", quietly = TRUE)) {
    stop("Please install GA: install.packages('GA')")
  }
  
  par_names <- c("delta_t", "lambda", "mu", "sd", "c0")
  
  # --- validate fixed names
  bad_fixed <- setdiff(names(fixed), par_names)
  if (length(bad_fixed) > 0) {
    stop("Unknown fixed params: ", paste(bad_fixed, collapse = ", "))
  }
  
  free_names <- setdiff(par_names, names(fixed))
  if (length(free_names) == 0) {
    stop("All parameters are fixed; nothing to optimize.")
  }
  
  # --- build vectors for GA bounds in correct order
  lower_vec <- unlist(lower[free_names], use.names = FALSE)
  upper_vec <- unlist(upper[free_names], use.names = FALSE)
  
  # --- objective for GA
  fitness_fun <- function(x_free) {
    theta <- start
    theta[names(fixed)] <- fixed
    theta[free_names] <- as.list(x_free)
    
    nll <- ll_function_VarC2(
      delta_t = theta$delta_t,
      lambda = theta$lambda,
      mu = theta$mu,
      sd = theta$sd,
      c0 = theta$c0,
      t = t,
      SF = SF,
      min_delta_t = min_delta_t,
      epsilon = epsilon,
      min_c0 = min_c0
    )
    
    -nll
  }
  
  set.seed(seed)
  
  ga_fit <- GA::ga(
    type = "real-valued",
    fitness = fitness_fun,
    lower = lower_vec,
    upper = upper_vec,
    popSize = popSize,
    maxiter = maxiter,
    run = run,
    parallel = parallel,
    seed = seed,
    monitor = TRUE
  )
  
  # --- extract best
  best_free <- as.numeric(ga_fit@solution[1, ])
  
  theta_hat <- start
  theta_hat[names(fixed)] <- fixed
  theta_hat[free_names] <- as.list(best_free)
  
  # --- compute final NLL
  nll_hat <- ll_function_VarC2(
    delta_t = theta_hat$delta_t,
    lambda = theta_hat$lambda,
    mu = theta_hat$mu,
    sd = theta_hat$sd,
    c0 = theta_hat$c0,
    t = t,
    SF = SF,
    min_delta_t = min_delta_t,
    epsilon = epsilon,
    min_c0 = min_c0
  )
  
  # --- natural scale
  natural <- list(
    delta_t = min_delta_t + exp(theta_hat$delta_t),
    lambda = exp(theta_hat$lambda),
    mu = exp(theta_hat$mu),
    sd = theta_hat$sd,
    c0 = min_c0 + exp(theta_hat$c0)
  )
  
  list(
    ga = ga_fit,
    nll = nll_hat,
    par_opt_scale = theta_hat,
    par_natural = natural,
    fixed = fixed,
    free_names = free_names
  )
}