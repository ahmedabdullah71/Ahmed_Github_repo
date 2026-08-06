
# Survival function: S0(tau) = (lambda * exp(-mu * tau)) / (lambda + mu * tau), tau >= 0
d_time_kill_VarC_0 <- function(tau, lambda, mu) {
  if (any(lambda < 0, na.rm = TRUE) || any(mu < 0, na.rm = TRUE)) {
    print(paste('lambda',lambda, sep=':'))
    print(paste('mu',mu, sep=':'))
    
    stop("lambda and mu must be > 0.")
  }
  if (any(tau < 0, na.rm = TRUE)) {
    stop("tau must be >= 0.")
  }
  (lambda * exp(-mu * tau)) / (lambda + mu * tau)
}

# Interval survival prediction with delay d:
# for tau < d -> 1
# else -> S0(tau+delta_t-d) - S0(tau-d)
d_time_kill_VarC <- function(tau, lambda, mu, delta_t, d) {
  # if (any(lambda < 0, na.rm = TRUE) || any(mu < 0, na.rm = TRUE)) {
  #   stop("lambda and mu must be > 0.")
  # }
  # if (any(tau < 0, na.rm = TRUE)) {
  #   stop("tau must be >= 0.")
  # }
  
  out <- rep(1, length(tau))  # default for tau < d
  
  idx <- tau >= d
  if (any(idx)) {
    out[idx] <- d_time_kill_VarC_0(tau[idx] - d, lambda, mu) - d_time_kill_VarC_0(tau[idx] + delta_t - d, lambda, mu) 
    
  }
  
  out
}




uniformly_divide_0To1 <- function(NN){
  seq(1,1/NN, length = NN)-1/(2*NN)
}


ll_function_VarC <- function(delta_t, lambda, mu, sd, d, t, SF,
                             min_delta_t = 0, epsilon = 1e-12 ,min_d =0) {
  d <- min_d+exp(d)
  delta_t <- min_delta_t + exp(delta_t)
  lambda <- exp(lambda)
  mu <- exp(mu)
  
  SF_pred <- d_time_kill_VarC(tau = t, delta_t = delta_t,
                              lambda = lambda, mu = mu, d = d)
  
  SF_pred <- pmax(pmin(SF_pred, 1 - epsilon), epsilon)
  
  y <- -sum(dnorm(log(SF),
                  mean = log(SF_pred + epsilon),
                  sd = sd,
                  log = TRUE))
  y
}


fit_VarC_GA <- function(t, SF,
                        fixed = list(),          # e.g. list(delta_t=log(fixed_delta_t), d=log(0))
                        lower = list(delta_t=log(1e-6), lambda=log(1e-6), mu=log(1e-6), sd=1e-6, d=log(1e-6)),
                        upper = list(delta_t=log(1e6),  lambda=log(1e6),  mu=log(1e2),  sd=10,   d=log(1e6)),
                        start = list(delta_t=log(50), lambda=log(10), mu=log(0.01), sd=1.4, d=log(91)),
                        min_delta_t = 0,
                        epsilon = 1e-12,
                        min_d =0, 
                        popSize = 100,
                        maxiter = 200,
                        run = 50,
                        seed = 1,
                        parallel = FALSE) {
  
  if (!requireNamespace("GA", quietly = TRUE)) {
    stop("Please install GA: install.packages('GA')")
  }
  
  par_names <- c("delta_t", "lambda", "mu", "sd", "d")
  
  # --- validate fixed names
  bad_fixed <- setdiff(names(fixed), par_names)
  if (length(bad_fixed) > 0) stop("Unknown fixed params: ", paste(bad_fixed, collapse = ", "))
  
  free_names <- setdiff(par_names, names(fixed))
  if (length(free_names) == 0) stop("All parameters are fixed; nothing to optimize.")
  
  # --- build vectors for GA bounds in correct order
  lower_vec <- unlist(lower[free_names], use.names = FALSE)
  upper_vec <- unlist(upper[free_names], use.names = FALSE)
  
  # --- objective for GA (GA maximizes fitness, so use -NLL)
  fitness_fun <- function(x_free) {
    theta <- start
    theta[names(fixed)] <- fixed
    theta[free_names] <- as.list(x_free)
    
    nll <- ll_function_VarC(delta_t = theta$delta_t,
                            lambda  = theta$lambda,
                            mu      = theta$mu,
                            sd      = theta$sd,
                            d       = theta$d,
                            t = t, SF = SF,
                            min_delta_t = min_delta_t,
                            epsilon = epsilon,
                            min_d = min_d)
    
    # GA maximizes
    -nll
  }
  
  set.seed(seed)
  
  ga_fit <- GA::ga(type = "real-valued",
                   fitness = fitness_fun,
                   lower = lower_vec,
                   upper = upper_vec,
                   popSize = popSize,
                   maxiter = maxiter,
                   run = run,
                   parallel = parallel,
                   seed = seed,
                   monitor = TRUE)
  
  # --- extract best
  
  best_free <- as.numeric(ga_fit@solution[1, ])
  # names(ga_sol) <- c("mu", "lambda", "t0", "sd")
  # 
  # best_free <- as.numeric(GA::ga_summary(ga_fit)$solution) # may be 1 row or multiple
  # if (is.matrix(ga_fit@solution)) best_free <- as.numeric(ga_fit@solution[1, ])
  # 
  theta_hat <- start
  theta_hat[names(fixed)] <- fixed
  theta_hat[free_names] <- as.list(best_free)
  
  # --- compute final NLL at best
  nll_hat <- ll_function_VarC(delta_t = theta_hat$delta_t,
                              lambda  = theta_hat$lambda,
                              mu      = theta_hat$mu,
                              sd      = theta_hat$sd,
                              d       = theta_hat$d,
                              t = t, SF = SF,
                              min_delta_t = min_delta_t,
                              epsilon = epsilon,
                              min_d = min_d)
  
  # --- transform to natural scale (matches your ll_function)
  natural <- list(
    delta_t = min_delta_t + exp(theta_hat$delta_t),
    lambda  = exp(theta_hat$lambda),
    mu      = exp(theta_hat$mu),
    sd      = theta_hat$sd,
    d       = min_d +exp(theta_hat$d)
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

