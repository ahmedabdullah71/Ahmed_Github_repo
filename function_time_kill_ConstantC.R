d_time_kill = Vectorize(function(t, lambda_c){
  K_t = 1-exp(-lambda_c/t)
  return(K_t)
},'t')
d_time_kill_min = Vectorize(function(t, lambda,c){
  K_t = 1-exp(-lambda*c/t)*exp(lambda)
  return(K_t)
},'t')
# d_time_kill_interval = Vectorize(function(t, delta_t, lambda_c){
#   K_t = exp(-lambda_c/(t+delta_t)) - exp(-lambda_c/t)
#   return(K_t)
# },'t')

# d_time_kill_interval = Vectorize(function(t, delta_t, lambda_c,log=F, n = 1, d=20){
#   if (t<d){
#     if (log) return(0)
#     return(1)
#   }
#   K_t = exp(-lambda_c/(t+delta_t-d)^n) - exp(-lambda_c/(t-d)^n)
#   if (log){
#     return(log(K_t))
#   }
#   return(K_t)
# },'t')
# 
# d_time_kill_interval = Vectorize(function(t, delta_t, lambda, c,log=F, n = 1, d=20){
#   if (t<d){
#     if (log) return(0)
#     return(1)
#   }
#   K_t = exp(-lambda*c/(t+delta_t-d)^n) - exp(-lambda*c/(t-d)^n)
#   if (log){
#     return(log(K_t))
#   }
#   return(K_t)
# },'t')
# 
# # d_time_kill_interval_min = Vectorize(function(t, delta_t, lambda, c,log=F, n = 1, d=20){
# #   if (t<d){
# #     if (log) return(0)
# #     return(1)
# #   }
# #   
# #   K_t = exp(lambda)*(exp(-lambda*c/(t+delta_t-d)^n) - exp(-lambda*c/(t-d)^n))
# #   #K_t = (exp(-lambda*c/(t+delta_t-d)^n) - exp(-lambda*c/(t-d)^n))
# #   if (log){
# #     return(log(K_t))
# #   }
# #   return(K_t)
# # },'t')
# 
# d_time_kill_interval_min = Vectorize(function(t, delta_t, lambda, c,log=F, n = 1, d=20, trunc = F){
#   if (t<d){
#     if (log) return(0)
#     return(1)
#   }
#   if (t >= c){
#     if (log) return(-Inf)
#     return(0)
#   }
#   K_t = exp(lambda)*(exp(-lambda*c/(t+delta_t-d)^n) - exp(-lambda*c/(t-d)^n))
#   #K_t = (exp(-lambda*c/(t+delta_t-d)^n) - exp(-lambda*c/(t-d)^n))
#   if (log){
#     return(log(K_t))
#   }
#   return(K_t)
# },'t')

# Base (non-truncated) window probability: K(t) = F(t+Δt) - F(t)
d_time_kill_interval = Vectorize(function(t, delta_t, lambda, c, log = FALSE, n = 1, d = 20){
  if (t < d){
    if (log) return(0)
    return(1)
  }
  K_t = exp(-lambda*c/(t + delta_t - d)^n) - exp(-lambda*c/(t - d)^n)
  if (log){
    return(log(K_t))
  }
  return(K_t)
}, 't')


# "Min-change" wrapper: optionally use truncated version by delegating to d_time_kill_interval
d_time_kill_interval_min = Vectorize(function(t, delta_t, lambda, c, log = FALSE, n = 1, d = 20, trunc = FALSE){
  
  # If trunc==TRUE, just use the base function (as you requested)
  if (isTRUE(trunc)){
    return(d_time_kill_interval(t = t, delta_t = delta_t, lambda = lambda, c = c, log = log, n = n, d = d))
  }
  
  # Otherwise use your scaled form with exp(lambda)
  if (t < d){
    if (log) return(0)
    return(1)
  }
  
  # (your suggested cutoff behavior when not truncating)
  if (t >= c){
    if (log) return(-Inf)
    return(0)
  }
  
  K_t = exp(lambda) * (exp(-lambda*c/(t + delta_t - d)^n) - exp(-lambda*c/(t - d)^n))
  # K_t = (exp(-lambda*c/(t+delta_t-d)^n) - exp(-lambda*c/(t-d)^n))
  
  if (log){
    return(log(K_t))
  }
  return(K_t)
  
}, 't')

uniformly_divide_0To1 <- function(NN){
  seq(1,1/NN, length = NN)-1/(2*NN)
}
fit_cons_GA <- function(t, SF,
                        ll_fun = ll_function_min,
                        trunc = FALSE,
                        fixed = list(),  # e.g. list(delta_t = log(fixed_delta_t))
                        lower = list(delta_t=log(1e-6), lambda=log(1e-6), c=log(1e-6), sd=1e-6, d=log(1e-6)),
                        upper = list(delta_t=log(1e6),  lambda=log(1e6),  c=log(1e6),  sd=10,   d=log(1e6)),
                        start = list(delta_t=log(50),   lambda=log(0.05), c=log(100),  sd=1.4227693, d=log(45)),
                        min_delta_t = 0,
                        min_d = 0,
                        popSize = 100,
                        maxiter = 200,
                        run = 50,
                        seed = 1,
                        parallel = FALSE,
                        monitor = TRUE) {
  
  if (!requireNamespace("GA", quietly = TRUE)) {
    stop("Please install GA: install.packages('GA')")
  }
  
  par_names <- c("delta_t", "lambda", "c", "sd", "d")
  
  bad_fixed <- setdiff(names(fixed), par_names)
  if (length(bad_fixed) > 0) stop("Unknown fixed params: ", paste(bad_fixed, collapse = ", "))
  
  free_names <- setdiff(par_names, names(fixed))
  if (length(free_names) == 0) stop("All parameters are fixed; nothing to optimize.")
  
  lower_vec <- unlist(lower[free_names], use.names = FALSE)
  upper_vec <- unlist(upper[free_names], use.names = FALSE)
  
  environment(ll_fun)$t  <- t
  environment(ll_fun)$SF <- SF
  
  fitness_fun <- function(x_free) {
    theta <- start
    theta[names(fixed)] <- fixed
    theta[free_names] <- as.list(x_free)
    
    # Call ll_fun with ONLY the parameters it likely expects
    nll <- ll_fun(delta_t = theta$delta_t,
                  lambda  = theta$lambda,
                  c       = theta$c,
                  sd      = theta$sd,
                  d       = theta$d,
                  trunc = trunc)
    
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
                   monitor = monitor)
  
  sol <- ga_fit@solution
  best_free <- if (is.null(dim(sol))) as.numeric(sol) else as.numeric(sol[1, ])
  
  theta_hat <- start
  theta_hat[names(fixed)] <- fixed
  theta_hat[free_names] <- as.list(best_free)
  
  nll_hat <- ll_fun(delta_t = theta_hat$delta_t,
                    lambda  = theta_hat$lambda,
                    c       = theta_hat$c,
                    sd      = theta_hat$sd,
                    d       = theta_hat$d,
                    trunc = trunc)
  
  natural <- list(
    delta_t = min_delta_t + exp(theta_hat$delta_t),
    lambda  = exp(theta_hat$lambda),
    c       = exp(theta_hat$c),
    sd      = theta_hat$sd,
    d       = min_d + exp(theta_hat$d)
  )
  
  list(
    ga = ga_fit,
    nll = nll_hat,
    par_opt_scale = theta_hat,
    par_natural = natural,
    fixed = fixed,
    free_names = free_names,
    trunc = trunc
  )
}



# analytic
# Fit at logscale 
ll_function <- function(delta_t, lambda, sd, d) {
  # n = exp(n)
  d= exp(d) 
  # if (d > t[2]){
  #   d =  t[2]
  # }
  delta_t = min_delta_t+exp(delta_t)
  lambda = exp(lambda)
  SF_pred <- d_time_kill_interval(t, delta_t, lambda, d =d)
  SF_pred <- pmax(pmin(SF_pred, 1 - epsilon), epsilon)  
  y = -sum(dnorm(log(SF), mean = log(SF_pred+epsilon), sd = sd, log = TRUE))
  return(y)
}

ll_function_min <- function(delta_t, lambda, c, sd, d, trunc = FALSE) {
  # n = exp(n)
  d= exp(d) 
  c= exp(c)
  
  delta_t = min_delta_t+exp(delta_t)
  lambda = exp(lambda)
  SF_pred <- d_time_kill_interval_min(t, delta_t=delta_t, lambda=lambda, c=c, d =d, trunc = trunc)
  SF_pred <- pmax(pmin(SF_pred, 1 - epsilon), epsilon)  
  y = -sum(dnorm(log(SF), mean = log(SF_pred+epsilon), sd = sd, log = TRUE))
  return(y)
}