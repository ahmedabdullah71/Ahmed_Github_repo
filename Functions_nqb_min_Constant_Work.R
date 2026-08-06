# Metabolic rate = 1

dcolony1 <- function(t, c,lambda, log_cond = FALSE, bin = 30){
  
  #f_t=c/t^2*lambda*exp(lambda)*exp(-lambda*c/t)
  f_t=log(c)+log(lambda)+lambda-2*log(t)-(lambda*c)/(t)
  f_t[t<=0|t>c] = -Inf
  if(log_cond) return(f_t)
  return(exp(f_t))
  
}

qcolony1 = function(pp, c, lambda){
  
  qq=c*lambda/(lambda-log(pp))
  return(qq)
}

compound_fun1 = function(t, Observed_T, c, lambda,t0, sd){
  h_t = dcolony1(t, c, lambda=lambda)*dnorm(Observed_T-t,t0,sd=sd)
  #print(h_t)
  return(h_t)
}

# For a given set of parameters the probability of the appearance time (Observed_T)
dcompound1 = Vectorize(function(c, lambda,t0, Observed_T, sd, log_cond = TRUE){
  
  b=integral_ahmed(fun = compound_fun1, 0, Observed_T, Observed_T=Observed_T, sd = sd, c=c, lambda=lambda, t0=t0,no_intervals=1/sd*10) 
  if (log_cond) return(log(b))
  return(b)
},'Observed_T')

# # Double Integral 
# pcompound1 = Vectorize(function(c, lambda,t0, sd, t_cdf,  log_cond = TRUE){
#   
#   b=integral_ahmed(fun = dcompound1, 0, t_cdf, c=c, lambda=lambda, t0=t0, sd = sd, no_intervals=1/sd*10, log_cond = F) 
#   if (log_cond) return(log(b))
#   return(b)
# },'t_cdf')


pcompound1 <- Vectorize(function(c, lambda, t0, sd, t, log_cond = FALSE){
  h <- function(u) dcolony1(u, c, lambda = lambda) *
    pnorm(t - u, mean = t0, sd = sd)
  hi <- t + 8*sd
  hi <- max(hi, 10*sd)
  val <- integral_ahmed(fun = h, 0,  hi, no_intervals=1/sd*10)
  if (log_cond) return(log(val))
  return(val)
}, "t")




fit_energy_density_min1_GA <- function(
    x, trunc = NULL, upper_limit, lower_limit, bin = bin, 
    ga_control = list(popSize = 120, maxiter = 400, run = 80, parallel = T), seed_i,
    ...  # passed to mle2 (e.g., method = "Nelder-Mead")
){
  # deps
  
  
  dots <- list(...)
  if (!"method" %in% names(dots)) dots$method <- "Nelder-Mead"
  
  # keep your start handling & transforms (mle2 uses transformed params)
  
  upper_limit <- upper_limit[c("c", "lambda", "t0", "sd")]
  upper_limit$lambda <- log(upper_limit$lambda)
  upper_limit$sd     <- log(upper_limit$sd)
  upper_limit$c <- log(upper_limit$c)
  print(upper_limit)
  
  lower_limit <- lower_limit[c("c", "lambda", "t0", "sd")]
  lower_limit$lambda <- log(lower_limit$lambda)
  lower_limit$sd     <- log(lower_limit$sd)
  lower_limit$c <- log(lower_limit$c)
  print(lower_limit)
  
  # collapse to counts like you do
  ctdf <- aggregate(1:length(x), by = data.frame(abundance = x), FUN = length)
  names(ctdf)[2] <- "n"
  print(ctdf)
  cat(sprintf("Fitting SAD on %d unique abundance/copy number pairs\n", nrow(ctdf)))
  
  # -----------------------------
  # negative log-likelihood (as in your code)
  # params are in transformed space: c_logit, log_lambda, t0, log_sd
  PENALTY <- 1e50
  log_eps <- -745  # ~ log(1e-323); safe floor to avoid -Inf in doubles
  

  LL_safeguarded <- function(c, lambda, t0, sd){
    # transforms (same as yours)
    c      <- exp(c)
    lambda <- exp(lambda)
    sd     <- exp(sd) + 1
    # print(paste(c, lambda, t0, sd))
    
    # print(paste(c, lambda, t0, sd))
    # quick domain checks
    if (!is.finite(c) || !is.finite(lambda) || !is.finite(t0) || !is.finite(sd) ||
        c <= 0 || lambda <= 0 || sd < 1) return(PENALTY)
    
    # model call: EXPECTS log-probabilities
    lps <- try(
      dcompound1(Observed_T = ctdf$abundance, c = c, lambda = lambda, t0 = t0, sd = sd),
      silent = TRUE
    )
    if (inherits(lps, "try-error") || any(!is.finite(lps))) return(PENALTY)
    
    
    # NLL: weight by counts; DO NOT take log() again
    nll <- -sum(lps * ctdf$n)
    # print(ctdf$abundance)
    if (!is.finite(nll)) nll <- PENALTY
    nll
  }
  # vectorized for GA
  LL_vec <- function(theta) LL_safeguarded(theta[1], theta[2], theta[3], theta[4])
  fitness_fun <- function(theta){
    nll <- LL_vec(theta)
    if (!is.finite(nll)) return(-Inf)
    -nll
  }
  
  # -----------------------------
  # GA bounds (tweak as needed)

  
  lower <- c(lower_limit$c, lower_limit$lambda, lower_limit$t0, lower_limit$sd)
  upper <- c(upper_limit$c, upper_limit$lambda, upper_limit$t0, upper_limit$sd)
  print(lower)
  print(upper)
  # GA maximizes; we want to minimize NLL → fitness = -NLL
  fitness_fun <- function(theta) {
    val <- -LL_vec(theta)
    if (!is.finite(val)) val <- -Inf
    val
  }
  
  if (!is.null(seed_i)){
    set.seed(seed_i)
  }
  
  
  ga_fit <- ga(
    type = "real-valued",
    lower = lower, upper = upper,
    fitness = fitness_fun,
    popSize = ga_control$popSize %||% 120,
    maxiter = ga_control$maxiter %||% 400,
    run     = ga_control$run %||% 80,
    keepBest = ga_control$keepBest %||% F,
    #parallel= ga_control$parallel %||% FALSE,
    parallel= ga_control$parallel %||% FALSE,
    
    optim   = T ,   # keep GA pure; we'll polish with mle2
    monitor  = FALSE
  )
  
  ga_sol <- as.numeric(ga_fit@solution[1, ])
  names(ga_sol) <- c("c", "lambda", "t0", "sd")
  print(ga_sol)
  if (ga_fit@iter < ga_fit@maxiter) {
    message("Converged before reaching maxiter.")
  } else {
    message("May not have converged — reached max iterations.")
  }
  
  
  # 
  # -----------------------------
  # polish with mle2 starting at GA solution (transformed)
  mle_fit <- mle2(
    minuslogl = function(c, lambda, t0, sd) LL_safeguarded(c, lambda, t0, sd),
    start     = as.list(ga_sol),
    control = list(maxit = 1000), 
    method    = dots$method
  )
  
  #unpack your post-processing (same as your code)
  #cat("Raw mle2 @fullcoef:\n"); print(mle_fit@fullcoef['c'])
  
  mle_fit@fullcoef['c']     <- exp(mle_fit@fullcoef['c'])
  mle_fit@fullcoef['lambda']<- exp(mle_fit@fullcoef['lambda'])
  mle_fit@fullcoef['sd']    <- exp(mle_fit@fullcoef['sd']) + 1
  
  #GA
  
  best_theta <- ga_fit@solution[1, ]
  best_NLL   <- LL_vec(best_theta)     # evaluate your minus log-likelihood directly
  k          <- length(best_theta)     # number of parameters
  AIC_GA     <- 2 * k + 2 * best_NLL
  
  # Assuming GA already ran and LL_vec(theta) gives minus log-likelihood
  pop_mat <-  as.data.frame(ga_fit@population)    # matrix: each row = individual, each col = parameter
  k <- ncol(pop_mat)              # number of parameters
  
  colnames(pop_mat) = c('c', 'lambda', 't0', 'sd')
  pop_trans = pop_mat
  
  pop_trans[, 1] <- exp(pop_mat[, 1])   # c
  pop_trans[, 2] <- exp(pop_mat[, 2])                                 # lambda
  # pop_trans[, 3] unchanged                                           # t0
  pop_trans[, 4] <- exp(pop_mat[, 4])  
  
  pop_trans$AIC = 2*k - 2*ga_fit@fitness
  
  
  
  ga_sol['c']     <- exp(ga_sol['c'])
  ga_sol['lambda']<- exp(ga_sol['lambda'])
  ga_sol['sd']    <- exp(ga_sol['sd']) + 1
  cat("GA solution (transformed space):\n"); print(ga_sol)
  #attr(mle_fit, "GA") <- ga_fit   # keep the GA object for inspection
  return(list(ga_fit=ga_fit, parameter_AIC_df = pop_trans, ga_coefs = ga_sol,AIC = AIC(mle_fit), mle_fit = mle_fit))
}



fit_energy_density_min1 <- function(
    x, trunc = NULL, start.value, bin = bin, 
    ...  # passed to mle2 (e.g., method = "Nelder-Mead")
){
  # deps
  
  
  dots <- list(...)
  if (!"method" %in% names(dots)) dots$method <- "Nelder-Mead"
  
  # keep your start handling & transforms (mle2 uses transformed params)
  
  start.value <- start.value[c("c", "lambda", "t0", "sd")]
  start.value$lambda <- log(start.value$lambda)
  start.value$sd     <- log(start.value$sd)
  start.value$c <- log(start.value$c)
  print(start.value)

  # collapse to counts like you do
  ctdf <- aggregate(1:length(x), by = data.frame(abundance = x), FUN = length)
  names(ctdf)[2] <- "n"
  print(ctdf)
  cat(sprintf("Fitting SAD on %d unique abundance/copy number pairs\n", nrow(ctdf)))
  
  # -----------------------------
  # negative log-likelihood (as in your code)
  # params are in transformed space: c_logit, log_lambda, t0, log_sd
  PENALTY <- 1e50
  log_eps <- -745  # ~ log(1e-323); safe floor to avoid -Inf in doubles
  
  
  LL_safeguarded <- function(c, lambda, t0, sd){
    # transforms (same as yours)
    c      <- exp(c)
    lambda <- exp(lambda)
    sd     <- exp(sd) + 1
    # print(paste(c, lambda, t0, sd))
    
    # print(paste(c, lambda, t0, sd))
    # quick domain checks
    if (!is.finite(c) || !is.finite(lambda) || !is.finite(t0) || !is.finite(sd) ||
        c <= 0 || lambda <= 0 || sd < 1) return(PENALTY)
    
    # model call: EXPECTS log-probabilities
    lps <- try(
      dcompound1(Observed_T = ctdf$abundance, c = c, lambda = lambda, t0 = t0, sd = sd),
      silent = TRUE
    )
    if (inherits(lps, "try-error") || any(!is.finite(lps))) return(PENALTY)
    
    
    # NLL: weight by counts; DO NOT take log() again
    nll <- -sum(lps * ctdf$n)
    # print(ctdf$abundance)
    if (!is.finite(nll)) nll <- PENALTY
    nll
  }
 
  mle_fit <- mle2(
    minuslogl = function(c, lambda, t0, sd) LL_safeguarded(c, lambda, t0, sd),
    start     = start.value,
    method    = dots$method
  )
  
  mle_fit@fullcoef['c']     <- exp(mle_fit@fullcoef['c'])
  mle_fit@fullcoef['lambda']<- exp(mle_fit@fullcoef['lambda'])
  mle_fit@fullcoef['sd']    <- exp(mle_fit@fullcoef['sd']) + 1

  return(mle_fit)
}


fit_energy_density_min1_trunc_GA <- function(
    x, trunc = NULL, upper_limit,lower_limit, bin = bin, 
    ga_control = list(popSize = 120, maxiter = 400, run = 80, parallel = T), seed_i,
    ...  # passed to mle2 (e.g., method = "Nelder-Mead")
){
  # deps
  
  
  dots <- list(...)
  if (!"method" %in% names(dots)) dots$method <- "Nelder-Mead"
  
  # keep your start handling & transforms (mle2 uses transformed params)
  
  upper_limit <- upper_limit[c("c", "lambda", "t0", "sd")]
  upper_limit$lambda <- log(upper_limit$lambda)
  upper_limit$sd     <- log(upper_limit$sd)
  upper_limit$c <- log(upper_limit$c)
  print(upper_limit)
  
  lower_limit <- lower_limit[c("c", "lambda", "t0", "sd")]
  lower_limit$lambda <- log(lower_limit$lambda)
  lower_limit$sd     <- log(lower_limit$sd)
  lower_limit$c <- log(lower_limit$c)
  # collapse to counts like you do
  ctdf <- aggregate(1:length(x), by = data.frame(abundance = x), FUN = length)
  t_max = max(ctdf$abundance)
  names(ctdf)[2] <- "n"
  print(ctdf)
  cat(sprintf("Fitting SAD on %d unique abundance/copy number pairs\n", nrow(ctdf)))
  
  # -----------------------------
  # negative log-likelihood (as in your code)
  # params are in transformed space: c_logit, log_lambda, t0, log_sd
  PENALTY <- 1e50
  log_eps <- -745  # ~ log(1e-323); safe floor to avoid -Inf in doubles
  

  LL_safeguarded <- function(c, lambda, t0, sd){
    # transforms (same as yours)
    c      <- exp(c)
    lambda <- exp(lambda)
    sd     <- exp(sd) + 1
    #print(paste(c, lambda, t0, sd))
    # quick domain checks
    if (!is.finite(c) || !is.finite(lambda) || !is.finite(t0) || !is.finite(sd) ||
        c <= 0 || lambda <= 0 || sd < 1) return(PENALTY)
    
    # model call: EXPECTS log-probabilities
    lps <- try(
      dcompound1(Observed_T = ctdf$abundance, c = c, lambda = lambda, t0 = t0, sd = sd) - pcompound1(c=c, lambda=lambda,t0=t0,sd=sd, t = t_max, log_cond = T) ,
      silent = TRUE
    )
    if (inherits(lps, "try-error") || any(!is.finite(lps))) return(PENALTY)
    
    # NLL: weight by counts; DO NOT take log() again
    nll <- -sum(lps * ctdf$n)
    if (!is.finite(nll)) nll <- PENALTY
    nll
  }
  # vectorized for GA
  LL_vec <- function(theta) LL_safeguarded(theta[1], theta[2], theta[3], theta[4])
  fitness_fun <- function(theta){
    nll <- LL_vec(theta)
    if (!is.finite(nll)) return(-Inf)
    -nll
  }
  
  # -----------------------------
  # GA bounds (tweak as needed)
  # c_logit in [-10, 10] ~ c in ~[~0, ~max_c]
  lower_log_c <- lower_limit$c
  upper_log_c <-  upper_limit$c
  
  
  # lambda on log scale (your Brent bounds suggested these)
  lower_log_lambda <- lower_limit$lambda
  upper_log_lambda <-   upper_limit$lambda
  
  # t0 domain: non-negative; cap by observed scale
  # choose [0, max_c] (adjust if you know a tighter prior range)
  lower_t0 <- lower_limit$t0
  upper_t0 <- upper_limit$t0
  
  # sd: sd = exp(logsd) + 1 >= 1; let sd-1 in [1e-6, 100]  -> logsd in [log(1e-6), log(100)]
  lower_logsd <- lower_limit$sd
  upper_logsd <- upper_limit$sd
  
  lower <- c(lower_log_c, lower_log_lambda, lower_t0, lower_logsd)
  upper <- c(upper_log_c, upper_log_lambda, upper_t0, upper_logsd)
  print(lower)
  print(upper)
  # GA maximizes; we want to minimize NLL → fitness = -NLL
  fitness_fun <- function(theta) {
    val <- -LL_vec(theta)
    if (!is.finite(val)) val <- -Inf
    val
  }
  
  if (!is.null(seed_i)){
    set.seed(seed_i)
  }
  
  
  ga_fit <- ga(
    type = "real-valued",
    lower = lower, upper = upper,
    fitness = fitness_fun,
    popSize = ga_control$popSize %||% 120,
    maxiter = ga_control$maxiter %||% 400,
    run     = ga_control$run %||% 80,
    keepBest = ga_control$keepBest %||% F,
    #parallel= ga_control$parallel %||% FALSE,
    parallel= ga_control$parallel %||% FALSE,
    
    optim   = T ,   # keep GA pure; we'll polish with mle2
    monitor  = FALSE
  )
  
  ga_sol <- as.numeric(ga_fit@solution[1, ])
  names(ga_sol) <- c("c", "lambda", "t0", "sd")
  
  if (ga_fit@iter < ga_fit@maxiter) {
    message("Converged before reaching maxiter.")
  } else {
    message("May not have converged — reached max iterations.")
  }
  
  
  # 
  # -----------------------------
  # polish with mle2 starting at GA solution (transformed)
  mle_fit <- mle2(
    minuslogl = function(c, lambda, t0, sd) LL_safeguarded(c, lambda, t0, sd),
    start     = as.list(ga_sol),
    method    = dots$method
  )
  
  #unpack your post-processing (same as your code)
  #cat("Raw mle2 @fullcoef:\n"); print(mle_fit@fullcoef['c'])
  
  mle_fit@fullcoef['c']     <- exp(mle_fit@fullcoef['c'])
  mle_fit@fullcoef['lambda']<- exp(mle_fit@fullcoef['lambda'])
  mle_fit@fullcoef['sd']    <- exp(mle_fit@fullcoef['sd']) + 1
  
  #GA
  
  best_theta <- ga_fit@solution[1, ]
  best_NLL   <- LL_vec(best_theta)     # evaluate your minus log-likelihood directly
  k          <- length(best_theta)     # number of parameters
  AIC_GA     <- 2 * k + 2 * best_NLL
  
  # Assuming GA already ran and LL_vec(theta) gives minus log-likelihood
  pop_mat <-  as.data.frame(ga_fit@population)    # matrix: each row = individual, each col = parameter
  k <- ncol(pop_mat)              # number of parameters
  
  colnames(pop_mat) = c('c', 'lambda', 't0', 'sd')
  pop_trans = pop_mat
  
  pop_trans[, 1] <- exp(pop_mat[, 1])   # c
  pop_trans[, 2] <- exp(pop_mat[, 2])                                 # lambda
  # pop_trans[, 3] unchanged                                           # t0
  pop_trans[, 4] <- exp(pop_mat[, 4])  
  
  pop_trans$AIC = 2*k - 2*ga_fit@fitness
  
  
  
  ga_sol['c']     <- exp(ga_sol['c'])
  ga_sol['lambda']<- exp(ga_sol['lambda'])
  ga_sol['sd']    <- exp(ga_sol['sd']) + 1
  cat("GA solution (transformed space):\n"); print(ga_sol)
  #attr(mle_fit, "GA") <- ga_fit   # keep the GA object for inspection
  return(list(ga_fit=ga_fit, parameter_AIC_df = pop_trans, ga_coefs = ga_sol,AIC = AIC(mle_fit), mle_fit = mle_fit))
}



fit_energy_density_min1_GA_constrained <- function(
    x, trunc = NULL, upper_limit, bin = bin, 
    ga_control = list(popSize = 120, maxiter = 400, run = 80, parallel = T), seed_i,
    ...  # passed to mle2 (e.g., method = "Nelder-Mead")
){
  # deps
  
  
  dots <- list(...)
  if (!"method" %in% names(dots)) dots$method <- "Nelder-Mead"
  
  # keep your start handling & transforms (mle2 uses transformed params)
  
  upper_limit <- upper_limit[c("c", "lambda", "t0", "sd")]
  upper_limit$lambda <- log(upper_limit$lambda)
  upper_limit$sd     <- log(upper_limit$sd)
  upper_limit$c <- log(upper_limit$c)
  print(upper_limit)
  
  # collapse to counts like you do
  ctdf1 <- aggregate(1:length(x[[1]]), by = data.frame(abundance = x[[1]]), FUN = length)
  ctdf2 <- aggregate(1:length(x[[2]]), by = data.frame(abundance = x[[2]]), FUN = length)
  
  names(ctdf1)[2] <- "n"
  names(ctdf2)[2] <- "n"
  print(ctdf1)
  cat(sprintf("Fitting SAD on %d unique abundance/copy number pairs\n", nrow(ctdf1)))
  print(ctdf2)
  cat(sprintf("Fitting SAD on %d unique abundance/copy number pairs\n", nrow(ctdf2)))
  # -----------------------------
  # negative log-likelihood (as in your code)
  # params are in transformed space: c_logit, log_lambda, t0, log_sd
  PENALTY <- 1e50
  log_eps <- -745  # ~ log(1e-323); safe floor to avoid -Inf in doubles
  
  
  LL_safeguarded <- function(c, lambda, t0, sd){
    # transforms (same as yours)
    c      <- exp(c)
    lambda <- exp(lambda)
    sd     <- exp(sd) + 1
    #print(paste(c, lambda, t0, sd))
    # quick domain checks
    if (!is.finite(c) || !is.finite(lambda) || !is.finite(t0) || !is.finite(sd) ||
        c <= 0 || lambda <= 0 || sd < 1) return(PENALTY)
    
    # model call: EXPECTS log-probabilities
    lps1 <- try(
      dcompound1(Observed_T = ctdf1$abundance, c = c, lambda = lambda, t0 = t0, sd = sd),
      silent = TRUE
    )
    if (inherits(lps1, "try-error") || any(!is.finite(lps1))) return(PENALTY)
    
    # NLL: weight by counts; DO NOT take log() again
    nll1 <- -sum(lps1 * ctdf1$n)
    
    
    # model call: EXPECTS log-probabilities
    lps2 <- try(
      dcompound1(Observed_T = ctdf2$abundance, c = c, lambda = lambda, t0 = t0, sd = sd),
      silent = TRUE
    )
    if (inherits(lps2, "try-error") || any(!is.finite(lps2))) return(PENALTY)
    
    # NLL: weight by counts; DO NOT take log() again
    nll2 <- -sum(lps2 * ctdf2$n)
    nll = nll1+nll2
    if (!is.finite(nll)) nll <- PENALTY
    nll
  }
  # vectorized for GA
  LL_vec <- function(theta) LL_safeguarded(theta[1], theta[2], theta[3], theta[4])
  fitness_fun <- function(theta){
    nll <- LL_vec(theta)
    if (!is.finite(nll)) return(-Inf)
    -nll
  }
  
  # -----------------------------
  # GA bounds (tweak as needed)
  # c_logit in [-10, 10] ~ c in ~[~0, ~max_c]
  lower_log_c <- -23
  upper_log_c <-  upper_limit$c
  
  
  # lambda on log scale (your Brent bounds suggested these)
  lower_log_lambda <- -23
  upper_log_lambda <-   upper_limit$lambda
  
  # t0 domain: non-negative; cap by observed scale
  # choose [0, max_c] (adjust if you know a tighter prior range)
  lower_t0 <- 0
  upper_t0 <- upper_limit$t0
  
  # sd: sd = exp(logsd) + 1 >= 1; let sd-1 in [1e-6, 100]  -> logsd in [log(1e-6), log(100)]
  lower_logsd <- log(1e-6)
  upper_logsd <- upper_limit$sd
  
  lower <- c(lower_log_c, lower_log_lambda, lower_t0, lower_logsd)
  upper <- c(upper_log_c, upper_log_lambda, upper_t0, upper_logsd)
  print(lower)
  print(upper)
  # GA maximizes; we want to minimize NLL → fitness = -NLL
  fitness_fun <- function(theta) {
    val <- -LL_vec(theta)
    if (!is.finite(val)) val <- -Inf
    val
  }
  
  if (!is.null(seed_i)){
    set.seed(seed_i)
  }
  
  
  ga_fit <- ga(
    type = "real-valued",
    lower = lower, upper = upper,
    fitness = fitness_fun,
    popSize = ga_control$popSize %||% 120,
    maxiter = ga_control$maxiter %||% 400,
    run     = ga_control$run %||% 80,
    keepBest = ga_control$keepBest %||% F,
    #parallel= ga_control$parallel %||% FALSE,
    parallel= ga_control$parallel %||% FALSE,
    
    optim   = T ,   # keep GA pure; we'll polish with mle2
    monitor  = FALSE
  )
  
  ga_sol <- as.numeric(ga_fit@solution[1, ])
  names(ga_sol) <- c("c", "lambda", "t0", "sd")
  
  if (ga_fit@iter < ga_fit@maxiter) {
    message("Converged before reaching maxiter.")
  } else {
    message("May not have converged — reached max iterations.")
  }
  
  
  # 
  # -----------------------------
  # polish with mle2 starting at GA solution (transformed)
  mle_fit <- mle2(
    minuslogl = function(c, lambda, t0, sd) LL_safeguarded(c, lambda, t0, sd),
    start     = as.list(ga_sol),
    method    = dots$method
  )
  
  #unpack your post-processing (same as your code)
  #cat("Raw mle2 @fullcoef:\n"); print(mle_fit@fullcoef['c'])
  
  mle_fit@fullcoef['c']     <- exp(mle_fit@fullcoef['c'])
  mle_fit@fullcoef['lambda']<- exp(mle_fit@fullcoef['lambda'])
  mle_fit@fullcoef['sd']    <- exp(mle_fit@fullcoef['sd']) + 1
  
  #GA
  
  best_theta <- ga_fit@solution[1, ]
  best_NLL   <- LL_vec(best_theta)     # evaluate your minus log-likelihood directly
  k          <- length(best_theta)     # number of parameters
  AIC_GA     <- 2 * k + 2 * best_NLL
  
  # Assuming GA already ran and LL_vec(theta) gives minus log-likelihood
  pop_mat <-  as.data.frame(ga_fit@population)    # matrix: each row = individual, each col = parameter
  k <- ncol(pop_mat)              # number of parameters
  
  colnames(pop_mat) = c('c', 'lambda', 't0', 'sd')
  pop_trans = pop_mat
  
  pop_trans[, 1] <- exp(pop_mat[, 1])   # c
  pop_trans[, 2] <- exp(pop_mat[, 2])                                 # lambda
  # pop_trans[, 3] unchanged                                           # t0
  pop_trans[, 4] <- exp(pop_mat[, 4])  
  
  pop_trans$AIC = 2*k - 2*ga_fit@fitness
  
  
  
  ga_sol['c']     <- exp(ga_sol['c'])
  ga_sol['lambda']<- exp(ga_sol['lambda'])
  ga_sol['sd']    <- exp(ga_sol['sd']) + 1
  cat("GA solution (transformed space):\n"); print(ga_sol)
  #attr(mle_fit, "GA") <- ga_fit   # keep the GA object for inspection
  return(list(ga_fit=ga_fit, parameter_AIC_df = pop_trans, ga_coefs = ga_sol,AIC = AIC(mle_fit), mle_fit = mle_fit))
}



## Constrained-c joint fit for two conditions
## - x: list of length 2; x[[1]] and x[[2]] are the Observed_T vectors for the two conditions
## - start.value must name: c, lambda1, t01, sd1, lambda2, t02, sd2
## - c is shared across conditions; (lambda,t0,sd) are condition-specific
## - Uses same parameter transforms as your original code:
##     lambda := exp(lambda), sd := exp(sd)+1, c := my_logit_transform(max_abundance, c)
## - Requires: my_logit_transform(), my_inverse_logit_transform(), dcompound1(), max_abundance in scope
##
fit_energy_density_min1_constrained_c <- function(x, trunc=NULL, start.value, bin = NULL, ...) {
  dots <- list(...)
  
  if (!is.list(x) || length(x) != 2) stop("x must be a list of length 2: x[[1]] (cond1), x[[2]] (cond2).")
  
  # Aggregate counts per condition (like your original ctdf)
  ctdf1 <- aggregate(seq_along(x[[1]]), by = data.frame(abundance = x[[1]]), FUN = length)
  ctdf2 <- aggregate(seq_along(x[[2]]), by = data.frame(abundance = x[[2]]), FUN = length)
  
  if (!"method" %in% names(dots)) dots$method <- "Nelder-Mead"
  
  # Reorder and transform starts
  start.value <- start.value[c("c","lambda1","t01","sd1","lambda2","t02","sd2")]
  # map starting values like your original
  start.value$lambda1 <- log(start.value$lambda1)
  start.value$lambda2 <- log(start.value$lambda2)
  start.value$sd1     <- log(start.value$sd1)
  start.value$sd2     <- log(start.value$sd2)
  # For c you used the 'inverse' transform at start, then transformed inside LL each call
  start.value$c <- my_inverse_logit_transform(m = max_abundance, s = start.value$c)
  
  message("Start values (transformed):")
  print(start.value)
  message(sprintf("Fitting on %d (cond1) + %d (cond2) unique abundance rows",
                  nrow(ctdf1), nrow(ctdf2)))
  
  # Negative log-likelihood combining both conditions with shared c
  LL <- function(c, lambda1, t01, sd1, lambda2, t02, sd2) {
    # back-transform params
    c      <- my_logit_transform(max_abundance, c)
    lambda1 <- exp(lambda1)
    lambda2 <- exp(lambda2)
    sd1     <- exp(sd1) + 1
    sd2     <- exp(sd2) + 1
    
    # per-condition contributions, weighted by counts
    # dcompound1 should return the pdf at Observed_T for given params
    pps1 <- dcompound1(Observed_T = ctdf1$abundance, c = c,
                       lambda = lambda1, t0 = t01, sd = sd1) * ctdf1$x
    pps2 <- dcompound1(Observed_T = ctdf2$abundance, c = c,
                       lambda = lambda2, t0 = t02, sd = sd2) * ctdf2$x
    
    nll <- -(sum(pps1) + sum(pps2))
    
    if (is.infinite(nll) || !is.finite(nll)) nll <- .Machine$double.xmax
    nll
  }
  
  # Fit with mle2
  args <- c(list(minuslogl = LL, start = start.value, data = list(x = x), control = list(maxit = 1000)), dots)
  result <- do.call(bbmle::mle2, args)
  
  # Back-transform reported coefficients for convenience (like your original)
  fc <- result@fullcoef
  # c back to original scale
  fc["c"]      <- my_logit_transform(max_abundance, fc["c"])
  # lambdas
  fc["lambda1"] <- exp(fc["lambda1"])
  fc["lambda2"] <- exp(fc["lambda2"])
  # sds
  fc["sd1"]     <- exp(fc["sd1"]) + 1
  fc["sd2"]     <- exp(fc["sd2"]) + 1
  
  message("Transformed coefficients (for reporting):")
  print(fc[c("c","lambda1","t01","sd1","lambda2","t02","sd2")])
  
  # Return fitted object; attach a convenience field with back-transformed coefs
  result@fullcoef <- fc
  return(result)
}
