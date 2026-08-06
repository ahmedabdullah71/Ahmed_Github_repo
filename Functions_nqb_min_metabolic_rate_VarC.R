# Number of simulations
library(tidyr)
library(dplyr)
### Load packages ###
library(egg)
library(ggplot2)
#change pnorm to dnorom
library(sads)
library("FRACTION")
#cdf diff t0 is normally distributed
#cdf diff
#probability is integration
#integration of density function
library(pracma)
library(grid)
library(gridExtra)

library(lamW)

library(scales)

source("/home/ahmed/nqb_R_script/quadgk51_23.R") #put file path of quadgk51")



data_load <- function(k){
  df1 = read.csv('/home/ahmed/Dropbox/Ahmed Abdullah/AMR/Balaban_data/lagscan1.csv')
  df2 = read.csv('/home/ahmed/Dropbox/Ahmed Abdullah/AMR/Balaban_data/lagscan2.csv')
  df3 = read.csv('/home/ahmed/Dropbox/Ahmed Abdullah/AMR/Balaban_data/lagscan3.csv')
  df4 = read.csv('/home/ahmed/Dropbox/Ahmed Abdullah/AMR/Balaban_data/lagscan4.csv')
  
  df_list = list(df1,df2,df3,df4)
  
  df = df_list[[k]]
  return(df)
}





getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}



uniformly_divide_0ToP <- function(NN,P){
  seq(P,P/NN, length = NN)-P/(2*NN)
}

# x is appearance time data, n is number of divisions, bins
counts <-  function(x,n) {
  xs = cut(x, breaks=seq(min(x),max(x)+1, length.out = n+1), right = FALSE)
  ys = as.vector(table(xs))
  return(ys)
}

uniformly_divide_0To1 <- function(NN){
  seq(1,1/NN, length = NN)-1/(2*NN)
}

#number has to be positive
binned_appearance_fun = function(myVec, bin=10){
  #myMode = getmode(myVec)
  myMode = 0
  myVec1 = c(-myMode, myVec) #add point
  n = (max(myVec1)-min(myVec1)) / bin 
  
  #appearance_count=counts(x,n)
  appearance_count = counts(myVec1, n)
  
  #interval= (max(x)-min(x))/n - Don't need to use these variables
  #intervalN = (max(myVec) - min(myVec))/n 
  #interval_density=appearance_count/(sum(appearance_count)*interval)
  
  #Ap_time=sort(min(x)+uniformly_divide_0To1(n)*(max(x)-min(x)), decreasing = F)
  Ap_time =sort(min(myVec1)+uniformly_divide_0To1(n)*(max(myVec1)-min(myVec1)), decreasing = F)
  
  # delete point
  Ap_time = Ap_time[-1]
  appearance_count = appearance_count[-1]
  
  # Plot if needed
  plot(Ap_time, appearance_count, type = "l")
  
  # Create CSV file 
  #filename = file.path(directory_path, 'lagscan.csv' )
  lagscan.df <- cbind(round(Ap_time), appearance_count) 
  colnames(lagscan.df) <- c("TimeAxis", "TotalDist")
  return(as.data.frame(lagscan.df))
  
}



### Fitting ###
freq_gen <- function(df, bin = bin , last_min =0){
  
  x=df$TimeAxis
  y=df$TotalDist
  # x=x[-1]
  # y=y[-1]
  TotalBact = sum(y)
  
  y = y/(TotalBact*bin)
  
  
  if (last_min == 0){
    Y = list(x=x,y=y)
    return(Y) 
  } else{
    pre_last_min = df$TimeAxis[length(df$TimeAxis)]
    
    x_tail = seq(pre_last_min, last_min, by=bin)
    y_tail = rep(0,length(x_tail))
    Y = list(x=c(x,x_tail), y = c(y, y_tail))
    return(Y)
  }
  
}

freqToAbundace = function(df){
  freqToAbundace_zero = function(df){
    if (df['TotalDist']==0){
      return(0)
    }
    else{
      n_rep = rep(df['TimeAxis'],df['TotalDist'])
      return(n_rep)
    }
    
  }
  v = apply(df,1, freqToAbundace_zero)
  v = unlist(v)
  names(v) = NULL
  time_abundance = v[v>0]  #truncate at 0
  time_abundance = time_abundance[order(time_abundance)] # order # order function produce the subscripts
  return(time_abundance)
}

time_abundance_for_datasets = function(lst){
  
  time_abundance = freqToAbundace(df)
  return(list_time_abundance)
}

# PDF of division time T = c/r
dcolony1_varC <- function(t, lambda, mu) {
  if (any(t < 0)) stop("t must be non-negative")
  lambda * mu * exp(-mu * t) * (1 / (lambda + mu * t) + 1 / (lambda + mu * t)^2)
}

pcolony1_varC <- function(t, lambda, mu) {
  if (any(t < 0)) stop("t must be non-negative")
  1 - (lambda * exp(-mu * t)) / (lambda + mu * t)
}



# qcolony1_varC <- function(pp, lambda, mu) {
#   t <- (-lambda + lambertW0( - (exp(lambda) * lambda) / (-1 + pp) )) / mu
#   return(t)
# }

# qcolony1_varC <- Vectorize(function(pp, lambda, mu) {
#   t <- uniroot(f = function(t) pcolony1_varC(t, lambda, mu) -pp, lower = 0, upper = 10^6)$root
#   return(t)
# }, 'pp')

# qcolony1_varC <- function(pp, lambda, mu) {
#   if (any(!is.finite(lambda)) || any(!is.finite(mu))) {
#     stop("lambda and mu must be finite numeric values.")
#   }
#   if (any(mu < 0)) {
#     stop("mu must be >= 0.")
#   }
#   
#   if (!requireNamespace("lamW", quietly = TRUE)) {
#     stop("Package 'lamW' is required. Install it with: install.packages('lamW')")
#   }
#   
#   z <- exp(lambda) * lambda/(1-pp)
#   w <- lamW::lambertW0(z)  # principal branch (ProductLog)
#   (-lambda + w) / mu
# }

lambertW_large_z <- function(lambda, pp) {
  if (all(pp >= 1)) stop("pp must be < 1")
  if (all(lambda < 0)) stop("lambda must be >= 0")
  
  # log(z) without overflow: log( exp(lambda) * lambda/(1-pp) )
  L1 <- lambda + log(lambda) - log1p(-pp)
  L2 <- log(L1)
  
  # asymptotic Lambert W (principal branch)
  L1 - L2 + L2 / L1
}

qcolony1_varC <- function(pp, lambda, mu, lambda_cutoff = 650) {
  if (any(!is.finite(pp)) || any(!is.finite(lambda)) || any(!is.finite(mu))) {
    stop("pp, lambda, mu must be finite numeric values.")
  }
  if (any(pp <= 0) || any(pp >= 1)) {
    stop("pp must be in (0, 1).")
  }
  if (any(lambda < 0)) {
    stop("lambda must be >= 0.")
  }
  if (any(mu < 0)) {
    stop("mu must be >= 0.")
  }
  
  # Decide per-element whether to use exact W or asymptotic W
  use_asymp <- lambda >= lambda_cutoff
  
  
  
  if (any(!use_asymp)) {
    if (!requireNamespace("lamW", quietly = TRUE)) {
      stop("Package 'lamW' is required for small/moderate lambda. Install it with: install.packages('lamW')")
    }
    idx <- which(!use_asymp)
    z <- exp(lambda) * lambda / (1 - pp)
    w <- lamW::lambertW0(z)
  }
  
  if (any(use_asymp)) {
    idx <- which(use_asymp)
    w <- lambertW_large_z(lambda, pp)
  }
  
  f_t = (-lambda + w) / mu
  f_t[f_t < 0] <- 0
  
  return(f_t)
}

# 0 energy
compound_fun0 = function(t, Observed_T, lambda, mu, t0, sd){
  h_t = dcolony1_varC(t,lambda=lambda, mu = mu)*dnorm(Observed_T-t,t0,sd=sd)
  #print(h_t)
  return(h_t)
}

dcompound = Vectorize(function(lambda, mu,t0, Observed_T, sd, log_cond = TRUE){
  if (sd <1){
    no_intervals = 1/sd*10
  } else{
    no_intervals =8
  }
  b=integral_ahmed(fun = compound_fun0, 0, Observed_T, Observed_T=Observed_T, sd = sd, lambda=lambda,mu = mu,  t0=t0,no_intervals=no_intervals) 
  if (log_cond) return(log(b))
  return(b)
},'Observed_T')





fit_energy_density_varC_GA <- function(
    x, trunc = NULL, upper_limit, lower_limit= list(mu = 1e-7, lambda = 1e-7, t0= 0, sd = 1e-6), bin = bin, 
    ga_control = list(popSize = 120, maxiter = 400, run = 80, parallel = T), seed_i,
    ...  # passed to mle2 (e.g., method = "Nelder-Mead")
){
  # deps
  
  
  dots <- list(...)
  if (!"method" %in% names(dots)) dots$method <- "Nelder-Mead"
  
  # keep your start handling & transforms (mle2 uses transformed params)
  
  upper_limit <- upper_limit[c("mu", "lambda", "t0", "sd")]
  upper_limit$lambda <- log(upper_limit$lambda)
  upper_limit$sd     <- log(upper_limit$sd)
  upper_limit$mu <- log(upper_limit$mu)
  print(upper_limit)
  
  lower_limit <- lower_limit[c("mu", "lambda", "t0", "sd")]
  lower_limit$lambda <- log(lower_limit$lambda)
  lower_limit$sd     <- log(lower_limit$sd)
  lower_limit$mu <- log(lower_limit$mu)
  
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
  
  # LL_safeguarded <- function(c, lambda, t0, sd){
  #   # transforms (same as yours)
  #   c      <- max_c * exp(c) / (1 + exp(c))
  #   lambda <- exp(lambda)
  #   sd     <- exp(sd) + 1
  #   #print(paste(c, lambda, t0, sd))
  #   # quick domain checks
  #   if (!is.finite(c) || !is.finite(lambda) || !is.finite(t0) || !is.finite(sd) ||
  #       c <= 0 || lambda <= 0 || sd < 1) return(PENALTY)
  #   
  #   # model call: EXPECTS log-probabilities
  #   lps <- try(
  #     dcompound1(Observed_T = ctdf$abundance, c = c, lambda = lambda, t0 = t0, sd = sd),
  #     silent = TRUE
  #   )
  #   if (inherits(lps, "try-error") || any(!is.finite(lps))) return(PENALTY)
  #   
  #   # lps are log-probs: handle -Inf / very small values robustly
  #   bad <- (lps <= log_eps) | !is.finite(lps)
  #   bad_frac <- mean(bad)
  #   if (bad_frac > 0.05) return(PENALTY)   # too many zeros -> implausible region
  #   lps[bad] <- log_eps                    # soft floor a few numerical zeros
  #   
  #   # NLL: weight by counts; DO NOT take log() again
  #   nll <- -sum(lps * ctdf$n)
  #   if (!is.finite(nll)) nll <- PENALTY
  #   nll
  # }
  LL_safeguarded <- function(mu, lambda, t0, sd){
    # transforms (same as yours)
    
    mu      <- exp(mu)
    lambda <- exp(lambda)
    sd     <- exp(sd) + 1
    #print(paste(mu, lambda, t0, sd))
    # quick domain checks
    if (!is.finite(mu) || !is.finite(lambda) || !is.finite(t0) || !is.finite(sd) ||
        mu <= 0 || lambda <= 0 || sd < 1) return(PENALTY)
    # 
    # model call: EXPECTS log-probabilities
    lps <- try(
      dcompound(Observed_T = ctdf$abundance, mu =mu, lambda = lambda, t0 = t0, sd = sd),
      silent = TRUE
    )
    if (inherits(lps, "try-error") || any(!is.finite(lps))) return(PENALTY)
    # 
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
  

  # sd: sd = exp(logsd) + 1 >= 1; let sd-1 in [1e-6, 100]  -> logsd in [log(1e-6), log(100)]
  
  
  
  lower <- c(lower_limit$mu, lower_limit$lambda, lower_limit$t0, lower_limit$sd)
  upper <- c(upper_limit$mu, upper_limit$lambda, upper_limit$t0, upper_limit$sd)
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
    parallel= T,
    optim   = T,    # keep GA pure; we'll polish with mle2
    monitor  = FALSE
  )
  
  ga_sol <- as.numeric(ga_fit@solution[1, ])
  names(ga_sol) <- c("mu", "lambda", "t0", "sd")
  
  # 
  # -----------------------------
  # polish with mle2 starting at GA solution (transformed)
  mle_fit <- mle2(
    minuslogl = function(mu, lambda, t0, sd) LL_safeguarded(mu, lambda, t0, sd),
    start     = as.list(ga_sol),
    method    = dots$method
  )
  
  #unpack your post-processing (same as your code)
  cat("Raw mle2 @fullcoef:\n"); print(mle_fit@fullcoef['mu'])
  
  mle_fit@fullcoef['mu']     <- exp(mle_fit@fullcoef['mu'])
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
  
  colnames(pop_mat) = c('mu', 'lambda', 't0', 'sd')
  pop_trans = pop_mat
  
  pop_trans[, 1] <- exp(pop_mat[, 1])
  pop_trans[, 2] <- exp(pop_mat[, 2])                                 # lambda
  # pop_trans[, 3] unchanged                                           # t0
  pop_trans[, 4] <- exp(pop_mat[, 4])  
  
  pop_trans$AIC = 2*k - 2*ga_fit@fitness
  
  
  
  ga_sol['mu']     <- exp(ga_sol['mu'])
  ga_sol['lambda']<- exp(ga_sol['lambda'])
  ga_sol['sd']    <- exp(ga_sol['sd']) + 1
  cat("GA solution (transformed space):\n"); print(ga_sol)
  #attr(mle_fit, "GA") <- ga_fit   # keep the GA object for inspection
  return(list(ga_fit=ga_fit, parameter_AIC_df = pop_trans, ga_coefs = ga_sol,AIC = AIC(mle_fit), mle_fit = mle_fit))
}



get_appearance_time = function(directory_path){
  # Get a list of all CSV files in the directory (ColoniesAppearance files)
  csv_files <- list.files(directory_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Use lapply to read all CSV files into a list of data frames
  data_frames <- lapply(csv_files, read.csv)
  
  # Only look at appearance time (disregard first column)
  ColoniesAppearance <- lapply(data_frames, function(df) df[[2]])
  GrowthTimes <- lapply(data_frames, function(df) df[[3]])
  # combine into one df
  combined_ap_df <- data.frame(TimeAxis = unlist(ColoniesAppearance),GrowthTimes = unlist(GrowthTimes))
  # 
  if (!is.null(GrowthTime_upper)){
    combined_ap_df = combined_ap_df[combined_ap_df$GrowthTimes<GrowthTime_upper,] 
  }
  
  combined_ap_df = combined_ap_df[combined_ap_df$GrowthTimes>GrowthTime_threshold,]
  
  col_ap <- combined_ap_df$TimeAxis
  # unlist the file 
  myVec = unlist(col_ap)
  return(myVec)
}