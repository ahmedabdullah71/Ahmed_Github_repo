library(microSAD)
library(bbmle) 
library(ggplot2)

exp_dist = function(t, k){
  k * exp(-k * t)
}

exp_free = function(t, A1, k){
  A1 * exp(-k * t)
}
dcolony1 <- function(t, c,lambda, log_cond = FALSE){
  
  #f_t=c/t^2*lambda*exp(lambda)*exp(-lambda*c/t)
  f_t=log(c)+log(lambda)+lambda-2*log(t)-(lambda*c)/(t)
  f_t[t<=0|t>c] = -Inf
  if(log_cond) return(f_t)
  return(exp(f_t))
}

pcolony1 <- function(t, c, lambda, log_cond = FALSE) {
  # initialize output
  f_t <- lambda - lambda * c / t   # works elementwise if t is a vector
  #print(f_t)
  # enforce truncation: for t > c, set to -Inf (log case) or 0 (prob case)
  if (log_cond) {
    f_t <- ifelse(t > c, 0, f_t)    # log probability
    return(f_t)
  } else {
    out <- exp(f_t)
    out <- ifelse(t > c, 1, out)       # probability
    return(out)
  }
}

qcolony1 = function(pp, c, lambda){
  
  qq=c*lambda/(lambda-log(pp))
  return(qq)
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


lag_fun <- function(t, A1, k, A2, b, tau0 = 93) {
  ifelse(t < tau0,
         A1 * exp(-k * t),
         A2 * t^-b)
}

lag_fun_norm <- function(t, A1, k, A2, b, tau0=93) {
  Z <- (A1 / k) * (1 - exp(-k * tau0)) + (A2 / (-(-b+1))) * tau0^(-b+1)
  f <- ifelse(t < tau0,
              A1 * exp(-k * t),
              A2 * t^-b)
  f / Z
}

lag_fun_cdf <- Vectorize(function(t, A1, k, A2, b, tau0 = 93) {
  f = pracma::integral(lag_fun_norm, 0, t , A1=A1, k=k, A2=A2, b=b, tau0=tau0)
  return(f)
},'t')

# 
# #Kill Curve
# lag_cdf_fun <- Vectorize(function(tau, A1, k, A2, b, tau0, d=100){
#   if (tau < d){
#     return(1)
#   }
#   f = 1 - pracma::integral(lag_fun, 0, tau - d , A1=A1, k=k, A2=A2, b=b, tau0=tau0)
#   return(f)
# },'tau')



# k = 0.063
# b = 2.1
# tau0 = 93
# A1 = 0.0622
# A2 = 2.42



fit_maxEnt_lagTime_cdf <- function(df , start.value, debug_flag = F){
  
  t = df$t
  counts = df$count
  width = df$width
  
  make_binned_nll <- function(t, counts, width, 
                              open_left = FALSE, open_right = FALSE,
                              eps = 1e-16 , min_c = 10) {
    stopifnot(length(t) == length(counts))
    if (length(width) == 1) width <- rep(width, length(t))
    stopifnot(length(width) == length(t))
    
    
    # compute bin edges from midpoints
    a <- c(0, head(t, -1))
    b <- t
    if (open_left)  a[1]           <- -Inf
    if (open_right) b[length(b)]   <-  Inf
    
    print(paste('a', a, 'b' , b , sep=':'))
    # Return minus log-likelihood function parameterized by (c, lambda)
    function(c, lambda) {
      c=max(t)+exp(c)
      lambda = exp(lambda)
      
      Fb <- pcolony1(b, c, lambda, log_cond = FALSE)
      Fa <- pcolony1(a, c, lambda, log_cond = FALSE)
      Z <- pcolony1(max(b), c, lambda, log_cond = FALSE)
      F_z <- (Fb - Fa)/Z
      # Bin probabilities
      p  <- pmax(F_z, eps)   # guard underflow/zeros
      #p <- pmin(p, Inf)
      #print(paste('p ',p))
      # Multinomial log-lik up to constant: sum n_i log p_i
      # Return NEGATIVE log-likelihood for optim/mle2
      pp = -sum(counts * log(p))
      if (debug_flag){
        print(paste('c ', c, 'lambda', lambda, sep=' '))
        print(paste('pp ', pp))
      }
      
      return(pp)
      
    }
  }
  
  minuslogl <- make_binned_nll(
    t       = df$t,
    counts     = df$count,
    width      = df$width,
    open_left  = FALSE,      # set as needed
    open_right = FALSE
  )
  
  # MLE with simple box constraints (adjust starts/bounds)
  fit <- mle2(
    minuslogl = minuslogl,
    start     = start.value,
    method    = "L-BFGS-B",
    
  )
  
  fit@fullcoef["lambda"] = exp(fit@fullcoef["lambda"])
  fit@fullcoef["c"] = max(t)+exp(fit@fullcoef["c"])
  
  return(fit)
  
}

fit_maxEnt_lagTime <- function(df , start.value, debug_flag = F){
  
  mids = df$mid
  counts = df$count
  width = df$width
  
  #added
  a <- mids - width/2
  b <- mids + width/2
  
  
  make_binned_nll <- function(mids, counts, width, 
                              open_left = FALSE, open_right = FALSE,
                              eps = 1e-16 , min_c = 10) {
    stopifnot(length(mids) == length(counts))
    if (length(width) == 1) width <- rep(width, length(mids))
    stopifnot(length(width) == length(mids))
    
    
    # compute bin edges from midpoints
    a <- mids - width/2
    b <- mids + width/2
    if (open_left)  a[1]           <- -Inf
    if (open_right) b[length(b)]   <-  Inf
    
    print(paste('a', a, 'b' , b , sep=':'))
    # Return minus log-likelihood function parameterized by (c, lambda)
    function(c, lambda) {
      c=max(b)+exp(c)
      lambda = exp(lambda)
      
      Fb <- pcolony1(b, c, lambda, log_cond = FALSE)
      Fa <- pcolony1(a, c, lambda, log_cond = FALSE)
      Z <- pcolony1(max(b), c, lambda, log_cond = FALSE)
      F_z <- (Fb - Fa)/Z
      # Bin probabilities
      p  <- pmax(F_z, eps)   # guard underflow/zeros
      #p <- pmin(p, Inf)
      #print(paste('p ',p))
      # Multinomial log-lik up to constant: sum n_i log p_i
      # Return NEGATIVE log-likelihood for optim/mle2
      pp = -sum(counts * log(p))
      if (debug_flag){
        print(paste('c ', c, 'lambda', lambda, sep=' '))
        print(paste('pp ', pp))
      }
      
      return(pp)
      
    }
  }
  
  minuslogl <- make_binned_nll(
    mids       = df$mid,
    counts     = df$count,
    width      = df$width,
    open_left  = FALSE,      # set as needed
    open_right = FALSE
  )
  
  # MLE with simple box constraints (adjust starts/bounds)
  fit <- mle2(
    minuslogl = minuslogl,
    start     = start.value,
    method    = "L-BFGS-B",
    
  )
  
  fit@fullcoef["lambda"] = exp(fit@fullcoef["lambda"])
  fit@fullcoef["c"] = max(b)+exp(fit@fullcoef["c"])
  
  return(fit)
  
}


# This is for fitting vector data c(100, 101, 104.2, ...)
fit_maxEnt_lagTime_dfun <- function(x , start.value, trunc = NULL, debug_flag = F){
 
  start.value =  lapply(start.value , log)
  t_max =max(x)

  make_nll <- function(x, eps = 1e-16 , min_c = 10) {
    # Return minus log-likelihood function parameterized by (c, lambda)
    
    
    function(c, lambda) {
      c=t_max+exp(c)
      lambda = exp(lambda)
      
      if (!is.null(trunc)){
        F_tmax = pcolony1(t_max, c, lambda, log_cond = T)
      } else{
        F_tmax = 0
      }
      
      
      Fx <- dcolony1(x, c, lambda, log_cond = T)
      Fz <- Fx - F_tmax
      # Bin probabilities
      p  <- pmax(Fz, log(eps))   # guard underflow/zeros
      #p <- pmin(p, Inf)
      #print(paste('p ',p))
      # Multinomial log-lik up to constant: sum n_i log p_i
      # Return NEGATIVE log-likelihood for optim/mle2
      pp = -sum(p)
      if (debug_flag){
        print(paste('c ', c, 'lambda', lambda, sep=' '))
        print(paste('pp ', pp))
      }
      
      return(pp)
      
    }
  }
  
  minuslogl <- make_nll(x=x)
  
  # MLE with simple box constraints (adjust starts/bounds)
  fit <- mle2(
    minuslogl = minuslogl,
    start     = start.value,
    method    = "L-BFGS-B",
    
  )
  
  fit@fullcoef["lambda"] = exp(fit@fullcoef["lambda"])
  fit@fullcoef["c"] = t_max+exp(fit@fullcoef["c"])
  
  return(fit)
  
}


fit_SimSek_lagTime_cdf <- function(df , start.value, debug_flag = F){
  
  t = df$t
  counts = df$count
  width = df$width
  
  make_binned_nll <- function(t, counts, width, 
                              open_left = FALSE, open_right = FALSE,
                              eps = 1e-16 , min_c = 10) {
    stopifnot(length(t) == length(counts))
    if (length(width) == 1) width <- rep(width, length(t))
    stopifnot(length(width) == length(t))
    
    
    # compute bin edges from midpoints
    aa <- c(0, head(t, -1))
    bb <- t
    print(paste('aa', aa, 'bb', bb, sep = ':'))
    if (open_left)  aa[1]           <- -Inf
    if (open_right) bb[length(bb)]   <-  Inf
    
    # Return minus log-likelihood function parameterized by (c, lambda)
    function(A1, k, A2, b) {
      k = exp(k)
      b = exp(b)
      A1 = exp(A1)
      A2 = exp(A2)
      
      
      Fb <- lag_fun_cdf(bb, A1, k, A2, b )
      Fa <- lag_fun_cdf(aa, A1, k, A2, b)
      
      Z <- lag_fun_cdf(max(bb), A1, k, A2, b)
      F_z <- (Fb - Fa)/Z
      # Bin probabilities
      p  <- pmax(F_z, eps)   # guard underflow/zeros
      #p <- pmin(p, Inf)
      #print(paste('p ',p))
      # Multinomial log-lik up to constant: sum n_i log p_i
      # Return NEGATIVE log-likelihood for optim/mle2
      pp = -sum(counts * log(p))
      if (debug_flag){
        #print(paste('Fa ',Fa,'Fb ',Fb, 'z', Z, 'A1 ', A1, 'k', k, 'A2', A2, 'b', b, sep=' '))
        print(paste('A1 ', A1, 'k', k, 'A2', A2, 'b', b, sep=' '))
        print(paste('pp ', pp))
      }
      
      return(pp)
      
    }
  }
  
  minuslogl <- make_binned_nll(
    t       = df$t,
    counts     = df$count,
    width      = df$width,
    open_left  = FALSE,      # set as needed
    open_right = FALSE
  )
  
  # MLE with simple box constraints (adjust starts/bounds)
  fit <- mle2(
    minuslogl = minuslogl,
    start     = start.value,
    method    = "L-BFGS-B",
    
  )
  
  fit@fullcoef["k"] = exp(fit@fullcoef["k"])
  fit@fullcoef["b"] = exp(fit@fullcoef["b"])
  fit@fullcoef["A1"] = exp(fit@fullcoef["A1"])
  fit@fullcoef["A2"] = exp(fit@fullcoef["A2"])
  
  return(fit)
  
}


fit_SimSek_lagTime <- function(df , start.value, debug_flag = F){
  
  mids = df$mid
  counts = df$count
  width = df$width
  
  make_binned_nll <- function(mids, counts, width, 
                              open_left = FALSE, open_right = FALSE,
                              eps = 1e-16 , min_c = 10) {
    stopifnot(length(mids) == length(counts))
    if (length(width) == 1) width <- rep(width, length(mids))
    stopifnot(length(width) == length(mids))
    
    
    # compute bin edges from midpoints
    aa <- mids - width/2
    bb <- mids + width/2
    if (open_left)  aa[1]           <- -Inf
    if (open_right) bb[length(bb)]   <-  Inf
    
    print(paste('aa', aa, 'bb', bb, sep = ':'))
    
    # Return minus log-likelihood function parameterized by (c, lambda)
    function(A1, k, A2, b) {
      k = exp(k)
      b = exp(b)
      A1 = exp(A1)
      A2 = exp(A2)
      
      
      Fb <- lag_fun_cdf(bb, A1, k, A2, b )
      Fa <- lag_fun_cdf(aa, A1, k, A2, b)
      
      Z <- lag_fun_cdf(max(bb), A1, k, A2, b)
      F_z <- (Fb - Fa)/Z
      # Bin probabilities
      p  <- pmax(F_z, eps)   # guard underflow/zeros
      #p <- pmin(p, Inf)
      #print(paste('p ',p))
      # Multinomial log-lik up to constant: sum n_i log p_i
      # Return NEGATIVE log-likelihood for optim/mle2
      pp = -sum(counts * log(p))
      if (debug_flag){
        #print(paste('Fa ',Fa,'Fb ',Fb, 'z', Z, 'A1 ', A1, 'k', k, 'A2', A2, 'b', b, sep=' '))
        print(paste('A1 ', A1, 'k', k, 'A2', A2, 'b', b, sep=' '))
        print(paste('pp ', pp))
      }
      
      return(pp)
      
    }
  }
  
  minuslogl <- make_binned_nll(
    mids       = df$mid,
    counts     = df$count,
    width      = df$width,
    open_left  = FALSE,      # set as needed
    open_right = FALSE
  )
  
  # MLE with simple box constraints (adjust starts/bounds)
  fit <- mle2(
    minuslogl = minuslogl,
    start     = start.value,
    method    = "L-BFGS-B",
    
  )
  
  fit@fullcoef["k"] = exp(fit@fullcoef["k"])
  fit@fullcoef["b"] = exp(fit@fullcoef["b"])
  fit@fullcoef["A1"] = exp(fit@fullcoef["A1"])
  fit@fullcoef["A2"] = exp(fit@fullcoef["A2"])
  
  return(fit)
  
}


fit_PowerLaw_lagTime <- function(df , start.value, debug_flag = F){
  
  mids = df$mid
  counts = df$count
  width = df$width
  print(width)
  make_binned_nll <- function(mids, counts, width, 
                              open_left = FALSE, open_right = FALSE,
                              eps = 1e-16 , min_c = 10) {
    stopifnot(length(mids) == length(counts))
    if (length(width) == 1) width <- rep(width, length(mids))
    stopifnot(length(width) == length(mids))
    
    
    # compute bin edges from midpoints
    a <- mids - width/2
    b <- mids + width/2
    if (open_left)  a[1]           <- -Inf
    if (open_right) b[length(b)]   <-  Inf
    
    # Return minus log-likelihood function parameterized by (c, lambda)
    function(s) {
      s=1+exp(s)
      # print(s)
      
      Fb <- ppower(b, s)
      Fa <- ppower(a, s)
      Z <- ppower(max(b), s)
      F_z <- (Fb - Fa)/Z
      # Bin probabilities
      p  <- pmax(F_z, eps)   # guard underflow/zeros
      #p <- pmin(p, Inf)
      #print(paste('p ',p))
      # Multinomial log-lik up to constant: sum n_i log p_i
      # Return NEGATIVE log-likelihood for optim/mle2
      pp = -sum(counts * log(p))
      if (debug_flag){
        print(paste('s', s))
        print(paste('Fb ', Fb, 'Fa', Fa, 'Z', Z ))
        print(paste('pp ', pp))
      }
      
      return(pp)
      
    }
  }
  
  minuslogl <- make_binned_nll(
    mids       = df$mid,
    counts     = df$count,
    width      = df$width,
    open_left  = FALSE,      # set as needed
    open_right = FALSE
  )
  
  # MLE with simple box constraints (adjust starts/bounds)
  fit <- mle2(
    minuslogl = minuslogl,
    start     = start.value,
    method    = "L-BFGS-B",
    
  )
  
  fit@fullcoef["s"] = 1+exp(fit@fullcoef["s"])
  
  
  return(fit)
  
}

## pcolony1 must be a CDF: pcolony1(x, c, lambda, log_cond = FALSE)


Simsek_pdf_calc = function(f_name, data_label){
  df <- read.csv(f_name, stringsAsFactors = FALSE)
  
  # Sort by bin edge (just in case)
  df <- df[order(df$TimeAxis), ]
  
  # Bin width and midpoint
  dt  <- c(df$TimeAxis[1]-0, diff(df$TimeAxis))                         # (T_i - T_{i-1})
  tau <- (df$TimeAxis + c(0, df$TimeAxis[-nrow(df)])) / 2 # (T_i + T_{i-1})/2
  
  # Probability density: f(tau) = (Ni / sum(Ni)) / width
  Ni       <- df$TotalDist
  probBin  <- Ni / sum(Ni, na.rm = TRUE)
  f_tau    <- probBin / dt
  
  df_out = data.frame(mid=tau , count = Ni, width = dt, frq=f_tau, Data = data_label)
}


Simsek_cdf_calc = function(f_name, data_label, df = NULL){
  if (is.null(df)){
    df <- read.csv(f_name, stringsAsFactors = FALSE)
  }

  
  # Sort by bin edge (just in case)
  df <- df[order(df[[1]]), ]
  
  # Bin width and midpoint
  dt  <- c(df[[1]][1]-0, diff(df[[1]]))                         # (T_i - T_{i-1})
  #t <- (df$t + c(0, df$t[-nrow(df)])) / 2 # (T_i + T_{i-1})/2
  
  # Probability density: f(tau) = (Ni / sum(Ni)) / width
  Ni       <- df[[2]]
  probBin  <- Ni / sum(Ni, na.rm = TRUE)
  f_tau    <- probBin 

  f_cdf <- cumsum(f_tau)

  df_out = data.frame(t=df[[1]] , count = Ni, width = dt, cdf_data=f_cdf, Data = data_label)
}




# calculate_mean_bin_density = function(df1,lagFunX, lag_fit ){
#   mids = df1$mid
#   counts = df1$count
#   width = df1$width
#   aa <- mids - width/2
#   bb <- mids + width/2
#   
#   Fb <- do.call(lagFunX, c(list(t=bb) , (lag_fit@fullcoef)))
#   Fa <- do.call(lagFunX, c(list(t=aa) , (lag_fit@fullcoef)))
#   
#   Z <- do.call(lagFunX, c(list(t=max(bb)) , (lag_fit@fullcoef)))
#   F_z <- (Fb - Fa)/Z
#   
#   pred = F_z/df1$width
#   return(pred)
# }

calculate_mean_bin_density = function(df1,lagFunX, fit_coefs ){
  mids = df1$mid
  counts = df1$count
  width = df1$width
  aa <- mids - width/2
  bb <- mids + width/2
  
  Fb <- do.call(lagFunX, c(list(t=bb) , fit_coefs))
  Fa <- do.call(lagFunX, c(list(t=aa) , fit_coefs))
  
  Z <- do.call(lagFunX, c(list(t=max(bb)) , fit_coefs))
  F_z <- (Fb - Fa)/Z
  
  pred = F_z/df1$width
  return(pred)
}

calculate_mean_bin_probability = function(df1,lagFunX, fit_coefs ){
  t = df1$t
  aa <- c(0, head(t, -1))
  bb <- t
  
  Fb <- do.call(lagFunX, c(list(t=bb) , fit_coefs))
  Fa <- do.call(lagFunX, c(list(t=aa) , fit_coefs))
  
  Z <- do.call(lagFunX, c(list(t=max(bb)) , fit_coefs))
  F_z <- (Fb - Fa)/Z
  
  pred = F_z
  return(pred)
}


calculate_mean_bin_density_temp = function(df1,lagFunX, fit_coefs ){
  mids = df1$mid
  counts = df1$count
  width = df1$width
  aa <- mids - width/2
  bb <- mids + width/2
  
  Fb <- do.call(lagFunX, c(list(q=bb) , fit_coefs))
  Fa <- do.call(lagFunX, c(list(q=aa) , fit_coefs))
  
  Z <- do.call(lagFunX, c(list(q=max(bb)) , fit_coefs))
  F_z <- (Fb - Fa)/Z
  
  pred = F_z/df1$width
  return(pred)
}

# 
# MaxEnt_mean_bin_density = function(df1,c, lambda ){
#   mids = df1$mid
#   counts = df1$count
#   width = df1$width
#   aa <- mids - width/2
#   bb <- mids + width/2
#   
#   Fb <- pcolony1(bb, c = c, lambda = lambda)
#   Fa <- pcolony1(aa, c = c, lambda = lambda)
#   
#   Z <- pcolony1(max(bb), c = c, lambda = lambda)
#   F_z <- (Fb - Fa)/Z
#   
#   pred = F_z/df1$width
#   return(pred)
# }
MaxEnt_mean_bin_density = function(t,c, lambda ){

  Fb <- pcolony1(bb, c = c, lambda = lambda)
  Fa <- pcolony1(aa, c = c, lambda = lambda)
  
  Z <- pcolony1(max(bb), c = c, lambda = lambda)
  F_z <- (Fb - Fa)/Z
  
  pred = F_z/df1$width
  return(pred)
}


# 
# SimSek_mean_bin_density = function(df1,A1, k, A2, b ){
#   mids = df1$mid
#   counts = df1$count
#   width = df1$width
#   aa <- mids - width/2
#   bb <- mids + width/2
#   Fb <- lag_fun_cdf(bb, A1=A1, k=k, A2=A2, b=b)
#   Fa <- lag_fun_cdf(aa, A1=A1, k=k, A2=A2, b=b)
#   
#   Z <- lag_fun_cdf(max(bb), A1=A1, k=k, A2=A2, b=b)
#   F_z <- (Fb - Fa)/Z
#   
#   pred = F_z/df1$width
#   return(pred)
# }


SimSek_mean_bin_density = function(t,A1, k, A2, b ){

  Fb <- lag_fun_cdf(bb, A1=A1, k=k, A2=A2, b=b)
  Fa <- lag_fun_cdf(aa, A1=A1, k=k, A2=A2, b=b)

  Z <- lag_fun_cdf(max(bb), A1=A1, k=k, A2=A2, b=b)
  F_z <- (Fb - Fa)/Z

  pred = F_z/df1$width
  return(pred)
}





calculate_mean_bin_density_nls = function(df1,lagFunX, lag_fit ){
  mids = df1$mid
  counts = df1$count
  width = df1$width
  aa <- mids - width/2
  bb <- mids + width/2
  
  Fb <- do.call(lagFunX, c(list(t=bb) , (lag_fit@fullcoef)))
  Fa <- do.call(lagFunX, c(list(t=aa) , (lag_fit@fullcoef)))
  
  Z <- do.call(lagFunX, c(list(t=max(bb)) , (lag_fit@fullcoef)))
  F_z <- (Fb - Fa)/Z
  
  pred = F_z/df1$width
  return(pred)
}


# it's actually PDF, But PDF calculated from cdf to deal with binning issue
# this fits a given  dataset two different models (MaxEnt) and Semsik biphasic model
fun_obs_pred_pdf = function(f_1, Data , width = NULL, Total_cells = NULL){
  if (is.null(width)){
    df1 = Simsek_pdf_calc(f_1,Data)
  } else{
    df1 = read.csv(f_1)
    frq = df1$frq
    mid = df1$t
    df1 = data.frame(mid=mid , count =round(frq*width*Total_cells),  width = width, frq=frq, Data = Data)
    
  }
  
  lag_fit1 = fit_maxEnt_lagTime(df1, start.value = list(c = log(800), lambda = log(0.1)))
  lag_fit2 = fit_SimSek_lagTime(df1, start.value = list(k = log(0.06), b= log(2), A1=log(0.03), A2 =log(5)))
  # summary(lag_fit)
  # coef(lag_fit)
  params1=coef(lag_fit1)
  print(params1)
  print(sprintf('AIC our model %s', AIC(lag_fit1)))
  
  params2=coef(lag_fit2)
  print(params2)
  print(sprintf('AIC Simsek model %s', AIC(lag_fit2)))
  
  pred1 = calculate_mean_bin_density(df1, pcolony1, lag_fit1@fullcoef )
  pred2 = calculate_mean_bin_density(df1, lag_fun_cdf, lag_fit2@fullcoef )
  
  
  df_1 = data.frame(t=df1$mid,  frq = (df1$frq), width = df1$width, Data = Data, line = 0)
  df_pred1 = data.frame(t=df1$mid,  frq = (pred1), width =df1$width, Data = sprintf('MaxEnt model - %s', Data), line =1)
  df_pred2 = data.frame(t=df1$mid,  frq = (pred2), width =df1$width, Data = sprintf('Exp + power-law (Şimşek et al. 2019) - %s', Data), line =2)
  
  return(list(df_obs = df_1, df_pred1 = df_pred1 , df_pred2 = df_pred2))
}


cross_validation_maxent = function(training_pdf, test_pdf,  width = NULL, Total_cells = NULL){
  
  lag_fit1 = fit_maxEnt_lagTime(training_pdf, start.value = list(c = log(800), lambda = log(0.1)))
   # summary(lag_fit)
  # coef(lag_fit)
  params1=coef(lag_fit1)
  print(params1)
  print(sprintf('AIC MaxEnt model %s', AIC(lag_fit1)))
  pred_training = calculate_mean_bin_density(training_pdf, pcolony1, lag_fit1@fullcoef )
  r2_maxent_training = r2modified(training_pdf$frq,pred_training, log = T)
  # test data 
  pred_test = calculate_mean_bin_density(test_pdf, pcolony1, lag_fit1@fullcoef )
  r2_maxent_test = r2modified(test_pdf$frq,pred_test, log = T)
  
  df_pred_training = data.frame(t=training_pdf$mid,  frq = (pred_training), width =training_pdf$width)
  
  df_pred_test = data.frame(t=test_pdf$mid,  frq = (pred_test), width =test_pdf$width)
  
   return(list(df_prediction = list(training = df_pred_training, test = df_pred_test), R2 = list(test = r2_maxent_test, training = r2_maxent_training)))
}


cross_validation_simsek = function(training_pdf, test_pdf,  width = NULL, Total_cells = NULL){
  
  lag_fit1 = fit_SimSek_lagTime(training_pdf, start.value = list(k = log(0.06), b= log(2), A1=log(0.03), A2 =log(5)))
  # summary(lag_fit)
  # coef(lag_fit)
  params1=coef(lag_fit1)
  print(params1)
  print(sprintf('AIC Simsek model %s', AIC(lag_fit1)))
  pred_training = calculate_mean_bin_density(training_pdf, lag_fun_cdf, lag_fit1@fullcoef )
  r2_simsek_training = r2modified(training_pdf$frq,pred_training, log = T)
  # test data 
  pred_test = calculate_mean_bin_density(test_pdf, lag_fun_cdf, lag_fit1@fullcoef )
  r2_simsek_test = r2modified(test_pdf$frq,pred_test, log = T)
  
  df_pred_training = data.frame(t=training_pdf$mid,  frq = (pred_training), width =training_pdf$width)
  df_pred_test = data.frame(t=test_pdf$mid,  frq = (pred_test), width =test_pdf$width)
  
  return(list(df_prediction = list(training = df_pred_training, test = df_pred_test), R2 = list(test = r2_simsek_test, training = r2_simsek_training)))
}


# it's actually PDF, But PDF calculated from cdf to deal with binning issue
fun_obs_pred_cdf2 = function(f_1, Data , width = NULL, Total_cells = NULL){
  if (is.null(width)){
    df1 = Simsek_cdf_calc(f_1,Data)
  } else{
    # Assumption is first column is time and second column is frequency (density)
    df1 = read.csv(f_1)
    
    t = df1[[1]]
    average_density = df1[[2]]
    
    Ni = round(average_density*width*Total_cells)
    probBin  <- Ni / sum(Ni, na.rm = TRUE)
    f_tau    <- probBin 
    
    f_cdf <- cumsum(f_tau)
    
    df1 = data.frame(t=t , count =Ni,  width = width, cdf_data=f_cdf, Data = Data)
    
  }
  
  lag_fit1 = fit_maxEnt_lagTime_cdf(df1, start.value = list(c = log(800), lambda = log(0.1)))
  lag_fit2 = fit_SimSek_lagTime_cdf(df1, start.value = list(k = log(0.06), b= log(2), A1=log(0.03), A2 =log(5)))
  # summary(lag_fit)
  # coef(lag_fit)
  params1=coef(lag_fit1)
  print(params1)
  print(sprintf('AIC our model %s', AIC(lag_fit1)))
  
  params2=coef(lag_fit2)
  print(params2)
  print(sprintf('AIC Simsek model %s', AIC(lag_fit2)))
  #print(df1)
  pred1 = calculate_mean_bin_probability(df1, pcolony1, lag_fit1@fullcoef )
  pred2 = calculate_mean_bin_probability(df1, lag_fun_cdf, lag_fit2@fullcoef )
  
  
  pred1 = cumsum(pred1)
  pred2 = cumsum(pred2)
  
  
  df_1 = data.frame(t=df1$t,  frq = (df1$cdf_data), width = df1$width, Data = Data, line = 0)
  df_pred1 = data.frame(t=df1$t,  frq = (pred1), width =df1$width, Data = sprintf('MaxEnt MLE fit (use cdf) Prediction %s', Data), line =1)
  df_pred2 = data.frame(t=df1$t,  frq = (pred2), width =df1$width, Data = sprintf('Exponential + power-law (Şimşek et al. 2019) %s', Data), line =2)
  
  return(list(df_obs = df_1, df_pred1 = df_pred1 , df_pred2 = df_pred2))
}

