

## Likelihood ``
neg_loglik_fun = function(x, params, f){
  # summarize abundance counts
  
  params = c(params, list(log_cond = T))
  
  ctdf = aggregate(1:length(x), by = data.frame(abundance = x), FUN = length)
  
  # get all formal argument names of f
  all_args = names(formals(f))
  
  # find which argument is not supplied in params
  missing_args = setdiff(all_args, names(params))
  
  # if there is exactly one missing argument, assume that’s the data input
  if (length(missing_args) != 1) {
    stop("Cannot determine which argument to fill — check params or f.")
  }
  
  data_arg_name = missing_args[1]
  
  # build argument list for do.call
  arg_list = c(setNames(list(ctdf$abundance), data_arg_name), params)
  
  # call f with both abundance and params
  pps = do.call(f, arg_list)
  
  # compute negative log-likelihood
  pp = -sum(pps * ctdf$x)
  
  return(pp)
}

constant_work_ll = function(x, lambda , c, t0, sd){
  ctdf = aggregate(1:length(x),by=data.frame(abundance=x),FUN=length)
  pps = dcompound1(Observed_T = ctdf$abundance, lambda = lambda, c = c, t0 = t0, sd =sd, log_cond = T)
  pp =-sum(pps*ctdf$x)
  return(pp)
}
