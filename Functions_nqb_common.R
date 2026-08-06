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
library(GA)

library(scales)
source("/home/ahmed/Dropbox/Ahmed Abdullah/NQBMatlab/Supporting_scripts/quadgk51_23.R") #put file path of quadgk51")



data_load <- function(k){
  df1 = read.csv('/home/ahmed/Education/UVA/Rotation/Martin Wu/colony_experiment/lagscan1.csv')
  df2 = read.csv('/home/ahmed/Education/UVA/Rotation/Martin Wu/colony_experiment/lagscan2.csv')
  df3 = read.csv('/home/ahmed/Education/UVA/Rotation/Martin Wu/colony_experiment/lagscan3.csv')
  df4 = read.csv('/home/ahmed/Education/UVA/Rotation/Martin Wu/colony_experiment/lagscan4.csv')
  
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



my_logit_transform = function(m, s){
  m*exp(s)/(1+exp(s))
}


my_inverse_logit_transform = function(m,s){
  log(s/(m-s))
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
  
  combined_ap_df = combined_ap_df[combined_ap_df$GrowthTimes>GrowthTime_threshold,]
  if (!is.null(GrowthTime_upper)){
    combined_ap_df = combined_ap_df[combined_ap_df$GrowthTimes<GrowthTime_upper,] 
  }
  
  col_ap <- combined_ap_df$TimeAxis
  # unlist the file 
  myVec = unlist(col_ap)
  return(myVec)
}



dcolony = Vectorize(function(t,lamda_c){
  f_t=1/t^2*lamda_c*exp(-lamda_c/t)
  return(f_t)
},'t')

pcolony = Vectorize(function(t, lamda_c){
  F_t=exp(-lamda_c/t)
  return(F_t)
}, 't')

qcolony = function(pp, lamda_c){
  qq=lamda_c/log(1/pp)
  return(qq)
}