# Read counts

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
  
  df_out = data.frame(t=tau , frq=f_tau, Data = data_label)
}


# 
# # Drop the first row (no left edge)
# out <- data.frame(tau = tau[-1], width = dt[-1], f_tau = f_tau[-1])
# 
# # Save
# write.csv(out, "Emrah_pdf.csv", row.names = FALSE)


data_label1 = 'Dataset 1'
data_label2 = 'Dataset 2'
data_label3 = 'Dataset 3'
df_obs1 = Simsek_pdf_calc('/home/ahmed/Dropbox/Ahmed Abdullah/AMR/Emrah.csv',data_label1)
df_obs2 = Simsek_pdf_calc('/home/ahmed/Dropbox/Ahmed Abdullah/AMR/Emrah2.csv',data_label2)
df_obs3 = Simsek_pdf_calc('/home/ahmed/Dropbox/Ahmed Abdullah/AMR/Triangles.csv',data_label3)


df_MaxEnt1 = df_MaxEnt_Prediction(df_obs1,  start_c=5000, start_lambda =0.11 ,data_label)

df_obs = rbind(df_obs1, df_obs2, df_obs3)
t =df_obs$t

#Prediction Emrah
y_pred=lag_fun(t, A1, k, A2, b, tau0)
df_Simsek_Pred = data.frame(t=t,frq=y_pred)
df_Simsek_Pred$Data = 'Exponential + power-law (Şimşek et al. 2019)'



df_MaxEnt = df_MaxEnt_Prediction(df_obs,  start_c=5000, start_lambda =0.11 ,'MaxEnt pred')
df_all = rbind(df_obs1,df_obs2,df_obs3, df_MaxEnt,df_Simsek_Pred)


#df_all = rbind(df_obs1,df_obs2,df_MaxEnt1,df_MaxEnt2,df_Simsek_Pred)

ggplot(data=df_all) + geom_line(aes(x=t,y=(frq), col = Data))+ geom_point(aes(x=t,y=(frq), col = Data)) + scale_y_continuous(trans='log10')+xlab('Lag time (min)')+ylab('Probalility density (log scale)')+
  theme(
    axis.title.x   = element_text(size = 16),
    axis.title.y   = element_text(size = 16),
    axis.text.x    = element_text(size = 14),
    axis.text.y    = element_text(size = 14),
    legend.title   = element_text(size = 14),
    legend.text    = element_text(size = 13),
    plot.title     = element_text(size = 18, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line      = element_line(),
    axis.ticks     = element_line(),
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white")
  )




r2_maxent = cor(df_obs$frq, df_MaxEnt$frq)^2
r2_simsek = cor(df_obs$frq, df_Simsek_Pred$frq )^2







