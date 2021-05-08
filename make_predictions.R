
rm(list=ls())

{library(dplyr)
  library(lubridate)
  library(forecast)
  library(forecastHybrid)
  #library(ForecastCombinations)
  #library(ForecastComb)
  library(ggplot2)
  library(stringr)
  #library(urca)
  library(tidyverse) 
  #library(dlm)
  library(skimr)
  library(xts)
  library(zoo)
  #require(nortest)
  #library(MARSS)
  #library(bsts)
  library(prophet)
  #library(fUnitRoots)
  #library(TSA)
  #library(caret)  
  #library(vars)  
  #library(BVAR)  
  #library(fpp)
  #library(forecTheta)
}

setwd("C:\\Users\\Alexander\\Documents\\Arbeit\\Siemens\\ts_stuff\\data\\")

## Save data as matrix and as multivar. ts obj.:
#-------------------------------------------------
data_cleaned = readRDS("cleaned_data.rds")
data = readRDS("raw_data.rds")
#ts_all = readRDS("data_mts_78.rds")
#mts_data = readRDS("data_mat_78.rds")
ts_all = readRDS("data_mts_42.rds")
mts_data = readRDS("data_mat_42.rds")
#--------------------------------------------------
dim(mts_data)

setwd("C:\\Users\\Alexander\\Documents\\Arbeit\\Siemens\\ts_stuff\\")

source("candidate_models.R")
source("BayesOpt_utils.R")
#source("hyper_opt.R")
source("hyper_opt2.R")
source("model_selection.R")

dev.off()
par(mfrow=c(2,2))

#---------------------
n.ahead = 6
(TT = ncol(mts_data))
hyper_tune = F
#-----------------------

set.seed(139)

setwd("C:\\Users\\Alexander\\Documents\\Arbeit\\Siemens\\ts_stuff\\data\\")

#------------------------- Loop over data set ----------------------------------------
model_labels = NULL
for(i in 1:nrow(mts_data)){
  
  cat(paste0("\nSeries (Nr.",i," of ",nrow(mts_data),"): ", colnames(ts_all)[i]),"\n")
  
  (y = ts(mts_data[i,], start = start(ts_all), end = end(ts_all), frequency = frequency(ts_all)))
  (y_obs = subset(y, start = 1, end = TT-n.ahead))   # observed data
  #(y_unseen = subset(y, start = TT-n.ahead+1, end = TT))
  
  lambda <- BoxCox.lambda(y_obs)
  
  (y_obs = BoxCox(y_obs, lambda = lambda))
  
  #(y_obs = log(1 + y_obs))
  #(y_obs = tsclean(y_obs))
  (y_obs = na.interp(y_obs))

  #png("descr.png")  
  par(mfrow=c(2,2))
  plot.ts(y_obs, main= paste("Series ", colnames(ts_all)[i]), col = "black", lty = 1, ylab="demand") ;
  points(y_obs, col="red")
  hist(y_obs,nclass=30, main="", col="grey", xlab="demand")
  acf(as.numeric(y_obs), lag.max = 50, main="ACF")
  pacf(as.numeric(y_obs), lag.max = 50, main = "PACF")
  #graphics.off()
  
  # Run:
  #-------
  out = model_selection(y_obs, p_oob = 0.8, n.ahead = n.ahead)

  # Plot RMSE as a function of forecast horizon:
  #------------------------------------------------
  dat <- data.frame(pro = cumsum(out$nstep_rmse$pro_RMSE), 
                    ets = cumsum(out$nstep_rmse$ets.RMSE),
                    ens = cumsum(out$nstep_rmse$ens.RMSE),
                    fnn = cumsum(out$nstep_rmse$fnn.RMSE),
                    rwd = cumsum(out$nstep_rmse$rw.RMSE))
  
  #png("RMSEs.png")
  matplot(dat, type = c("b"),pch=1,col = 1:ncol(dat), ylab="RMSE", xlab="forecast horizon", 
          main = paste("Winner:", out$winner)) 
  legend("topleft", legend = c('Prophet', 'ETS', 'ENSEM.','FNN', 'RwD'), col=1:ncol(dat), pch=1, bty="n")
  #dev.off()
  #-------------------------------------------------------------------------------------
  
  if(out$winner == 'PROPHET'){
    if(hyper_tune){
      
        best_params_ba = hyper_optim(y_obs, p_oob=0.85, n.ahead=6, winner = out$winner) 
        forec = run_prophet(y.window = y_obs, n.ahead = n.ahead, 
                            growth = 'linear',       
                            daily.seasonality = F,
                            weekly.seasonality = F,
                            seasonality.prior.scale = best_params_ba[['seasonality_prior_scale']], 
                            changepoint.prior.scale = best_params_ba[['changepoint_prior_scale']],
                            n.changepoints = best_params_ba[['n_changepoints']])
     } else {
      forec = run_prophet(y.window = y_obs, n.ahead = n.ahead, 
                           growth = 'linear',
                           daily.seasonality = F,
                           weekly.seasonality = F,
                           n.changepoints = 5) 
     }
    
    yhat_pro = forec$forecast$yhat
    (y_hat6 = yhat_pro[(length(yhat_pro)-n.ahead+1):length(yhat_pro)])
   } 
  #--------------------------------------------------------------------------
  if(out$winner == 'ARFIMA'){

    tryCatch({    # cannot cope with NAs!!!
      fore_arfima = run_arfima(y.window = y_obs, n.ahead = n.ahead, estim = "ls")
      (y_hat6 = fore_arfima$forecast$mean)
    }, warning = function(w) {
      print(w)
    }, error = function(e) {
      print(e)
    })
  }
  #-----------------------------------------------------------------
  if(out$winner == 'STLM-ETS'){
      stlm = run_stlm_ets(y.window = y_obs, n.ahead, modelfunction=ar, 
                          allow.multiplicative.trend = T)
      (y_hat6  = as.numeric(stlm$forecast$mean))
  }
  #-----------------------------------------------------------------
  if(out$winner == 'ETS'){
      fore_ets = run_ets(y.window = y_obs, n.ahead)
      (y_hat6 = as.numeric(fore_ets$forecast$mean))
  }  
  #-----------------------------------------------------------------
  if(out$winner == 'ENSEM'){
    
    fore_ensem = run_ensem(y.window = y_obs, 
                      n.ahead,models = "efn", 
                      weights="equal", errorMethod = "RMSE")
    (y_hat6 = as.numeric(fore_ensem$forecast$mean))
  } 
  #-----------------------------------------------------------------
  if(out$winner == 'TBATS'){
    
    fore_tbats = run_tbats(y.window = y_obs, n.ahead, 
                           seasonal.periods = frequency(y_obs))
    (y_hat6 = as.numeric(fore_tbats$forecast$mean))
  } 
  #--------------------------------------------------------------
  if(out$winner == 'RWD'){
    
    tryCatch({
      rwd = run_rwd(y.window = y_obs, n.ahead = n.ahead)
      (y_hat6 = rwd$forecast$mean)  
    }, warning = function(w) {
      print(w)
    }, error = function(e) {
      print(e)
    }#, finally = {print("Out")}
    )
  }
  #-----------------------------------------------------------------
  if(out$winner == 'NNAR'){
    
    if(hyper_tune){
        best_params_ba = hyper_optim(y_obs, p_oob = 0.85, n.ahead = 6, winner = out$winner) 
        
        fnn = run_fnn(y.window = y_obs, n.ahead = n.ahead, 
                      p = best_params_ba[['p_nonsea']], size = best_params_ba[['nof_nodes']])
    } else{
      fnn = run_fnn(y.window = y_obs, n.ahead = n.ahead)
    }
    (y_hat6 = as.numeric(fnn$forecast$mean))
  }
  #--------------------------------------------------------------------

  # Save series name and winner:
  model_labels = rbind(model_labels, c(series = colnames(ts_all)[i], model = out$winner))

  #(yhat6_retransf = exp(y_hat6)-1)

  # Reverse Box-Cox transformation:
  (yhat6_retransf = InvBoxCox(y_hat6, lambda = lambda))
  
  stopifnot(length(yhat6_retransf) == n.ahead)
  #up95 = fore_ets$forecast$upper[,'95%']
  #low95 = fore_ets$forecast$lower[,'95%']
  
  # Assign and retransform:
  mts_data[i,(ncol(mts_data)-n.ahead+1):ncol(mts_data)] = yhat6_retransf
  
  #dev.off()
  #png("series_and_fore.png")
  par(mfrow=c(1,1), ask=F)
  plot.ts(mts_data[i,], lty=2,main= paste("Series ", colnames(ts_all)[i]), ylab="demand")
  lines(c(ncol(mts_data)-n.ahead+1):ncol(mts_data), yhat6_retransf, col="red")
  #points(c(ncol(mts_data)-n.ahead+1):ncol(mts_data), yhat6_retransf, col="red")
  abline(v = ncol(mts_data)-n.ahead+1, lty = 3, col="grey")
  #lines(c(ncol(mts_data)-n.ahead+1):ncol(mts_data), exp(up95)-1, col="red", lty=3)
  #lines(c(ncol(mts_data)-n.ahead+1):ncol(mts_data), exp(low95)-1, col="red", lty=3)
  #dev.off() 
  
}; print("FINISHED!")
#-----------------------------------------------------------------------------------  

graphics.off()

model_labels

## Save results:
#-------------------------
setwd("C:\\Users\\Alexander\\Documents\\Arbeit\\Siemens\\ts_stuff\\data\\")
#saveRDS(mts_data, "mts_data_final78.rds")
#dput(model_labels, 'model_labels_78.R')
saveRDS(mts_data, "mts_data_final42.rds")
dput(model_labels, 'model_labels_42.R')


colnames(mts_data) = format(as.yearmon(time(ts_all)), "%Y-%m")
row.names(mts_data) = str_replace(row.names(mts_data),'-',' ')

dim(mts_data)

output = data.frame(model_labels[,'model'], date=mts_data)
colnames(output) = c("model", colnames(mts_data))

View(output)

# Save results:
#----------------
write.csv(output, file = "output78.csv")
#write.csv(output, file = "output42.csv")

output42 = read.csv2("output42.csv")
output78 = read.csv2("output78.csv")

View(output78)

head(output42)

dim(output42)
dim(output78)

ts_all78 = readRDS("data_mts_78.rds")
mts_data78 = readRDS("data_mat_78.rds")
ts_all42 = readRDS("data_mts_42.rds")
mts_data42 = readRDS("data_mat_42.rds")

colnames(output78)= c('Country/Product','model',format(as.yearmon(time(ts_all78)), "%Y-%m"))
colnames(output42)= c('Country/Product','model',format(as.yearmon(time(ts_all42)), "%Y-%m"))

dput(output78, 'output78.R') 
dput(output42, 'output42.R')  

i = 1
output78_form = apply(output78[,-c(1,2)], 1,function(x) as.numeric(x))


as.character(output78[i,'Country/Product'])

(y = ts(output78_form[i,], start = start(ts_all78), end = end(ts_all78), frequency = frequency(ts_all42)))


plot.ts(mts_data[i,], lty=2,main= paste("Series ", colnames(ts_all)[i]), ylab="demand")
lines(c(ncol(mts_data)-n.ahead+1):ncol(mts_data), yhat6_retransf, col="red")
#points(c(ncol(mts_data)-n.ahead+1):ncol(mts_data), yhat6_retransf, col="red")
abline(v = ncol(mts_data)-n.ahead+1, lty = 3, col="grey")


y78 = output78 %>% 
  filter(`Country/Product`=='Country_0 Product_1') %>% 
  dplyr::select(-one_of('Country/Product', 'model')) ; y78

y42 = output42 %>% 
        filter(`Country/Product`=='Country_0 Product_4') %>% 
        dplyr::select(-one_of('Country/Product', 'model')); y42

glimpse(y42)


# Prepare results for presentation:
#-------------------------------------------
model_labels_78 = dget('model_labels_78.R')
model_labels_42 = dget('model_labels_42.R')
mts_data_final42 = readRDS("mts_data_final42.rds")
mts_data_final78 = readRDS("mts_data_final78.rds")

dim(model_labels_78)
dim(model_labels_42)
dim(mts_data_final42)
dim(mts_data_final78)

colnames(mts_data_final78) = format(as.yearmon(time(ts_all78)), "%Y-%m")
row.names(mts_data_final78) = str_replace(row.names(mts_data78),'-',' ')
colnames(mts_data_final42) = format(as.yearmon(time(ts_all42)), "%Y-%m")
row.names(mts_data_final42) = str_replace(row.names(mts_data42),'-',' ')


output42 = data.frame(model_labels_42[,'model'], date=mts_data_final42)
colnames(output42) = c("model", colnames(mts_data_final42))
output78 = data.frame(model_labels_78[,'model'], date=mts_data_final78)
colnames(output78) = c("model", colnames(mts_data_final78))

View(output42)
View(output78)












