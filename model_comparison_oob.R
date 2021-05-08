
rm(list=ls())

{library(dplyr)
  library(lubridate)
  #library(prophet)
  library(forecast)
  library(forecastHybrid)
  #library(ForecastCombinations)
  library(ForecastComb)
  library(ggplot2)
  #library(urca)
  #library(dlm)
  library(skimr)
  library(xts)
  library(zoo)
  #require(nortest)
  library(MARSS)
  library(bsts)
  library(prophet)
  #library(fUnitRoots)
  #library(TSA)
  #library(caret)  
  #library(vars)  
  library(BVAR)  
  library(fpp)
  library(forecTheta)
}

setwd("C:\\Users\\Alexander\\Documents\\Arbeit\\Siemens\\ts_stuff\\data\\")
setwd("C:\\Users\\Alexander\\Documents\\Arbeit\\Siemens\\ts_stuff\\")

source("candidate_models.R")
source("BayesOpt_utils.R")

#------------------------------------------------------------------------------------
model_selection = function(y_obs, p_oob = 0.85, n.ahead = 6)
{
  (T_obs = length(y_obs))
  (R = floor(T_obs * p_oob))                       #the higher n, the more accurate the forecasts! 
  
  # Rolling window approach to compute out-of-bag error:
  #-------------------------------------------------------------------
  j = 0 ; all_metrics = arfima.MAE = NULL
  while(j<=(T_obs-n.ahead-R)) {
    
    train_index = 1:(R+j)
    valid_index = (R+j+1):(R+j+n.ahead) 
    (y.window = ts(y_obs[train_index], fre=frequency(y_obs),start=start(y_obs)))       #used subsample
    (y.true = y_obs[valid_index])                             #actual values for comparison 
    
    # Impute missing values and smooth outliers:
    #--------------------------------------------
    (y.window = tsclean(y.window))
    
    cat("\nTest error run nr.",j," (out of ",T_obs-n.ahead-R,") [Training size R: ",length(y.window),"]\n------------------------------------------------------------------\n",sep=""); 

    # Exponential smoothing state space model:
    #-------------------------------------------
    #print("Fitting ETS")
    #fore_ets = y.window %>% ets() %>% forecast(h = n.ahead) #%>% autoplot()
    #(yhat_ets = fore_ets$mean)
    
    # Automatic ARIMA forecasts:
    #-----------------------------
    # print("Fitting ARIMA")
    # fore_arima = y.window %>% auto.arima(lambda = "auto", stationary = T) %>% 
    #   forecast(h = n.ahead) 
    # (yhat_arima = fore_arima$mean)
    
    # Automatic ARFIMA forecasts:
    #-----------------------------
    #print("Fitting ARFIMA")
    #fore_arfima = y.window %>% arfima() %>% forecast(h = n.ahead) #%>% autoplot()
    #(yhat_arfima = fore_arfima$mean)
    #tsdisplay(residuals(fit))
    
    tryCatch({
      fore_arfima = run_arfima(y.window, n.ahead, estim = c("ls"))
      (yhat_arfima = fore_arfima$forecast$mean)
    }, warning = function(w) {
      print(w)
    }, error = function(e) {
      print(e)
    }#, finally = {print("Out")}
    )
    
    err_arf = y.true - yhat_arfima
    (arfima.MAE = rbind(arfima.MAE, abs(err_arf))); 
    
    
    # Forecasting with STL
    #------------------------
    print("Fitting STLM")
    fore_stlm = y.window %>% stlm(modelfunction=ar) %>% forecast(h = n.ahead) #%>% autoplot()
    (yhat_stlm = fore_stlm$mean)
    
    #print("Fitting STLF")
    #fore_stlf = y.window %>% stlf(lambda=0) %>% forecast(h = n.ahead) #%>% autoplot()
    #(yhat_stlf = fore_stlf$mean)
    
    #print("Fitting STL")
    #fore_stl = y.window %>% stl(s.window='periodic') %>% forecast(h = n.ahead) #%>% autoplot()
    #(yhat_stl = fore_stl$mean)
    
    # TBATS : Exponential smoothing state space model 
    # with Box-Cox transformation, ARMA errors, 
    # Trend and Seasonal componentsforecasts
    #-----------------------------------------------------
    #print("Fitting TBATS")
    #fore_tbats = y.window %>% tbats() %>% forecast(h = n.ahead) #%>% autoplot()
    #(yhat_tbats = fore_tbats$mean)
    
    # Feed-forward neural networks with a single hidden layer and lagged inputs
    #----------------------------------------------------------------------------
    fnn = run_fnn(y.window, n.ahead)
    (yhat_nnetar = fnn$forecast)
    
    # Linear time series model:
    #---------------------------
    fit_lm <- tslm(y.window ~ trend + season)
    fore_lm = forecast(fit_lm, h = n.ahead)
    (yhat_lm = fore_lm$mean)
    
    # Bsts model
    ----------------------------
    # print("Fitting BSTS")
    # ss <- AddLocalLinearTrend(list(), y.window)
    # ss <- AddSeasonal(ss, y.window, nseasons = 12)
    # bsts.model <- bsts(y.window, state.specification = ss, niter = 500, ping=0, seed=2016)
    # burn <- SuggestBurn(0.1, bsts.model)    ### Get a suggested number of burn-ins
    # p <- predict.bsts(bsts.model, horizon = n.ahead, burn = burn, quantiles = c(.025, .975))
    # (fore_bsts = p$mean)
    
    # Naive RW
    #----------------------------------------------------------------------------
    tryCatch({
      #fore_rw = y.window %>% naive() %>% forecast(h = n.ahead) #%>% autoplot()
      #fore_rw = y.window %>% rwf(drift = T) %>% forecast(h = n.ahead)
      #(yhat_rw = fore_rw$mean)
      rwd = run_rwd(y.window, n.ahead)
      (yhat_rwd = rwd$forecast$mean)  
    }, warning = function(w) {
      print(w)
    }, error = function(e) {
      print(e)
    }#, finally = {print("Out")}
    )
    
    # Naive SRW
    #----------------------------------------------------------------------------
    #print("Fitting SRW")
    #fore_srw = y.window %>% snaive() %>% forecast(h = n.ahead) 
    #(yhat_srw = fore_srw$mean)
    
    # Model ensemble:
    #------------------
    ensem = run_ensem(y.window, n.ahead,models = "afnst", weights="equal", errorMethod = "MAE")
    (yhat_ensem = ensem$forecast$mean)
    
    # Dynamic optimized theta model:
    #---------------------------------
    #print("Fitting Optimized Theta model")
    #dynopt_theta <- dotm(y.window, h=n.ahead, s='additive')
    #yhat_dtheta = dynopt_theta$mean
    
    # Prophet:
    #--------------------
    forec = run_prophet(y.window, n.ahead, growth = 'linear', n.changepoints = 3)
    yhat_pro = forec$forecast[valid_index, 'yhat']
    
    # Random Forest: only 1-step ahead
    #-----------------------------------
    # lags = 1
    # train_set = data.frame(embed(y.window,d=lags+1)) ;#y0 = yy[,1] ; train_set = yy[,-1]
    # colnames(train_set) = c("y",paste0("x",1:lags))
    # if(j==0)
    #    lrn = makeLearner("regr.randomForest")     # Random Forests 
    # 
    # tail(train_set)
    # trainTask = makeRegrTask(data = train_set, target = "y")
    # model = train(lrn, trainTask); 
    # pred = predict(model, newdata = train_set)         # predict 'new' data, i.e. of the test set (or validation set)     
    
    #---------------
    # Performances:
    #---------------
    metrics = c(#ETS = accuracy(fore_ets, y.true)["Test set",c("MAPE")],
                #ARIMA = accuracy(fore_arima, y.true)["Test set",c("MAPE")],
                ARFIMA = accuracy(yhat_arfima, y.true)["Test set",c("MAPE")],
                `STLM-ETS` = accuracy(fore_stlm, y.true)["Test set",c( "MAPE")],
                #`STL-ETS` = accuracy(fore_stl, y.true)["Test set",c("MAPE")],
                #`STLF-ETS` = accuracy(fore_stlf, y.true)["Test set",c( "MAPE")],
                NNAR = accuracy(yhat_nnetar, y.true)["Test set",c( "MAPE")],
                #TBATS = accuracy(fore_tbats, y.true)["Test set",c( "MAPE")],
                RWD = accuracy(yhat_rwd, y.true)["Test set",c( "MAPE")],
                #SRW_B2 = accuracy(fore_srw, y.true)["Test set",c( "MAPE")],
                Naive_Avg = accuracy(meanf(y.window, h = n.ahead), y.true)["Test set",c( "MAPE")],
                LM = accuracy(fore_lm, y.true)["Test set",c( "MAPE")],
                #OptThe = accuracy(yhat_dtheta, y.true)["Test set",c( "MAPE")],
                ENSEM = accuracy(yhat_ensem , y.true)["Test set",c( "MAPE")],
                PROPHET = accuracy(yhat_pro, y.true)["Test set",c("MAPE")]
                #BSTS = accuracy(fore_bsts, y.true)["Test set",c("MAPE")]
    )
    all_metrics = rbind(all_metrics, metrics)
    j=j+1
  }
  (mean_acc = colMeans(all_metrics)) 
  
  print("Finished!")
  cat("Best model:", names(mean_acc)[which.min(mean_acc)])
  rownames(all_metrics) = paste0("i=",1:nrow(all_metrics))
  
  return(list(all_metrics = all_metrics, mean_acc = mean_acc, 
              winner = names(mean_acc)[which.min(mean_acc)]))
}
#------------------------------------------------------------------------------


i = 25
(TT = ncol(mts_data))
(y = ts(mts_data[i,], start = start(ts_all), end = end(ts_all), frequency = frequency(ts_all)))
(y_obs = subset(y, start = 1, end = TT-6))   # observed data
(y_unseen = subset(y, start = TT-5, end = TT))

# Plot:
#----------------
par(mfrow=c(2,2))
plot.ts(y_obs, main= paste("Series ", colnames(ts_all)[i]), col = "black", lty = 1, ylab="demand") ;
points(y_obs, col="red")
hist(y_obs,nclass=30, main="", col="grey", xlab="demand")
acf(ts(na.omit(y_obs)), lag.max = 25, main="ACF")
pacf(ts(na.omit(y_obs)), lag.max = 25, main = "PACF")

#(y_obs = tsclean(y_obs))
#(y_obs = diff(log(1 + y_obs)))
#(y_obs = log(1 + y_obs))

# Plot:
#---------
plot.ts(y_obs, main= paste("Series ", colnames(ts_all)[i]), col = "black", lty = 1, ylab="demand") ;
points(y_obs, col="red")
hist(y_obs,nclass=30, main="", col="grey", xlab="demand")
acf(ts(na.omit(y_obs)), lag.max = 25, main="ACF")
pacf(ts(na.omit(y_obs)), lag.max = 25, main = "PACF")


# Run:
#-------
out = model_selection(y_obs, p_oob = 0.9, n.ahead = 6)

out$winner




