
rm(list=ls())

{library(dplyr)
library(lubridate)
#library(prophet)
library(forecast)
#library(forecastHybrid)
#library(ForecastCombinations)
library(ggplot2)
library(urca)
library(dlm)
library(skimr)
library(xts)
library(zoo)
#require(nortest)
library(MARSS)
library(prophet)
library(bsts)
library(KFAS)
#library(mlr)
}

setwd("C:\\Users\\Alexander\\Documents\\Arbeit\\Siemens\\ts_stuff\\data\\")


## Save data as matrix and as multivar. ts obj.:
#-------------------------------------------------
data_cleaned = readRDS("cleaned_data.rds")
data = readRDS("raw_data.rds")
ts_all = readRDS("data_mts.rds")
mts_data = readRDS("data_mat.rds")
#--------------------------------------------------
dim(mts_data)


# Filter out missings (to be predicted):
#---------------------------------------
#data_compl = mts_data[,-c(73:78)]      # numeric matrix
#data_compl_new = ts_all[-c(73:78),]   # ts matrix object
#dim(data_compl)       # M * T , with t = 1...T and m=1...M time series
#dim(data_compl_new)     # T * M

dim(data_compl)

# Pearson correlations: 
cm = cor(t(data_compl))
dim(cm)
cm


i = 3

n.ahead = 6

(TT = ncol(mts_data))

(y = ts(mts_data[i,], start = start(ts_all), end = end(ts_all), frequency = frequency(ts_all)))
(y_obs = subset(y, start = 1, end = TT-n.ahead))   # observed data

#(y_obs = log(1 + y_obs))
#y_obs = na.interp(y_obs)

plot(y_obs)

(T_obs = length(y_obs))
(R = floor(T_obs * .85))                       #the higher n, the more accurate the forecasts! 

#monthdays()
#easter(tsclean(y_obs))

#y = tsclean(y_obs, lambda="auto")
#y = y_obs
#fit <- tslm(y ~ trend + season, biasadj = T)

# plot(forecast(fit, h = n.ahead, robust = T))
# summary(fit)
# CV(fit)
# 
# autoplot(y, series="Data") +
#   autolayer(fitted(fit), series="Fitted") +
#   xlab("Year") + ylab("Megalitres") +
#   ggtitle("Quarterly Beer Production")


#model <- SSModel(Nile ~ SSMtrend(1,Q=1469), distribution="poisson")
#pred <- predict(model,n.ahead=10)


#far2 <- function(x, h){forecast(ets(y_obs), h=h)}
#e <- tsCV(lynx, far2, h=6)

# Rolling window approach to compute out-of-bag error:
#-------------------------------------------------------------------
j = 0 ; all_metrics = ets.RMSE = arima.RMSE = arfima.RMSE = fnn.RMSE = rw.RMSE = pro.RMSE = ens.RMSE = NULL
while(j<(T_obs-n.ahead-R)) 
{
  train_index = 1:(R+j)
  valid_index = (R+j+1):(R+j+n.ahead) 
  (y.window = ts(y_obs[train_index], fre=frequency(y_obs),start=start(y_obs)))       #used subsample
  (y.true = y_obs[valid_index])                             #actual values for comparison 
  
  if(any(is.na(y.true))){
     warning('NAs in y.true!!!')
  }
  
  # Impute missing values and smooth outliers:
  #--------------------------------------------
  #(y.window = tsclean(y.window))
  (y.window = na.interp(y.window))
  
  cat("\nTest error run nr.",j+1," (out of ",T_obs-n.ahead-R,") [Training size R: ",length(y.window),"]\n------------------------------------------------------------------\n",sep=""); 
  cat("\nTraining batch (index) ", min(train_index), max(train_index), "\n")
  cat("Validation batch (index) ", min(valid_index), max(valid_index),"\n\n")
  
  #length(y.window)
  #length(y.true)
  #y_obs
  
  # Exponential smoothing state space model:
  #-------------------------------------------
  print("Fitting ETS")
  fore_ets = y.window %>% ets(lambda = "auto") %>% forecast(h = n.ahead) #%>% autoplot()
  (yhat_ets = fore_ets$mean)
  
  bias_ets = yhat_ets - y.true
  (ets.RMSE = rbind(ets.RMSE, bias_ets^2)); 
  
  # Automatic ARIMA forecasts:
  #-----------------------------
  print("Fitting ARIMA")
  #xreg = fourier(y.window, K = 1)    # fourier pairs; maximum here: K = (seasonal period)/2 = 6
  #xreg_new = fourier(y.window, K = 2, h = n.ahead)    # fourier pairs future
  
  fore_arima = y.window %>% auto.arima(lambda = "auto") %>% 
                            forecast(h = n.ahead) #%>% autoplot()
  (yhat_arima = fore_arima$mean)
  
  bias_arima = yhat_arima - y.true
  (arima.RMSE = rbind(arima.RMSE, bias_arima^2)); 
  
  # Automatic ARFIMA forecasts:
  #-----------------------------
  print("Fitting ARFIMA")
  fore_arfima = y.window %>% arfima(lambda = "auto") %>% forecast(h = n.ahead) #%>% autoplot()
  (yhat_arfima = fore_arfima$mean)
  
  bias_arfima = yhat_arfima - y.true
  (arfima.RMSE = rbind(arfima.RMSE, bias_arfima^2)); 
  
  # Forecasting with STL
  #------------------------
  print("Fitting STLM")
  fore_stlm = y.window %>% stlm(modelfunction=ar, lambda = "auto") %>% forecast(h = n.ahead) #%>% autoplot()
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
  print("Fitting TBATS")
  fore_tbats = y.window %>% tbats(lambda = "auto") %>% forecast(h = n.ahead) #%>% autoplot()
  (yhat_tbats = fore_tbats$mean)
  
  # Feed-forward neural networks with a single hidden layer and lagged inputs
  #----------------------------------------------------------------------------
  print("Fitting NNETAR")
  xreg = fourier(y.window, K = 2)    # fourier pairs
  xreg_new = fourier(y.window, K = 2, h = n.ahead)    # fourier pairs future
  
  fore_nnetar = y.window %>% nnetar(p=2, P=1, xreg = xreg, lambda = "auto") %>% 
                          forecast(h = n.ahead, xreg = xreg_new) #%>% autoplot()
  (yhat_nnetar = fore_nnetar$mean)
  
  bias_fnn = yhat_nnetar - y.true
  (fnn.RMSE = rbind(fnn.RMSE, bias_fnn^2)); 
  
  #fore_baggedETS = y.window %>% baggedETS() %>% forecast(h = n.ahead) #%>% autoplot()
  #(yhat_baggedETS = fore_baggedETS$mean)
  
  # Linear time series model:
  #---------------------------
  fit_lm <- tslm(y.window ~ trend + season, biasadj = T, lambda = "auto")
  fore_lm = forecast(fit_lm, h = n.ahead)
  (yhat_lm = fore_lm$mean)
  
  # Bsts model
  #----------------------------
  # print("Fitting BSTS")
  # ss <- AddLocalLinearTrend(list(), y.window)
  # ss <- AddSeasonal(ss, y.window, nseasons = 12)
  # bsts.model <- bsts(y.window, state.specification = ss, niter = 500, ping=0, seed=2016)
  # burn <- SuggestBurn(0.1, bsts.model)    ### Get a suggested number of burn-ins
  # p <- predict.bsts(bsts.model, horizon = n.ahead, burn = burn, quantiles = c(.025, .975))
  # (fore_bsts = p$mean)
  
  # Naive RW
  #----------------------------------------------------------------------------
  print("Fitting RW")
  fore_rw = y.window %>% naive() %>% forecast(h = n.ahead) #%>% autoplot()
  (yhat_rw = fore_rw$mean)
  
  bias_rw = yhat_rw - y.true
  (rw.RMSE = rbind(rw.RMSE, bias_rw^2)); 
  
  # Naive SRW
  #----------------------------------------------------------------------------
  print("Fitting SRW")
  fore_srw = y.window %>% snaive() %>% forecast(h = n.ahead) #%>% autoplot()
  (yhat_srw = fore_srw$mean)
  
  # Model ensemble 1:
  #------------------
  # print("Fitting ENSEMBLE1")
  # fore_ensem <- y.window %>% 
  #                 hybridModel(models = "fnstz", weights="equal", errorMethod = "MAE") %>% 
  #                 forecast(h = n.ahead) 
  # (yhat_ensem = fore_ensem$mean)
  # 
  # bias_ens = yhat_ensem - y.true
  # (ens.RMSE = rbind(ens.RMSE, bias_ens^2)); 

  # Dynamic optimized theta model:
  #---------------------------------
  #print("Fitting Optimized Theta model")
  #dynopt_theta <- dotm(y.window, h=n.ahead, s='additive')
  #yhat_dtheta = dynopt_theta$mean

  
  # Prophet:
  #--------------------
  print("Fitting PROPHET")
  df = data.frame(ds = as.yearmon(time(y.window)), y = y.window)
  m = prophet(df, growth = 'linear', n.changepoints = 3)
  future = make_future_dataframe(m, periods = n.ahead, freq = 'month')
  forecast = predict(m, future)
  fitt_train = forecast[train_index, 'yhat']
  (fore_pro = forecast[valid_index, 'yhat'])
  
  bias_pro = fore_pro - y.true
  (pro.RMSE = rbind(pro.RMSE, bias_pro^2)); 
  
  
  # Model Ensemble 2:
  #-------------------
  print("Fitting ENSEMBLE2")
  fore_comb = cbind(a = fitt_train,
                    b=fitted(fore_nnetar), 
                    c=fitted(fore_tbats), 
                    e = fitted(fore_arfima))
  
  new_preds = cbind(a = fore_pro, 
    b = yhat_nnetar, 
    c = yhat_tbats, 
    #d = yhat_dtheta,
    e = yhat_arfima)

  input_data <- foreccomb(observed_vector = y.window, 
                          prediction_matrix = fore_comb,
                          newobs = y.true, newpreds = new_preds)
  #comb <- comb_OLS(input_data)
  comb <- rolling_combine(input_data, "comb_OLS")
  yhat_ensem2 = comb$Forecasts_Test
  #plot(comb)
  #colMeans(comb$Weights)
  

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
  metrics = c(ETS = accuracy(fore_ets, y.true)["Test set",c("MAPE")],
    ARIMA = accuracy(fore_arima, y.true)["Test set",c("MAPE")],
    ARFIMA = accuracy(fore_arfima, y.true)["Test set",c("MAPE")],
    `STLM-ETS` = accuracy(fore_stlm, y.true)["Test set",c( "MAPE")],
    #`STL-ETS` = accuracy(fore_stl, y.true)["Test set",c("MAPE")],
    #`STLF-ETS` = accuracy(fore_stlf, y.true)["Test set",c( "MAPE")],
    NNAR = accuracy(fore_nnetar, y.true)["Test set",c( "MAPE")],
    TBATS = accuracy(fore_tbats, y.true)["Test set",c( "MAPE")],
    RW_B1 = accuracy(fore_rw, y.true)["Test set",c( "MAPE")],
    SRW_B2 = accuracy(fore_srw, y.true)["Test set",c( "MAPE")],
    Naive_Avg_B3 = accuracy(meanf(y.window, h = n.ahead), y.true)["Test set",c( "MAPE")],
    LM = accuracy(fore_lm, y.true)["Test set",c( "MAPE")],
    #OptThe = accuracy(yhat_dtheta, y.true)["Test set",c( "MAPE")],
    #ENSEM = accuracy(fore_ensem, y.true)["Test set",c( "MAPE")],
    ENSEM2 = accuracy(yhat_ensem2, y.true)["Test set",c( "MAPE")],
    PROPHET = accuracy(fore_pro, y.true)["Test set",c("MAPE")]
    #BSTS = accuracy(fore_bsts, y.true)["Test set",c("MAPE")]
    )
  
  all_metrics = rbind(all_metrics, metrics)
  j=j+1
}
(mean_acc = colMeans(all_metrics)) 
#---------------------------------------------------------------------------------------

all_metrics

(pro_RMSE = sqrt(colMeans(pro.RMSE))) 



plot(fore_ensem)

autoplot(y_comp) +
  autolayer(fore_ets, series="ETS", PI=FALSE) +
  autolayer(fore_arima, series="ARIMA", PI=FALSE) +
  autolayer(fore_stl, series="STL", PI=FALSE) +
  autolayer(fore_nnetar, series="NNAR", PI=FALSE) +
  autolayer(fore_tbats, series="TBATS", PI=FALSE) +
  autolayer(fore_ensem, series="Ensemble") +
  xlab("months") + ylab("demand") +
  ggtitle("Demand forecast next six months")


# par(mfrow=c(1,2))
# 
# k = 1:n.ahead; ylim=c(min(sarma.rmse,bpar.rmse,sdum.rmse),max(sarma.rmse,bpar.rmse,sdum.rmse));
# par(las=0,cex.main=1.3,font.main=11)
# plot(k,sarma.rmse,type="l",xlim=c(1.5,n.ahead),ylim=ylim,ylab="RMSE",xlab="Forecast horizon (in months)",lty=2)
# #plot(k,colMeans(pmse12),type="l",main="",xlab="forecast horizon 'k'",ylab="",sub=expression(S==12))     #single PMSE's for each single horizon k
# #points(k,colMeans(pmse12),type="h",lty=3,col="gray");
# #points(k,cumsum(colMeans(pmse12)),type="l",lty=1,pch=20,col="red4")
# #title("Single PMSEs of a monthly PAR(1) model",sub=expression(S==12))
# lines(k,bpar.rmse,lty=1)
# lines(k,sdum.rmse,lty=3)
# lines(k,bpar1.RMSE_cond,lty=4)
# legend("bottomright",legend=c("BMA","SARMA","Dummy","BPAR(1)"),lty=c(1,2,3,4),bty="n",lwd=1)
# 
# graphics.off() 




# Partition into training and test set:
#---------------------------------------
TT_comp = length(y_comp)
nTest = hstep  
nTrain = TT_comp - nTest

(y_train = subset(y_comp, start = 1, end = TT_comp - nTest))

#-----Plot univariate time series----------------------------------------------
dev.off()
par(mfrow=c(2,2))
plot.ts(y_comp, main= paste("Series ", rownames(mts_data)[i] ), col = "black", lty = 1, ylab="demand") ;
#points(y_comp, col="red")
hist(y_comp,nclass=30, main="", col="grey", xlab="demand")
acf(y_comp, lag.max = 50, main="ACF")
pacf(y_comp, lag.max = 50, main = "PACF")
#------------------ Differences--------------------------------------------------
plot(ts(diff(y_comp)), main= paste("First diff. ", rownames(mts_data)[i] ), col = "black", lty = 1, ylab="demand") ;
points(diff(y_comp), col="black", bg="orange", pch=21)
#plot(density(data_compl[i,],bw=.3),main="Density est.",lwd=2)
hist(diff(y_comp),nclass=30, main="", col="grey", xlab="demand")
acf(ts(diff(y_comp)), lag.max = 50, main="ACF")
pacf(ts(diff(y_comp)), lag.max = 50, main = "PACF")
#------------------ Differences logs--------------------------------------------------
plot(log(1 + y_comp), main= paste("Log(1 + x): ", rownames(mts_data)[i] ), col = "black", lty = 1, ylab="demand") ;
points(log(1 + y_comp), col="black", bg="orange", pch=21)
#plot(density(data_compl[i,],bw=.3),main="Density est.",lwd=2)
hist(log(1 + y_comp),nclass=30, main="", col="grey", xlab="demand")
acf(log( 1+ y_comp), lag.max = 50, main="ACF")
pacf(log(1 + y_comp), lag.max = 50, main = "PACF")
#----------------------------------------------------------------------------






