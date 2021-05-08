
print("Read in model_selection.")


model_selection = function(y_obs, p_oob = 0.85, n.ahead = 6)
{
  (T_obs = length(y_obs))
  (R = floor(T_obs * p_oob))                       #the higher n, the more accurate the forecasts! 
  
  # Rolling window approach to compute out-of-bag error:
  #-------------------------------------------------------------------
  j = 0 ; all_metrics = ets.RMSE = arima.RMSE = arfima.RMSE = tbats.RMSE = NULL;
  stlm.RMSE = fnn.RMSE = rw.RMSE = pro.RMSE = ens.RMSE = NULL
  
  while(j<(T_obs-n.ahead-R)) {
    
    train_index = 1:(R+j)
    valid_index = (R+j+1):(R+j+n.ahead) 
    (y.window = ts(y_obs[train_index], fre=frequency(y_obs),start=start(y_obs)))       #used subsample
    (y.true = y_obs[valid_index])                             #actual values for comparison 
    
    if(any(is.na(y.true))){
      warning('NAs in y.true')
    }
    
    # Impute missing values and smooth outliers:
    #--------------------------------------------
    #(y.window = tsclean(y.window))
    (y.window = na.interp(y.window))
    
    cat("\nTest error run nr.",j+1," (out of ",T_obs-n.ahead-R,") [Training size R: ",length(y.window),"]\n------------------------------------------------------------------\n",sep=""); 
    cat("\nTraining batch (index) ", min(train_index), max(train_index), "\n")
    cat("Validation batch (index) ", min(valid_index), max(valid_index),"\n\n")
    
    # Exponential smoothing state space model:
    #-------------------------------------------
    fore_ets = run_ets(y.window, n.ahead, nmse = n.ahead)
    (yhat_ets = as.numeric(fore_ets$forecast$mean))
    
    bias_ets = yhat_ets - y.true
    (ets.RMSE = rbind(ets.RMSE, bias_ets^2)); 
    
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
    
    
    # tryCatch({
    #   fore_arfima = run_arfima(y.window, n.ahead, estim = c("ls"))
    #   (yhat_arfima = as.numeric(fore_arfima$forecast$mean))
    #   
    #   bias_arfima = yhat_arfima - y.true
    #   (arfima.RMSE = rbind(arfima.RMSE, bias_arfima^2)); 
    #   
    # }, warning = function(w) {
    #   print(w)
    # }, error = function(e) {
    #   print(e)}
    # )
    
    # Forecasting with STL
    #------------------------
    # stlm = run_stlm_ets(y.window, n.ahead, modelfunction=ar, allow.multiplicative.trend = T)
    # (yhat_stlm  = as.numeric(stlm$forecast$mean))
    # 
    # bias_stlm = yhat_stlm - y.true
    # (stlm.RMSE = rbind(stlm.RMSE, bias_stlm^2));

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
    #Tbats = run_tbats(y.window, n.ahead, seasonal.periods = frequency(y.window))
    #(yhat_tbats = as.numeric(Tbats$forecast$mean))

    #bias_tbats = yhat_tbats - y.true
    #(tbats.RMSE = rbind(tbats.RMSE, bias_tbats^2));
    
    # Feed-forward neural networks with a single hidden layer and lagged inputs
    #----------------------------------------------------------------------------
    fnn = run_fnn(y.window, n.ahead)
    (yhat_nnetar = as.numeric(fnn$forecast$mean))
    
    bias_fnn = yhat_nnetar - y.true
    (fnn.RMSE = rbind(fnn.RMSE, bias_fnn^2));
    
    # Linear time series model:
    #---------------------------
    #fit_lm <- tslm(y.window ~ trend + season)
    #fore_lm = forecast(fit_lm, h = n.ahead)
    #(yhat_lm = fore_lm$mean)
    
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
    tryCatch({
      rwd = run_rwd(y.window, n.ahead)
      (yhat_rwd = rwd$forecast$mean)  
    }, warning = function(w) {
      print(w)
    }, error = function(e) {
      print(e)
    })
    
    bias_rw = yhat_rwd - y.true
    (rw.RMSE = rbind(rw.RMSE, bias_rw^2)); 
    
    # Naive SRW
    #----------------------------------------------------------------------------
    #print("Fitting SRW")
    #fore_srw = y.window %>% snaive() %>% forecast(h = n.ahead) 
    #(yhat_srw = fore_srw$mean)
    
    # Model ensemble 1:
    #------------------
    ensem = run_ensem(y.window, n.ahead,models = "efn", weights="equal", errorMethod = "RMSE")
    (yhat_ensem = ensem$forecast$mean)
    
    bias_ens = yhat_ensem - y.true
    (ens.RMSE = rbind(ens.RMSE, bias_ens^2)); 
    
    # Dynamic optimized theta model:
    #---------------------------------
    #print("Fitting Optimized Theta model")
    #dynopt_theta <- dotm(y.window, h=n.ahead, s='additive')
    #yhat_dtheta = dynopt_theta$mean
    
    # Prophet:
    #--------------------
    forec = run_prophet(y.window, n.ahead, 
                        growth = 'linear',
                        daily.seasonality = F,
                        weekly.seasonality = F,
                        n.changepoints = 5)
    
    fitt_train = forec$forecast[train_index, 'yhat']
    yhat_pro = forec$forecast[valid_index, 'yhat']
    
    bias_pro = yhat_pro - y.true
    (pro.RMSE = rbind(pro.RMSE, bias_pro^2)); 
    

    # Random Forest: only 1-step ahead
    #-----------------------------------
    # lags = 6
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
    metrics = c(ETS = accuracy(yhat_ets, y.true)["Test set",c("MAPE")],
      #ARIMA = accuracy(fore_arima, y.true)["Test set",c("MAPE")],
      #ARFIMA = accuracy(yhat_arfima, y.true)["Test set",c("MAPE")],
      #`STLM-ETS` = accuracy(yhat_stlm, y.true)["Test set",c( "MAPE")],
      
      #`STL-ETS` = accuracy(fore_stl, y.true)["Test set",c("MAPE")],
      #`STLF-ETS` = accuracy(fore_stlf, y.true)["Test set",c( "MAPE")],
      NNAR = accuracy(yhat_nnetar, y.true)["Test set",c( "MAPE")],
      #TBATS = accuracy(yhat_tbats, y.true)["Test set",c( "MAPE")],
      RWD = accuracy(yhat_rwd, y.true)["Test set",c( "MAPE")],
      #SRW_B2 = accuracy(fore_srw, y.true)["Test set",c( "MAPE")],
      #Naive_Avg = accuracy(meanf(y.window, h = n.ahead), y.true)["Test set",c( "MAPE")],
      #LM = accuracy(fore_lm, y.true)["Test set",c( "MAPE")],
      #OptThe = accuracy(yhat_dtheta, y.true)["Test set",c( "MAPE")],
      ENSEM = accuracy(yhat_ensem , y.true)["Test set",c( "MAPE")],
      PROPHET = accuracy(yhat_pro, y.true)["Test set",c("MAPE")]
      #BSTS = accuracy(fore_bsts, y.true)["Test set",c("MAPE")]
    )
    all_metrics = rbind(all_metrics, metrics)
    j=j+1
  }
  
  (mean_acc = colMeans(all_metrics)) 
  
  nstep_rmse = list(pro_RMSE = sqrt(colMeans(pro.RMSE)),
                  ets.RMSE = sqrt(colMeans(ets.RMSE)),
                  #tbats.RMSE = sqrt(colMeans(tbats.RMSE)),
                  #stlm.RMSE = sqrt(colMeans(stlm.RMSE)),
                  fnn.RMSE = sqrt(colMeans(fnn.RMSE)),
                  rw.RMSE = sqrt(colMeans(rw.RMSE)),
                  ens.RMSE = sqrt(colMeans(ens.RMSE))
                  )
  
  print("Finished!")
  cat("Best model:", names(mean_acc)[which.min(mean_acc)],"\n")
  
  rownames(all_metrics) = paste0("i=",1:nrow(all_metrics))
  
  return(list(all_metrics = all_metrics, 
              nstep_rmse = nstep_rmse, 
              mean_acc = mean_acc, 
              winner = names(mean_acc)[which.min(mean_acc)]))
}
#------------------------------------------------------------------------------





