
print("Read in Bayesian Optimization.")

# FNN
#-------------------------------------------------------------------------------
bo_optim_fnn = function(p_nonsea, nof_nodes) {
  
  n.ahead = 6
  
  # Rolling window approach to compute out-of-bag error:
  #--------------------------------------------------------------------------
  j = 0 ; error = c()
  while(j<(T_obs-n.ahead-R)) 
  {
    train_index = 1:(R+j)
    valid_index = (R+j+1):(R+j+n.ahead) 
    (y.window = ts(y_obs[train_index],fre=frequency(y_obs),start=start(y_obs)))       #used subsample
    (y.true = y_obs[valid_index])                             #actual values for comparison 
    
    cat("\nTest error run nr.",j+1," (out of ",T_obs-n.ahead-R,") [Training size R: ",length(y.window),"]\n------------------------------------------------------------------\n",sep=""); 
    
    # FNN
    #---------------------------------------------------------------------
    #xreg = fourier(y.window, K = 2)    # fourier pairs
    #xreg_new = fourier(y.window, K = 2, h = n.ahead)    # fourier pairs future
    fore_nnetar = y.window %>% nnetar(p=p_nonsea, P=1, size=nof_nodes) %>% 
                                    forecast(h = n.ahead) 
    #fore_nnetar = y.window %>%
    #               nnetar(p=1, P=1, size=nof_nodes) %>%
    #               forecast(h = n.ahead)
    
    (yhat_nnetar = fore_nnetar$mean)
    error_d = forecast::accuracy(yhat_nnetar, y.true)["Test set",c("MAPE")]
    #------------------------------------------------------------------------
    error = c(error, error_d)
    j=j+1
  }
  return(list(Score = -mean(error), Pred = 0))
}
#----------------------------------------------------------------------------------


# Prophet - scoring/objective function to be maximized
#-------------------------------------------------------------------------------
bo_optim_prophet = function(changepoint_prior_scale, 
                            seasonality_prior_scale, 
                            n_changepoints) 
{
  assign('y_obs', y_obs, envir = global_env())
  
  (T_obs = length(y_obs))
  (R = floor(T_obs * .85))                       #the higher n, the more accurate the forecasts! 
  n.ahead = 6
  
  # Rolling window approach to compute out-of-bag error:
  #--------------------------------------------------------------------------
  j = 0 ; error = c()
  while(j<=(T_obs-n.ahead-R)) 
  {
    train_index = 1:(R+j)
    valid_index = (R+j+1):(R+j+n.ahead) 
    (y.window = ts(y_obs[train_index],fre=frequency(y_obs),start=start(y_obs)))       #used subsample
    (y.true = y_obs[valid_index])                             #actual values for comparison 
    
    cat("\nTest error run nr.",j," (out of ",T_obs-n.ahead-R,") [Training size R: ",length(y.window),"]\n------------------------------------------------------------------\n",sep=""); 
    
    # Prophet
    #--------------------
    df = data.frame(ds = as.yearmon(time(y.window)), y = y.window)
    m = prophet(df, growth = 'linear',
                seasonality.prior.scale = seasonality_prior_scale,
                changepoint.prior.scale = changepoint_prior_scale,
                n.changepoints = n_changepoints,
                weekly.seasonality = F,
                daily.seasonality = F)

    future = make_future_dataframe(m, periods = n.ahead, freq = 'month')
    forecast = predict(m, future)
    fore_pro = forecast[valid_index, 'yhat']
    error_d = forecast::accuracy(fore_pro, y.true)["Test set",c("MAPE")]
    #---------------------------------------------------------------------
    error = c(error, error_d)
    j=j+1
  }
  return(list(Score = -mean(error)))
}
#----------------------------------------------------------------------------------


