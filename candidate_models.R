
print("Read in models.")

# Prophet:
run_prophet = function(y.window, n.ahead,...){
  print("Fitting PROPHET")
  df = data.frame(ds = as.yearmon(time(y.window)), y = y.window)
  m = prophet(df, ...)
  future = make_future_dataframe(m, periods = n.ahead, freq = 'month')
  forecast = predict(m, future)
  return(list(fitted_model = m, data_new = future, forecast = forecast))
}


# FNN:
run_fnn = function(y.window, n.ahead,...){
  print("Fitting NNETAR")
  #xreg = fourier(y.window, K = 2)    # fourier pairs
  #xreg_new = fourier(y.window, K = 2, h = n.ahead)    # fourier pairs future
  fit_nnetar = y.window %>% nnetar(...) 
  fore_nnetar = fit_nnetar %>% forecast(h = n.ahead) 
  return(list(fitted_model = fit_nnetar, data_new = NULL, forecast = fore_nnetar))
}

# Ensemble:
run_ensem = function(y.window, n.ahead,...){
    print("Fitting ENSEMBLE")
    fit_ensem = y.window %>% hybridModel(...)
    fore_ensem = fit_ensem %>% forecast(h = n.ahead) 
    return(list(fitted_model = fit_ensem, data_new = NULL, forecast = fore_ensem))
}

  
# ARFIMA:
run_arfima = function(y.window, n.ahead,...){
  print("Fitting ARFIMA")
  #xreg = fourier(y.window, K = 2)    # fourier pairs
  #xreg_new = fourier(y.window, K = 2, h = n.ahead)    # fourier pairs future
  fit_arfima = y.window %>% arfima(...) 
  fore_arfima = fit_arfima %>% forecast(h = n.ahead) 
  return(list(fitted_model = fit_arfima, data_new = NULL, forecast = fore_arfima))
}


# ETS
run_ets = function(y.window, n.ahead,...){
  print("Fitting ETS")
  fit_ets = y.window %>% ets(...) 
  fore_ets = fit_ets %>% forecast(h = n.ahead) #%>% autoplot()
  return(list(fitted_model = fit_ets, data_new = NULL, forecast = fore_ets))
}


# Forecasting with STLM-ETS
run_stlm_ets = function(y.window, n.ahead,...){
  print("Fitting STLM-ETS")
  fit_stlm = y.window %>% stlm(...) 
  fore_stlm = fit_stlm %>% forecast(h = n.ahead) #%>% autoplot()
  return(list(fitted_model = fit_stlm, data_new = NULL, forecast = fore_stlm))
}


# TBATS : Exponential smoothing state space model 
# with Box-Cox transformation, ARMA errors, 
# Trend and Seasonal componentsforecasts
run_tbats = function(y.window, n.ahead,...){
  print("Fitting TBATS")
  fit_tbats = y.window %>% tbats(...) 
  fore_tbats = fit_tbats %>% forecast(h = n.ahead)
  return(list(fitted_model = fit_tbats, data_new = NULL, forecast = fore_tbats))
}



# RwD:
run_rwd = function(y.window, n.ahead,...){
  print("Fitting RW")
  #fore_rw = y.window %>% naive() %>% forecast(h = n.ahead) #%>% autoplot()
  fit_rwd = y.window %>% rwf(drift = T) 
  fore_rwd = fit_rwd %>% forecast(h = n.ahead) 
  return(list(fitted_model = fit_rwd, data_new = NULL, forecast = fore_rwd))
}


