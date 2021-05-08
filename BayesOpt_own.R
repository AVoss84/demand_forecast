
library(rBayesianOptimization)
library(ggplot2)

i = 40

n.ahead = 6
(TT = ncol(mts_data))

(y = ts(mts_data[i,], start = start(ts_all), end = end(ts_all), frequency = frequency(ts_all)))
(y_obs = subset(y, start = 1, end = TT-n.ahead))   # observed data

(y_obs = log(1 + y_obs))
y_obs = na.interp(y_obs)

plot(y_obs)

(T_obs = length(y_obs))
(R = floor(T_obs * .90))                       #the higher n, the more accurate the forecasts! 

#-------------------------------------------------------------------------------
#prophet_fit_bayes = function(changepoint_prior_scale, seasonality_prior_scale, n_changepoints) {
fit_bayes = function(p_nonsea, nof_nodes) {
  
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
    
    # FNN
    #---------------------------------------------------------------------
    xreg = fourier(y.window, K = 2)    # fourier pairs
    xreg_new = fourier(y.window, K = 2, h = n.ahead)    # fourier pairs future
    fore_nnetar = y.window %>% nnetar(p=p_nonsea, P=1, size=nof_nodes, xreg = xreg) %>% 
                            forecast(h = n.ahead, xreg = xreg_new) 
    #fore_nnetar = y.window %>%
    #               nnetar(p=1, P=1, size=nof_nodes) %>%
    #               forecast(h = n.ahead)
    (yhat_nnetar = fore_nnetar$mean)
    error_d = forecast::accuracy(yhat_nnetar, y.true)["Test set",c("MAPE")]
    #------------------------------------------------------------------------
    
    # Prophet
    #--------------------
    # df = data.frame(ds = as.yearmon(time(y.window)), y = y.window)
    # m = prophet(df, growth = 'linear',
    #             seasonality.prior.scale = seasonality_prior_scale,
    #             changepoint.prior.scale = changepoint_prior_scale,
    #             n.changepoints = n_changepoints,
    #             weekly.seasonality = F,
    #             daily.seasonality = F)
    # 
    # future = make_future_dataframe(m, periods = n.ahead, freq = 'month')
    # forecast = predict(m, future)
    # fore_pro = forecast[valid_index, 'yhat']
    # error_d = forecast::accuracy(fore_pro, y.true)["Test set",c("MAPE")]
    #---------------------------------------------------------------------
    
    error = c(error, error_d)
    j=j+1
  }
return(list(Score = -mean(error), Pred = 0))
}
#----------------------------------------------------------------------------------

#out = fit_bayes(nof_nodes = 100)


grid_len = 10

rand_search_grid =  data.frame( 
  nof_nodes          = sample(20:60, grid_len, replace = F),
  p_nonsea          = sample(1:11, grid_len, replace = F),
  Value                   = rep(0, grid_len)
  ) ; rand_search_grid

bayesian_search_bounds = list(nof_nodes = range(rand_search_grid$nof_nodes),
                              p_nonsea = range(rand_search_grid$p_nonsea))

# Run:
ba_search = BayesianOptimization(fit_bayes,
                                 bounds = bayesian_search_bounds,
                                 init_grid_dt = rand_search_grid, 
                                 init_points = 1, 
                                 n_iter = 15,
                                 acq = 'ucb', # 'ei'
                                 kappa = 1, 
                                 eps = 0,
                                 verbose = TRUE)

#---------------------------- Prophet ---------------------------------------------------

rand_search_grid =  data.frame(
          changepoint_prior_scale = sort(runif(10, 0.01, 0.1)),
          seasonality_prior_scale = c(sort(sample(c(runif(5, 0.01, 0.05), runif(5, 1, 10)), 5, replace = F)),
                                      sort(sample(c(runif(5, 0.01, 0.05), runif(5, 1, 10)), 5, replace = F))),
          n_changepoints          = sample(1:2, 2, replace = F),
          Value                   = rep(0, 10))

changepoint_bounds    = range(rand_search_grid$changepoint_prior_scale)
n_changepoint_bounds  = as.integer(range(rand_search_grid$n_changepoints))
seasonality_bounds    = range(rand_search_grid$seasonality_prior_scale)

bayesian_search_bounds = list(changepoint_prior_scale = changepoint_bounds,
                              seasonality_prior_scale = seasonality_bounds,
                              n_changepoints = as.integer(n_changepoint_bounds))

# Run:
ba_search = BayesianOptimization(prophet_fit_bayes,
                                 bounds = bayesian_search_bounds,
                                 init_grid_dt = rand_search_grid, 
                                 init_points = 1, 
                                 n_iter = 15,
                                 acq = 'ucb', 
                                 kappa = 1, 
                                 eps = 0,
                                 verbose = TRUE)

best_params_ba  = c(ba_search$Best_Par, Value = -1*ba_search$Best_Value)

df = data.frame(ds = as.yearmon(time(y.window)), y = y.window)

# Conditional on optimized hyper parameter - retrain the model:
m = prophet(df, growth = 'linear',
            seasonality.prior.scale = best_params_ba[['seasonality_prior_scale']], 
            changepoint.prior.scale = best_params_ba[['changepoint_prior_scale']],
            n.changepoints = best_params_ba[['n_changepoints']])


future = make_future_dataframe(m, periods = n.ahead, freq = 'month')

forecast = predict(m, future)
#forecast$ds = as.Date(forecast$ds)
forecast[valid_index, 'yhat']

# p = ggplot() + 
#   geom_point(data = subset(cv_set, ds >= origin_week - days(7*52)), aes(x = as.Date(ds), y = y), size = 1) +
#   geom_line(data = subset(forecast, ds >= origin_week - days(7*52)), aes(x = as.Date(ds), y = yhat), color = "#0072B2", size = 1) +
#   geom_ribbon(data = subset(forecast, ds >= origin_week - days(7*52)), aes(x = as.Date(ds), ymin = yhat_lower, ymax = yhat_upper), fill = "#0072B2", alpha = 0.3) +
#   geom_point(data = test, aes(x = as.Date(ds), y = y), size = 1, color = '#4daf4a') ;p

