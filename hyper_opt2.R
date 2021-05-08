

print("Read in hyper_optim.")

require(rBayesianOptimization)

#----------------------------------------------------------------------------------------
hyper_optim = function(y_obs, p_oob = 0.8, n.ahead = 6, winner)
{
  #==============================================================
  
  (T_obs = length(y_obs))
  (R = floor(T_obs * p_oob))                       #the higher n, the more accurate the forecasts! 
  
  # FNN
  #-------------------------------------------------------------------------------
  bo_optim_fnn = function(p_nonsea, nof_nodes) {
    
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

      (yhat_nnetar = fore_nnetar$mean)
      error_d = forecast::accuracy(yhat_nnetar, y.true)["Test set",c("MAPE")]
      #------------------------------------------------------------------------
      error = c(error, error_d)
      j=j+1
    }
    return(list(Score = -mean(error)))
  }
  #----------------------------------------------------------------------------------
  
  
  # Prophet - scoring/objective function to be maximized
  #-------------------------------------------------------------------------------
  bo_optim_prophet = function(changepoint_prior_scale, 
                              seasonality_prior_scale, 
                              n_changepoints) 
  {
    #assign('y_obs', y_obs, envir = global_env())

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
  
  #==============================================================
  if(winner == 'PROPHET'){
    
    print("Optimize Prophet hyperparameter.")
    
    rand_search_grid =  data.frame(
      changepoint_prior_scale = sort(runif(10, 0.01, 0.1)),
      seasonality_prior_scale = c(sort(sample(c(runif(5, 0.01, 0.05), runif(5, 1, 10)), 5, replace = F)),
                                  sort(sample(c(runif(5, 0.01, 0.05), runif(5, 1, 10)), 5, replace = F))),
      n_changepoints          = sample(1:5, 5, replace = F),
      Value                   = rep(0, 10))
    
    changepoint_bounds    = range(rand_search_grid$changepoint_prior_scale)
    n_changepoint_bounds  = as.integer(range(rand_search_grid$n_changepoints))
    seasonality_bounds    = range(rand_search_grid$seasonality_prior_scale)
    
    bayesian_search_bounds = list(changepoint_prior_scale = changepoint_bounds,
                                  seasonality_prior_scale = seasonality_bounds,
                                  n_changepoints = as.integer(n_changepoint_bounds))
    
    #tryCatch({
    ba_search = BayesianOptimization(bo_optim_prophet,
                                     bounds = bayesian_search_bounds,
                                     init_grid_dt = rand_search_grid, 
                                     init_points = 5, 
                                     n_iter = 6,
                                     acq = 'ucb', # 'ei'
                                     kappa = 1, 
                                     #kern = "Gaussian",  # "Matern52"
                                     eps = 0,    # tuning parameter for aquisition function
                                     verbose = TRUE)
    
    best_params_ba  = c(ba_search$Best_Par, Value = -1*ba_search$Best_Value)
    
    # }, warning = function(w) {
    #   print(w)
    # }, error = function(e) {
    #   print(e)
    # })
    
  }
  #------------------------------------------------------------------------------------
  if(winner == 'NNAR'){
    
    print("Optimize FNN hyperparameter.")
    
    grid_len = 10
    rand_search_grid =  data.frame( 
      nof_nodes          = sample(20:40, grid_len, replace = F),
      p_nonsea          = sample(1:11, grid_len, replace = F),
      Value                   = rep(0, grid_len)
    ) 
    
    bayesian_search_bounds = list(nof_nodes = range(rand_search_grid$nof_nodes),
                                  p_nonsea = range(rand_search_grid$p_nonsea))
    
    # Run:
    ba_search = BayesianOptimization(bo_optim_fnn,
                                     bounds = bayesian_search_bounds,
                                     init_grid_dt = rand_search_grid, 
                                     init_points = 5, 
                                     n_iter = 6,
                                     acq = 'ucb', # 'ucb'
                                     kappa = 1, 
                                     #kern = "Gaussian",  # "Matern52"
                                     eps = 0,
                                     verbose = TRUE)
    
    best_params_ba  = c(ba_search$Best_Par, Value = -1*ba_search$Best_Value)
    
  }
  return(best_params_ba)
}
#--------------------------------------------------------------------------------------------



