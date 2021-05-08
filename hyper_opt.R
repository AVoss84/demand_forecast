
print("Read in hyper_optim.")

require(rBayesianOptimization)


#----------------------------------------------------------------------------------------
hyper_optim = function(winner)
{
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



