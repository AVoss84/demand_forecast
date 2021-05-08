
# https://rpubs.com/tdneumann/351073

library(RSocrata)
library(dplyr)
library(lubridate)

fbi_code = "'05'"
url = sprintf("https://data.cityofchicago.org/resource/6zsd-86xi.json?$select=*&$where=fbi_code=%s", fbi_code)
burglaries = read.socrata(url)

burglaries$date_clean = as.Date(as.POSIXct(burglaries$date, format = '%Y-%m-%d %H:%M:%S'))
burglaries$week       = as.Date(cut(burglaries$date_clean, 'week'))

burglary_counts       = burglaries %>%
  group_by(district, week) %>%
  summarise(COUNTCRIMES = length(district))

head(burglary_counts)

library(prophet)

### rolling cv period
set.seed(12345)
num_weeks   = 6
dist        = '004'
origin_week = as.Date(cut(Sys.Date(), 'week')) - days(7*num_weeks) 
end_week    = as.Date(cut(Sys.Date(), 'week')) - days(7) 
date_seq    = seq(origin_week, end_week, by = '7 days')
date_seq

cv_set  = burglary_counts %>%
  ungroup() %>%
  filter(week <= end_week & district == dist) %>%
  dplyr::select(week, COUNTCRIMES)

cv_set

plot(ts(cv_set$y))

names(cv_set) = c('ds', 'y')

rand_search_grid =  data.frame( 
  changepoint_prior_scale = sort(runif(10, 0.01, 0.1)),
  seasonality_prior_scale = c(sort(sample(c(runif(5, 0.01, 0.05), runif(5, 1, 10)), 5, replace = F)),
                              sort(sample(c(runif(5, 0.01, 0.05), runif(5, 1, 10)), 5, replace = F))),
  n_changepoints          = sample(5:25, 10, replace = F),
  Value                   = rep(0, 10)
)
rand_search_grid 

# Search best parameters
for (i in 1:nrow(rand_search_grid)) {
  parameters = rand_search_grid[i, ]
  error = c()
  for (d in date_seq) {
    train = subset(cv_set, ds < d)
    test  = subset(cv_set, ds == d)
    
    m = prophet(train, growth = 'linear',
                seasonality.prior.scale = parameters$seasonality_prior_scale, 
                changepoint.prior.scale = parameters$changepoint_prior_scale,
                n.changepoints = parameters$n_changepoints,
                weekly.seasonality = F,
                daily.seasonality = F)
    
    future = make_future_dataframe(m, periods = 1, freq = 'week')
    
    # NOTE: There's a problem in function names with library(caret)
    forecast = predict(m, future)
    forecast$ds = as.Date(forecast$ds)
    
    error_d = forecast::accuracy(forecast[forecast$ds %in% test$ds, 'yhat'], test$y)[ , 'MAPE']
    error = c(error, error_d)
  }
  mean_error = mean(error)
  print(mean_error)
  rand_search_grid$Value[i] = -mean_error
}

best_cv_value = -1*rand_search_grid[which.max(rand_search_grid$Value),'Value']
rand_search_grid = arrange(rand_search_grid, -Value)
head(rand_search_grid)



library(rBayesianOptimization)
library(ggplot2)

prophet_fit_bayes = function(changepoint_prior_scale, seasonality_prior_scale, n_changepoints) {
  
  error = c()
  for (d in date_seq) {
    train = subset(cv_set, ds < d)
    test  = subset(cv_set, ds == d)
    
    m = prophet(train, growth = 'linear',
                seasonality.prior.scale = seasonality_prior_scale, 
                changepoint.prior.scale = changepoint_prior_scale,
                n.changepoints = n_changepoints,
                weekly.seasonality = F,
                daily.seasonality = F)
    
    future = make_future_dataframe(m, periods = 1, freq = 'week')
    
    # NOTE: There's a problem in function names with library(caret)
    forecast = predict(m, future)
    forecast$ds = as.Date(forecast$ds)
    
    error_d = forecast::accuracy(forecast[forecast$ds %in% test$ds, 'yhat'], test$y)[ , 'MAPE']
    error = c(error, error_d)
  }
  
  ## The function wants to _maximize_ the outcome so we return 
  ## the negative of the resampled MAPE value. `Pred` can be used
  ## to return predicted values but we'll avoid that and use zero
  list(Score = -mean(error), Pred = 0)
}

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
                                 init_points = 1,  # do not set to zero as in post!
                                 n_iter = 15,
                                 acq = 'ucb', 
                                 kappa = 1, 
                                 eps = 0,
                                 verbose = TRUE)


best_params_ba  = c(ba_search$Best_Par, Value = -1*ba_search$Best_Value)

df = data.frame(ds = as.yearmon(time(y.window)), y = y.window)

m = prophet(df, growth = 'linear',
            seasonality.prior.scale = best_params_ba[['seasonality_prior_scale']], 
            changepoint.prior.scale = best_params_ba[['changepoint_prior_scale']],
            n.changepoints = best_params_ba[['n_changepoints']])


future = make_future_dataframe(m, periods = 1, freq = 'week')

# NOTE: There's a problem in function names with library(caret)
forecast = predict(m, future)
forecast$ds = as.Date(forecast$ds)

p = ggplot() + 
  geom_point(data = subset(cv_set, ds >= origin_week - days(7*52)), aes(x = as.Date(ds), y = y), size = 1) +
  geom_line(data = subset(forecast, ds >= origin_week - days(7*52)), aes(x = as.Date(ds), y = yhat), color = "#0072B2", size = 1) +
  geom_ribbon(data = subset(forecast, ds >= origin_week - days(7*52)), aes(x = as.Date(ds), ymin = yhat_lower, ymax = yhat_upper), fill = "#0072B2", alpha = 0.3) +
  geom_point(data = test, aes(x = as.Date(ds), y = y), size = 1, color = '#4daf4a') ;p



