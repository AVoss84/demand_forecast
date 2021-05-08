
library(forecast)
library(ggplot2)
# ARFIMA forecasts
library(fracdiff)

(train <- log(1 + y_train))
(train = y_train)

# ETS forecasts
train %>%
  naive() %>%
  forecast() %>%
  autoplot()

# ETS forecasts
train %>%
  ets() %>%
  forecast() %>%
  autoplot()

# Automatic ARIMA forecasts
train %>%
  auto.arima() %>%
  forecast(h = nTest) %>%
  autoplot()

arfima(train) %>%
  forecast(h = nTest) %>%
  autoplot()

# Forecasting with STL
train %>%
  stlm(modelfunction=ar) %>%
  forecast(h = nTest) %>%
  autoplot()

train %>%
  stlf(lambda=0) %>%
  autoplot()

train %>%
  stl(s.window='periodic') %>%
  forecast(h = nTest) %>%
  autoplot()

# TBATS forecasts
train %>%
  tbats() %>%
  forecast(h = nTest) %>%
  autoplot()

train %>%
  nnetar(lambda=0) %>%
  forecast(h = nTest) %>%
  autoplot()
  
plot(decompose(train))

fit <- nnetar(sunspotarea, lambda=0)
autoplot(forecast(fit,h=30))

#-----------------------------------------------------------
etsfc <- train %>% 
            ets() %>% 
            forecast(h = nTest)

baggedfc <- train %>% 
              baggedETS() %>% 
              forecast(h = nTest)

autoplot(train) +
  autolayer(baggedfc, series="BaggedETS", PI=FALSE) +
  autolayer(etsfc, series="ETS", PI=FALSE) +
  guides(colour=guide_legend(title="Forecasts"))


# Same by hand:
sim <- bld.mbb.bootstrap(train, 10) %>%
           as.data.frame() %>%
           ts(frequency=12, start=2000)

fc <- purrr::map(as.list(sim),
                 function(x){forecast(ets(x))[["mean"]]}) %>%
            as.data.frame() %>%
            ts(frequency=12, start=start)

autoplot(train) +
  autolayer(sim, colour=TRUE) +
  autolayer(fc, colour=TRUE) +
  autolayer(train, colour=FALSE) +
  ylab("Bootstrapped series") +
  guides(colour="none")

