
library(BVAR)

# Access a subset of the fred_qd dataset and transform it to be stationary
data("fred_qd")
data <- fred_qd[, c("CPIAUCSL", "UNRATE", "FEDFUNDS")]
head(data)

data[5:nrow(data), 1] <- diff(log(data[, 1]), lag = 4) * 100
  
data <- data[5:nrow(data), ]
dim(data)

data_compl = mts_data[,-c(73:78)]
dim(data_compl)
dat = t(data_compl)
dim(dat)

dat = apply(dat, 2, na.interp)

dat = apply(dat, 2, function(x) diff(log(1+x)))

head(dat, 20)

dat[,73]

#any(is.na(dat))

# Compute VAR using 2 lags and a ridiculously low number of draws
x <- bvar(data = dat, lags = 1,
  n_draw = 500, n_burn = 400, n_thin = 2, verbose = FALSE)

#x <- bvar(data, lags = 2, irf = bv_irf(), fcast = bv_fcast())

#bv_fcast(horizon = n.ahead, conditional = FALSE)

fitted(x, conf_bands = 0.10)

logLik(x)

# Only get the density of the marginal likelihood
density(x, vars = "ml")

# Check out some of the outputs generated
plot(x)

# Only plot the marginal likelihood's density
plot(x, "dens", "ml")

# Adjust, compute and store a longer forecast
x$fcast <- predict(x, horizon = n.ahead)

dim(x$fcast$fcast)
x$fcast$quants

# Lower draws, use `bv_fcast()` to set options and add confidence bands
predict(x, bv_fcast(n.ahead), n_thin = 10L, conf_bands = c(0.05, 0.16))

# Use new data to calculate a prediction
predict(x, newdata = matrix(rnorm(200), ncol = 2))

# Get a summary of the last saved forecast
summary(x)
# Limit the summary to variable #2
summary(x, vars = 2L)


plot(predict(x))

irf(x)
plot(irf(x))



