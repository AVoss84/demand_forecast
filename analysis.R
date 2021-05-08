
rm(list=ls())

{library(dplyr)
library(lubridate)
#library(prophet)
library(forecast)
library(forecastHybrid)
#library(ForecastCombinations)
library(ForecastComb)
library(ggplot2)
library(urca)
library(dlm)
library(skimr)
library(xts)
library(zoo)
#require(nortest)
library(MARSS)
library(bsts)
library(prophet)
#library(fUnitRoots)
#library(TSA)
#library(caret)  
#library(vars)  
#library(BVAR)  
#library(fpp)
#library(forecTheta)
 library(tidyverse) 
}
  
#install.packages('tidyverse', repos = 'https://cloud.r-project.org')

setwd("C:\\Users\\Alexander\\Documents\\Arbeit\\Siemens\\ts_stuff\\data\\")

#----------------
# Read dataset:
#----------------
data = read.csv("DataChallenge_2019-07.csv")

data$Date = as.Date(data$Date) 

View(data)

glimpse(data)

#uni_count = unique(data$Country)
#uni_prod = unique(data$Product)

data_cleaned = data %>% mutate(Demand_cleaned = factor(ifelse(Demand < 0 , NA, Demand)))

dgr = data_cleaned %>% group_by(Country, Product) %>% 
               summarise(n = n(), start_date = min(Date), end = max(Date)) %>% 
               arrange(Country, Product) ;dgr
#View(dgr)

unique(dgr$n)


# Use series with full observation period: T = 78
#---------------------------------------------------
fil78 = dgr %>% filter(n == 78) %>% 
          filter((Country != "Country_1") & (Product != "Product_15")) 
fil78
#View(fil78)

dim(fil78)

#fil78 = dgr %>% filter(n == 78) %>% filter(Product == "Product_10")

# Use series with shorter observation period: T = 42
#----------------------------------------------------
fil42= dgr %>% filter(n == 42) %>% 
           filter((Country != "Country_1") & (Product != "Product_15"))

fil42
dim(fil42)

#my_ts = data_cleaned %>% filter((Country == 'Country_0') & (Product == i)) %>% 
#            filter(complete.cases(.)) 

### PREPARE DATA:
###################

#------ Only Country/Product series ------------------------------------
#-------------------------------------------------------------------------
#mts_data = matrix(nr=nrow(fil78), ncol=78)
mts_data = matrix(nr=nrow(fil42), ncol=42)
my_cols = ts_all = NULL
for(i in 1:nrow(mts_data)){
  
  print(i) 
  #(coun = fil78$Country[i]) 
  #(pro = fil78$Product[i])
  
  coun = fil42$Country[i] 
  pro = fil42$Product[i]

  my_cols = c(my_cols, paste0(coun, "-", pro))
  #cat("Country: ",   print(coun), " and Product: ", print(pro))  
  (my_ts = data_cleaned %>% filter((Country == coun) & (Product == pro))) #%>% filter(complete.cases(.)) 
  #my_ts = xts(my_ts$Demand, order.by = my_ts$Date)
  
  # Caused problem!!!!!
  (mts_data[i,] = as.numeric(as.matrix(my_ts$Demand_cleaned)[,1]))
  
  (ts_all = ts.union(ts_all, ts1 = ts(mts_data[i,], start = c(year(min(my_ts$Date)), 
                                                             month(min(my_ts$Date))), 
                                     end = c(year(max(my_ts$Date)), 
                                             month(max(my_ts$Date))), deltat = 1/12)))
}

dim(mts_data)
dim(ts_all)

date_index = my_ts$Date
colnames(ts_all) = my_cols
rownames(mts_data) = my_cols
#--------------------------------------------------------------



######################################
## Given country - product combination
######################################

#layout(matrix(c(1,1,2,3), byrow = T, nr=2, nc=2))
dev.off()
par(mfrow=c(2,2))
par(ask=T)

data_compl = mts_data[,-c(37:42)]
data_compl_new = ts_all[-c(37:42),]   # without missings 
#data_compl = mts_data[,-c(73:78)]
#data_compl_new = ts_all[-c(73:78),]   # without missings 
for(i in 1:nrow(mts_data)){
  
  print(i)
  plot.ts(data_compl_new[,i], main= paste("Series ", colnames(data_compl_new)[i]), col = "black", lty = 1, ylab="demand") ;
  points(data_compl_new[,i], col="red")
  hist(data_compl_new[,i],nclass=30, main="", col="grey", xlab="demand")
  acf(na.omit(data_compl_new[,i]), lag.max = 50, main="ACF")
  pacf(na.omit(data_compl_new[,i]), lag.max = 50, main = "PACF")
  
  #-----------------------
  # Unit root testing:
  #----------------------
  #(TT = length(data_compl_new[,i]))
  #(lags = trunc(TT^(1/3)))					    #Vorschlag von Said/Dickey verwenden
  #adf = ur.df(data_compl_new[,i], type="trend",lags=lags)		#vergleiche tau
  #summary(ur.df(data_compl_new[,i], type="trend",lags=lags,select="BIC"))
  
  # Case 1
  #print(ifelse(adf@teststat[1] < adf@cval[1,2], yes = "Not reject H0! (Nonstat.)", no = "Reject H0!"))		#H0: "pi=0" nicht ablehnen
  
  # Case 4
  #print(ifelse(adf@teststat[3] > adf@cval[3,2], yes = "Not reject H0! (Nonstat.)", no = "Reject H0!"))		#H0: "pi=0" nicht ablehnen
  #------------------------------------------------------  
      
  #plot(ts(diff(data_compl_new[,i])), main= paste("First diff. ", rownames(data_compl)[i]), col = "black", lty = 1, ylab="demand") ;
  #points(diff(data_compl_new[,i]), col="black", bg="orange", pch=21)
  #hist(diff(data_compl_new[,i]),nclass=30, main="", col="grey", xlab="demand")
  #acf(ts(diff(data_compl_new[,i])), lag.max = 50, main="ACF")
  #pacf(ts(diff(data_compl_new[,i])), lag.max = 50, main = "PACF")

  plot(ts(log(1 + data_compl_new[,i])), main= paste("Log(1 + x): ", rownames(mts_data)[i] ), col = "black", lty = 1, ylab="demand") ;
  points(log(1 + data_compl_new[,i]), col="black", bg="orange", pch=21)
  #plot(density(data_compl[i,],bw=.3),main="Density est.",lwd=2)
  hist(log(1 + data_compl_new[,i]),nclass=30, main="", col="grey", xlab="demand")
  acf(na.omit(log( 1+ data_compl_new[,i])), lag.max = 50, main="ACF")
  pacf(na.omit(log(1 + data_compl_new[,i])), lag.max = 50, main = "PACF")
}
graphics.off() ;par(ask=F)


## Save data as matrix and as multivar. ts obj.:
#-------------------------------------------------
saveRDS(data, "raw_data.rds")
saveRDS(data_cleaned, "cleaned_data.rds")
saveRDS(ts_all, "data_mts_78.rds")
saveRDS(mts_data, "data_mat_78.rds")
#saveRDS(ts_all, "data_mts_42.rds")
#saveRDS(mts_data, "data_mat_42.rds")

#########################################################################################
######################### END ############################################################

# other code

Box.test(residuen, lag = 20, type="Ljung-Box")

tsdiag(model3)   				#Visualisierung Ergebnisse der Modelldiagnose

qqnorm(data.used)		#indiziert rechtsschiefe Verteilung!
qqline(data.used)

jarque.bera.test(data.used)   

lillie.test(data.used)			#KS-Test für Test auf NV
shapiro.test(data.used)				#also nicht annäherend nv!


(model3 = arima(data,order=c(0,1,1))) 

(pred.all = predict(model3, n.ahead = hstep)$pred)   			#h- Schritt Prognose

par(mfrow=c(2,1),lwd=1, lty=1, cex.main=1.5, font.main=4)

(result1 = plot(model3,n.ahead=hstep, ylab='',xlab='Year', pch=19,col="blue4"))    #siehe auch transform-Argument in plot.Arima!

data.frame(orig=data,ind=c(rep("#",length(data.used)),rep("!",hstep)),preds=c(data.used,result1$pred))

# Daten rücktransformieren -> Differenzenbildung!
# Yhat_t+k = Wt+k + Yhat_t+k-1  (vgl. Chan, S.209)
#------------------------------------------------#

(theta = as.double(model3$coef[1]))

(Wt = result1$pred)
(yt = data[used])				
(yhat = yt + theta*res.hat[length(res.hat)])		# Initialisierung

#(yhat2 = Wt[2]+ yhat)
#(yhat3 = Wt[3]+ yhat2)

for(k in 2:hstep){(yhat =  c(yhat, Wt[k] + yhat[k-1]))}
yhat

plot(1:(fcst+hh),c(data.used,yhat), type="l", xlab="", ylab="",ylim=c(1,6))
lines(1:fcst,data, type="l", xlab="", ylab="CPI",ylim=c(1,6))
lines((used+1):(fcst+hh),yhat,col="red3",lwd=2,lty=1)
abline(v=fcst,lty=3)
(outsample.f = (fcst+1):(fcst+hh))
points(x=outsample.f,y=yhat[outsample.f-100],col="blue4",pch=1)
title("h - Schritt Prognose")

# Prognosegüte checken:

(ergeb = data.frame(data=data[(used+1):fcst], diff.pred=pred.all[1:(hstep-hh)], diffback=yhat[1:(hstep-hh)]))

( MSFE = mean((ergeb[,"data"]-ergeb[,"diffback"])^2) )     		#Mean Square (Forecast) Error

# Running MSFE plotten:

Ns = 1:(hstep-hh)
run.mse = ((ergeb[,"data"]-ergeb[,"diffback"])^2)/Ns

dev.off()
plot(Ns,run.mse,col="red4",type="o",,xlab="Vorhersagehorizont",lwd=2,main="Laufender mittlerer\n quadratischer Vorhersagefehler (MSFE)",font.main=13)



