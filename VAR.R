
library(vars)

data_compl = mts_data[,-c(73:78)]
dim(data_compl)

X = t(data_compl)
dim(X)

X = apply(X, 2, na.interp)
X

X_z = X

X = apply(dat, 2, function(x) diff(log(1+x)))
X = apply(dat, 2, function(x) diff(x))

colnames(X) = paste0("V_",1:ncol(X))

head(X, 20)

X[,73]

X_z = as.data.frame(X) %>% mutate_at(.funs = funs(scale),.vars=colnames(X))

dim(X_z)

X_z[,1:5]

fit = VAR(y = as.matrix(X_z[,1:5]), type="const", lag.max = 1, ic = "SC")   #compare

#print(coef(fit));

BIC(fit)

## Checking the roots:
##---------------------##
(roots <- roots(fit, modulus=T))  #moduli of eigenvalues of the companion matrix A (=VAR(1)-form). ev < 1 ?


## Forecasting objects of class varest
##-------------------------------------##
# args(vars:::predict.varest)

(predictions <- predict(out, n.ahead = 6, ci = 0.95))

attributes(predictions)

## Plot of predictions for inc.adj
##----------------------------------##
plot(predictions)
#plot(predictions, names = "con1")

## or use Fanchart:
##-----------------##
fanchart(predictions)





## Testing serial correlation
##----------------------------##
## Residual Autocorrelations and multivariate Portmanteau-Test:
## H_{0}: "Residuals obey model's assumptions", i.e. are White Noise

dim(vardat)    #Sample size!! -> PT.adjusted for finite sample correction

(var2c.serial = serial.test(varsimest, lags.pt = 16,type = "PT.adjusted"))    #see LP, p.169,510

plot(var2c.serial, names = "inv1")
plot(var2c.serial, names = "inc1")
plot(var2c.serial, names = "con1")

## Testing heteroscedasticity:
##-----------------------------##
# args(arch.test)
(var2c.arch <- arch.test(varsimest, lags.multi = 5, multivariate.only = T))   #ARCH-LM test of H0: no arch effects

## Testing for normality (siehe Kapitel 3):
##------------------------------------------##
(var2c.norm <- normality.test(varsimest, multivariate.only = T))   #siehe Mardias Schiefe und Wölbungsmaß, Folie 9, Kapitel 3



## Testing serial correlation
##----------------------------##
## Residual Autocorrelations and multivariate Portmanteau-Test:
## H_{0}: "Residuals obey model's assumptions", i.e. are White Noise

dim(vardat)    #Sample size!! -> PT.adjusted for finite sample correction

(var2c.serial = serial.test(varsimest, lags.pt = 16,type = "PT.adjusted"))    #see LP, p.169,510

plot(var2c.serial, names = "inv1")
plot(var2c.serial, names = "inc1")
plot(var2c.serial, names = "con1")

## Testing heteroscedasticity:
##-----------------------------##
# args(arch.test)
(var2c.arch <- arch.test(varsimest, lags.multi = 5, multivariate.only = T))   #ARCH-LM test of H0: no arch effects

## Testing for normality (siehe Kapitel 3):
##------------------------------------------##
(var2c.norm <- normality.test(varsimest, multivariate.only = T))   #siehe Mardias Schiefe und Wölbungsmaß, Folie 9, Kapitel 3

# U = t(residuals(varsimest)); dim(U)
# N = ncol(U)
# (u.mean = U%*%rep(1,N)/N)                                  #residual mean vector
# (U.mean = matrix(rep(U.mean,N),byrow=F,nc=N,nr=nrow(U)))    #matrix of repeated means
# (S = (U-U.mean)%*%t(U-U.mean)/(N-1))             #sample covariance matrix of residuals, siehe Folie 5, Kapitel 2
# var(residuals(varsimest))                   #check

## class and methods for diganostic tests
# class(var2c.serial)
# class(var2c.arch)
# class(var2c.norm)
# methods(class = "varcheck")
# 
# ## Plot of objects "varcheck"
# args(vars:::plot.varcheck)
#----Check for structural change------------------------------------------------------#
#(reccusum <- stability(varsimest,type = "OLS-CUSUM"))
#(fluctuation <- stability(varsimest, type = "fluctuation"))
#-----------------------------------------------------------#

## Causality tests (Wald tests, p.102 ff.)
## Granger and Instantaneous causality:
##------------------------------------------##
(var.causal <- causality(varsimest, cause = c("inc1","con1"), boot=F, boot.runs=1000))       # Income/Cons. -> Invest., see p.104

# -> there is instantaneous causality between (Income/Cons.) and Invest. and vice versa (symmetry!). Feedback system!

#round(var(residuals(varsimest)),6)                   #estimated residual covariance matrix

## Forecasting objects of class varest
##-------------------------------------##
# args(vars:::predict.varest)

(predictions <- predict(varsimest, n.ahead = 8, ci = 0.95))

attributes(predictions)

## Plot of predictions for inc.adj
##----------------------------------##
plot(predictions)
#plot(predictions, names = "con1")

## or use Fanchart:
##-----------------##
fanchart(predictions)


## Transform predicted first differences back ('integrate') 
## to (log)level form:
##----------------------------------------------------------##

# Investments:
(fore.inv = diffinv(predictions$fcst$inv1[,"fcst"],xi=inv[length(inv)]))
(inv.low = diffinv(predictions$fcst$inv1[,"lower"], xi = inv[length(inv)] ))
(inv.up = diffinv(predictions$fcst$inv1[,"upper"], xi = inv[length(inv)] ))

# Income:
(fore.inc = diffinv(predictions$fcst$inc1[,"fcst"],xi=inc[length(inc)]))
(inc.low = diffinv(predictions$fcst$inc1[,"lower"], xi = inc[length(inc)] ))
(inc.up = diffinv(predictions$fcst$inc1[,"upper"], xi = inc[length(inc)] ))

# Consumption:
(fore.con = diffinv(predictions$fcst$con1[,"fcst"],xi=con[length(con)]))
(con.low = diffinv(predictions$fcst$con1[,"lower"], xi = con[length(con)] ))
(con.up = diffinv(predictions$fcst$con1[,"upper"], xi = con[length(con)] ))

graphics.off() ; par(mfrow=c(3,1),las=1,cex.main=1.5)

plot(ts(c(inv,fore.inv),start=start(inv),fre=4),type="l",lty=3,col="red",ylab="",main="Invest.",ylim=c(min(inv),max(inv.up))); 
lines(inv,lwd=2) ; 
lines(ts(c(inv,inv.low),start=start(inv),fre=4),lty=2)
lines(ts(c(inv,inv.up),start=start(inv),fre=4),lty=2)

plot(ts(c(inc,fore.inc),start=start(inc),fre=4),type="l",lty=3,col="red",ylab="",main="Income",ylim=c(min(inc),max(inc.up))); 
lines(inc,lwd=2)
lines(ts(c(inc,inc.low),start=start(inc),fre=4),lty=2)
lines(ts(c(inc,inc.up),start=start(inc),fre=4),lty=2)

plot(ts(c(con,fore.con),start=start(con),fre=4),type="l",lty=3,col="red",ylab="",main="Consump.",ylim=c(min(con),max(con.up))); 
lines(con,lwd=2)
lines(ts(c(con,con.low),start=start(con),fre=4),lty=2)
lines(ts(c(con,con.up),start=start(con),fre=4),lty=2)

# class(predictions)
# args(vars:::plot.varprd)

## Impulse response analysis:
##----------------------------##

#round(var(residuals(varsimest)),6)                   #estimated residual covariance matrix, see also p.75 (3.2.19)

(ir.fct <- irf(varsimest, impulse = c("inc1","con1"), response = "inv1", n.ahead = 10,
               ortho = T, cumulative = F))       #cumulative: accumulated effect over several periods of a shock in one variable, see p.55 below

plot(ir.fct)                #impulse resp. are zero if one of the variables does not granger cause the other variables taken as a group, p.54

(ir.fct <- irf(varsimest, n.ahead = 10,ortho = T, cumulative = F))

plot(ir.fct)
