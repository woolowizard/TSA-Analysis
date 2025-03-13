######################################### BEGIN IMPORT LIB & FUNCTIONS #####################################
file_path = "/Users/andreacommodi/Desktop/Commodi_Urzi/" # Master folder with all resources needed 
source(paste0(file_path, "Functions/TSA-Predict-Student-Functions.R"))
source(paste0(file_path, "Functions/CalendarEffects-Student-Functions.R"))
source(paste0(file_path, "Functions/TSA-Useful-Functions.R"))
library(forecast)
library(lmtest)       ## For better model print
library(tsoutliers)   ## For outliers
library(urca)         ## For UR tests
library(uroot)        ## For UR test
library(FinTS)        ## For ArchTest (from RForge)
library(lubridate) 
######################################### END IMPORT LIB & FUNCTIONS #####################################

######################################### BEGIN SCRIPT #####################################
#### BEGIN DATA PREPARATION #####################################
file.data = paste0(file_path, "df/NGRCUUS.csv")
name <- "Prezzo del gas naturale per i consumatori residenziali"
unit <- "Dollari"

#### Read data
data <- read.table(file = file.data, header = TRUE, sep = ",", 
                   na.strings = "NA")
data = data[!grepl("13$", as.character(data$YYYYMM)), ]
data = data[as.numeric(substr(data$YYYYMM, 1, 4)) > 1984, ]
rownames(data) = NULL
date <- as.Date(paste0(data$YYYYMM, "01"), format = "%Y%m%d") # Aggiungo "01" come giorno accanto ad ogni YYYY-MM
g.transf = "log"
y = data$Value  
y_log = log(as.numeric(y)) # Cast y from double --> int
y = as.numeric(y)
y = ts(data=y, frequency = 12) # Creazione oggetto TS; frequency è s --> stagionalità
yor = y
y_log = ts(data=y_log, frequency = 12)
#### END DATA PREPARATION #####################################

#### Ts plot, acf, pacf of the original series #######
x=y

#### Lag con la regola di Jenkins
t=dim(data)[1]
h=round(sqrt(t)+45) 

plot(x=date, y=x, type = "l", main=name, ylab=unit)
Acf(x = x, type = "correlation", na.action = na.pass, lag.max = h, main = name)
Acf(x = x, type = "partial", na.action = na.pass, lag.max = h, main = name)

#### Ts plot, acf, pacf of the diff(12) series
x_1 <- diff(y, lag=1)
x_fill_1 <- c(rep(NA, 1), x_1)

x_12 = diff(y, lag=12)
x_fill_12 <- c(rep(NA, 12), x_12)

## Plot differenze prime
plot(x = date, y = x_fill_1, type = "l", main = name, ylab = unit)
Acf(x = x_1, type = "correlation", na.action = na.pass, lag.max = h, main = name)
Acf(x = x_1, type = "partial",     na.action = na.pass, lag.max = h, main = name)

## Plot differenze dodicesime
plot(x = date, y = x_fill_12, type = "l", main = name, ylab = unit)
Acf(x = x_12, type = "correlation", na.action = na.pass, lag.max = h, main = name)
Acf(x = x_12, type = "partial",     na.action = na.pass, lag.max = h, main = name)

## Model selection with x = y (no trasformazioni)
xreg <- NULL
model_base <- Arima(
  y,
  order = c(1, 0, 1), 
  seasonal = list(order = c(1, 1, 2), period = 12), 
  xreg = xreg, 
  include.constant = FALSE
)
.print.arima(x = model_base)

res1  <- residuals(model_base)    

## Plot residui + ArchTest
main <- "residuals^2"
x1 <- res1^2
plot(x1, type = "l", main = main, xlab = "", ylab = "")
Acf(x = x1, type = "correlation", lag.max = 60, na.action = na.pass, main = main)
ArchTest (x = res1^2, lags = 12, demean = FALSE) 

main <- "|residuals|"
x1 <- abs(res1)
plot(x1, type = "l", main = main, xlab = "", ylab = "")
Acf(x = x1, type = "correlation", lag.max = 60, na.action = na.pass, main = main)
ArchTest (x = sqrt(abs(res1)), lags = 12, demean = FALSE)


###
.loglik(model_base, g.transf)
###### SERIE LOGARTMICA PLOT ########
#### Ts plot, acf, pacf of the original series
x=y_log

#### Lag con la regola di Jenkins
t=dim(data)[1]
h=round(sqrt(t)+45) 

plot(x=date, y=x, type = "l", main=name, ylab=unit)
Acf(x = x, type = "correlation", na.action = na.pass, lag.max = h, main = name)
Acf(x = x, type = "partial", na.action = na.pass, lag.max = h, main = name)

#### Ts plot, acf, pacf of the diff(12) series
x_1 <- diff(y_log, lag=1)
x_fill_1 <- c(rep(NA, 1), x_1)

x_12 = diff(y_log, lag=12)
x_fill_12 <- c(rep(NA, 12), x_12)

## Plot differenze prime
plot(x = date, y = x_fill_1, type = "l", main = name, ylab = unit)
Acf(x = x_1, type = "correlation", na.action = na.pass, lag.max = h, main = name)
Acf(x = x_1, type = "partial",     na.action = na.pass, lag.max = h, main = name)

## Plot differenze dodicesime
plot(x = date, y = x_fill_12, type = "l", main = name, ylab = unit)
Acf(x = x_12, type = "correlation", na.action = na.pass, lag.max = h, main = name)
Acf(x = x_12, type = "partial",     na.action = na.pass, lag.max = h, main = name)

########### Model selection #################
xreg <- NULL
model_base <- Arima(
  y_log,
  order = c(0, 0, 0), 
  seasonal = list(order = c(0, 0, 0), period = 12), 
  xreg = xreg, 
  include.constant = FALSE
)
.print.arima(x = model_base)
###
.loglik(model_base, g.transf)

########### Model selection #################
xreg <- NULL
model_base <- Arima(
  y_log,
  order = c(0, 0, 0), 
  seasonal = list(order = c(0, 1, 0), period = 12), 
  xreg = xreg, 
  include.constant = FALSE
)
.print.arima(x = model_base)
###
.loglik(model_base, g.transf)

########### Model selection #################
xreg <- NULL
model_base <- Arima(
  y_log,
  order = c(1, 0, 0), 
  seasonal = list(order = c(0, 1, 0), period = 12), 
  xreg = xreg, 
  include.constant = FALSE
)
.print.arima(x = model_base)
###
.loglik(model_base, g.transf)

########### Model selection #################
xreg <- NULL
model_base <- Arima(
  y_log,
  order = c(1, 0, 0), 
  seasonal = list(order = c(0, 1, 1), period = 12), 
  xreg = xreg, 
  include.constant = FALSE
)
.print.arima(x = model_base)
###
.loglik(model_base, g.transf)

########### Model selection #################
xreg <- NULL
model_base <- Arima(
  y_log,
  order = c(1, 0, 0), 
  seasonal = list(order = c(1, 1, 1), period = 12), 
  xreg = xreg, 
  include.constant = FALSE
)
.print.arima(x = model_base)
###
.loglik(model_base, g.transf)

########### Model selection #################
xreg <- NULL
model_base <- Arima(
  y_log,
  order = c(1, 0, 1), 
  seasonal = list(order = c(1, 1, 1), period = 12), 
  xreg = xreg, 
  include.constant = FALSE
)
.print.arima(x = model_base)
###
.loglik(model_base, g.transf)

########### Model selection #################
xreg <- NULL
model_base <- Arima(
  y_log,
  order = c(1, 0, 1), 
  seasonal = list(order = c(1, 1, 2), period = 12), 
  xreg = xreg, 
  include.constant = FALSE
)
.print.arima(x = model_base)
###
.loglik(model_base, g.transf)

########### Model selection #################
xreg <- NULL
model_base <- Arima(
  y_log,
  order = c(2, 1, 1), 
  seasonal = list(order = c(0, 0, 1), period = 12), 
  xreg = xreg, 
  include.constant = FALSE
)
.print.arima(x = model_base)
###
.loglik(model_base, g.transf)

########### Model selection #################
xreg <- NULL
model_base <- Arima(
  y_log,
  order = c(2, 0, 2), 
  seasonal = list(order = c(0, 1, 1), period = 12), 
  xreg = xreg, 
  include.constant = FALSE
)
.print.arima(x = model_base)
###
.loglik(model_base, g.transf)

########### Model selection #################
xreg <- NULL
model_base <- Arima(
  y_log,
  order = c(1, 1, 0), 
  seasonal = list(order = c(1, 1, 0), period = 12), 
  xreg = xreg, 
  include.constant = FALSE
)
.print.arima(x = model_base)
###
.loglik(model_base, g.transf)

########### Model selection #################
xreg <- NULL
model_base <- Arima(
  y_log,
  order = c(1, 1, 0), 
  seasonal = list(order = c(1, 1, 1), period = 12), 
  xreg = xreg, 
  include.constant = FALSE
)
.print.arima(x = model_base)
###
.loglik(model_base, g.transf)

########### Model selection #################
xreg <- NULL
model_base <- Arima(
  y_log,
  order = c(1, 1, 1), 
  seasonal = list(order = c(1, 1, 1), period = 12), 
  xreg = xreg, 
  include.constant = FALSE
)
.print.arima(x = model_base)
###
.loglik(model_base, g.transf)

########### Model selection #################
xreg <- NULL
model_base <- Arima(
  y_log,
  order = c(1, 1, 2), 
  seasonal = list(order = c(1, 1, 1), period = 12), 
  xreg = xreg, 
  include.constant = FALSE
)
.print.arima(x = model_base)
###
.loglik(model_base, g.transf)

########### Model selection #################
xreg <- NULL
model_base <- Arima(
  y_log,
  order = c(1, 1, 3), 
  seasonal = list(order = c(1, 1, 1), period = 12), 
  xreg = xreg, 
  include.constant = FALSE
)
.print.arima(x = model_base)
###
.loglik(model_base, g.transf)

########### Model selection Normale #################
xreg <- NULL
model_base <- Arima(
  y_log,
  order = c(2, 1, 3), 
  seasonal = list(order = c(1, 1, 1), period = 12), 
  xreg = xreg, 
  include.constant = FALSE
)
.print.arima(x = model_base)
###
.loglik(model_base, g.transf)

#### residuals #############################################

## Useful quantities
res1  <- residuals(model_base)                          
resst1 <- scale(res1)   
res_2 = scale(res1^2)

### Plot

main <- "residuals"
x1 <- res1
par(mfrow = c(3, 1))  # 2 righe, 2 colonne
plot(x1, type = "l", main = main, xlab = "", ylab = "")
Acf(x = x1, type = "correlation", lag.max = h, na.action = na.pass, main = main)
Acf(x = x1, type = "partial",     lag.max = h, na.action = na.pass, main = main)

########### Model selection CALENDAR #################

cal <- .calendarEffects(time = date, country = "us")
ind1 <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun", "sh", "lh")
ind2 <- c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun", "sh")
ind3 <- c("wi")
cal <- cal[, ind3, drop = FALSE]
cal <- as.matrix(cal)

xreg <- cal
model_calendar <- Arima(y_log,
                        order = c(1, 0, 1), seasonal = list(order = c(1, 1, 2)),
                        xreg = xreg, include.constant = FALSE)
.print.arima(x = model_calendar)
.loglik(model_calendar, g.transf)

##### Diagnostiche ########################

fit = model_base
res1 <- residuals(fit)                          ## Residuals
resst1 <- scale(res1)

main <- "residuals"
x1 <- res1
Acf(x = x1, type = "correlation", lag.max = h, na.action = na.pass, main = main)
Acf(x = x1, type = "partial",     lag.max = h, na.action = na.pass, main = main)

npar1 <- NROW(fit$coef)                       ## Number of parameters
fitdf1 <- 0                                   ## If we want to remove np from df
lag1  <- fitdf1+c(1, 2, 5, 10, 15, 20)           ## lag
lb <- mapply(FUN=Box.test, lag=lag1,
             MoreArgs = list(x = x1, type = "Ljung-Box", fitdf = 0))[1:3, , drop = FALSE]
rbind(lag = lag1, lb)

####### Eteroschedasticità + ARCH #######

main <- "|residuals|"
x1 <- abs(res1)
Acf(x = x1, type = "correlation", lag.max = h, na.action = na.pass, main = main)
ArchTest (x = sqrt(abs(res1)), lags = 12, demean = FALSE)

###### QQ-Plot e Shapiro###
shapiro.test(x = res1)
hist(x = resst1, breaks = 25, freq = FALSE, main = "residuals", xlab = "")
x1 <- seq(from = min(resst1), to = max(resst1)+1, length.out = 100)
lines(x = x1, y = dnorm(x = x1, mean = 0, sd = 1), col = "red")
qqnorm(y = resst1, main = "Normal Q-Q Plot",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE)
abline(a = 0, b = 1, col = "red")

###### Previsione Ex-Post #####
col.list <- c("brown", "violet", "orange")
model.list <- c("ARIMA", "ARIMA + calendar", "Naive") 

y=y_log

J <- 12 # Numero di obs che tengo fuori
H <- 1# Orizzonte della previsione

t1 <- .predict.t1(nobs = NROW(y), J = J, n.ahead = H)

# Arima Only
xreg <- NULL
pred1.1 <- .predict(object = model_base, n.ahead = H, t = t1, y = y, xreg = xreg,
                    fixed.n.ahead = TRUE)

# Arima + Calendar
xreg <- cal
pred2.1 <- .predict(object = model_calendar, n.ahead = H, t = t1, y = y, xreg = xreg,
                    fixed.n.ahead = TRUE)

# Naive
predn.1 <- .predict.naive(fit = model_base, J = J, n.ahead = H, g = g.transf)

x1 <- .pred.bands(pred = pred1.1, alpha = 0.05, g = g.transf)
x2 <- .pred.bands(pred = pred2.1, alpha = 0.05, g = g.transf)

em1 <- .ErrorMeasures(y = yor, fit = x1$mean, naive = predn.1)
em2 <- .ErrorMeasures(y = yor, fit = x2$mean, naive = predn.1)
emn <- .ErrorMeasures(y = yor, fit = predn.1, naive = predn.1)

data.frame(
  model = model.list, 
  rbind(em1, em2, emn))

ind <- x1$t
time1 <- date[ind]
ylim <- range(x1$lower, x1$upper, x2$lower, x2$upper, predn.1)
plot(x = time1, y = yor[ind], ylim = ylim, xlab = "", ylab = "Price", main = "Ex-Post")
lines(x = time1, y = x1$lower, col = col.list[1], lty = 2)
lines(x = time1, y = x1$mean, col = col.list[1])
lines(x = time1, y = x1$upper, col = col.list[1], lty = 2)
lines(x = time1, y = x2$mean, col = col.list[2])
lines(x = time1, y = predn.1, col = col.list[3])
legend(x = "bottomleft", fill = NULL, legend = model.list, 
       col = col.list, lty = 1, border = FALSE, text.col = col.list)

#### Previsioni Ex-Ante ###########
y=y_log

H <- 12            ## horizon
t1 <- NROW(y)      ## origin

# Arima Only
pred1 <- .predict(object = model_base, n.ahead = H, t = t1, y = y, xreg = NULL, 
                  fixed.n.ahead = FALSE)

# ARIMA + Calendar
time1 <- .extend.time(x = date, n.ahead = H, by = "month")
cal1 <- .calendarEffects(time = time1, country = "it")
xreg <- rbind(cal, cal1[, colnames(cal), drop = FALSE])
xreg <- as.matrix(xreg)
pred2 <- .predict(object = model_calendar, n.ahead = H, t = t1, y = y, xreg = xreg,
                  fixed.n.ahead = FALSE)

# Naive
predn <- .predict.naive(fit = model_base, J = 0, n.ahead = H, g = g.transf)

x1 <- .pred.bands(pred = pred1, alpha = 0.05, g = g.transf)
x2 <- .pred.bands(pred = pred2, alpha = 0.05, g = g.transf)

ylim <- range(x1$lower, x1$upper, x2$lower, x2$upper, predn)
plot(x = time1, y = x1$mean, type = "l", col = col.list[1], ylim = ylim, xlab = "", ylab = "Price", main = "Ex-Ante")
lines(x = time1, y = x1$lower, col = col.list[1], lty = 2)
lines(x = time1, y = x1$upper, col = col.list[1], lty = 2)
lines(x = time1, y = x2$mean, col = col.list[2])
lines(x = time1, y = predn, col = col.list[3])
par(mfrow = c(1,1), mar = c(2, 4, 2, 0.5))
legend(x = "bottomleft", fill = NULL, legend = model.list, 
       col = col.list, lty = 1, border = FALSE, text.col = col.list)

### Outliers

outlier <- tso(y_log, xreg = NULL, cval = 6, delta = 0.7, 
               types = c("AO", "LS", "TC"), 
               maxit = 10, maxit.iloop = 100, maxit.oloop = 10,
               tsmethod = "arima", 
               args.tsmethod = list(order = c(1, 0, 1), seasonal = list(order = c(1, 1, 2))))
.plot.tso(outlier)            ## Use .plot.tso(fit) if plot(fit) stops with an error 

### Test Dicky Fuller
adf.1 <- ur.df(y_log, type = "trend", lags = 24, selectlags = "AIC")
data.frame(t(adf.1@teststat), adf.1@cval, check.names = FALSE)

adf.2 <- ur.df(y_log, type = "drift", lags = 24, selectlags = "AIC")
data.frame(t(adf.2@teststat), adf.2@cval, check.names = FALSE)

adf.3 <- ur.df(y_log, type = "none", lags = 24, selectlags = "AIC")
data.frame(t(adf.3@teststat), adf.3@cval, check.names = FALSE)

x1_diff <- diff(y_log, 1)
## Start from drift: in diff() for sure there is no trend
adf.1 <- ur.df(y = x1_diff, type = "drift", lags = 24, selectlags = "AIC")
data.frame(t(adf.1@teststat), adf.1@cval, check.names = FALSE)

## Start from drift: in diff() for sure there is no trend
x12_diff <- diff(y_log, 12)
adf.1 <- ur.df(y = x12_diff, type = "drift", lags = 24, selectlags = "AIC")
data.frame(t(adf.1@teststat), adf.1@cval, check.names = FALSE)