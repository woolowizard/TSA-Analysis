######################################### BEGIN IMPORT LIB & FUNCTIONS #####################################
file_path = "/Users/andreacommodi/Desktop/Commodi_Urzi/" # Master folder with all resources needed 
source(paste0(file_path, "Functions/TSA-Finance-Functions.R"))
source(paste0(file_path, "Functions/TSA-Predict-Student-Functions.R"))

library(tseries)  
library(sandwich)
library(lmtest)
library(urca)     ## For unit root
library(rugarch)  ## For GARCH models
library(FinTS)    ## For ArchTest (download from RForge)
library(car)
library(forecast) 
library(quantmod) ## For downloading data
library(xts)

### Function specific for this script
.MincerZarnowitz.local <- function(y, fit)
{
  #### Make the test
  x1 <- .MincerZarnowitz(y = y, fit = fit)
  ####
  x2 <- data.frame(
    NA, 
    NA, x1$x$F[2], x1$x$`Pr(>F)`[2], 
    NA, x1$xHC$F[2], x1$xHC$`Pr(>F)`[2], 
    NA, x1$xHAC$F[2], x1$xHAC$`Pr(>F)`[2])
  colnames(x2) <- colnames(x1$coef)
  ####
  coef <- rbind(x1$coef, Ftest = x2)
  #### Select
  ind <- c("estimate", "HAC.s.e.", "HAC.tstat", "HAC.pvalue")
  coef <- coef[, ind, drop = FALSE]
  colnames(coef) <- c("estimate", "HAC.s.e.", "HAC.stat", "HAC.pvalue")
  data.frame(name = rownames(coef), coef, check.names = FALSE, 
             fix.empty.names = FALSE)
}

.DieboldMariano.local <- function(y, f1, f2, h, loss)
{
  #### Make the test
  x1 <- .DieboldMariano(y = y, f1 = f1, f2 = f2, h = h, loss = loss)
  #### Readable output
  data.frame(
    horizon = as.numeric(x1$parameter["Forecast horizon"]),
    loss = loss,
    statistic = as.numeric(x1$statistic))
}
# ------------------------------------------------------------------------------
######################################### END IMPORT LIB & FUNCTIONS #####################################

######################################### BEGIN SCRIPT #####################################
#### BEGIN DATA PREPARATION #####################################
# Get data from local source
file.data = paste0(file_path, "df/CRM20240116.csv")
data <- read.table(file = file.data, header = TRUE, sep = "\t", check.names = FALSE, comment.char = "", na.strings = ".")
head(data)
tail(data)

# Extract period
data$Date = as.Date(x=data$Date)
ind   <- as.Date(x = "2018-06-01") <= as.Date(x = data$Date)
data  <- data[ind, , drop = FALSE]

# Add variables
data <- data.frame(data,
                   cc.ret = c(NA, diff(log(data$Adjusted))),
                   gkVol = .garmanklass(data = data, sd = TRUE),
                   check.names = TRUE)
data <- data[-1, , drop = FALSE]

# Split data
data.tot <- data
ind <- data.tot$Date <= as.Date(x = "2023-12-31")
data <- data.tot[ind, , drop = FALSE]
data.out <- data.tot[!ind, , drop = FALSE]

time  <- as.Date(data$Date, format = "%Y-%m-%d")
yc    <- data$Close
yclog <- log(yc)
y     <- data$Adjusted
ylog  <- log(y)
nobs <- NROW(y)

par(mfrow = c(2,1), mar = c(2.5,2.5,2,0.5))
#plot(x = time, y = yc,    main = "Close",        xlab = "", ylab = "", type = "l")
#plot(x = time, y = yclog, main = "Ln(close)",    xlab = "", ylab = "", type = "l")
plot(x = time, y = y,     main = "AdjClose",     xlab = "", ylab = "", type = "l")
plot(x = time, y = ylog,  main = "Ln(AdjClose)", xlab = "", ylab = "", type = "l")

Acf(x = y, lag.max = 75, type = "correlation", main = "Price")
Acf(x = y, lag.max = 75, type = "partial", main = "Price")
#### END DATA PREPARATION #####################################

### BEGIN UR TEST #####################################
adf.1 <- ur.df(y = y, type = "drift", lags = 5, selectlags = "AIC")
data.frame(t(adf.1@teststat), adf.1@cval, check.names = FALSE)

adf.2 <- ur.df(y = y, type = "none", lags = 5, selectlags = "AIC")
data.frame(t(adf.2@teststat), adf.2@cval, check.names = FALSE)
### END UR TEST #####################################

#### Percentage log-returns
yret = xts(x = 100 * data$cc.ret, order.by = time)
par(mfrow = c(1,1))
plot(x = time, y = yret, main = "Returns", xlab = "", ylab = "", type = "l") 
abline(h = mean(data$cc.ret, na.rm = TRUE), col = "red")
par(mfrow = c(1,1))
Acf(x = yret, lag.max = 75, type = "correlation", main = "Returns")
Acf(x = abs(yret), lag.max = 75, type = "correlation", main = "|Returns|")
Acf(x = yret^2, lag.max = 75, type = "correlation", main = expression(Returns^2))

#### BEGIN LJ-BOX TEST SUI RENDIMENTI #####################################
npar <- 0
lag <- c(2, 5, 10, 15, 20, 30, 50) + npar
lb <- mapply(FUN = Box.test, lag = lag, 
             MoreArgs = list(x = yret, type = "Ljung-Box", fitdf = npar))[1:3,]
rbind(lag = lag, lb)

lag <- c(4, 8, 12, 16)
at <- mapply(FUN = ArchTest, lags = lag, 
             MoreArgs = list(x = yret, demean = TRUE))
print(rbind(lag = lag, at[1:3,]))
#### END LJ-BOX TEST SUI RENDIMENTI #####################################

##### BEGIN PLOT DISTRIBUZIONE DEI RENDIMENTI E QQ-PLOT #######################
.hist(x = yret, xlim = c(-10, 10), n = 200, breaks = 200, main = "Returns")
qqnorm(y = scale(yret))
abline(a = 0, b = 1, col = "red")
##### END PLOT DISTRIBUZIONE DEI RENDIMENTI E QQ-PLOT #######################

############################## BEGIN MODEL SELECTION #######################################
###CONST + Simple Garch v(1,1)####
spec1 <- ugarchspec(
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE,
                    external.regressors = NULL),
  variance.model = list(model = "sGARCH", garchOrder = c(1,1), 
                        submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
  distribution.model = "std")
#### Estimation
fit <- ugarchfit(spec = spec1, data = yret)
#### Copy
fit1 <- fit

#### Information criteria
infocriteria(fit)
#### Estimated coefficients
fit@fit$robust.matcoef
#### Standardized residuals
zres <- fit@fit$z

###CONST + GJR-GARCH(1,1)####
spec2 <- ugarchspec(
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE,
                    external.regressors = NULL),
  variance.model = list(model = "gjrGARCH", garchOrder = c(1,1), 
                        submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
  distribution.model = "std")
#### Estimation
fit <- ugarchfit(spec = spec2, data = yret)
#### Copy
fit2 <- fit

infocriteria(fit2)
#### Estimated coefficients
fit@fit$robust.matcoef
#### Standardized residuals
zres <- fit@fit$z

###CONST + GJR-GARCH(1,1) + VT####
spec2_vt <- ugarchspec(
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE,
                    external.regressors = NULL),
  variance.model = list(model = "gjrGARCH", garchOrder = c(1,1), 
                        submodel = NULL, external.regressors = NULL, variance.targeting = TRUE),
  distribution.model = "std")
#### Estimation
fit_spec2vt <- ugarchfit(spec = spec2_vt, data = yret)
#### Copy
fit2_vc <- fit_spec2vt

infocriteria(fit2_vc)
#### Estimated coefficients
fit2_vc@fit$robust.matcoef
#### Standardized residuals
zres <- fit@fit$z

####CONST + T-Garch(1,1)####
spec3 <- ugarchspec(
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE,
                    external.regressors = NULL),
  variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
                        submodel = "TGARCH", external.regressors = NULL, variance.targeting = FALSE),
  distribution.model = "std")
## Include the line below to estimate a TGARCH without asymmetric effect
# fixed.pars = list(eta11 = 0))
#### Estimation
fit <- ugarchfit(spec = spec3, data = yret)
#### Copy
fit3 <- fit
fit3c <- .fgarch.2.gjr(fit = fit)

infocriteria(fit3)
#### Estimated coefficients
fit@fit$robust.matcoef ## Do not inlude in the report
fit3c$robust.matcoef
#### Standardized residuals
zres <- fit@fit$z

#### CONST + TGARCH(1,1) + VT ####
spec3_vt <- ugarchspec(
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE,
                    external.regressors = NULL),
  variance.model = list(model = "fGARCH", garchOrder = c(1,1), 
                        submodel = "TGARCH", external.regressors = NULL, variance.targeting = TRUE),
  distribution.model = "std")
## Include the line below to estimate a TGARCH without asymmetric effect
# fixed.pars = list(eta11 = 0))
#### Estimation
fit <- ugarchfit(spec = spec3_vt, data = yret)
#### Copy
fit3vt <- fit
fit3_vt <- .fgarch.2.gjr(fit = fit)

infocriteria(fit3vt)
#### Estimated coefficients
fit3_vt$robust.matcoef
#### Standardized residuals
zres <- fit@fit$z

#### Constant + IGARCH(1,1) ####
spec4 <- ugarchspec(
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE,
                    external.regressors = NULL),
  variance.model = list(model = "iGARCH", garchOrder = c(1,1), 
                        submodel = NULL, external.regressors = NULL, variance.targeting = FALSE),
  distribution.model = "std")
#### Estimation
fit <- ugarchfit(spec = spec4, data = yret)
#### Copy
fit4 <- fit

infocriteria(fit4)
#### Estimated coefficients
fit@fit$robust.matcoef
#### Standardized residuals
zres <- fit@fit$z

############################## END MODEL SELECTION #######################################

#### ACF sui residui, residui^2 e |residui| ####
Acf(x = zres, lag.max = 75, type = "correlation", main = "Residuals")
Acf(x = zres^2, lag.max = 75, type = "correlation", main = "Residuals^2")
Acf(x = abs(zres), lag.max = 75, type = "correlation", main = "|Residuals|")

#### Distribuzione e QQ-Plot residui del modello migliore ####
par(mfrow = c(1,2), mar = c(4,4,3,0.5))
.hist.fit(fit = fit3, xlim = c(-5, 5), n = 200, breaks = 100, main = "Residuals")
.qqplot.fit(fit = fit3)

### Test Nyblom sul modello migliore ###
nyblom(fit3)

#### Sign test ####
signbias(fit3)

##### BEGIN EX POST PREDICTION ####

yret.tot <- xts(x = data.tot$cc.ret * 100, order.by = data.tot$Date)
y.tot <- data.tot$gkVol * 100
nobs.tot <- NROW(y.tot)
y <- data$gkVol * 100
y.out <- data.out$gkVol * 100
nobs.out <- NROW(y.out)
H <- 1
n.ahead <- H                     ## The usual h
out.sample <- nobs.out - 1 + H   ## Corresponding to J
n.roll <- nobs.out - 1           ## How many forecasts to do in addition to 1st

# sGarch
fit <- fit1
specx <- getspec(fit)
setfixed(specx) <- as.list(coef(fit))
forc <- ugarchforecast(fitORspec = specx, n.ahead = n.ahead,  
                       data = yret.tot, out.sample = out.sample, n.roll = n.roll)
## Examine forc@forecast
forc1 <- forc@forecast$sigmaFor[H,]

#### gjrGARCH
fit <- fit2
specx <- getspec(fit)
setfixed(specx) <- as.list(coef(fit))
forc <- ugarchforecast(fitORspec = specx, n.ahead = n.ahead,  
                       data = yret.tot, out.sample = out.sample, n.roll = n.roll)
forc2 <- forc@forecast$sigmaFor[H,]

#### TGARCH
fit <- fit3
specx <- getspec(fit)
setfixed(specx) <- as.list(coef(fit))
forc <- ugarchforecast(fitORspec = specx, n.ahead = n.ahead,  
                       data = yret.tot, out.sample = out.sample, n.roll = n.roll)
forc3 <- forc@forecast$sigmaFor[H,]

#### IGARCH
fit <- fit4
specx <- getspec(fit)
setfixed(specx) <- as.list(coef(fit))[-4]
forc <- ugarchforecast(fitORspec = specx, n.ahead = n.ahead,  
                       data = yret.tot, out.sample = out.sample, n.roll = n.roll)
forc4 <- forc@forecast$sigmaFor[H,]

#### TGARCH + VT
fit <- fit3vt
specx <- getspec(fit)
setfixed(specx) <- as.list(coef(fit))
forc <- ugarchforecast(fitORspec = specx, n.ahead = n.ahead,  
                       data = yret.tot, out.sample = out.sample, n.roll = n.roll)
forc5 <- forc@forecast$sigmaFor[H,]

#### GJR + VT
fit <- fit2_vc
specx <- getspec(fit)
setfixed(specx) <- as.list(coef(fit))
forc <- ugarchforecast(fitORspec = specx, n.ahead = n.ahead,  
                       data = yret.tot, out.sample = out.sample, n.roll = n.roll)
forc6 <- forc@forecast$sigmaFor[H,]

##### END EX POST PREDICTION ####

##### BEGIN ERROR MEASURES #######################
naive.vol <- sd(yret) 
naive.var <- naive.vol^2 

ErrorMeas <- data.frame(
  measure = c("Volatility", "Volatility", "Volatility", "Volatility", "Volatility", "Variance", "Variance",
              "Variance", "Variance", "Variance"), 
  model = c("GARCH", "GJR-GARCH", "T-GARCH", "IGARCH", "naive",
            "GARCH", "GJR-GARCH", "T-GARCH", "IGARCH", "naive"), 
  rbind( 
    .ErrorMeasures(y = y.out,   fit = forc1,   naive = naive.vol), 
    .ErrorMeasures(y = y.out,   fit = forc2,   naive = naive.vol), 
    .ErrorMeasures(y = y.out,   fit = forc3,   naive = naive.vol), 
    .ErrorMeasures(y = y.out,   fit = forc4,   naive = naive.vol), 
    .ErrorMeasures(y = y.out,   fit = rep(naive.vol, nobs.out),  naive = naive.vol), 
    .ErrorMeasures(y = y.out^2, fit = forc1^2, naive = naive.var), 
    .ErrorMeasures(y = y.out^2, fit = forc2^2, naive = naive.var), 
    .ErrorMeasures(y = y.out^2, fit = forc3^2, naive = naive.var), 
    .ErrorMeasures(y = y.out^2, fit = forc4^2, naive = naive.var),
    .ErrorMeasures(y = y.out^2, fit = rep(naive.var, nobs.out),  naive = naive.var) ) ) 
print(ErrorMeas, row.names = FALSE)
##### END ERROR MEASURES #######################

#### BEGIN Garman e Klass ####
par(mfrow = c(1,1), lwd = 1)
legend <- c("GK volatility", "GJR-GARCH", "T-GARCH")
col  <- c("black", "red", "blue")
plot(x = time, y = y, type = "l", col = col[1], ylab = "Garman-Klass volatility")
lines(x = time, y = fit2@fit$sigma, col = col[3])
lines(x = time, y = fit3@fit$sigma, col = col[2]) 
legend(x = "topright", y = NULL, legend = legend, border = FALSE, col = col, 
       lty = 1, text.col = col)
#### END Garman e Klass ####

############### BEGIN Mincer Zarnowitz test #######################
MincerZarnowitz <- rbind(
  data.frame(model = "GARCH", .MincerZarnowitz.local(y = y.out, fit = forc1)),
  data.frame(model = "GJR-GARCH", .MincerZarnowitz.local(y = y.out, fit = forc2)),
  data.frame(model = "T-GARCH", .MincerZarnowitz.local(y = y.out, fit = forc3)),
  data.frame(model = "IGARCH", .MincerZarnowitz.local(y = y.out, fit = forc4)),
  data.frame(model = "T-GARCH + VT", .MincerZarnowitz.local(y = y.out, fit = forc5)),
  data.frame(model = "GJR-GARCH + VT", .MincerZarnowitz.local(y = y.out, fit = forc6)) )
print(MincerZarnowitz, row.names = FALSE)
forc3 <- forc@forecast$sigmaFor[H,]
############### END Mincer Zarnowitz test #######################

##### BEGIN Diebold Mariano test ###############
h <- 1
DieboldMariano <- data.frame(
  measure = c(
    "Volatility",
    "Variance"),
  model = c(
    "GJR-GARCH VS T-GARCH",
    "GJR-GARCH VS T-GARCH"),
  rbind(
    .DieboldMariano.local(y = y.out, f1 = forc2, f2 = forc3, h = h, loss ="LLE"),
    .DieboldMariano.local(y = y.out^2, f1 = forc2^2, f2 = forc3^2, h = h, loss ="LLE") ),
  p_value = c(
    2*(1-pnorm(abs(.DieboldMariano.local(y = y.out, f1 = forc2, f2 = forc3, h =h, loss = "LLE")$statistic))),
    2*(1-pnorm(abs(.DieboldMariano.local(y = y.out^2, f1 = forc2^2, f2 =forc3^2, h = h, loss = "LLE")$statistic)))
  )
)


print(DieboldMariano, row.names = FALSE)
##### END Diebold Mariano test ###############

######################### END EX-POST FORECAST ###########################

######################### BEGIN EX-ANTE FORECAST ###########################
specx <- getspec(object = fit3vt)
fit <- ugarchfit(spec = specx, data = yret.tot)
H <- 10
####  ex-ante, h = 1:H
forc <- ugarchforecast(fitORspec = fit, n.ahead = H, 
                       data = NULL, out.sample = 0, n.roll = 0)
forc@forecast
par(mfrow = c(2,1))
plot(forc, which = 1)
plot(forc, which = 3)

######################### END EX-ANTE FORECAST ###########################
######################################### END SCRIPT #####################################



