
############################################
#Install required packages
# Loading packages
###########################################

# Install necessary packages
install.packages(c("quantmod", "rugarch", "PerformanceAnalytics", "forecast"))
install.packages("urca")

# Load libraries
library(quantmod) # To get the data
library(rugarch) # To estimate the GARCH model
library(PerformanceAnalytics) # For evaluation
library(forecast) # For forecasting 
library(urca) # for unit root
library(ggplot2) # For ploting 
library(tseries) # for unit root
library(MTS) # for ARCH test
library(FinTS) #for function `ArchTest()`

#########################################################
# Downloading data and cleaning the data
########################################################


# Download TCS stock data
getSymbols("SBIN.NS", src = "yahoo", from = "1996-01-01", to = "2024-11-01")

# Taking the concerned data 
wdata <- Cl(SBIN.NS)

# examining the data
head(wdata)
start(wdata)
end(wdata)
summary(wdata)

#check for missing values

sum(is.na(wdata))

wdata <- na.omit(wdata)


# Plotting series including ACF and PACF

ggtsdisplay(wdata, plot.type = "partial",smooth = TRUE, theme=theme_classic())

#Check for stationarity (ADF-test)


summary(ur.df(wdata, type = "drift", selectlags = "AIC"))
#Value of t-statistics = -0.523 > -3.43 = Critical value 
# => H0 hypothesis of non-stationarity is failed to be rejected


#Take log of the prices
rtcs = diff(log(wdata))
rtcs <- rtcs[-1,] # dropping first observation for NA
#Plot the graph
plot(rtcs, type = "l", xlab = "")


#Check for stationarity again for transformed series(ADF-test)
summary(ur.df(rtcs, type = "drift", selectlags = "AIC"))
#Value of t-statistics = -35.6766 < -3.43 = Critical value 
# => H0 hypothesis of non-stationarity of the return is rejected


# Visualize returns
ts.plot(rtcs, main = "Daily Log Returns", col = "blue")
ggtsdisplay(rtcs, plot.type = "partial",smooth = TRUE, theme=theme_classic())
# Looks pretty stationary

##################################################################
# GARCH Modelling
##################################################################
# Check for the presence of ARCH effect 

#ARCH Test
rtcsArchTest <- ArchTest(rtcs, lags=1, demean=TRUE)
rtcsArchTest


# Presence of GARCH effect

# Plot squared returns
srtcs <- rtcs^2
plot(srtcs, main = "Squared Returns", col = "red")

#If volatility clustering exists, squared returns will show bursts of high values clustered together.

# ACF of squared returns
acf(srtcs, main = "ACF of Squared Returns")

#Significant autocorrelations in the squared returns suggest the presence of volatility clustering.

# Perform Ljung-Box test on squared returns
Box.test(srtcs, lag = 12, type = "Ljung-Box")

#the p-value is small, reject the null hypothesis, indicates potential GARCH

###############################################
#GARCH Specification and estimation
###############################################

# Specify a GARCH(1,1) model
garch_spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE)
)

# Fit the GARCH model
garch11 <- ugarchfit(spec = garch_spec,out.sample = 10, data = rtcs)

# Display the results of the model
print(garch11)


# Plot results
plot(garch11, which = "all") 

# Check for autocorrelation in residuals
acf(residuals(garch11, standardize = TRUE), main = "ACF of Standardized Residuals")


# Ljung-Box test for ARCH effects
Box.test(residuals(garch11, standardize = TRUE)^2, lag = 12, type = "Ljung-Box")

g11 <- infocriteria(garch11)[1]

# Check for GARCH(1,2) model
garch12 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,2)))
garch12 = ugarchfit(garch12, data = rtcs)
garch12

g12 <- infocriteria(garch12)[1]


garch12a <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,2)),
                       mean.model = list(armaOrder = c(0, 0), include.mean = TRUE))
garch12a = ugarchfit(garch12a, data = rtcs)
garch12a

g12a <- infocriteria(garch12)[1]


# GRACH-M

garchm_spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE, archm = 1, archpow = 2), # archm=1 for variance
    distribution.model = "norm" # Normal distribution for residuals
)

garchm <- garchm_fit <- ugarchfit(spec = garchm_spec, data = rtcs)
garchm 
gm <- infocriteria(garchm)[1]


# TGARCH

tgarch_spec <- ugarchspec(
    variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE)
)

tgarch <- ugarchfit(spec = tgarch_spec, data = rtcs)
print(tgarch)

tg <- infocriteria(tgarch)[1]

# EGARCH

egarch_spec <- ugarchspec(
    variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE)
)

egarch <- ugarchfit(spec = egarch_spec, data = rtcs)
print(egarch)

eg <- infocriteria(egarch)[1]


# print information criteria

# Extract information criteria for each model
garch_ic11 <- infocriteria(garch11)
garch_ic12 <- infocriteria(garch12)
garch_ic12a <- infocriteria(garch12a)
mgarch_ic <- infocriteria(garchm)
tgarch_ic <- infocriteria(tgarch)
egarch_ic <- infocriteria(egarch)



# Collecting the information criteria into a single data frame
info_criteria_df <- data.frame(
    Model = c("GARCH(1,1)", "GARCH(1,2)", "GARCH(1,2a)", "MGARCH", "TGARCH", "EGARCH"),
    AIC = c(garch_ic11[1], garch_ic12[1], garch_ic12a[1], mgarch_ic[1], tgarch_ic[1], egarch_ic[1]),
    BIC = c(garch_ic11[2], garch_ic12[2], garch_ic12a[2], mgarch_ic[2], tgarch_ic[2], egarch_ic[2]),
    HQIC = c(garch_ic11[3], garch_ic12[3], garch_ic12a[3], mgarch_ic[3], tgarch_ic[3], egarch_ic[3])
)

# Display the data frame
print(info_criteria_df)





# Reshape the data frame for plotting
library(reshape2)
info_criteria_long <- melt(info_criteria_df, id.vars = "Model")

# Plot the criteria
library(ggplot2)
ggplot(info_criteria_long, aes(x = Model, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Comparison of Information Criteria for GARCH Models",
         x = "Model", y = "Value", fill = "Criterion") +
    theme_minimal()


# Forecast 10-step ahead volatility
garch_forecast <- ugarchforecast(garch12, n.ahead = 10)

# Extract forecasted volatility
forecasted_volatility <- sigma(garch_forecast)
print(forecasted_volatility)

# Plot forecasted volatility
plot(garch_forecast, which = 3)


# Final Model: Re-estimate GARCH with out.sample
garch12r <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,2)),
                         mean.model = list(armaOrder = c(0, 0), include.mean = TRUE ))

garch12r = ugarchfit(garch12r, data = rtcs,out.sample = 10)

print(garch12r)

#Diagnostics

# Plot results
plot(garch11, which = "all") 

# Check for autocorrelation in residuals
acf(residuals(garch12r, standardize = TRUE), main = "ACF of Standardized Residuals")


# Ljung-Box test for ARCH effects
Box.test(residuals(garch12r, standardize = TRUE)^2, lag = 12, type = "Ljung-Box")

# Forecast
garch_forecastr <- ugarchforecast(garch12r, n.ahead = 10)
fpm(garch_forecastr)


# Extract conditional volatility
forecasted_volatility <- sigma(garch_forecastr)
print(forecasted_volatility)

# Plot conditional volatility
plot(garch_forecastr, which = 3)



