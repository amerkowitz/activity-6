ghg <- read.csv("/cloud/project/activity06/Deemer_GHG_Data.csv")

install.packages("dplyr")
install.packages("ggplot2")
install.packages("olsrr")
install.packages("PerformanceAnalytics")

library(dplyr)
library(ggplot2)
library(olsrr)
library(PerformanceAnalytics)
library(forecast)
library(lubridate)

# log transform methane fluxes
ghg$log.ch4 <- log(ghg$ch4+1)

ghg$log.age <- log(ghg$age)
ghg$log.DIP <- log(ghg$DIP+1)
ghg$log.precip <- log(ghg$precipitation)

unique(ghg$Region)

# binary variable for boreal region
ghg$BorealV <- ifelse(ghg$Region == "Boreal",1,0)
# binary variable for tropical region
ghg$TropicalV <- ifelse(ghg$Region == "Tropical",1,0)

# binary variable for alpine region
ghg$AlpineV <- ifelse(ghg$Alpine == "yes",1,0)

# binary variable for known hydropower
ghg$HydroV <- ifelse(ghg$hydropower == "yes",1,0)


# multiple regression
# creates a model object
mod.full <- lm(log.ch4 ~ airTemp+
                 log.age+mean.depth+
                 log.DIP+
                 log.precip+ BorealV, data=ghg) #uses the data argument to specify dataframe
summary(mod.full)

res.full <- rstandard(mod.full)
fit.full <- fitted.values(mod.full)

# qq plot
qqnorm(res.full, pch=19, col="grey50")
qqline(res.full)

# shapiro-wilks test
shapiro.test(res.full)

plot(fit.full,res.full, pch=19, col="grey50")
abline(h=0)

# isolate continuous model variables into data frame:

reg.data <- data.frame(ghg$airTemp,
                       ghg$log.age,ghg$mean.depth,
                       ghg$log.DIP,
                       ghg$log.precip)

# make a correlation matrix 
chart.Correlation(reg.data, histogram=TRUE, pch=19)

# run stepwise
full.step <- ols_step_forward_aic(mod.full)
# view table
full.step 

# plot AIC over time
plot(full.step )

# prediction with interval for predicting a point
predict.lm(mod.full, data.frame(airTemp=20,log.age=log(2),
                                mean.depth=15,log.DIP=3,
                                log.precip=6, BorealV=0),
           interval="prediction")

# look at prediction with 95% confidence interval of the mean

predict.lm(mod.full, data.frame(airTemp=20,log.age=log(2),
                                mean.depth=15,log.DIP=3,
                                log.precip=6, BorealV=0),
           interval="confidence")

#in class activity
ETdat <- read.csv("/cloud/project/activity06/ETdata.csv")
unique(ETdat$crop)

# average fields for each month for almonds
almond <- ETdat %>% # ET data
  filter(crop == "Almonds") %>% # only use almond fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

# visualize the data
ggplot(almond, aes(x=ymd(date),y=ET.in))+
  geom_point()+
  geom_line()+
  labs(x="year", y="Monthy evapotranspiration (in)")

# almond ET time series
almond_ts <- ts(almond$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit

# decompose almond ET time series
almond_dec <- decompose(almond_ts)
# plot decomposition
plot(almond_dec)

acf(na.omit(almond_ts), # remove missing data
    lag.max = 24) # look at 2 years (24 months)

pacf.plot <- pacf(na.omit(almond_ts))

almond_y <- na.omit(almond_ts)
model1 <- arima(almond_y , # data 
                order = c(1,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model1

model4 <- arima(almond_y , # data 
                order = c(4,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
model4

# calculate fit
AR_fit1 <- almond_y - residuals(model1) 
AR_fit4 <- almond_y - residuals(model4)
#plot data
plot(almond_y)
# plot fit
points(AR_fit1, type = "l", col = "tomato3", lty = 2, lwd=2)
points(AR_fit4, type = "l", col = "darkgoldenrod4", lty = 2, lwd=2)
legend("topleft", c("data","AR1","AR4"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3","darkgoldenrod4"),
       bty="n")

newAlmond <- forecast(model4)
newAlmond

#make dataframe for plotting
newAlmondF <- data.frame(newAlmond)

# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newAlmondF$dateF <- ymd(paste(years,"/",month,"/",1))

# make a plot with data and predictions including a prediction interval
ggplot() +
  geom_line(data = almond, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(almond$date[1]),newAlmondF$dateF[24])+  # Plotting original data
  geom_line(data = newAlmondF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newAlmondF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)")

#homework 6

ghg$transformco2 <- 1/(ghg$co2+1000)
ghg$log.transformco2 <- log(ghg$transformco2+1)


#check chlorophyll
plot(ghg$chlorophyll.a, # x data
     ghg$log.transformco2, # y data
     type = "p", 
     pch = 19, # symbol shape
     ylab = "CO2 Transformation", #y axis label
     xlab = "Chlorophyll") #x axis label

ghg$log.chlorophyll.a <- log(ghg$chlorophyll.a+1)

plot(ghg$log.chlorophyll.a, # x data
     ghg$log.transformco2, # y data
     type = "p", 
     pch = 19, # symbol shape
     ylab = "CO2 Transformation", #y axis label
     xlab = "Log Chlorophyll") #x axis label

#check DIP
plot(ghg$DIP, # x data
     ghg$log.transformco2, # y data
     type = "p", 
     pch = 19, # symbol shape
     ylab = "CO2 Transformation", #y axis label
     xlab = "DIP") #x axis label

ghg$log.DIP <- log(ghg$DIP+1)
plot(ghg$log.DIP, # x data
     ghg$log.transformco2, # y data
     type = "p", 
     pch = 19, # symbol shape
     ylab = "CO2 Transformation", #y axis label
     xlab = "DIP") #x axis label
#check precipitation

plot(ghg$precipitation, # x data
     ghg$log.transformco2, # y data
     type = "p", 
     pch = 19, # symbol shape
     ylab = "CO2 Transformation", #y axis label
     xlab = "Precipitation") #x axis label
#check runoff
plot(ghg$log.runoff, # x data
     ghg$transformco2, # y data
     type = "p",
     pch = 19, # symbol shape
     ylab = "CO2 Transformation", #y axis label
     xlab = "Runoff") #x axis label
ghg$log.runoff <- log(ghg$runoff)
plot(ghg$log.runoff, # x data
     ghg$log.transformco2, # y data
     type = "p", 
     pch = 19, # symbol shape
     ylab = "CO2 Transformation", #y axis label
     xlab = "Runoff") #x axis label
#check surface area
plot(ghg$surface.area, # x data
     ghg$log.transformco2, # y data
     type = "p", 
     pch = 19, # symbol shape
     ylab = "CO2 Transformation", #y axis label
     xlab = "Surface Area") #x axis label

ghg$log.surfacearea <- log(ghg$surface.area+1)
plot(ghg$log.surfacearea, # x data
     ghg$log.transformco2, # y data
     type = "p", 
     pch = 19, # symbol shape
     ylab = "CO2 Transformation", #y axis label
     xlab = "Log Surface Area") #x axis label
#check age
plot(ghg$age, # x data
     ghg$log.transformco2, # y data
     type = "p", 
     pch = 19, # symbol shape
     ylab = "CO2 Transformation", #y axis label
     xlab = "Age") #x axis label

#check volume
plot(ghg$volume, # x data
     ghg$log.transformco2, # y data
     type = "p", 
     pch = 19, # symbol shape
     ylab = "CO2 Transformation", #y axis label
     xlab = "Volume") #x axis label
ghg$log.volume <- log(ghg$volume+1)
plot(ghg$log.volume, # x data
     ghg$log.transformco2, # y data
     type = "p", 
     pch = 19, # symbol shape
     ylab = "CO2 Transformation", #y axis label
     xlab = "Log Volume") #x axis label
plot(ghg$mean.depth, # x data
     ghg$log.transformco2, # y data
     type = "p", 
     pch = 19, # symbol shape
     ylab = "CO2 Transformation", #y axis label
     xlab = "Mean Depth") #x axis label
ghg$log.mean.depth <- log(ghg$mean.depth+1)
plot(ghg$log.mean.depth, # x data
     ghg$log.transformco2, # y data
     type = "p", 
     pch = 19, # symbol shape
     ylab = "CO2 Transformation", #y axis label
     xlab = "Mean Depth") #x axis label

# multiple regression
mod.full2 <- lm(log.transformco2 ~ log.mean.depth+log.volume+age+
                 log.DIP+precipitation+
                 log.chlorophyll.a,data=ghg)
summary(mod.full2)

res.full2 <- rstandard(mod.full2)
fit.full2 <- fitted.values(mod.full2)

qqnorm(res.full2, pch=19, col="grey50")
qqline(res.full2)

# shapiro-wilks test
shapiro.test(res.full2)

# isolate continuous model variables into data frame:

reg.data2 <- data.frame(ghg$log.volume,
                       ghg$log.mean.depth, 
                       ghg$log.DIP, ghg$age, ghg$precipitation,
                       ghg$log.chlorophyll.a)

# make a correlation matrix 
chart.Correlation(reg.data2, histogram=TRUE, pch=19)

# run stepwise
full.step2 <- ols_step_forward_aic(mod.full2)
# view table
full.step2 




#question 2
# decompose almond ET time series
almond_dec <- decompose(almond_ts)
# plot decomposition
plot(almond_dec)


# decompose pistachio ET time series

pistachio <- ETdat %>% # ET data
  filter(crop == "Pistachios") %>% # only use pistachio fields
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

pistachio_ts <- ts(pistachio$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit
pistachio_dec <- decompose(pistachio_ts)
# plot decomposition
plot(pistachio_dec)

# decompose Fallow/Idle Cropland ET time series

fallow <- ETdat %>% # ET data
  filter(crop == "Fallow/Idle Cropland") %>%
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

fallow_ts <- ts(fallow$ET.in, # data
                   start = c(2016,1), #start year 2016, month 1
                   #first number is unit of time and second is observations within a unit
                   frequency= 12) # frequency of observations in a unit
# decompose fallow ET time series
fallow_dec <- decompose(fallow_ts)
# plot decomposition
plot(fallow_dec)

# decompose corn ET time series

corn <- ETdat %>% # ET data
  filter(crop == "Corn") %>%
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

corn <- ts(corn$ET.in, # data
               start = c(2016,1), #start year 2016, month 1
               #first number is unit of time and second is observations within a unit
               frequency= 12) # frequency of observations in a unit
corn <- decompose(corn)
# plot decomposition
plot(corn)

# decompose grapes ET time series

grape <- ETdat %>% # ET data
  filter(crop == "Grapes (Table/Raisin)") %>%
  group_by(date) %>% # calculate over each date
  summarise(ET.in = mean(Ensemble.ET, na.rm=TRUE)) # average fields

grape_ts <- ts(grape$ET.in, # data
                start = c(2016,1), #start year 2016, month 1
                #first number is unit of time and second is observations within a unit
                frequency= 12) # frequency of observations in a unit
grape_dec <- decompose(grape_ts)
# plot decomposition
plot(grape_dec)


#question 3

#pistachio
pacf.plot2 <- pacf(na.omit(pistachio_ts))

pistachio_y <- na.omit(pistachio_ts)
modelone <- arima(pistachio_y , # data 
                order = c(1,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
modelone

modelfour <- arima(pistachio_y , # data 
                order = c(4,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
modelfour

# calculate fit

#make dataframe for plotting
newpistachioF <- data.frame(newpistachio)

# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newpistachioF$dateF <- ymd(paste(years,"/",month,"/",1))

ggplot() +
  geom_line(data = pistachio, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(pistachio$date[1]),newpistachioF$dateF[24])+  # Plotting original data
  geom_line(data = newpistachioF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newpistachioF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)")



#fallow
pacf.plot3 <- pacf(na.omit(fallow_ts))

fallow_y <- na.omit(fallow_ts)
fmodel1 <- arima(fallow_y , # data 
                order = c(1,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
fmodel1
fmodel4 <- arima(fallow_y , # data 
                order = c(4,0,0)) # first number is AR order all other numbers get a 0 to keep AR format
fmodel4

# calculate fit
fAR_fit1 <- fallow_y - residuals(fmodel1) 
fAR_fit4 <- fallow_y - residuals(fmodel4)
#plot data
plot(fallow_y)
# plot fit
points(fAR_fit1, type = "l", col = "tomato3", lty = 2, lwd=2)
points(fAR_fit4, type = "l", col = "darkgoldenrod4", lty = 2, lwd=2)
legend("topleft", c("data","AR1","AR4"),
       lty=c(1,2,2), lwd=c(1,2,2), 
       col=c("black", "tomato3","darkgoldenrod4"),
       bty="n")
newfallow <- forecast(fmodel4)
newfallow

#make dataframe for plotting
newfallowF <- data.frame(newfallow)

# set up dates
years <- c(rep(2021,4),rep(2022,12), rep(2023,8))
month <- c(seq(9,12),seq(1,12), seq(1,8))
newfallowF$dateF <- ymd(paste(years,"/",month,"/",1))

# make a plot with data and predictions including a prediction interval
ggplot() +
  geom_line(data = fallow, aes(x = ymd(date), y = ET.in))+
  xlim(ymd(fallow$date[1]),newfallowF$dateF[24])+  # Plotting original data
  geom_line(data = newfallowF, aes(x = dateF, y = Point.Forecast),
            col="red") +  # Plotting model forecasts
  geom_ribbon(data=newfallowF, 
              aes(x=dateF,ymin=Lo.95,
                  ymax=Hi.95), fill=rgb(0.5,0.5,0.5,0.5))+ # uncertainty interval
  theme_classic()+
  labs(x="year", y="Evapotranspiration (in)")
