library(forecast)
library(beepr)
library(nnfor)
library(forecastHybrid)
library(dlm)
library(pracma)

cpiQuarterly = read.csv("quarterly cpi change for many countries.csv")
gdpQuarterly = read.csv("quarterly gdp change for many countries.csv")

largestEconomies = c("USA", "JPN", "DEU", "GBR", "FRA", "IND", "ITA", "CAN", "KOR", "AUS", "TUR")

#kalman function for cvts ####
myKalmanFilter = function(ts, modelType){
  if(modelType == "local trend"){
    buildFun <- function(parm) {
      dlmModPoly(2, dV = exp(parm[1]), dW = exp(parm[2:3]))
    }
    fit = dlmMLE(ts, rep(0,3), buildFun, method = "Nelder-Mead")
  }else if(modelType == "bsm"){
    buildFun <- function(p) {
      mod <- dlmModPoly() + dlmModSeas(4)
      V(mod) <- exp(p[1])
      diag(W(mod))[1:3] <- exp(p[2:4])
      return(mod)
    }
    fit = dlmMLE(ts, rep(0,4), buildFun, method = "Nelder-Mead")
  }else if(modelType == "local level"){
    buildFun <- function(p) {
      dlmModPoly(1, dV = exp(p[1]), dW = exp(p[2]))
    }
    fit = dlmMLE(ts, rep(0,2), buildFun, method = "Nelder-Mead")
  }else if(modelType == "regar"){
    buildFun <- function(p) {
      dlmModReg(1, dV = 0, dW = c(0,0)) +
        dlmModARMA(ar = c(p[1], p[2]), sigma = exp(p[3]))
    }
    fit = dlmMLE(ts, rep(0,3), buildFun, method = "Nelder-Mead")
  }
  
  # fit = dlmMLE(myts, rep(0,3), buildFun)
  fit$convergence
  myMod = buildFun(fit$par)
  filter = dlmFilter(ts, myMod)
  
  return(filter)
}

myKalmanForecast = function(filter, h){
  temp = dlmForecast(filter, nAhead = h) 
  kalmanForecasts = list()
  kalmanForecasts$mean = temp$f
  kalmanForecasts$fitted = c(filter$f, temp$f)
  return(kalmanForecasts)
}

# cpiAnnual = read.csv("cpiAnnualWorldBank.csv")
# cpiAnnual = cpiAnnual[which(cpiAnnual$Country.Code %in% largestEconomies),]
# gdpAnnual = read.csv("gdpAnnualWorldBank.csv")
# gdpAnnual = gdpAnnual[which(gdpAnnual$Country.Code %in% largestEconomies),]
# input = na.omit(input)

# listOfCountries = unique(input$Country.Name)
# ts(raw_cpi[which(raw_cpi$LOCATION == "AUS"),"Value"], freq = 4)



#disable scientific notation
options(scipen = 999)

cpiQuarterlyList = list()
cpiQuarterlyHurst = list()
gdpQuarterlyList = list()
gdpQuarterlyHurst = list()

bestCpiPredictions = list()
bestGdpPredictions = list()

#loop ###################################
for(i in 1:length(largestEconomies)){
  # country = listOfCountries[i]
  country = largestEconomies[i]
  #using quarterly data
  # myts = ts(cpiQuarterly[which(cpiQuarterly[,1] == country),"Value"], freq = 4)
  myts = ts(gdpQuarterly[which(gdpQuarterly[,1] == country),"Value"], freq = 4)
  
  #using annual data
  # myts = na.omit(as.numeric(ts(cpiAnnual[cpiAnnual$Country.Code == country,5:61])))
  # myts = na.omit(as.numeric(tCs(gdpAnnual[gdpAnnual$Country.Code == country,5:61])))
  # myts = na.omit(as.numeric(ts(cpiAnnual[cpiAnnual$Country.Code == country,5:61])))
  
  
  # myts = ts(as.numeric(input[i, 5:61]))

  #cvts ##################
  myts = ts(myts, frequency = 4)
  
  
  numObservations = NA
  numObservations = length(myts)
  numObservations
  
  arimaCvts = NA
  localTrendCvts = NA
  bsmCvts = NA
  localLevelCvts = NA
  arfimaCvts = NA
  
  forecastHorizon = 4
  testPeriodLength = forecastHorizon * 10
  initialTrainSize = numObservations - testPeriodLength
  
  arimaCvts = cvts(myts, FUN = auto.arima, rolling = TRUE, windowSize = initialTrainSize, maxHorizon = forecastHorizon, num.cores = 4)
  arfimaCvts = cvts(myts, FUN = arfima, rolling = TRUE, windowSize = initialTrainSize, maxHorizon = forecastHorizon, num.cores = 4)
  localLevelCvts = cvts(myts, FUN = myKalmanFilter, FCFUN = myKalmanForecast, rolling = TRUE, windowSize = initialTrainSize, maxHorizon = forecastHorizon, num.cores = 4, extraPackages = "dlm", modelType = "local level")
  localTrendCvts = cvts(myts, FUN = myKalmanFilter, FCFUN = myKalmanForecast, rolling = TRUE, windowSize = initialTrainSize, maxHorizon = forecastHorizon, num.cores = 4, extraPackages = "dlm", modelType = "local trend")
  bsmCvts = cvts(myts, FUN = myKalmanFilter, FCFUN = myKalmanForecast, rolling = TRUE,  windowSize = initialTrainSize, maxHorizon = forecastHorizon, num.cores = 4, extraPackages = "dlm", modelType = "bsm")
  
  multiStepMeanErrorMtx = matrix(nrow = 5, ncol = 4)
  colnames(multiStepMeanErrorMtx) = c("1stepTestErr", "2stepTestErr", "3stepTestErr", "4stepTestErr")
  rownames(multiStepMeanErrorMtx) = c("arima", "arfima", "localLevel", "localTrend", "bsm")
  
  multiStepMeanErrorMtx[1,] = colMeans(abs(arimaCvts$residuals))
  multiStepMeanErrorMtx[2,] = colMeans(abs(arfimaCvts$residuals))
  multiStepMeanErrorMtx[3,] = colMeans(abs(localLevelCvts$residuals))
  multiStepMeanErrorMtx[4,] = colMeans(abs(localTrendCvts$residuals))
  multiStepMeanErrorMtx[5,] = colMeans(abs(bsmCvts$residuals))
  multiStepMeanErrorMtx
  
  hurstExponent = hurstexp(myts)
  h = mean(unlist(hurstExponent))
  
  # cpiQuarterlyList[[i]] = multiStepMeanErrorMtx
  # cpiQuarterlyHurst[[i]] = h
  
  # gdpQuarterlyList[[i]] = multiStepMeanErrorMtx
  # gdpQuarterlyHurst[[i]] = h
  
  # par(xpd=FALSE)
  # plot((1:4),multiStepMeanErrorMtx["arima",], type = "l", col = 1,
  #      xlab = "num steps ahead", ylab = "mean absolute error", lwd = 2, xaxt="n",
  #      main = "US quarterly GDP percent change forecast mean absolute error for the last 40 quarters",
  #      ylim = c(min(multiStepMeanErrorMtx), max(multiStepMeanErrorMtx)))
  # axis(side = 1, at = 1:4, labels = c("1Q ahead", "2Q ahead", "3Q ahead", "4Q ahead"))
  # lines(multiStepMeanErrorMtx["arfima",], type = "l", col = 2, lwd = 2)
  # lines(multiStepMeanErrorMtx["localLevel",], type = "l", col = 3, lwd = 2)
  # lines(multiStepMeanErrorMtx["localTrend",], type = "l", col = 4, lwd = 2)
  # lines(multiStepMeanErrorMtx["bsm",], type = "l", col = 6, lwd = 2)
  # legend("topleft", 95, legend=c("arima", "arfima", "localLevel", "localTrend", "bsm"),
  #        col=c((1:4), 6), lty = 1, cex=0.9)
  
  
  #forecast for the year ahead ####################
  bestModelName = NA
  # bestModelName = names(which.min(testMAPEs))
  

  predictions = NA
  oneStepMeanPred = NA
  if(bestModelName == "arima"){
    tryCatch({arimaFit = auto.arima(myts, D = 1)}, error=function(e){force(do.next)})
    bestModelFit = arimaFit
    predictions = forecast(arimaFit, h = 4)
    res = residuals(arimaFit)
  }else if(bestModelName == "arfima"){
    tryCatch({arfimaFit = arfima(myts)}, error=function(e){force(do.next)})
    predictions = forecast(arfimaFit, h = 4)
    res = residuals(arfimaFit)
  }else if(bestModelName == "localTrend"){
    kalmanFit = myKalmanFilter(myts, "local trend")
    predictions = myKalmanForecast(kalmanFit, h = 4)
    res = na.omit(myts-kalmanFit$f)
  }else if(bestModelName == "bsm"){
    kalmanFit = myKalmanFilter(myts, "bsm")
    predictions = myKalmanForecast(kalmanFit, h = 4)
    res = na.omit(myts-kalmanFit$f)
  }else if(bestModelName == "localLevel"){
    kalmanFit = myKalmanFilter(myts, "local level")
    predictions = myKalmanForecast(kalmanFit, h = 4)
    res = na.omit(myts-kalmanFit$f)
  }

  oneStepMeanPred = predictions$mean[1]
  oneStepMeanPred
  
  plot(predictions)
  # cpiIncreaseNext4Quarters = ((predictions$mean[4] - predictions$mean[1]) / predictions$mean[1])*100
  # bestCpiPredictions[[i]] = cpiIncreaseNext4Quarters
  
  gdpIncreaseNext4Quarters = ((predictions$mean[4] - predictions$mean[1]) / predictions$mean[1])*100
  bestGdpPredictions[[i]] = gdpIncreaseNext4Quarters
  

  print(country)

  
  
}
beep()


#plotting and aggregation of results ##################################
cpiMeanMultiStepErr = lapply(cpiQuarterlyList, colMeans)
cpiMeanMultiStepErr = matrix(unlist(cpiMeanMultiStepErr), nrow=11, byrow=T)
rownames(cpiMeanMultiStepErr) = largestEconomies
cpiMeanMultiStepErr = data.frame(cpiMeanMultiStepErr)
cpiMeanMultiStepErr = cbind(cpiMeanMultiStepErr, largestEconomies)
colnames(cpiMeanMultiStepErr) = c("1Q-ahead", "2Q-ahead", "3Q-ahead", "4Q-ahead", "country")

gdpMeanMultiStepErr = lapply(gdpQuarterlyList, colMeans)
gdpMeanMultiStepErr = matrix(unlist(gdpMeanMultiStepErr), nrow=11, byrow=T)
rownames(gdpMeanMultiStepErr) = largestEconomies
gdpMeanMultiStepErr = data.frame(gdpMeanMultiStepErr)
gdpMeanMultiStepErr = cbind(gdpMeanMultiStepErr, largestEconomies)
colnames(gdpMeanMultiStepErr) = c("1Q-ahead", "2Q-ahead", "3Q-ahead", "4Q-ahead", "country")

library(ggplot2)
library(reshape2)
df = melt(gdpMeanMultiStepErr, variable.name = "steps")
df = melt(cpiMeanMultiStepErr, variable.name = "steps")

p = ggplot(df, aes(country, value, fill=steps)) + 
  geom_bar(stat = "identity", position = "dodge")

p + labs(title = "multi-step ahead quarterly GDP % change prediction errors") + ylab("absolute mean error")


temp = temp[order(factor(names(temp), levels = c('AUS', 'CAN', 'DEU', 'FRA', 'GBR',  
                                          "IND", "ITA", "JPN", "KOR", "TUR", "USA")))]




df <- data.frame(id=LETTERS[1:4], min=lower, max=upper)
library(ggplot2)
ggplot(df, aes(x=id))+
  geom_linerange(aes(ymin=min,ymax=max),linetype=2,color="blue")+
  geom_point(aes(y=min),size=3,color="red")+
  geom_point(aes(y=max),size=3,color="red")+
  theme_bw()


df <- data.frame(largestEconomies, cpi = unlist(bestCpiPredictions), gdp = unlist(bestGdpPredictions))
f = ggplot(df, aes(df$cpi, df$gdp))
f1 = f + geom_point(aes(col=largestEconomies), size = 4) + geom_text(aes(label=largestEconomies), size=6) 
f1 = f1 + geom_hline(yintercept=0, linetype="dashed", color = "red") + geom_vline(xintercept=0, linetype="dashed", color = "red")
f1 + labs(title = "QoQ CPI growth rate vs QoQ GDP growth rate") + ylab("Percent increase GDP QoQ % over 4 quarters") + xlab("Percent increase GDP QoQ % over 4 quarters")
