makeInputData = function(rawData,
                         species,
                         selectedStations,
                         plot = FALSE) {
  
  if ("dplyr" %in% installed.packages() == FALSE) {
    stop("You need the 'dplyr' package to run this function")
  }
  require(dplyr)
  if(is.null(species)) {
    stop("Please specify a species")
  } else if (length(species) != 1) {
    stop("Specify only one species at a time")
  }
  
  sumData = rawData %>% 
    group_by(YEAR, STATION, Genus_species) %>%
    filter(Genus_species %in% species) %>%
    filter(STATION %in% selectedStations) %>%
    summarise(logCPUE = mean(log(Total_NUM+2), na.rm = TRUE),
              temp    = mean(Temp_C, na.rm = TRUE),
              sal     = mean(SAL, na.rm = TRUE),
              maxTemp = mean(Daily_Max_Temp, na.rm = TRUE),
              daysAb33= mean(Days_Above_33_6, na.rm = TRUE),
              consHWds= mean(Consecutive_Heatwave_Days, na.rm = TRUE))
  
  # make STATION a character vector 
  sumData$STATION = paste0(sumData$STATION, ".", sep = "")
  
  splitData = split(sumData, sumData$STATION)
  
  # MAKE RESPONSE DATA
  Ydf = NULL
  for (i in 1:length(splitData)) {
    Ydf[i] = list(splitData[[i]]$logCPUE)
  }
  Ydf           = list2DF(lapply(Ydf, `length<-`, max(lengths(Ydf)))) # MISSING OBSERVATIONS ARE ASSIGNED NAs
  colnames(Ydf) = selectedStations
  Ymat          = as.matrix(Ydf)
  
  # MAKE COVARIATE DATA
  tempdf = NULL
  for (i in 1:length(splitData)) {
    tempdf[i] = list(splitData[[i]]$temp)
  }
  tempdf           = list2DF(lapply(tempdf, `length<-`, max(lengths(tempdf))))
  colnames(tempdf) = selectedStations
  tempmat          = as.matrix(tempdf)
  
  saldf = NULL
  for (i in 1:length(splitData)) {
    saldf[i] = list(splitData[[i]]$sal)
  }
  saldf           = list2DF(lapply(saldf, `length<-`, max(lengths(saldf))))
  colnames(saldf) = selectedStations
  salmat          = as.matrix(saldf)
  
  maxTempdf = NULL
  for (i in 1:length(splitData)) {
    maxTempdf[i] = list(splitData[[i]]$maxTemp)
  }
  maxTempdf           = list2DF(lapply(maxTempdf, `length<-`, max(lengths(maxTempdf))))
  colnames(maxTempdf) = selectedStations
  maxTempmat          = as.matrix(maxTempdf)
  
  daysAb33df = NULL
  for (i in 1:length(splitData)) {
    daysAb33df[i] = list(splitData[[i]]$daysAb33)
  }
  daysAb33df           = list2DF(lapply(daysAb33df, `length<-`, max(lengths(daysAb33df))))
  colnames(daysAb33df) = selectedStations
  daysAb33mat          = as.matrix(daysAb33df)
  
  consHWdsdf = NULL
  for (i in 1:length(splitData)) {
    consHWdsdf[i] = list(splitData[[i]]$consHWds)
  }
  consHWdsdf           = list2DF(lapply(consHWdsdf, `length<-`, max(lengths(consHWdsdf))))
  colnames(consHWdsdf) = selectedStations
  consHWdsmat          = as.matrix(consHWdsdf)
  
  nStations = length(selectedStations)
  X0 = NULL;
  for(i in 1:nStations) {X0[i] = t(Ymat)[,1][i]} # GET INITIAL STATES

  # HANDLING MISSING VALUES
  fulldf    = as.data.frame(cbind(Ymat, tempmat, salmat,
                                  maxTempmat, daysAb33mat,
                                  consHWdsmat))

  splitFullDF = split.default(fulldf,
                              ceiling(seq_along(fulldf)/length(selectedStations)))
  
  Y       = as.matrix(splitFullDF$`1`)
  for(j in 1:ncol(Y)) {
    for (i in 1:nrow(Y)) {
      Y[i,j] = ifelse(is.na(Y[i,j]), mean(Y[j]), Y[i,j])
    }
  }
  Temp    = as.matrix(splitFullDF$`2`)
  for(j in 1:ncol(Temp)) {
    for (i in 1:nrow(Temp)) {
      Temp[i,j] = ifelse(is.na(Temp[i,j]), mean(Temp[j]), Temp[i,j])
    }
  }
  Sal     = as.matrix(splitFullDF$`3`)
  for(j in 1:ncol(Sal)) {
    for (i in 1:nrow(Sal)) {
      Sal[i,j] = ifelse(is.na(Sal[i,j]), mean(Sal[j]), Sal[i,j])
    }
  }
  maxTemp = as.matrix(splitFullDF$`4`)
  for(j in 1:ncol(maxTemp)) {
    for (i in 1:nrow(Temp)) {
      maxTemp[i,j] = ifelse(is.na(maxTemp[i,j]), mean(maxTemp[j]), maxTemp[i,j])
    }
  }
  daysAb33= as.matrix(splitFullDF$`5`)
  for(j in 1:ncol(daysAb33)) {
    for (i in 1:nrow(daysAb33)) {
      daysAb33[i,j] = ifelse(is.na(daysAb33[i,j]), mean(daysAb33[j]), daysAb33[i,j])
    }
  }
  consHWds= as.matrix(splitFullDF$`6`)
  for(j in 1:ncol(consHWds)) {
    for (i in 1:nrow(consHWds)) {
      consHWds[i,j] = ifelse(is.na(consHWds[i,j]), mean(consHWds[j]), consHWds[i,j])
    }
  }
  
  output = list(
    Y    = t(Y),
    Sal  = Sal,
    Temp = Temp,
    maxTemp = maxTemp,
    daysAb33 = daysAb33,
    consHWds = consHWds,
    T    = dim(Y)[1],
    S    = dim(Y)[2],
    X0   = X0
    )
  
  return(output)  
}

#-------------------------------------------------------------------------------
plotAndGetStates = function(fit, nStations, stationIDs,
                            obData, startYr, endYr,
                            panelRows, panelCols, plot = TRUE) {
  
  require(tidyr)
  
  attach.jags(fit)
  means = apply(X, c(2,3), mean)
  upper = apply(X, c(2,3), quantile, 0.95)
  lower = apply(X, c(2,3), quantile, 0.05)
  par(mfrow = c(panelRows, panelCols))
  years  = startYr:endYr
  years  = as.vector(years)
  nYears = length(years)
  
  if (plot == TRUE) {

  for (i in 1:nrow(means)) {
    plot(means[i,], lwd = 2, ylim = c(min(lower),
                                      max(upper)),
         type = "n", main = stationIDs[i], ylab = "log CPUE",
         xlab = "Years")
    polygon(c(1:nYears, nYears:1,1), c(upper[i,], rev(lower[i,]),
                                       upper[i,1]),
            col = "lightblue", lty = 0)
    lines(means[i, ], lwd = 2)
    }
  }
  
  meandf = as.data.frame(t(means))
  colnames(meandf) = as.character(stationIDs)
  meandf = tidyr::pivot_longer(meandf,
                        cols = 1:ncol(meandf),
                        names_to = "station", values_to = "estState")
  
  upperdf = as.data.frame(t(upper))
  colnames(upperdf) = as.character(stationIDs)
  upperdf = tidyr::pivot_longer(upperdf,
                         cols = 1:ncol(upperdf),
                         names_to = "station", values_to = "upperEstState")
  
  lowerdf = as.data.frame(t(lower))
  colnames(lowerdf) = as.character(stationIDs)
  lowerdf = tidyr::pivot_longer(lowerdf,
                         cols = 1:ncol(lowerdf),
                         names_to = "station", values_to = "lowerEstState")
  
  observeddf = as.data.frame(t(obData))
  colnames(observeddf) = as.character(stationIDs)
  observeddf = tidyr::pivot_longer(observeddf,
                                   cols = 1:ncol(observeddf),
                                   names_to = "station", values_to = "observed")
  
 estStates = data.frame(
   gl(nYears, nStations),
   meandf$station,
   meandf$estState,
   upperdf$upperEstState,
   lowerdf$lowerEstState,
   observeddf$observed
 )
 names(estStates) = c("years", "station", "mean", "upper", "lower", "observations")
 return(estStates)
  
}

# ------------------------------------------------------------------------------

plotGammas = function (fit, sppName) {
  
  require(dplyr)
  
  attach.jags(fit)
  gammaSal  = as.data.frame(gammaSal)
  gammaSal  = gammaSal %>% mutate(Parameter = "gamma[Salinity]")
  names(gammaSal) = c("value", "Parameter")
  gammaTemp = as.data.frame(gammaTemp)
  gammaTemp  = gammaTemp %>% mutate(Parameter = "gamma[Temperature]")
  names(gammaTemp) = c("value", "Parameter")
  gammaMaxTemp = as.data.frame(gammaMaxTemp)
  gammaMaxTemp  = gammaMaxTemp %>% mutate(Parameter = "gamma[MaxTemperature]")
  names(gammaMaxTemp) = c("value", "Parameter")
  
  df = rbind(gammaSal, gammaTemp, gammaMaxTemp)
  dfSum = df %>%
    group_by(Parameter) %>%
    summarise(mean = median(value),
              sd = sd(value))
  
  ggplot(data = df, aes(x = Parameter, y = value)) +
    see::geom_violinhalf(fill = "skyblue3", alpha = 0.3) +
    geom_point(data = dfSum,
               aes(x = Parameter,
                   y = mean),
               pch = 20, size = 3) +
    coord_flip() +
    scale_x_discrete(labels = scales::label_parse()) +
    theme_linedraw() +
    ggtitle(label = sppName) +
    theme(plot.title = element_text(size = 9)) +
    xlab("") + 
    ylab("") +
    geom_hline(aes(yintercept = 0),
               linetype = "dashed",
               col = "red")
  
}















