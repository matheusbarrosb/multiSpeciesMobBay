require(rstan)
require(dplyr)
require(ggplot2)
require(ggmcmc)

rawData = FAMPtrawl

selectedSpecies = c("CALLINECTES SAPIDUS")
selectedStations = c(28)

sumData = rawData %>% 
  group_by(YEAR, STATION, Genus_species) %>%
  filter(Genus_species %in% selectedSpecies) %>%
  filter(STATION %in% selectedStations) %>%
  summarise(logCPUE = mean(log(Total_NUM+1), na.rm = TRUE),
            sdCPUE = sd(log(Total_NUM+1), na.rm = TRUE),
            seCPUE = sdCPUE/sqrt(n()),
            temp = mean(Temp_C, na.rm = TRUE),
            sdTemp = sd(Temp_C, na.rm = TRUE),
            seTemp = sdTemp/sqrt(n()),
            sal = mean(SAL, na.rm = TRUE),
            sdSal = sd(SAL, na.rm = TRUE),
            seSal = sdSal/sqrt(n()))

ggplot(sumData, aes(x = YEAR, y = logCPUE)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = logCPUE - 2*sdCPUE,
                  ymax = logCPUE + 2*sdCPUE), alpha = 0.3) +
  facet_wrap(Genus_species~STATION, scales = "free_y") +
  theme_linedraw()

# fitting a model
Y = as.matrix(sumData$logCPUE)
X = as.matrix(data.frame(sumData$temp/mean(sumData$temp),
                         sumData$sal/mean(sumData$sal)))
stations = as.factor(sumData$STATION)
st = as.numeric(stations)


dataList = list(Y = Y,
                T = dim(Y)[1],
                X = X,
                P = 1,
                G = dim(Y)[2],
                K = dim(X)[2],
                S = length(unique(stations)),
                station = st)

stanc("SPATIALMAR02.stan")

CALSAP_AR = stan(file = "SPATIALMAR02.stan",
                data = dataList,
                warmup = 300,
                iter = 3000,
                chains = 2,
                cores = 8,
                control = list(max_treedepth = 10),
                verbose = F)

stan_trace(CALSAP_AR, pars = "phi0")
stan_trace(CALSAP_AR, pars = "gamma")
stan_trace(CALSAP_AR, pars = "phi")
stan_dens(CALSAP_AR, pars = "SIGMA")


mcmc = ggs(CALSAP_AR)
mcmc = mcmc %>% 
  filter(grepl('mu', Parameter)) %>%
  group_by(Parameter) %>%
  summarise(mean = mean(value),
           sd = sd(value))


plot(dataList$Y[,1]/mean(dataList$Y[,1]), col = "red", pch = 20, ylim = c(0,2.5))
lines(mcmc$mean/mean(mcmc$mean))
lines(mcmc$mean/mean(mcmc$mean) - mcmc$sd/mean(mcmc$sd)/2, lty = "dashed")
lines(mcmc$mean/mean(mcmc$mean) + mcmc$sd/mean(mcmc$sd)/2, lty = "dashed")

plot(dataList$Y[,1], col = "red", pch = 20)
lines(mcmc$mean+2)
lines(mcmc$mean - 0.5, lty = "dashed")
lines(mcmc$mean + 0.5, lty = "dashed")

stanc("SpatialMAR.stan")






