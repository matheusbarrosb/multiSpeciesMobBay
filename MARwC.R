require(MARSS)

fulldat <- lakeWAplanktonTrans 
years <- fulldat[, "Year"] >= 1965 & fulldat[, "Year"] < 1975 
dat <- t(fulldat[years, c("Greens", "Bluegreens")]) 
the.mean <- apply(dat, 1, mean, na.rm = TRUE) 
the.sigma <- sqrt(apply(dat, 1, var, na.rm = TRUE)) 
dat <- (dat- the.mean) * (1 / the.sigma)

covariates <- rbind(Temp = fulldat[years, "Temp"],
                    TP = fulldat[years, "TP"] ) 
covariates <- zscore(covariates)

Y = t(dat)
X = t(covariates)
data = as.data.frame(cbind(Y,X))
data = na.exclude(data)
Y = as.matrix(data[,2])
X = as.matrix(data[,3:4])
station = as.vector(rbind(rep(1, length(Y[,1])/2), rep(1, length(Y[,1])/2)))
dataList = list(Y = Y,
                T = dim(Y)[1],
                X = X,
                G = dim(Y)[2],
                P = 1,
                K = dim(X)[2],
                S = length(unique(station)),
                station = station)

stanc("MARwC.stan")

MARwCfit = stan(file = "SPATIALMAR02.stan",
               data = dataList,
               warmup = 500,
               iter = 2000,
               chains = 2)

stan_dens(MARwCfit, pars = "SIGMA")


mcmc = ggs(MARwCfit)
mcmc = mcmc %>% 
  filter(grepl('mu', Parameter)) %>%
  group_by(Parameter) %>%
  summarise(mean = mean(value),
            sd = sd(value))


plot(Y[,1], col = "red", pch = 20, ylim = c(-3,3))
lines(mcmc$mean)
lines(mcmc$mean - 1.96*mcmc$sd, lty = "dashed")
lines(mcmc$mean + 1.96*mcmc$sd, lty = "dashed")
lines(dataList$X[,2], col = "blue")

