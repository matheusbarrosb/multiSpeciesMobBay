require(rstan)
require(ggmcmc)
stanc("AR01.stan")

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

Y = as.matrix(data[,1:2])
X = as.matrix(data[,3:4])

dataList = list(Y = Y[,2],
                X = X[,2],
                N = length(Y[,2]))

stanc("MARwC.stan")

AR01fit = stan(file = "AR01.stan",
               data = dataList,
               warmup = 500,
               iter = 2000,
               chains = 2)


mcmc = ggs(AR01fit)
mcmc = mcmc %>% 
  filter(grepl('Yhat', Parameter)) %>%
  group_by(Parameter) %>%
  summarise(mean = mean(value),
            sd = sd(value))


plot(Y[,2], col = "red", pch = 20)
lines(mcmc$mean)
lines(mcmc$mean - mcmc$sd, lty = "dashed")
lines(mcmc$mean + mcmc$sd, lty = "dashed")
lines(dataList$X, col = "blue")




