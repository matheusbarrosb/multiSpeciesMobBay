require(R2jags)
require(rjags)
require(coda)

modelString <- cat("
model {

# MULTIVARIATE AUTOREGRESSIVE STATE-SPACE MODEL FOR MODELLING RECRUITMENT TIME SERIES
# AUTHOR: MATHEUS DE BARROS
# DATE: 29/04/2024 (dd/mm/yyyy)

# S = NUMBER OF STATIONS
# T = LENGTH OF TIME SERIES
# K = NUMBER OF COVARIATES

# PRIOR DISTRIBUTIONS ----------------------------------------------------------
gammaTemp ~ dnorm(0, 0.01)        # EFFECT OF TEMPERATURE
gammaSal  ~ dnorm(0, 0.01)        # EFFECT OF SALINITY
for (s in 1:S) {
  phi[s]  ~ dnorm(0, 0.01)        # AUTOREGRESSIVE COEFFICIENT
  phi0[s] ~ dnorm(0, 0.01)        # INTERCEPT
  tauR[s] ~ dgamma(0.01,0.01)     # INVERSE OF OBSERVATION VARIANCE
  R[s] <- 1/tauR[s]               # OBSERVATION VARIANCE, I.E. VARIANCE-COVARIANCE MATRIX SINCE Y ~ MVN(X,R)
}
for (t in 2:T) {
  tauQ[t] ~ dgamma(0.01, 0.01)    # INVERSE OF PROCESS VARIANCE
  Q[t]    = 1/tauQ[t]             # PROCESS VARIANCE, TIME VARYING
}

# INITIAL STATE VECTOR T = 1 ---------------------------------------------------
for (s in 1:S) {
  X[s,1] ~ dnorm(X0[s], 10)
}

# AUTOREGRESSIVE PROCESS T + 1 : T ---------------------------------------------
for (t in 2:T) {
  for (s in 1:S) {
    Xhat[s,t] <- phi0[s] + phi[s]*(X[s,t-1] - phi0[s]) + gammaSal*Sal[t,s] + gammaTemp*Temp[t,s]
    X[s,t] ~ dnorm(Xhat[s,t], tauQ[t])
  }
}

# OBSERVATION COMPONENT LIKELIHOOD ---------------------------------------------
for (t in 1:T) {
  for (s in 1:S) {
    Y[s,t] ~ dmnorm(X[s,t], tauR[s])
  }
}

}",
file = "MARSS_jags.txt")

Y = as.matrix(data.frame(cbind(sumDataCALSAP03$logCPUE,
                               sumDataCALSAP04$logCPUE)))
tempmat = as.matrix(data.frame(
  sumDataCALSAP03$temp/mean(sumDataCALSAP03$temp),
  sumDataCALSAP04$temp/mean(sumDataCALSAP04$temp))
  )
colnames(tempmat) = c("1", "2")

salmat = as.matrix(data.frame(
  sumDataCALSAP03$sal/mean(sumDataCALSAP03$sal),
  sumDataCALSAP04$sal/mean(sumDataCALSAP04$sal))
)
colnames(salmat) = c("1", "2")

X0 = NULL; for (i in 1:dim(Y)[2]) {X0[i] = Y[,i][1]} # get initial states

salmat[33,2] = NA
salmat[34,2] = NA

dataList = list(
  Y = t(Y),
  Sal = salmat,
  Temp = tempmat,
  T = dim(Y)[1], # length of time series 
  S = dim(Y)[2],
  X0 = X0
)


params = c("X", "Q", "R", "gammaTemp", "gammaSal", "phi0", "phi")
model = "MARSS_jags.txt"
fitJags = jags(a, parameters.to.save = params, model.file = model,
               n.chains = 2, n.burnin = 3000, n.thin = 1, n.iter = 10000)

plot(fitJags, pars = c("gammaTemp", "gammaSal"))

attach.jags(fitJags)
means <- apply(X, c(2, 3), mean)
upperCI <- apply(X, c(2, 3), quantile, 0.975)
lowerCI <- apply(X, c(2, 3), quantile, 0.025)
par(mfrow = c(3, 3))
nYears <- 34
for (i in 1:nrow(means)) {
  plot(means[i, ], lwd = 3, ylim = c(0,6),
       type = "n",
       ylab = "log CPUE", 
  xlab = "time step")
  polygon(c(1:nYears, nYears:1, 1), c(upperCI[i, ], rev(lowerCI[i, 
  ]), upperCI[i, 1]), col = "skyblue", lty = 0)
  lines(means[i, ], lwd = 2)
  title(rownames(Y)[i])
}

par(mfrow = c(2, 2))
plot(Y[,1], col = "red", pch = 20)
lines(means[1,])
plot(Y[,2], col = "red", pch = 20)
lines(means[2,])
plot(Y[,3], col = "red", pch = 20)
lines(means[3,])


