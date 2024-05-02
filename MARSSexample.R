install.packages("MARSS")
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

Q <- "unconstrained" 
B <- "diagonal and unequal" 
A <- U <- x0 <- "zero" 
R <- "diagonal and equal" 
d <- covariates 
D <- "unconstrained" 
Z <- "identity"
y <- dat 
model.list <- list(B = B, U = U, Q = Q, Z = Z,
                   A = A, R = R, D = D, d = d, x0 = x0 ) 
control.list <- list(maxit = 1500) 
kem <- MARSS(y, model = model.list, control = control.list)

states = t(kem$states)

plot(states[,1], type = "l")
lines(t(y)[,1], col = "red")
plot(states[,2], type = "l")

plot(t(y)[,1], col = "red", pch = 20)
lines(states[,1])

kem$par
kem$coef
