data {
  int N;          // TIME SERIES LENGTH
  vector[N]   Y;  // TIME SERIES
  vector[N]   X;  // DESIGN MATRIX OF COVARIATES
}
parameters {
  real        beta;   // AUTOREGRESSIVE COEFFICIENT
  real       alpha;
  real           C;     // MATRIX OF COEFFICIENTS FOR COVARIATES
  real<lower=0> sigma;  // VARIANCE
}
model {
  alpha ~ normal(0,5);
  beta  ~ normal(0,5);
  sigma ~ cauchy(0,5);
  C     ~ normal(0,5);
  for (t in 2:N) {
      Y[t] ~ normal(alpha + beta*Y[t-1] + C*X[t], sigma);
  }
}
generated quantities {
  vector[N] Yhat;
  for (t in 2:N) {
    Yhat[t] = normal_rng(alpha + beta*Y[t-1] + C*X[t], sigma);
  }
}


























