// MULTIVARIATE AUTOREGRESSIVE MODEL AR(P) WITH COVARIATES //
// Author: Matheus de Barros //
// Date: 24/04/2024 //
data {
  int<lower=0> T;     // TIME SERIES LENGTH
  int<lower=0> G;     // NUMBER OF SPECIES
  int<lower=0> S;     // NUMBER OF STATIONS
  int<lower=1> P;     // LAG ORDER
  int<lower=1> K;     // NUMBER OF COVARIATES
  vector[G]    Y[T];  // OBSERVATIONS, ARRAY WITH G COLUMNS AND T ROWS
  vector[K]    X[T];  // COVARIATES, ARRAY WITH K COLUMNS AND T ROWS
  int    station[T];
}
transformed data {
  vector[G] Yobs[T-P];
  for (t in 1:(T-P)) {
    Yobs[t] = Y[t+P];
  }
}
parameters {
  vector[G]               phi0[S];    // INTERCEPT, ARRAY WITH G COLUMNS AND S ROWS         
  matrix[G,G]             phi[P,S];   // AUTOREGRESSIVE PARAMETER
  cholesky_factor_corr[G] rho;     // CORRELATION PARAMETER, ARRAY WITH G COLUMNS AND S ROWS 
  vector<lower=0>[G]      sigma;   // NOISE, ARRAY WITH G COLUMNS AND S ROWS 
  matrix[G,K]             gamma;      // COVARIATE EFFECTS, S-DIMENSIONAL ARRAY OF G x K MATRICES
}
transformed parameters {
  matrix[G,G] SIGMA;  // VARIANCE-COVARIANCE MATRIX
  SIGMA = diag_pre_multiply(sigma, rho);

  vector[G] mu[T-P,S];
  for (t in 1:(T-P)) {
    for (s in 1:S) {
      mu[t,s] = phi0[station[s]] + gamma*X[t+P];
      for (p in 1:P) {
        mu[t,s] = mu[t,s] + phi[p, station[s]]*Y[t+p-1];
      }
    }
  }
}
model {
  sigma ~ cauchy(0,1);
  rho   ~ lkj_corr_cholesky(2.0);
  to_vector(gamma) ~ normal(0,1);
  for (s in 1:S) {
    phi0[s]  ~ normal(0,1);
  }
  for (p in 1:P) {
    for (s in 1:S) {
      to_vector(phi[p,s]) ~ normal(0,1);
      Yobs[s] ~ multi_normal_cholesky(mu[s], SIGMA);
    }
  }
}
generated quantities {
  vector[G] Yhat[T-P,S];
  for (t in 1:(T-P)) {
    for (s in 1:S) {
      Yhat[s] = multi_normal_cholesky_rng(mu[s], SIGMA);
    }
  }
}











