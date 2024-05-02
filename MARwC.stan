// MULTIVARIATE AUTOREGRESSIVE MODEL AR(P) WITH COVARIATES //
// Author: Matheus de Barros //
// Date: 24/04/2024 //
data {
  int<lower=0> T;     // TIME SERIES LENGTH
  int<lower=0> G;     // NUMBER OF SPECIES
  int<lower=1> P;     // LAG ORDER
  int<lower=1> K;     // NUMBER OF COVARIATES
  vector[G]    Y[T];  // OBSERVATIONS, ARRAY WITH G COLUMNS AND T ROWS
  vector[K]    X[T];  // COVARIATES, ARRAY WITH K COLUMNS AND T ROWS
}
transformed data {
  vector[G] Yobs[T-P];
  for (t in 1:(T-P)) {
    Yobs[t] = Y[t+P];
  }
}
parameters {
  vector[G]               phi0;     // INTERCEPT           
  matrix[G,G]             phi[P];   // AUTOREGRESSIVE PARAMETER
  cholesky_factor_corr[G] rho;      // CORRELATION PARAMETER
  vector<lower=0>[G]      sigma;    // NOISE
  matrix[G,K]             gamma;    // COVARIATE EFFECTS
}
transformed parameters {
  matrix[G,G] SIGMA;  // VARIANCE-COVARIANCE MATRIX
  SIGMA = diag_pre_multiply(sigma, rho);
  vector[G] mu[T-P];
  for (t in 1:(T-P)) {
    mu[t] = phi0 + gamma*X[t+P];
    for (p in 1:P) {
      mu[t] = mu[t] + phi[p]*Y[t+p-1];
    }
  }
}
model {
  rho ~ lkj_corr_cholesky(2.0);
  sigma ~ cauchy(0,1);
  phi0 ~ normal(0,1);
  for (p in 1:P) {
    to_vector(phi[p]) ~ normal(0,1);
  }
  Yobs ~ multi_normal_cholesky(mu, SIGMA);
}
generated quantities {
  vector[G] Yhat[T-P];
  for (t in 1:(T-P)) {
   Yhat = multi_normal_cholesky_rng(mu, SIGMA);
  }
}











