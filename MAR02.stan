data {
  int<lower=0> T;
  int<lower=0> S;
  vector[S] Y[T];
}
parameters {
  vector[S] X0;
  real       U;
  vector[S] dev[T];
  real<lower=0> sigmaQ;
  real<lower=0> sigmaR[S];
}
transformed parameters {
  vector[S] X[T];
  for (s in 1:S) {
    X[1,s] = X0[s] + U + dev[1,s];
    for (t in 2:T) {
      X[t,s] = X[t-1,s] + U + dev[t,s];
    }
  }
}
model {
  sigmaQ ~ cauchy(0,1);
  for (s in 1:S) {
    X0[s] ~ normal(Y[s], 10);
    sigmaR[s] ~ cauchy(0,5);
    for (t in 1:T) {
      dev[t,s] ~ normal(0, sigmaQ);
    }
  }
  U ~ normal(0,2);
  for (t in 1:T) {
    for (s in 1:S) {
      Y[s,t] ~ normal(X[s,t], sigmaR[s]);
    }
  } 
}




















