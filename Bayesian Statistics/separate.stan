data {
  int<lower=0> N;
  int<lower=0> J;
  array[N] vector[J] y;
}

parameters {
  vector[J] mu;
  vector<lower=0>[J] sigma;
}

model {
  // priors
  mu ~ normal(100, 100);
  sigma ~ normal(0, 50);

  // likelihood
  for (j in 1:J)
    y[, j] ~ normal(mu[j], sigma[j]);
}

generated quantities {
  real ypred = normal_rng(mu[6], sigma[6]);
  vector[N * J] log_lik;

  for (j in 1:J)
    for (i in 1:N)
      log_lik[(j - 1) * N + i] = normal_lpdf(y[i, j] | mu[j], sigma[j]);
}
