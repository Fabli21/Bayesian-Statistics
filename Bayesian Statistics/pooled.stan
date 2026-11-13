data {
  int<lower=0> N;
  int<lower=0> J;
  array[N] vector[J] y;
}

parameters {
  real mu;
  real<lower=0> sigma;
}

model {
  // priors
  mu ~ normal(100, 100);
  sigma ~ normal(0, 50);

  // likelihood
  for (j in 1:J)
    y[, j] ~ normal(mu, sigma);
}

generated quantities {
  real ypred = normal_rng(mu, sigma);
  vector[N * J] log_lik;

  for (j in 1:J)
    for (i in 1:N)
      log_lik[(j - 1) * N + i] = normal_lpdf(y[i, j] | mu, sigma);
}
