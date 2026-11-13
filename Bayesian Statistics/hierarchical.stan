data {
  int<lower=0> N;
  int<lower=0> J;
  array[N] vector[J] y;
}

parameters {
  real theta;
  real<lower=0> tau;
  vector[J] mu;
  real<lower=0> sigma;
}

model {
  // priors
  theta ~ normal(100, 100);
  tau ~ normal(0, 50);
  sigma ~ normal(0, 50);
  mu ~ normal(theta, tau);

  // likelihood
  for (j in 1:J)
    y[, j] ~ normal(mu[j], sigma);
}

generated quantities {
  real ypred = normal_rng(mu[6], sigma);
  real mu7 = normal_rng(theta, tau);
  vector[N * J] log_lik;

  for (j in 1:J)
    for (i in 1:N)
      log_lik[(j - 1) * N + i] = normal_lpdf(y[i, j] | mu[j], sigma);
}
