data {
  int<lower=0> N;                    // Number of dose groups
  vector[N] x;                       // Dose (log g/ml)
  array[N] int<lower=0> n;           // Number of animals per group
  array[N] int<lower=0> y;           // Number of deaths per group
  vector[2] mu0;                     // Prior mean
  matrix[2, 2] Sigma0;               // Prior covariance
}
parameters {
  vector[2] theta;                   // theta[1] = alpha, theta[2] = beta
}
model {
  theta ~ multi_normal(mu0, Sigma0);  // Prior
  y ~ binomial_logit(n, theta[1] + theta[2] * x); // Likelihood
}