data {
  int<lower=0> N;       // Number of data points
  vector[N] x;          // Observation year
  vector[N] y;          // Observation number of drowned
  real xpred;           // prediction year
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] mu = alpha + beta * x;
}
model {
  alpha ~ normal(138, 50); // Weakly-informative prior for the intercept
  beta ~ normal(0, 26.787);   // Weakly-informative prior for the slope
  y ~ normal(mu, sigma);
}
generated quantities {
  real ypred = normal_rng(alpha + beta * xpred, sigma);
}