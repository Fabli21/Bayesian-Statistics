data {
  int<lower=0> N;
  int<lower=0> J;
  array[N] vector[J] y;
  int i_out;   
  int j_out;  
}

parameters {
  vector[J] mu;
  vector<lower=0>[J] sigma;
}

model {
  mu ~ normal(100, 100);
  sigma ~ normal(0, 50);

  for (j in 1:J)
    for (i in 1:N)
      if (i != i_out || j != j_out)         
        y[i, j] ~ normal(mu[j], sigma[j]);
}

generated quantities {
  real log_lik_out = normal_lpdf(y[i_out, j_out] | mu[j_out], sigma[j_out]);
}
