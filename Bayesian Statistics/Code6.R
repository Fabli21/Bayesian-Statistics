library(bsda)
library(rstan)
library(MASS)
data("bioassay")
# Prior
mu0 <- c(0, 10)
Sigma0 <- matrix(c(2^2, 12,
                   12, 10^2), nrow = 2, byrow = TRUE)

# 5 chains 
n_chains <- 5
set.seed(47000)
# Inflate Sigma for the starting points
inflate <- 4
starts <- mvrnorm(n_chains, mu = mu0, Sigma = inflate * Sigma0)
colnames(starts) <- c("alpha", "beta")

# Create vector for the starting points
inits <- lapply(1:n_chains, function(k) list(theta = as.numeric(starts[k, ])))

# Data for the model
stan_data <- list(
  N = nrow(bioassay),
  x = bioassay$x,
  n = as.integer(bioassay$n),
  y = as.integer(bioassay$y),
  mu0 = mu0,
  Sigma0 = Sigma0
)
# Fit the model
fit <- stan(file = "/Users/fabianlindh/Desktop/Bayes/HA6/mod.stan", data = stan_data,
            chains = n_chains, iter = 10000, warmup = 5000,
            init = inits, seed = 47000)

# Rhat
draws <- rstan::extract(fit, permuted = FALSE)
Rhat(draws[, , "theta[1]"])   # alpha
Rhat(draws[, , "theta[2]"])   # beta


# Scatter plot
post  <- rstan::extract(fit, pars = "theta")
alpha <- post$theta[, 1]
beta  <- post$theta[, 2]

plot(alpha, beta,
     pch = 19, cex = 0.2, col = rgb(0, 0, 0, 0.2),
     xlab = expression(alpha), ylab = expression(beta),
     main = " ")


# 3
# i)
# Prepare data for the model to predict drownings for 2020
stan_data_drown <- list(
  N = nrow(drowning),
  x = drowning$year- 1980, # To shift the data to make the interpretation easier
  y = drowning$drownings,
  xpred = 2020-1980
)

# Model
fit_drown <- stan(
  file   = "/Users/fabianlindh/Desktop/Bayes/HA6/drowning.stan",
  data   = stan_data_drown,
  chains = 4, iter = 2000, warmup = 1000, seed = 47000
)

# Print fit for beta
print(fit_drown, "beta",
      probs = c(0.025, 0.5, 0.975))

# Extract for histogram
post_d <- rstan::extract(fit_drown, pars = c("beta", "ypred"))

# Plot beta
hist(post_d$beta, breaks = 40,
     main = " ",
     xlab = expression(beta), 
     col = "aquamarine")

# Print prediction for 2020
print(fit_drown, pars = c("ypred"),
      probs = c(0.025, 0.5, 0.975))

# Plot prediction 2020
hist(post_d$ypred, breaks = 40, 
     main = " ", 
     xlab = "Drownings", 
     col = "aquamarine")


# a)
# Debug in test debug.stab

# b)
qnorm(0.995)
sigma_beta <- 69 / qnorm(0.995)
sigma_beta

# Controll
pnorm(69, 0, sigma_beta) - pnorm(-69, 0, sigma_beta)

# c) 
pnorm(0, 138, 50, lower.tail = FALSE)


#   alpha ~ normal(138, 50); // Weakly-informative prior for the intercept
# beta ~ normal(0, 26.787);   // Weakly-informative prior for the slope

View(drowning)
nrow(drowning)
