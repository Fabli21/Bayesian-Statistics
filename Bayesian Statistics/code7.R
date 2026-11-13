library(bsda)
library(posterior)
data("factory")

set.seed(47000)
stan_data <- list(
  y = factory,
  N = nrow(factory),
  J = ncol(factory)
)

# 2 b)

mod_sep <- rstan::stan_model(file = "separate.stan")
mod_pool <- rstan::stan_model(file = "pooled.stan")
mod_hier <- rstan::stan_model(file = "hierarchical.stan")

fit_sep <- rstan::sampling(mod_sep, data = stan_data)
fit_pool <- rstan::sampling(mod_pool, data = stan_data)
fit_hier <- rstan::sampling(mod_hier, data = stan_data, control = list(adapt_delta = 0.9999))


print(fit_sep, pars = "sigma", probs = c(0.025, 0.5, 0.975))
print(fit_pool, pars = "sigma", probs = c(0.025, 0.5, 0.975))
print(fit_hier, pars = "sigma", probs = c(0.025, 0.5, 0.975))







# 2 c) i)
draws_sep  <- rstan::extract(fit_sep)
draws_pool <- rstan::extract(fit_pool)
draws_hier <- rstan::extract(fit_hier)

mu6_sep  <- draws_sep$mu[, 6]
mu6_pool <- draws_pool$mu
mu6_hier <- draws_hier$mu[, 6]

# Function to obtain the MCSE 
est_mcse <- function(fit, par) {
  x <- as.array(fit, pars = par)[, , 1]   
  c(mean(x),
    posterior::mcse_mean(x),
    quantile(x, 0.05),
    posterior::mcse_quantile(x, probs = 0.05),
    quantile(x, 0.95),
    posterior::mcse_quantile(x, probs = 0.95))
}
# Table
tab_mu6 <- rbind(
  Separate     = est_mcse(fit_sep,  "mu[6]"),
  Pooled       = est_mcse(fit_pool, "mu"),
  Hierarchical = est_mcse(fit_hier, "mu[6]")
)
colnames(tab_mu6) <- c("Mean", "MCSE Mean", "5%", "MCSE 5%", "95%", "MCSE 95%")
round(tab_mu6, 2)
# Plots

br <- seq(min(mu6_sep, mu6_pool, mu6_hier),
          max(mu6_sep, mu6_pool, mu6_hier),
          length.out = 40)

par(mfrow = c(1, 3))
hist(mu6_sep,  breaks = br, main = "Separate",     
     xlab = expression(mu[6]), col = "aquamarine")
hist(mu6_pool, breaks = br, main = "Pooled",       
     xlab = expression(mu[6]), col = "aquamarine")
hist(mu6_hier, breaks = br, main = "Hierarchical", 
     xlab = expression(mu[6]), col = "aquamarine")
par(mfrow = c(1, 1))


# 2 c) ii)

ypred_sep  <- draws_sep$ypred
ypred_pool <- draws_pool$ypred
ypred_hier <- draws_hier$ypred

tab_ypred <- rbind(
  Separate     = est_mcse(fit_sep,  "ypred"),
  Pooled       = est_mcse(fit_pool, "ypred"),
  Hierarchical = est_mcse(fit_hier, "ypred")
)
colnames(tab_ypred) <- c("Mean", "MCSE Mean", "5%", "MCSE 5%", "95%", "MCSE 95%")
round(tab_ypred, 2)

br <- seq(min(ypred_sep, ypred_pool, ypred_hier),
          max(ypred_sep, ypred_pool, ypred_hier),
          length.out = 40)

par(mfrow = c(1, 3))
hist(ypred_sep,  breaks = br, main = "Separate",     xlab = expression(tilde(y)), col = "aquamarine")
hist(ypred_pool, breaks = br, main = "Pooled",       xlab = expression(tilde(y)), col = "aquamarine")
hist(ypred_hier, breaks = br, main = "Hierarchical", xlab = expression(tilde(y)), col = "aquamarine")
par(mfrow = c(1, 1))


# 2 c) iii)
mu7_pool <- draws_pool$mu
mu7_hier <- draws_hier$mu7

tab_mu7 <- rbind(
  Pooled       = est_mcse(fit_pool, "mu"),
  Hierarchical = est_mcse(fit_hier, "mu7")
)
colnames(tab_mu7) <- c("Mean", "MCSE Mean", "5%", "MCSE 5%", "95%", "MCSE 95%")
round(tab_mu7, 2)

br <- seq(min(mu7_pool, mu7_hier), max(mu7_pool, mu7_hier), length.out = 60)

par(mfrow = c(1, 2))
hist(mu7_pool, breaks = br, main = "Pooled",       xlab = expression(mu[7]), col = "aquamarine")
hist(mu7_hier, breaks = br, main = "Hierarchical", xlab = expression(mu[7]), col = "aquamarine")
par(mfrow = c(1, 1))

# 2 d) 
theta   <- draws_hier$theta
logtau  <- log(draws_hier$tau)


# plot
par(mfrow = c(1, 2))
hist(theta,  breaks = 40, main =  " ",
     xlab = expression(theta), col = "aquamarine")
hist(logtau, breaks = 40, main =  " ",
     xlab = expression(log(tau)), col = "aquamarine")
par(mfrow = c(1, 1))

# 2 e)

mu3    <- draws_hier$mu[, 3]

plot(mu3, logtau,
     pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.3), # opacity
     xlab = expression(mu[3]), ylab = expression(log(tau)),
     main =  " ")
abline(v = mean(theta), col = "red", lty = 2)


# 2 f)

mod_sep_g  <- rstan::stan_model(file = "separate.stan")
mod_pool_g <- rstan::stan_model(file = "pooled.stan")
mod_hier_g <- rstan::stan_model(file = "hierarchical.stan")

fit_sep_g  <- rstan::sampling(mod_sep_g,  data = stan_data)
fit_pool_g <- rstan::sampling(mod_pool_g, data = stan_data)
fit_hier_g <- rstan::sampling(mod_hier_g, data = stan_data, 
                              control = list(adapt_delta = 0.9999))


d_sep_g  <- rstan::extract(fit_sep_g)
d_pool_g <- rstan::extract(fit_pool_g)
d_hier_g <- rstan::extract(fit_hier_g)

d_sep_g  <- rstan::extract(fit_sep_g)
d_pool_g <- rstan::extract(fit_pool_g)
d_hier_g <- rstan::extract(fit_hier_g)

mu6_sep_g  <- d_sep_g$mu[, 6]
mu6_pool_g <- d_pool_g$mu
mu6_hier_g <- d_hier_g$mu[, 6]

tab_f <- rbind(
  Separate     = est_mcse(fit_sep_g,  "mu[6]"),
  Pooled       = est_mcse(fit_pool_g, "mu"),
  Hierarchical = est_mcse(fit_hier_g, "mu[6]")
)
colnames(tab_f) <- c("Mean", "MCSE Mean", "5%", "MCSE 5%", "95%", "MCSE 95%")
round(tab_f, 2)

2 * pnorm(5, 0, 50, lower.tail = FALSE)      # N+(0, 50)
pgamma(5, shape = 1, rate = 1, lower.tail = FALSE)  # Gamma(1, 1)


