library(bsda)
library(rstan)
library(loo)

data("factory")

# 2 1)
set.seed(47000)
stan_data <- list(
  y = factory,
  N = nrow(factory),
  J = ncol(factory)
)

# Models
mod_sep  <- stan_model(file = "separate.stan")
mod_pool <- stan_model(file = "pooled.stan")
mod_hier <- stan_model(file = "hierarchical.stan")

fit_sep  <- sampling(mod_sep,  data = stan_data, seed = 47000)
fit_pool <- sampling(mod_pool, data = stan_data, seed = 47000)
fit_hier <- sampling(mod_hier, data = stan_data, seed = 47000,
                     control = list(adapt_delta = 0.9999))

# Log-likelihood: S x N matrix. S = 4000 draws N = 30 observations
log_lik_sep  <- extract_log_lik(fit_sep,  merge_chains = FALSE)
log_lik_pool <- extract_log_lik(fit_pool, merge_chains = FALSE)
log_lik_hier <- extract_log_lik(fit_hier, merge_chains = FALSE)

# 2 2)
# Relative effective sample size
r_eff_sep  <- relative_eff(exp(log_lik_sep))
r_eff_pool <- relative_eff(exp(log_lik_pool))
r_eff_hier <- relative_eff(exp(log_lik_hier))

# PSIS-LOO
loo_sep  <- loo(log_lik_sep,  r_eff = r_eff_sep)
loo_pool <- loo(log_lik_pool, r_eff = r_eff_pool)
loo_hier <- loo(log_lik_hier, r_eff = r_eff_hier)

# elpd values
print(loo_sep)
print(loo_pool)
print(loo_hier)

# k-hat values
par(mfrow = c(1, 3))
plot(loo_sep,  main = "Separate")
plot(loo_pool, main = "Pooled")
plot(loo_hier, main = "Hierarchical")
par(mfrow = c(1, 1))



# 2 3)
# lppd (BDA3 eq. 7.5)
lpd <- function(fit) {
  ll <- extract_log_lik(fit)                        
  sum(apply(ll, 2, function(x) matrixStats::logSumExp(x) - log(length(x))))
}

lpd_sep  <- lpd(fit_sep)
lpd_pool <- lpd(fit_pool)
lpd_hier <- lpd(fit_hier)

# p_eff (BDA3 eq. 7.15): p_loo = lppd - lppd_{loo-cv}
elpd_sep <- loo_sep$estimates["elpd_loo", "Estimate"]
elpd_pool <- loo_pool$estimates["elpd_loo", "Estimate"]
elpd_hier <- loo_hier$estimates["elpd_loo", "Estimate"]

p_eff_sep  <- lpd_sep  - elpd_sep
p_eff_pool <- lpd_pool - elpd_pool
p_eff_hier <- lpd_hier - elpd_hier

# Table: p_eff compared with the number of estimated parameters
tab_peff <- data.frame(
  lpd      = c(lpd_sep, lpd_pool, lpd_hier),
  elpd_loo = c(elpd_sep, elpd_pool, elpd_hier),
  p_eff    = c(p_eff_sep, p_eff_pool, p_eff_hier),
  n_params = c(get_num_upars(fit_sep),  
              get_num_upars(fit_pool),
              get_num_upars(fit_hier)),
  row.names = c("Separate", "Pooled", "Hierarchical")
)

round(tab_peff, 2)


# 2 5)
comp <- loo_compare(list(Separate = loo_sep,
                         Pooled = loo_pool,
                         Hierarchical = loo_hier))
print(comp)



# 3)

mod_loo <- stan_model(file = "separate_loo.stan")
# Store the valeus 
elpd_i <- c()

for (j in 1:6) {
  for (i in 1:5) {
    stan_data$i_out <- i
    stan_data$j_out <- j
    fit <- sampling(mod_loo, data = stan_data, seed = 47000, refresh = 0)
    ll <- rstan::extract(fit)$log_lik_out
    elpd_i <- c(elpd_i, log(mean(exp(ll))))
  }
}
sum(elpd_i)


