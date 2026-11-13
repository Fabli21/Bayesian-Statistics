# =====================================================================
# 02_fit_stan.R  -  fit the model with rstan and check diagnostics
#
# Data prep takes a couple of minutes (26 600 optimisations), so cache it.
# The cache is rebuilt automatically whenever 01_prepare_data.R changes.

file.remove("prep_PL.rds")

# ---------------------------------------------------------------------
# ANSVARSFÖRDELNING I DENNA FIL
#   FABIAN  = setup, sampling och diagnostik (huvudansvar)
#   GIBB    = LOO-delen i slutet
# ---------------------------------------------------------------------
# =====================================================================
# =====================================================================
# >>> FABIAN: setup <<<
# =====================================================================
setwd("/Users/fabianlindh/Desktop/Bayes/Proj/New")


library(rstan)
library(ggplot2)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("01_prepare_data.R")

# Must be the same league as in 03_compare_with_paper.R
league <- "PL"                                        # "PL", "BundesLiga", "Liga", "serieA"

# Data prep takes a couple of minutes (26 600 optimisations), so cache it
prep_file <- paste0("prep_", league, ".rds")
if (file.exists(prep_file)) {
  prep <- readRDS(prep_file)
} else {
  prep <- prepare_league(league, data_dir = ".")
  saveRDS(prep, prep_file)
}

stan_data <- prep$stan_data
# old prep files use N / N_test; the Stan model now uses M / M_test (as in the paper)
if (is.null(stan_data$M)) {
  stan_data$M      <- stan_data$N
  stan_data$M_test <- stan_data$N_test
}

# =====================================================================
# >>> FABIAN: kompilering och sampling <<<
# =====================================================================
mod <- stan_model("egidi_model.stan")

# The JAGS run used 3 chains, 1000 adapt + 1000 burn-in + 2500 samples, thin 2.
# HMC needs far fewer iterations.
fit <- sampling(
  mod, data = stan_data,
  chains = 3, iter = 5000, warmup = 1000, seed = 2025, # 3 chains, 5000 iter, 1000 burnin
  init_r = 0.5,            # smaller random initial values (default 2) to avoid extreme starting rates
  save_warmup = FALSE,     # warmup draws are discarded anyway
  # Sampler settings (kept from the original Stan implementation, not from
  # the paper): stricter than Stan's defaults (adapt_delta 0.8,
  # max_treedepth 10) to reduce the risk of divergent transitions.
  control = list(adapt_delta = 0.9, max_treedepth = 12),
  pars = c("z_att", "z_def", "rate1", "rate2", "theta1", "theta2"),
  include = FALSE          # don't save these, to save memory
)
saveRDS(fit, paste0("fit_", league, ".rds"))

# =====================================================================
# >>> FABIAN: diagnostik <<<
# =====================================================================
# ---- HMC diagnostics: divergences, tree depth, E-BFMI -----------------
check_hmc_diagnostics(fit)
cat("Divergenser:", get_num_divergent(fit), "\n")
cat("Max treedepth nådd:", get_num_max_treedepth(fit), "\n")

# ---- global parameters: posterior summary and diagnostics -------------
main <- c("mu", "home", "mu_att", "mu_def", "sigma2_att", "sigma2_def", "tau1", "tau2")
print(fit, pars = main, probs = c(0.025, 0.5, 0.975), digits_summary = 5)
library(posterior)
dr <- as_draws_df(as.array(fit, pars = main))
summarise_draws(dr, sd, mcse_sd)

# NOTE: mu_att and mu_def cancel in the zero-sum constraint, so they are not
# identified by the likelihood: their posterior equals their N(0, 100) prior.
# This is a property of the original model, not a sampling problem.

# Table for the report. Posterior means are rounded to the decimal position
# of the first significant digit of their Monte Carlo standard error (MCSE).
options(scipen = 999)   # avoid scientific notation in printouts
sm <- summary(fit, pars = main)$summary
digits <- pmax(-floor(log10(sm[, "se_mean"])), 0)
tab <- data.frame(
  mean = mapply(function(x, d) format(round(x, d), nsmall = d), sm[, "mean"], digits),
  sd   = format(signif(sm[, "sd"], 2), scientific = FALSE, drop0trailing = TRUE),
  MCSE = format(signif(sm[, "se_mean"], 2), scientific = FALSE, drop0trailing = TRUE),
  ESS  = round(sm[, "n_eff"]),
  Rhat = sprintf("%.2f", sm[, "Rhat"])
)
print(tab, right = TRUE)

# ---- all parameters: R-hat and ESS summary ----------------------------
# generated quantities (simulated values) are left out
s  <- summary(fit)$summary
gq <- "^(y1_rep|y2_rep|log_lik|y1_test|y2_test|theta1_test|theta2_test)\\[|^lp__$"
s  <- s[!grepl(gq, rownames(s)), ]
rh <- s[, "Rhat"]

cat("\nAntal parametrar:", sum(!is.na(rh)), "\n")
cat("Antal med Rhat > 1.01:", sum(rh > 1.01, na.rm = TRUE), "\n")
cat("Andel med Rhat > 1.01:", round(mean(rh > 1.01, na.rm = TRUE) * 100, 2), "%\n")
cat("Antal med Rhat > 1.1:", sum(rh > 1.1, na.rm = TRUE), "\n")
cat("Max Rhat:", round(max(rh, na.rm = TRUE), 3), "\n")
cat("Minsta n_eff:", round(min(s[, "n_eff"], na.rm = TRUE)), "\n")

bad <- names(rh[rh > 1.01 & !is.na(rh)])
cat("\nParametergrupper med Rhat > 1.01:\n"); print(table(sub("\\[.*", "", bad)))
cat("\nLargest Rhat:\n");   print(head(sort(rh, decreasing = TRUE), 10))
cat("\nSmallest n_eff:\n"); print(head(sort(s[, "n_eff"]), 10))

# ---- trace plots (post-warmup draws) ----------------------------------
traceplot(fit, pars = main) + labs(x = "Iteration")       # global parameters
worst <- names(head(sort(rh, decreasing = TRUE), 6))
traceplot(fit, pars = worst) + labs(x = "Iteration")      # the 6 worst R-hat


traceplot(fit, pars = main) +
  labs(x = "Iteration") +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.83, 0.11),   # lägre än 0.17
    legend.justification   = c(0.5, 0.5),
    legend.title           = element_text(size = 18),
    legend.text            = element_text(size = 18),
    legend.key.height      = unit(1.2, "cm"), # lite tätare rader
    legend.key.width       = unit(1.9, "cm")
  ) +
  guides(colour = guide_legend(override.aes = list(linewidth = 2)))



# =====================================================================
# >>> GIBB: LOO <<<
# =====================================================================
# ---- LOO for the score likelihood (useful if you compare variants) ----
library(loo)
ll <- extract_log_lik(fit, merge_chains = FALSE)
loo_fit <- loo(ll, r_eff = relative_eff(exp(ll)))
print(loo_fit)


packageVersion("rstan"); packageVersion("StanHeaders"); R.version.string; Sys.info()[["machine"]]
fit@stan_args[[1]][c("seed", "iter", "warmup", "chains", "thin", "init")]


p_med <- apply(p_draws, 2, median)
low <- order(p_med)[1:280]
table(paste(d$y1[low], d$y2[low], sep = "-"))

high <- order(p_med, decreasing = TRUE)[1:100]
summary(d$y1[high] + d$y2[high])   # expected: many goals





stan_data <- prep$stan_data
stan_data$M      <- stan_data$N
stan_data$M_test <- stan_data$N_test

fields <- c("M","M_test","nteams","T","S","team1","team2","season","y1","y2",
            "team1_test","team2_test","season_test",
            "theta_bm1","theta_bm2","theta_bm1_test","theta_bm2_test",
            "alpha1","alpha2","alpha1_test","alpha2_test")
pdb_data <- stan_data[fields]

jsonlite::write_json(pdb_data, "egidi_epl_0717.json",
                     auto_unbox = TRUE, digits = NA)



names(stan_data)


prep <- readRDS("prep_PL.rds")
sd <- prep$stan_data
sd$M <- sd$N; sd$M_test <- sd$N_test
sd$N <- NULL; sd$N_test <- NULL
prep$stan_data <- sd[c("M","M_test", setdiff(names(sd), c("M","M_test")))]
saveRDS(prep, "prep_PL.rds")

names(prep$stan_data)
chk <- jsonlite::read_json("egidi_epl_0717.json", simplifyVector = TRUE)
dim(chk$theta_bm1)        # ska vara 3420 x 7
dim(chk$theta_bm1_test)   # ska vara 380 x 7




main <- c("mu", "home", "mu_att", "mu_def", "sigma2_att", "sigma2_def", "tau1", "tau2")
library(posterior)

dr <- as_draws_df(as.array(fit, pars = main))

tab <- summarise_draws(
  dr,
  "mean", "mcse_mean", "sd", "mcse_sd",
  ~quantile(.x, probs = c(0.025, 0.5, 0.975)),
  "ess_bulk", "rhat"
)

tab <- as.data.frame(tab)
rownames(tab) <- tab$variable
tab$variable <- NULL
names(tab) <- c("mean", "se_mean", "sd", "se_sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat")

print(round(tab, 10))


ex <- rstan::extract(fit)
ex$lp__ <- NULL
dims <- lapply(ex, function(a) {
  d <- dim(a)[-1]
  if (length(d) == 0) 1 else d
})
str(dims)
