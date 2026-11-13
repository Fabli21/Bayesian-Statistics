library(MASS)
library(bsda)
data("bioassay")

# 1a) # Function compute the density ratio function 11.1 in BDA
density_ratio <- function(alpha_propose, alpha_previous,
                          beta_propose, beta_previous,
                          x, y, n) {
  # Set Prior parameters
  mu0    <- c(0, 10)
  Sigma0 <- matrix(c(2^2, 12,
                     12,  10^2), nrow = 2, byrow = TRUE)
  
  # Log unnormalized posterior = log-likelihood + log-prior
  log_post_propose <- bioassaylp(alpha_propose, beta_propose, x, y, n) +
    dmvnorm(c(alpha_propose, beta_propose), mean = mu0, sigma = Sigma0, log = TRUE)
  
  log_post_previous <- bioassaylp(alpha_previous, beta_previous, x, y, n) +
    dmvnorm(c(alpha_previous, beta_previous), mean = mu0, sigma = Sigma0, log = TRUE)
  
  exp(log_post_propose - log_post_previous)
}

# 1b) # Create function for Metropolis algorithm
Metropolis_bioassay <- function(starting_point, n_iter, sigma, x, y, n) {
  draws <- matrix(NA, nrow = n_iter, ncol = 2,
                  dimnames = list(NULL, c("alpha", "beta")))
  draws[1, ] <- starting_point
  n_accepted <- 0
  
  for (t in 2:n_iter) {
    previous <- draws[t - 1, ]
    
    # Proposal normal distribution with prevoius values
    proposal <- rnorm(2, mean = previous, sd = sigma)
    
    # Density ratio from a)
    r <- density_ratio(alpha_propose  = proposal[1],
                       alpha_previous = previous[1],
                       beta_propose   = proposal[2],
                       beta_previous  = previous[2],
                       x = x, y = y, n = n)
    
    # Accept with probability min(1, r)
    if (runif(1) < r) {
      draws[t, ] <- proposal
      n_accepted <- n_accepted + 1
    } else {
      draws[t, ] <- previous
    }
  }
  
  list(draws = draws, acceptance_rate = n_accepted / (n_iter - 1))
}

set.seed(47000)

# 2c) Starting points for the chains
# Starting points for the chains with the prior with infalted Sigma 
mu0    <- c(0, 10)
Sigma0 <- matrix(c(2^2, 12,
                   12,  10^2), nrow = 2, byrow = TRUE)

n_chains <- 5
set.seed(47000)
# Sample from the norm for the starting points
inflate <- 4
starts <- mvrnorm(n_chains, mu = mu0, Sigma = inflate * Sigma0)
colnames(starts) <- c("alpha", "beta")
starts

# 2d) number of iterations for each chain
n_iter <- 10000

# 2e) the warm-up length which will be discarded
warmup <- 5000

# 2f) # Number of chains
n_chains <- nrow(starts)     


# 2b) Test several proposals, jump scales and compare acceptance rate

# Function that runs for all chains with help of the 
#metropolis_bioassay function
run_chains <- function(sigma) {
  lapply(1:n_chains, function(k)
    Metropolis_bioassay(starts[k, ], n_iter = n_iter, sigma = sigma,
                        x = bioassay$x, y = bioassay$y, n = bioassay$n))
}

for (sig in list(c(0.2, 1), c(0.8, 4), c(1, 5), c(1.2, 6), c(5, 25))) {
  res <- run_chains(sig)
  cat("sigma =", sig, " acceptance rate =",
      round(mean(sapply(res, `[[`, "acceptance_rate")), 2), "\n")
}

# Choose the most reasonable sigma (sigma = (1.2, 6) acceptance rate = 0.37)
sigma_final <- c(1.2, 6)
chains <- run_chains(sigma_final)


# collect all chains in matrices
alpha_all <- sapply(chains, function(ch) ch$draws[, "alpha"])
beta_all  <- sapply(chains, function(ch) ch$draws[, "beta"])

# 2g) plot all chains for alpha add warm-up as jagged line
matplot(alpha_all, type = "l", lty = 1, col = 1:n_chains,
        xlab = "Iteration", ylab = expression(alpha), main = " ")
abline(v = warmup, lty = 2)  

# 2h) Same for beta
matplot(beta_all, type = "l", lty = 1, col = 1:n_chains,
        xlab = "Iteration", ylab = expression(beta), main = " ")
abline(v = warmup, lty = 2)




library(rstan)


alpha_post <- alpha_all[(warmup + 1):nrow(alpha_all), ]
beta_post  <- beta_all[(warmup + 1):nrow(beta_all), ]

# 3) R-hat for alpha and beta
Rhat(alpha_post)
Rhat(beta_post)
# 4) plot
plot(alpha_post, beta_post,
     pch  = 19,                      
     cex  = 0.2,
     col  = rgb(0, 0, 0, alpha = 0.2), # Opacity
     xlab = expression(alpha),
     ylab = expression(beta))

