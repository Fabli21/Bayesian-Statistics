library(markmyassignment)
assignment_path <-
  paste("https://github.com/MansMeg/BSDA/",
        "blob/main/assignments/tests/assignment3.yml", sep="")
set_assignment(assignment_path)
# To check your code/functions, just run
mark_my_assignment()


# Install and load the BSDA package
if (!requireNamespace("BSDA", quietly = TRUE)) {
  install.packages("BSDA")
}
library(BSDA)

# Load dataset
data("windshieldy1")



# Bayesian point estimate: posterior mean E(mu | y) = y_bar
mu_point_est <- function(data) {
  mean(data)
}

mu_point_est(windshieldy1)

# Posterior credible interval for mu using t-distribution
mu_interval <- function(data, prob) {
  n <- length(data)
  df <- n-1
  y_bar <- mean(data)
  s <- sd(data)
  se <- s / sqrt(n)
  
  # Alpha tail probability
  alpha <- 1 - prob
  t_crit <- qt(1 - alpha / 2, df)
  
  c(y_bar - t_crit * se, y_bar + t_crit * se)
}

mu_interval(windshieldy1, 0.95)


# Parameters
data <- windshieldy1
n <- length(data)
df <- n - 1
y_bar <- mean(data)
se <- sd(data) / sqrt(n)

# Simulate values
mu_vals <- seq((y_bar- 4*se), (y_bar+ 4*se), length.out = 1000)

# Create the density using dtnew from bsda
density_vals <- dtnew(mu_vals, df = df, mean = y_bar, scale = se)

# Plot density
plot(mu_vals, density_vals, type = "l", lwd = 2,
     xlab = expression(mu), ylab = " ")

# Point estimate for new observation
mu_pred_point_est <- function(data) {
  mean(data)
}
mu_pred_point_est(windshieldy1)


# Predictioninterval
mu_pred_interval <- function(data, prob) {
  # Parameters
  n <- length(data)
  df <- n - 1
  y_bar <- mean(data)
  s <- sd(data)
  pred_se <- s * sqrt(1 + 1 / n)
  
  # Formula
  alpha <- 1 - prob
  t_crit <- qt(1 - alpha / 2, df = df)
  
  c(y_bar - t_crit * pred_se, y_bar + t_crit * pred_se)
}

mu_pred_interval(windshieldy1, 0.95)

# Scale for new observations
pred_se <- sd(windshieldy1) * sqrt(1 + 1 / n)

# Simulate values
y_new_vals <- seq(y_bar - 4 * pred_se, y_bar + 4 * pred_se, length.out = 1000)

# Simulate values
y_new_vals <- seq(y_bar -  pred_se^2, y_bar +  pred_se^2, length.out = 1000)


# Create the density using dtnew from bsda
pred_dens <- dtnew(y_new_vals, df = df, mean = y_bar, scale = pred_se)
# Plot
plot(y_new_vals, pred_dens, type = "l", lwd = 2,
     xlab = expression(tilde(y)), ylab = " ",
     main = " ")


# Function to compute oddsratio
posterior_odds_ratio_point_est <- function(p0, p1) {
  odds_ratio <- (p1 / (1 - p1)) / (p0 / (1 - p0))
  mean(odds_ratio)
}

# Define function for odds ratio interval
posterior_odds_ratio_interval <- function(p0, p1, prob) {
  odds_ratio <- (p1 / (1 - p1)) / (p0 / (1 - p0))
  quantile(odds_ratio, c((1 - prob) / 2, 1 - (1 - prob) / 2))
}

y0 <- 43
n0 <- 674

y1 <- 22
n1 <- 680

alpha0 <- 90
beta0 <- 630

alpha1 <- 14
beta1 <- 686

alpha0_post <- alpha0 + y0
beta0_post <- beta0 + n0 - y0

alpha1_post <- alpha1 + y1
beta1_post <- beta1 + n1 - y1

# Set seed and draw 100 000 samples from each distribution
set.seed(4711)
p0 <- rbeta(100000, alpha0_post, beta0_post)
p1 <- rbeta(100000, alpha1_post, beta1_post)

posterior_odds_ratio_point_est(p0, p1)


posterior_odds_ratio_interval(p0, p1, 0.95)

# Compute the odds ratio
odds_ratio <- (p1 / (1 - p1)) / (p0 / (1 - p0))

# mean and quantiles
or_mean <- posterior_odds_ratio_point_est(p0, p1)
or_int   <- posterior_odds_ratio_interval(p0, p1, 0.95)

# Compute the odds ratio
odds_ratio <- (p1 / (1 - p1)) / (p0 / (1 - p0))

# Plot
hist(odds_ratio, 
     breaks = 60, 
     col = "lightblue", 
     border = "white",
     main = " ",
     xlab = "Odds Ratio (OR)",
     ylab = "Frequency")

# Add expected value and 95% credible interval
abline(v = or_mean, col = "red", lwd = 2)
abline(v = or_int, col = "darkblue", lwd = 2, lty = 2)



# obtain parameters
n1 <- length(windshieldy1)
ybar1 <- mean(windshieldy1)
s1 <- sd(windshieldy1)
n2 <- length(windshieldy2)
ybar2 <- mean(windshieldy2)
s2 <- sd(windshieldy2)

# Simulate from each posterior
S <- 100000
mu1_sim <- ybar1 + (s1 / sqrt(n1)) * rt(S, df = n1 - 1)
mu2_sim <- ybar2 + (s2 / sqrt(n2)) * rt(S, df = n2 - 1)

# Compute the difference
mu_d <- mu1_sim - mu2_sim

# Mean point estiamte
point_estimate <- mean(mu_d)
point_estimate

# 95% posterior interval
ci_95 <- quantile(mu_d, probs = c(0.025, 0.975))
ci_95

# Plot histogram
hist(mu_d, breaks = 50, col = "skyblue", border = "white",
     main = " ",
     xlab = expression(mu[d]), 
     xlim = range(mu_d))
abline(v = point_estimate, col = "red", lwd = 2)
abline(v = ci_95, col = "darkblue",lwd = 2 lty = 2)


