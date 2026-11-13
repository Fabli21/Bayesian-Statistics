library(markmyassignment)
assignment_path <-
  paste("https://github.com/MansMeg/BSDA/",
        "blob/main/assignments/tests/assignment2.yml", sep="")
set_assignment(assignment_path)
# To check your code/
mark_my_assignment()
install.packages("remotes")
remotes::install_github("MansMeg/BSDA", subdir = "rpackage", upgrade = "never")

install.packages("bsda")
library(bsda)
data("algae")
head(algae)

round(qbeta(c(0.025, 0.975), 1.5, 0.5),4)

round(pbeta(0.5, 3/2, 1/2, lower.tail = FALSE),3)

round(qbeta(c(0.025, 0.975), 5/2, 3/2),4)

round(pbeta(0.5, 5/2, 3/2, lower.tail = FALSE),3)



beta_point_est <- function(prior_alpha, prior_beta, data){
 # Obtain y and n
  y <- sum(data)
  n <- length(data)
  
  # Posteroir parameters
  post_alpha <- prior_alpha + y
  post_beta <- prior_beta + n -y
  
  # Point estimate (mean for the beta distribution)
  point_est <- post_alpha/(post_alpha + post_beta)
  
  return(point_est)
}

beta_point_est(2, 10, algae)

beta_interval <- function(prior_alpha, prior_beta, data, prob){
  # Obtain y and n
  y <- sum(data)
  n <- length(data)
  
  # Posteroir parameters
  post_alpha <- prior_alpha + y
  post_beta <- prior_beta + n -y
  
  # Compute the tail probabilities  
  lower_tail <- (1- prob)/2
  upper_tail <- 1 - lower_tail
  
  # Find the upper and lower bound of the credible interval
  lower_bound <- qbeta(lower_tail, post_alpha, post_beta)
  upper_bound <- qbeta(upper_tail, post_alpha, post_beta)
  
  return(c(lower_bound, upper_bound))
}

beta_interval(2, 10, algae, 0.9)


beta_low <- function(prior_alpha, prior_beta, data, pi_0){
  # Obtain y and n
  y <- sum(data)
  n <- length(data)
  
  # Posteroir parameters
  post_alpha <- prior_alpha + y
  post_beta <- prior_beta + n -y
  
  # Compute the probability pi is less than pi_0
  prob <- pbeta(pi_0, post_alpha, post_beta)
  
  return(prob)
}

beta_low(2, 10, algae, 0.2)


y <- sum(algae)
n <- length(algae)

# Create sequence of pi values between 0 and 1
pi <- seq(0, 1, length.out = 1000)

# Create three different posteriors Beta(a + y, b + n - y)

# Posterior 1, with Prior Beta (2,10) the original as a reference
posterior_1 <- dbeta(pi, 2 + y, 10 + n - y)

# Posterior 2, with Prior Beta(1,1)
posterior_2 <- dbeta(pi, 1 + y, 1 + n - y)

# Posterior 3, with Prior Beta(8,4)
posterior_3 <- dbeta(pi, 8 + y, 4 + n - y)

# Find max to adjust height of plot
max_y <- max(c(posterior_1, posterior_2, posterior_3))

plot(pi, posterior_1, type = "l", col = "blue", lwd = 2, 
     ylab = "Density", xlab = expression(pi), 
     xlim = c(0.1, 0.25),
     ylim = c(0, max_y * 1.4))

# Add the other posteriors
lines(pi, posterior_2, col = "red", lwd = 2, lty = 2)
lines(pi, posterior_3, col = "green", lwd = 2, lty = 3)

# Explanation (legend)
legend("topright",
       legend = c("Prior Beta(2,10)", "Prior Beta(1,1)", "Prior Beta(8,4)"),
       col = c("blue", "red", "green"),
       lwd = 2, lty = c(1, 2, 3), 
       bty = "n", cex = 0.9)





library(datasets)
data("warpbreaks")
head(warpbreaks$breaks)

# y and n
y_bar <- mean(warpbreaks$breaks)
n <- length(warpbreaks$breaks)
# sum of y
sum_y <- n*y_bar

# posterior parameters
alpha_post <- 1.5 + sum_y
beta_post <- 0.1 + n

lower_bound <- qgamma(0.025, shape = alpha_post, rate = beta_post)
upper_bound <- qgamma(0.975, shape = alpha_post, rate = beta_post)

round(c(lower_bound, upper_bound), 3)


# Compute prob
prob <- beta_post/(beta_post + 1)

# Predict y
y_pred <- rnbinom(n, alpha_post, prob)

# Change settings
par(mfrow = c(1, 2)) 

# Plot histograms for the simulated data and the actual data
hist(y_pred, xlim = c(5, 75), main = "Posterior Predictive", 
     xlab = " ", col = "aquamarine")

hist(y, xlim = c(5, 75), main = "Observed", 
     xlab = " ", col = "coral")


var(warpbreaks$breaks)
var(y_pred)

var(warpbreaks$breaks)/var(y_pred)


0.9*14
