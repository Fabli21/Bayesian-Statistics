install.packages("markmyassignment")
library(markmyassignment)
assignment_path <-
  paste("https://github.com/MansMeg/BSDA/",
        "blob/main/assignments/tests/assignment1.yml", sep="")
set_assignment(assignment_path)

mu <- 0.2
sigma2 <- 0.01

# Calculate alpha and beta
a <- mu * ((mu * (1 - mu)) / sigma2 - 1)
b <- (a * (1 - mu)) / mu

# Create a sequence of x values from 0 to 1
x <- seq(0, 1, length.out = 1000)

# Calculate the Beta density
y <- dbeta(x, a, b)

# Plot
plot(x, y, type = "l", lwd = 2)


# Draw 1000 random samples
samples <- rbeta(1000, a, b)

# Plot histogram
hist(samples, freq = FALSE)

# Add the true density line from part a
lines(x, y, lwd = 2)

# Sample mean and sample variance
round(mean(samples),3)

round(var(samples),4)

# 95% probability interval
quantile(samples, probs = c(0.025, 0.975))


boxes <- matrix(c(2,4,1,5,1,3), ncol = 2,
                dimnames = list(c("A", "B", "C"), c("red", "white")))

p_red <- function(boxes) {
  p_box_selection <- c(0.4, 0.1, 0.5)
  
  
  p_red_given_box <- boxes[, "red"] / rowSums(boxes)
  
  # Law of Total Probability
  total_prob <- sum(p_red_given_box * p_box_selection)
  
  return(total_prob)
}

p_red(boxes)

p_box <- function(boxes) {
  p_box_selection <- c(0.4, 0.1, 0.5)
  
  #Probability of picking red from each box
  p_red_given_box <- boxes[, 1] / rowSums(boxes)
  
  #Probability picking a red ball all scenarios
  joint_probs <- p_red_given_box * p_box_selection
  
  #Probability of picking red
  p_red_total <- sum(joint_probs)
  
  # Bayes theorem
  posterior_probs <- joint_probs / p_red_total
  
  return(posterior_probs)
}

p_box(boxes)

p_identical_twin <- function(fraternal_prob, identical_prob)
  {
  p_bb_i <- 0.5
  p_bb_f <- 0.25
  
  p_bb <- p_bb_i* identical_prob + p_bb_f*fraternal_prob
  

  p_bb_i_posterior <- (p_bb_i * identical_prob)/p_bb
  
  
  return(p_bb_i_posterior)
}

p_identical_twin(1/150, 1/400)  


mark_my_assignment()



