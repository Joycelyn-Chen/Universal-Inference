# Required Libraries
library(MASS)    # For mvrnorm, if needed
library(mclust)  # For GMM fitting via EM

# Function to simulate data from a GMM with k components
simulate_gmm <- function(n, k, mu, sigma, pi) {
  comps <- sample(1:k, size = n, replace = TRUE, prob = pi)
  rnorm(n, mean = mu[comps], sd = sigma[comps])
}

# Function to fit a GMM model and return log-likelihood
fit_gmm <- function(data, k) {
  model <- Mclust(data, G = k, warn = FALSE)
  list(
    log_likelihood = model$loglik
  )
}

# Sieve-selection helper
run_sieve <- function(data, alpha, max_k) {
  lr_vec <- numeric(max_k - 1)
  prev_fit <- fit_gmm(data, 1)
  selected <- 1
  
  for (k in 2:max_k) {
    curr_fit <- fit_gmm(data, k)
    lr <- exp(curr_fit$log_likelihood - prev_fit$log_likelihood)
    lr_vec[k - 1] <- lr
    
    # If we fail to reject H_{k-1}, pick k-1
    if (lr <= 1/alpha) {
      selected <- k - 1
      break
    }
    
    prev_fit <- curr_fit
    selected <- k
  }
  
  list(likelihood_ratios = lr_vec, selected = selected)
}

# Simulation parameters
n      <- 200
alpha  <- 0.1
max_k  <- 10

# TRUE-K = 3 scenario
true_k3   <- 3
mu3       <- c(0, 2, 5)
sigma3    <- c(1, 1, 1)
pi3       <- c(0.3, 0.4, 0.3)
set.seed(42)
data3     <- simulate_gmm(n, true_k3, mu3, sigma3, pi3)
res3      <- run_sieve(data3, alpha, max_k)
cat("Selected model for true_k = 3:", res3$selected, "components\n")

# TRUE-K = 4 scenario
true_k4   <- 4
mu4       <- c(0, 2, 5, 8)
sigma4    <- rep(1, 4)
pi4       <- rep(1/4, 4)
set.seed(123)
data4     <- simulate_gmm(n, true_k4, mu4, sigma4, pi4)
res4      <- run_sieve(data4, alpha, max_k)
cat("Selected model for true_k = 4:", res4$selected, "components\n")

# Combined visualization
x_vals <- 2:max_k
plot(x_vals, res3$likelihood_ratios, type="b", pch=16,
     xlab="Number of Components", ylab="Likelihood Ratio",
     main="Sieve Model Selection: true_k=3 vs true_k=4",
     ylim = range(c(res3$likelihood_ratios, res4$likelihood_ratios, 1/alpha)))
lines(x_vals, res4$likelihood_ratios, type="b", pch=17, col="blue")
abline(h = 1/alpha, col = "red", lty = 2)
legend("topright",
       legend = c("true_k = 3", "true_k = 4", "Threshold 1/α"),
       col    = c("black",      "blue",      "red"),
       pch    = c(16,           17,          NA),
       lty    = c(1,            1,           2))


# # Required Libraries
# library(MASS)    # For mvrnorm, if you ever need multivariate normals

# # Install the package: install.packages("mclust")
# library(mclust)  # For GMM fitting via EM

# # Function to simulate data from a GMM with k components
# simulate_gmm <- function(n, k, mu, sigma, pi) {
#   comps <- sample(1:k, size = n, replace = TRUE, prob = pi)
#   rnorm(n, mean = mu[comps], sd = sigma[comps])
# }

# # Function to fit a GMM model and return log-likelihood and parameters
# fit_gmm <- function(data, k) {
#   model <- Mclust(data, G = k, warn = FALSE)
#   list(
#     log_likelihood = model$loglik,
#     mu     = model$parameters$mean,
#     sigma  = sqrt(model$parameters$variance$sigmasq),
#     pi     = model$parameters$pro
#   )
# }

# # Simulation parameters
# n        <- 200
# true_k   <- 3
# mu_true  <- c(0, 2, 5)
# sigma_true <- c(1, 1, 1)
# pi_true  <- c(0.3, 0.4, 0.3)

# # Simulate data
# set.seed(42)
# data <- simulate_gmm(n, true_k, mu_true, sigma_true, pi_true)

# # Sieve parameters
# alpha  <- 0.1
# max_k  <- 5

# # Storage
# likelihood_ratios <- numeric(max_k - 1)
# model_selected   <- NA

# # Fit the first model (k = 1) by default
# prev <- fit_gmm(data, 1)
# model_selected <- 1

# # Loop through k = 2 to max_k
# for (k in 2:max_k) {
#   curr <- fit_gmm(data, k)
  
#   # Compute likelihood ratio: L_k / L_{k-1}
#   lr <- exp(curr$log_likelihood - prev$log_likelihood)
#   likelihood_ratios[k - 1] <- lr
  
#   # Test H_{k-1}: reject if lr > 1/alpha
#   if (lr <= 1/alpha) {
#     model_selected <- k - 1
#     break
#   }
  
#   prev <- curr
#   model_selected <- k
# }

# # Print selected model
# cat("Selected model has", model_selected, "components\n")

# # Visualization
# x_vals <- 2:max_k

# plot(x_vals, likelihood_ratios,
#      type = "b", pch = 16,
#      xlab = "Number of Components",
#      ylab = "Likelihood Ratio",
#      main = "Model Selection Using Sieves",
#      ylim = range(c(likelihood_ratios, 1/alpha))  # ensure threshold line fits
# )

# abline(h = 1/alpha, col = "red", lty = 2)

# # Highlight the selected model (if > 1)
# if (model_selected > 1) {
#   points(model_selected, likelihood_ratios[model_selected - 1],
#          col = "blue", pch = 19, cex = 1.5)
# }

# legend("topright",
#        legend = c("Likelihood Ratios", "Threshold 1/α", "Selected Model"),
#        col    = c("black",      "red",        "blue"),
#        pch    = c(16,           NA,           19),
#        lty    = c(1,            2,            NA))
