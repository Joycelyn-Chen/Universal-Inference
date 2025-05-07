# Number of simulations
N_sims <- 1000

# Parameters for simulation
n <- 100  # Sample size
mu_vals <- seq(0, 3, by = 0.1)  # Varying mean separation values (from 0 to 3)
alpha <- 0.1  # Significance level

# Power for each mu value
power_universal <- rep(0, length(mu_vals))
power_bootstrap <- rep(0, length(mu_vals))

# Function to compute the log likelihood
log_L <- function(mu, x) {
  sum(log(1/2*dnorm(x+0.5*mu) + 1/2*dnorm(x-0.5*mu)))
}

# Run simulations for each mu value
for (mu_index in 1:length(mu_vals)) {
  mu <- mu_vals[mu_index]
  
  # For each simulation
  ratio_arr <- rep(0, N_sims)
  power_bootstrap_sim <- 0
  
  for (i in 1:N_sims) {
    # Simulate the mixture data
    B <- rbinom(2 * n, 1, 0.5)
    component <- 1 - 2 * B
    y <- rnorm(2 * n, mean = 0.5*mu * component)
    
    y1 <- y[1:n]
    y0 <- y[-(1:n)]
    
    # Maximum likelihood # Number of simulations
N_sims <- 1000

# Parameters for simulation
n <- 100  # Sample size
mu_vals <- seq(0, 3, by = 0.1)  # Varying mean separation values (from 0 to 3)
alpha <- 0.1  # Significance level

# Power for each mu value
power_universal <- rep(0, length(mu_vals))
power_bootstrap <- rep(0, length(mu_vals))

# Function to compute the log likelihood
log_L <- function(mu, x) {
  sum(log(1/2*dnorm(x+0.5*mu) + 1/2*dnorm(x-0.5*mu)))
}

# Run simulations for each mu value
for (mu_index in 1:length(mu_vals)) {
  mu <- mu_vals[mu_index]
  
  # For each simulation
  ratio_arr <- rep(0, N_sims)
  power_bootstrap_sim <- 0
  
  for (i in 1:N_sims) {
    # Simulate the mixture data
    B <- rbinom(2 * n, 1, 0.5)
    component <- 1 - 2 * B
    y <- rnorm(2 * n, mean = 0.5*mu * component)
    
    y1 <- y[1:n]
    y0 <- y[-(1:n)]
    
    # Maximum likelihood estimation
    y1_L <- function(mu) {
      -log_L(mu, y1)
    }
    mu1_hat <- abs(optim(mu, y1_L, method = 'Brent', lower = min(y1), upper = max(y1))$par)
    mu0_hat <- 0
    
    
    # Store the log-likelihood ratio
    log_ratio_arr <- log_L(mu1_hat, y0) - log_L(mu0_hat, y0)
    ratio_arr[i] <- exp(log_ratio_arr)  # Convert back to likelihood ratio

    # Bootstrap test
    bootstrap_sample <- sample(1:(2 * n), n, replace = TRUE)
    y_bootstrap <- y[bootstrap_sample]
    
    # Estimate mu for the bootstrap sample
    mu_bootstrap_hat <- abs(optim(mu, y1_L, method = 'Brent', lower = min(y_bootstrap), upper = max(y_bootstrap))$par)
    
    # Compute the bootstrap likelihood ratio statistic
    bootstrap_statistic <- exp(log_L(mu_bootstrap_hat, y) - log_L(mu0_hat, y))
    
    # Bootstrap power estimation
    if (bootstrap_statistic > 1 / alpha) {
      power_bootstrap_sim <- power_bootstrap_sim + 1
    }
  }
  
  # Compute the power of the universal test and the bootstrap test
  power_universal[mu_index] <- mean(ratio_arr > 1 / alpha)
  power_bootstrap[mu_index] <- power_bootstrap_sim / N_sims
}


# Plot the results (power vs mean separation)
plot(mu_vals, power_universal, type = "l", col = "black", lwd = 2, 
     xlab = expression(mu), ylab = "Power", main = "Power of Universal and Bootstrap Tests")
lines(mu_vals, power_bootstrap, col = "red", lwd = 2)
legend("bottomright", legend = c("Universal Test", "Bootstrap Test"), 
       col = c("black", "red"), lwd = 2)
# lines(mu_vals, xx, col='blue')
estimation
    y1_L <- function(mu) {
      -log_L(mu, y1)
    }
    mu1_hat <- abs(optim(mu, y1_L, method = 'Brent', lower = min(y1), upper = max(y1))$par)
    mu0_hat <- 0
    
    
    # Store the log-likelihood ratio
    log_ratio_arr <- log_L(mu1_hat, y0) - log_L(mu0_hat, y0)
    ratio_arr[i] <- exp(log_ratio_arr)  # Convert back to likelihood ratio

    # Bootstrap test
    bootstrap_sample <- sample(1:(2 * n), n, replace = TRUE)
    y_bootstrap <- y[bootstrap_sample]
    
    # Estimate mu for the bootstrap sample
    mu_bootstrap_hat <- abs(optim(mu, y1_L, method = 'Brent', lower = min(y_bootstrap), upper = max(y_bootstrap))$par)
    
    # Compute the bootstrap likelihood ratio statistic
    bootstrap_statistic <- exp(log_L(mu_bootstrap_hat, y) - log_L(mu0_hat, y))
    
    # Bootstrap power estimation
    if (bootstrap_statistic > 1 / alpha) {
      power_bootstrap_sim <- power_bootstrap_sim + 1
    }
  }
  
  # Compute the power of the universal test and the bootstrap test
  power_universal[mu_index] <- mean(ratio_arr > 1 / alpha)
  power_bootstrap[mu_index] <- power_bootstrap_sim / N_sims
}


# Plot the results (power vs mean separation)
plot(mu_vals, power_universal, type = "l", col = "black", lwd = 2, 
     xlab = expression(mu), ylab = "Power", main = "Power of Universal and Bootstrap Tests")
lines(mu_vals, power_bootstrap, col = "red", lwd = 2)
legend("bottomright", legend = c("Universal Test", "Bootstrap Test"), 
       col = c("black", "red"), lwd = 2)
# lines(mu_vals, xx, col='blue')
