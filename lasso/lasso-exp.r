## run_experiments.R

# — Install/load dependencies
if (!require(glmnet))   install.packages("glmnet",   repos="https://cloud.r-project.org")
if (!require(doParallel)) install.packages("doParallel", repos="https://cloud.r-project.org")
if (!require(foreach))    install.packages("foreach",    repos="https://cloud.r-project.org")
if (!require(ggplot2))    install.packages("ggplot2",    repos="https://cloud.r-project.org")

library(glmnet)
library(doParallel)
library(foreach)
library(ggplot2)

# — User‐tuneable settings
repCount <- 30        # how many reps per combo
sigma_true <- 10      # fixed noise sd
alpha <- 0.005        # test level (not used in LC)
ncores <- 38          # number of cores to use

# — Parameter grid
paramGrid <- expand.grid(
  n      = c(100,   1000, 10000),
  k      = c(10,    50,   100,  500,  2000, 5000),
  beta1  = c(0.1,   1.0,  10),
  beta2  = c(0.1,   1.0,  10),
  beta3  = c(0.1,   1.0,  10),
  stringsAsFactors = FALSE
)

# — parallel setup
cl <- makeCluster(ncores)
registerDoParallel(cl)

# — run simulations
results <- foreach(pi = iter(paramGrid, by="row"), 
                   .combine = rbind, 
                   .packages = c("glmnet")) %dopar% {
  local_out <- data.frame()
  
  for (rep in seq_len(repCount)) {
    set.seed(1000 + rep)  # reproducibility
    
    # unpack params
    n  <- pi$n
    k  <- pi$k
    p  <- k + 1
    b1 <- pi$beta1
    b2 <- pi$beta2
    b3 <- pi$beta3
    
    # simulate design + response
    X <- cbind(1, matrix(rnorm(n * k), nrow = n))
    beta_true <- numeric(p)
    beta_true[1] <- b1
    beta_true[2] <- b2
    beta_true[3] <- b3
    y <- as.vector(X %*% beta_true + rnorm(n, sd = sigma_true))
    
    # split train/test
    n2 <- n/2
    train_idx <- (n2+1):n
    test_idx  <- 1:n2
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    X_test  <- X[test_idx,  , drop = FALSE]
    y_test  <- y[test_idx]
    
    # null log-likelihood (intercept only)
    beta0_hat <- c(mean(y_train), rep(0, k))
    sigma0_hat <- sqrt(sum((y_train - X_train %*% beta0_hat)^2) / n2)
    ll0 <- sum(dnorm(y_test - X_test %*% beta0_hat, sd = sigma0_hat, log = TRUE))
    
    # LASSO fit (alt)
    cvfit <- cv.glmnet(x = X_train[,-1], y = y_train,
                       alpha = 1, intercept = TRUE, standardize = TRUE)
    lam <- cvfit$lambda.1se
    pred_train <- predict(cvfit, newx = X_train[,-1], s = lam)
    sigma1_hat <- sqrt(sum((y_train - pred_train)^2) / n2)
    pred_test <- predict(cvfit, newx = X_test[,-1], s = lam)
    ll1 <- sum(dnorm(y_test - pred_test, sd = sigma1_hat, log = TRUE))
    
    # LC = exp(ll0 - ll1)
    LCval <- exp(ll0 - ll1)
    LRval <- exp(ll1 - ll0)
    
    # collect
    local_out <- rbind(
      local_out,
      data.frame(n = n, k = k,
                 beta1 = b1, beta2 = b2, beta3 = b3,
                 rep   = rep,    LC = LCval, LR = LRval)
    )
  }
  
  local_out
}

stopCluster(cl)

# — save raw results
write.csv(results, "results.csv", row.names = FALSE)

# — prepare for plotting
results$k     <- factor(results$k)
results$beta1 <- factor(results$beta1)
results$beta2 <- factor(results$beta2)
results$beta3 <- factor(results$beta3)

# — example visualization strategy:
#    boxplots of LC vs k, colored by beta1, with small multiples for beta2×beta3 at each sample size
results$beta23 <- paste0("b2=", results$beta2, ", b3=", results$beta3)

gg <- ggplot(results, aes(x = k, y = LC, fill = beta1)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  facet_grid(beta23 ~ n, scales = "free_x", labeller = label_value) +
  labs(
    x    = "Number of predictors (k)",
    y    = "LC (p-value)",
    fill = expression(beta[1]),
    title = "Distribution of LC across k, β₂/β₃ settings and sample sizes"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8)
  )

ggsave("LC_boxplot.png", gg, width = 12, height = 9)
