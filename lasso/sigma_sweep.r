############################################################
# sigma_sweep.r
#
#   Sweep over sigma_true ∈ [4,6] (by 0.05),
#   for n ∈ {100,1000,10000}, two β‐settings,
#   1000 simulations each, then plot P(LC<th) vs sigma_true.
#
############################################################

# — Load/install dependencies
if (!require(glmnet))   install.packages("glmnet", repos="https://cloud.r-project.org")
if (!require(ggplot2))  install.packages("ggplot2", repos="https://cloud.r-project.org")
if (!require(reshape2)) install.packages("reshape2", repos="https://cloud.r-project.org")

library(glmnet)
library(ggplot2)
library(reshape2)

# — Parameters
sigma_values <- seq(4.00, 6.00, by = 0.05)   # 4.00, 4.05, 4.10, …, 5.95, 6.00
ns           <- c(100, 1000, 10000)         # three sample sizes
reps         <- 1000                        # # of sims per combination
k            <- 100                         # number of predictors (dropping intercept)
p            <- k + 1                       # total columns in X (including intercept)
alpha        <- 0.005                       # (not used directly for P(LC<…))
n_beta_sets  <- 2

# Two β‐settings (only β₂ and β₃ vary; β₁ is fixed = 2)
beta_settings <- list(
  list(b2 = 1.5, b3 = 0.5, label = "β₂=1.5,β₃=0.5", linetype = "solid"),
  list(b2 = 0.0, b3 = 0.0, label = "β₂=0,  β₃=0",   linetype = "dashed")
)

# Prepare a data.frame to collect results
results <- data.frame(
  n         = integer(0),
  beta_lbl  = character(0),
  sigma_true= numeric(0),
  prop_005  = numeric(0),
  prop_01   = numeric(0),
  prop_05   = numeric(0),
  stringsAsFactors = FALSE
)

set.seed(2048)  # for reproducibility of random draws

for (n in ns) {
  n2 <- n / 2
  
  for (bs in beta_settings) {
    b2       <- bs$b2
    b3       <- bs$b3
    beta_lbl <- bs$label
    linetype <- bs$linetype
    
    cat("======== n =", n, ", setting:", beta_lbl, "========\n")
    
    for (sigma_true in sigma_values) {
      # Pre‐allocate vectors for LC values in this block of 1000 sims
      LCs <- numeric(reps)
      LRs <- numeric(reps)
      # (We compute LC only; we don’t actually need to store LRs for the final plot.)
      
      for (i in seq_len(reps)) {
        # — 1) Simulate X and y under true model
        #    X has an intercept column of 1's, then k columns N(0,1)
        X <- cbind(1, matrix(rnorm(n * (p - 1)), nrow = n))
        
        # True beta: only components 1,2,3 are nonzero; rest = 0
        beta_true <- numeric(p)
        beta_true[1] <- 2.0
        beta_true[2] <- b2
        beta_true[3] <- b3
        
        # Generate y = X β + ε, ε ~ N(0, σₜᵣᵤₑ²)
        y <- as.vector(X %*% beta_true + rnorm(n, sd = sigma_true))
        
        # — 2) Split into train / test (half‐half)
        X_test  <- X[       1:n2,     , drop = FALSE]
        y_test  <- y[       1:n2]
        X_train <- X[-(1:n2),        , drop = FALSE]
        y_train <- y[-(1:n2)]
        
        # — 3) Null model on TRAIN (intercept only)
        beta0_hat  <- c(mean(y_train), rep(0, k))
        resid0_tr  <- y_train - (X_train %*% beta0_hat)
        sigma0_hat <- sqrt(sum(resid0_tr^2) / n2)
        resid0_te  <- y_test  - (X_test  %*% beta0_hat)
        loglik0    <- sum(dnorm(resid0_te, mean = 0, sd = sigma0_hat, log = TRUE))
        
        # — 4) Alternative = LASSO on TRAIN (glmnet CV)
        cvfit       <- cv.glmnet(
                          x           = X_train[ , -1, drop = FALSE],  # drop intercept col
                          y           = y_train,
                          alpha       = 1,
                          intercept   = TRUE,
                          standardize = TRUE
                       )
        lam         <- cvfit$lambda.1se
        fitted_tr   <- predict(cvfit, newx = X_train[ , -1, drop = FALSE], s = lam)
        sigma1_hat  <- sqrt(sum((y_train - fitted_tr)^2) / n2)
        fitted_te   <- predict(cvfit, newx = X_test[ , -1, drop = FALSE],     s = lam)
        loglik1     <- sum(dnorm(y_test - fitted_te, mean = 0, sd = sigma1_hat, log = TRUE))
        
        # — 5) Compute LC = exp(loglik0 – loglik1)
        LCs[i] <- exp(loglik0 - loglik1)
        LRs[i] <- exp(loglik1 - loglik0)
      }
      
      # — 6) Summarize: proportions of LC < thresholds
      prop_005 <- mean(LCs < 0.005)
      prop_01  <- mean(LCs < 0.01)
      prop_05  <- mean(LCs < 0.05)
      
      # Append to results
      results <- rbind(results, data.frame(
        n          = n,
        beta_lbl   = beta_lbl,
        sigma_true = sigma_true,
        prop_005   = prop_005,
        prop_01    = prop_01,
        prop_05    = prop_05,
        stringsAsFactors = FALSE
      ))
      
      # Print progress
      cat(sprintf(
        "  σ=%.2f → P(LC<0.005)=%.3f, P(LC<0.01)=%.3f, P(LC<0.05)=%.3f\n",
        sigma_true, prop_005, prop_01, prop_05
      ))
      
    } # end for sigma_true
  }   # end for each beta_setting
}     # end for each n

############################################################
# 7) Melt results for plotting
############################################################

#   We want three separate plots: 
#     – Plot A: threshold = “LC < 0.05”  (prop_05)
#     – Plot B: threshold = “LC < 0.01”  (prop_01)
#     – Plot C: threshold = “LC < 0.005” (prop_005)

# Melt into long form:
df_long <- melt(
  results,
  id.vars       = c("n", "beta_lbl", "sigma_true"),
  measure.vars  = c("prop_05", "prop_01", "prop_005"),
  variable.name = "threshold",
  value.name    = "proportion"
)

# Rename thresholds for prettier axis/legends
df_long$threshold <- factor(
  df_long$threshold,
  levels = c("prop_05", "prop_01", "prop_005"),
  labels = c("LC < 0.05", "LC < 0.01", "LC < 0.005")
)

############################################################
# 8) Plotting
############################################################

# We’ll draw three separate PNG files, one per threshold.
for (thr in levels(df_long$threshold)) {
  df_sub <- subset(df_long, threshold == thr)
  
  p <- ggplot(df_sub,
              aes(x = sigma_true,
                  y = proportion,
                  color = factor(n),
                  linetype = factor(beta_lbl)
              )
         ) +
       geom_line(size = 1) +
       geom_point(size = 1.5) +
       scale_color_brewer(palette = "Dark2") +
       scale_linetype_manual(
         values = c("solid", "dashed"),
         breaks = c("β₂=1.5,β₃=0.5", "β₂=0, β₃=0"),
         labels = c("β₂=1.5,β₃=0.5", "β₂=0, β₃=0")
       ) +
       theme_minimal(base_size = 12) +
       labs(
         title = paste0("Proportion of LC < ", sub("LC <", "", thr)),
         x     = expression(sigma_true),
         y     = "Proportion of Simulations",
         color = "Sample size (n)",
         linetype = "β‐setting"
       ) +
       theme(
         plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
         axis.text.x  = element_text(angle = 45, hjust = 1),
         legend.position = "right"
       )
  
  fname <- paste0("Proportion_vs_sigma_", gsub("[ <]", "", thr), ".png")
  ggsave(filename = fname, plot = p,
         width = 7, height = 5, dpi = 150)
  cat("Saved plot for", thr, "→", fname, "\n")
}

