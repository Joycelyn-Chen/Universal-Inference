#— Install/load required package
if(!require(glmnet)) install.packages("glmnet")
library(glmnet)

set.seed(2025)

### 1) Data construction and split
n  <- 200          # total sample size (must be even)
k <- 100
p  <- k + 1      # intercept + 100 predictors
alpha <- 0.1       # test level
sigma_true <- 1    # true noise SD

# True data‐generating betas: intercept = 2, rest = 0
beta_true <- numeric(p)
beta_true[1] <- 2

# Simulate predictors (columns 2..p) and add intercept column
X <- cbind(1, matrix(rnorm(n*(p-1)), nrow=n))
y <- as.vector(X %*% beta_true + rnorm(n, 0, sigma_true))

# Split into two equal halves
n2 <- n/2
X_test  <- X[     1:n2, ]
y_test  <- y[     1:n2   ]
X_train <- X[(n2+1):n, ]
y_train <- y[(n2+1):n   ]

### 2) Fit null model (intercept only) on TRAINING
#   Under H0: beta = (beta0, 0, 0, …, 0).
beta0_hat       <- rep(0, p)
beta0_hat[1]    <- mean(y_train)
resid0_train    <- y_train - X_train %*% beta0_hat
sigma0_hat      <- sqrt(sum(resid0_train^2) / n2)

# Evaluate Null log-likelihood on TEST
resid0_test  <- y_test  - X_test  %*% beta0_hat
loglik0      <- -n2/2 * log(2*pi*sigma0_hat^2) -
                sum(resid0_test^2)/(2*sigma0_hat^2)

### 3) Fit alternative via LASSO on TRAINING
#    We let glmnet pick lambda by 5-fold CV.
cvfit   <- cv.glmnet(
               x      = X_train[,-1],   # drop intercept col
               y      = y_train,
               alpha  = 1,              # LASSO penalty
               intercept = TRUE,
               standardize = TRUE
           )
lambda_star  <- cvfit$lambda.1se   # conservative choice

# Extract full p-vector of coefficients (including intercept)
bmat         <- coef(cvfit, s=lambda_star)
beta1_hat    <- as.numeric(bmat)   # length p

# Noise variance under alternative (on TRAINING residuals)
fitted1_train <- predict(cvfit, newx=X_train[,-1], s=lambda_star)
sigma1_hat    <- sqrt(sum((y_train - fitted1_train)^2) / n2)

### 4) Evaluate alternative log-likelihood on TEST
fitted1_test <- predict(cvfit, newx=X_test[,-1], s=lambda_star)
loglik1      <- -n2/2 * log(2*pi*sigma1_hat^2) -
                sum((y_test - fitted1_test)^2)/(2*sigma1_hat^2)

### 5) Likelihood‐ratio test
LR <- exp(loglik1 - loglik0)

cat(sprintf("Null log-likelihood:      %.3f\n", loglik0))
cat(sprintf("Alt  log-likelihood:      %.3f\n", loglik1))
cat(sprintf("Likelihood ratio = LR:    %.3g\n", LR))
if(LR > 1/alpha) {
  cat("=> Reject H0 at level alpha\n")
} else {
  cat("=> Fail to reject H0\n")
}
