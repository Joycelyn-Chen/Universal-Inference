#— dependencies
if(!require(glmnet)) install.packages("glmnet")
library(glmnet)

# set.seed(2048)

### fixed params
n           <- 100    # increase to 3000
k           <- 100
p           <- k + 1
alpha       <- 0.005
sigma_true  <- 2
n2          <- n/2

### simulate one design matrix once
X <- cbind(1, matrix(rnorm(n*(p-1)), nrow=n))

# 1) set up true betas & simulate y
beta_true    <- numeric(p)
beta_true[1] <- 2
beta_true[2] <- 1.5
beta_true[3] <- 0.5
y            <- as.vector(X %*% beta_true + rnorm(n, 0, sigma_true))

# split into test / train
y_test      <- y[    1:n2 ]
y_train     <- y[ -(1:n2) ]
X_test      <- X[    1:n2,]
X_train     <- X[ -(1:n2),]

LCs <- numeric(30)
LRs <- numeric(30)

for (i in 1:30){
  # 2) null fit on TRAIN (intercept only)
  beta0_hat       <- rep(0, p)
  beta0_hat[1]    <- mean(y_train)
  resid0_train    <- y_train - X_train %*% beta0_hat
  sigma0_hat      <- sqrt(sum(resid0_train^2) / n2)

  # null LL on TEST
  resid0_test <- y_test - X_test %*% beta0_hat
  loglik0     <- sum(dnorm(resid0_test, mean=0, sd=sigma0_hat, log=TRUE))

  # 3) alternative = LASSO on TRAIN
  cvfit      <- cv.glmnet(
                    x           = X_train[,-1],  # drop intercept col
                    y           = y_train,
                    alpha       = 1,
                    intercept   = TRUE,
                    standardize = TRUE
                )
  lambda_star <- cvfit$lambda.1se

  coef_star <- coef(cvfit, s = lambda_star)
  beta1_hat <- as.vector(coef_star)

  # re-calc noise var on TRAIN
  fitted1_train <- predict(cvfit, newx = X_train[,-1], s = lambda_star)
  sigma1_hat    <- sqrt(sum((y_train - fitted1_train)^2) / n2)

  # 4) alt log-LL on TEST
  fitted1_test <- predict(cvfit, newx = X_test[,-1], s = lambda_star)
  loglik1      <- sum(dnorm(y_test - fitted1_test, mean=0, sd=sigma1_hat, log=TRUE))
    
    # 5) split-LRT “p-value”
  LCs[i] <- exp(loglik0 - loglik1)
  LRs[i] <- exp(loglik1 - loglik0)
}

cat(sprintf("Null log-likelihood:      %.3f\n", loglik0))
cat(sprintf("Alt  log-likelihood:      %.3f\n", loglik1))
cat(sprintf("Likelihood ratio = LR:    %.3g\n", LR))
cat(sprintf("p-value = LC:    %.20f\n", LC))
if(LR > (1/alpha)) {
  cat("=> Reject H0 at level alpha\n")
} else {
  cat("=> Fail to reject H0\n")
}

boxplot(LCs)