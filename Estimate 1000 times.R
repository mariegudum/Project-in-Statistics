set.seed(123)

# Parameters
T <- 1
N <- 500
dt <- T / N
t <- seq(0, T, length.out = N + 1)

# True parameters
true_params <- list(
  beta0 = 0.04,
  beta1 = -0.1,
  beta2 = 0.6,
  theta1 = 0.43,
  theta2 = 0.21,
  lambda10 = 80,
  lambda11 = 10,
  lambda20 = 40,
  lambda21 = 3
)

# Negative log-likelihood function (same as your current version)
neg_log_likelihood <- function(params, X) {
  beta0 <- params[1]
  beta1 <- -0.1
  beta2 <- 0.6
  theta1 <- 0.43
  theta2 <- 0.21
  lambda10 <- 80
  lambda11 <- 10
  lambda20 <- 40
  lambda21 <- 3
  ll <- 0
  for (k in 1:N) {
    S_k <- X[k, 1]
    R1tk <- X[k, 2]
    R10k <- X[k, 3]
    R2tk <- X[k, 4]
    R20k <- X[k, 5]
    R11k <- (R1tk - (1 - theta1) * R10k) / theta1
    sqrt_term <- sqrt(pmax(R2tk, 0))
    beta_term <- beta0 + beta1 * R1tk + beta2 * sqrt_term
    Fk <- c(
      0,
      - (1 - theta1) * lambda10 * R10k - theta1 * lambda11 * R11k,
      - lambda10 * R10k
    )
    Sigmak <- matrix(c(
      beta_term * S_k,
      ((1 - theta1) * lambda10 + theta1 * lambda11) * beta_term,
      lambda10 * beta_term
    ), ncol = 1)
    diff <- X[k + 1, 1:3] - X[k, 1:3] - dt * Fk
    SigmaSigmaT <- dt * Sigmak %*% t(Sigmak)
    epsilon <- 1e-8 
    SigmaSigmaT_reg <- SigmaSigmaT + epsilon * diag(nrow(SigmaSigmaT))
    inv_Sigma <- solve(SigmaSigmaT_reg)
    det_val <- det(SigmaSigmaT_reg)
    ll <- ll + t(diff) %*% inv_Sigma %*% diff + log(det_val)
  }
  return(-as.numeric(ll))
}

# Simulation and estimation function
simulate_and_estimate <- function() {
  S <- numeric(N + 1)
  R10 <- numeric(N + 1)
  R11 <- numeric(N + 1)
  R20 <- numeric(N + 1)
  R21 <- numeric(N + 1)
  sigma_t <- numeric(N)
  
  S[1] <- 1
  R10[1] <- -0.10392006
  R11[1] <- 0.05214729
  R20[1] <- 0.00655017
  R21[1] <- 0.01580152
  
  for (i in 1:N) {
    R1t <- (1 - true_params$theta1) * R10[i] + true_params$theta1 * R11[i]
    R2t <- (1 - true_params$theta2) * R20[i] + true_params$theta2 * R21[i]
    sigma_t[i] <- true_params$beta0 + true_params$beta1 * R1t + true_params$beta2 * sqrt(max(R2t, 0))
    dW <- rnorm(1, mean = 0, sd = sqrt(dt))
    S[i + 1] <- S[i] * exp(-0.5 * sigma_t[i]^2 * dt + sigma_t[i] * dW)
    R10[i + 1] <- R10[i] + true_params$lambda10 * (-R10[i]) * dt + true_params$lambda10 * sigma_t[i] * dW
    R11[i + 1] <- R11[i] + true_params$lambda11 * (-R11[i]) * dt + true_params$lambda11 * sigma_t[i] * dW
    R20[i + 1] <- R20[i] + true_params$lambda20 * (sigma_t[i]^2 - R20[i]) * dt
    R21[i + 1] <- R21[i] + true_params$lambda21 * (sigma_t[i]^2 - R21[i]) * dt
  }
  
  R1t <- (1 - true_params$theta1) * R10 + true_params$theta1 * R11
  R2t <- (1 - true_params$theta2) * R20 + true_params$theta2 * R21
  X <- cbind(S, R1t, R10, R2t, R20)
  
  init_guess <- c(0.04)
  lower_bounds <- c(0.03)
  upper_bounds <- c(0.05)
  
  fit <- optim(par = init_guess, fn = neg_log_likelihood, method = "L-BFGS-B",
               lower = lower_bounds, upper = upper_bounds, X = X)
  
  return(fit$par)
}

# Run the estimation 1000 times
n_iter <- 100
estimates <- numeric(n_iter)

for (i in 1:n_iter) {
  estimates[i] <- simulate_and_estimate()
  if (i %% 50 == 0) cat("Iteration", i, "done\n")
}

hist(estimates, breaks = 40, main = "Distribution of Estimated beta0",
     xlab = expression(hat(beta)[0]), col = "lightblue", border = "gray")
abline(v = true_params$beta0, col = "red", lwd = 2, lty = 2) 
abline(v = ci_lower, col = "darkgreen", lwd = 2, lty = 3)
abline(v = ci_upper, col = "darkgreen", lwd = 2, lty = 3)
legend("topright", legend = c("True beta0", "95% CI"),
       col = c("red", "darkgreen"), lty = c(2, 3), lwd = 2)

# Compute summary statistics
mean_beta0 <- mean(estimates)
sd_beta0 <- sd(estimates)

# 95% Confidence interval (using empirical quantiles)
ci_lower <- quantile(estimates, 0.025)
ci_upper <- quantile(estimates, 0.975)

cat("Mean estimate of beta0:", mean_beta0, "\n")
cat("Standard deviation of estimates:", sd_beta0, "\n")
cat("95% confidence interval (empirical): [", ci_lower, ",", ci_upper, "]\n")
