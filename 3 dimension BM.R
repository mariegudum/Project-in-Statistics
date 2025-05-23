set.seed(123)

T <- 1
N <- 500
dt <- T / N
t <- seq(0, T, length.out = N + 1)

# True parameters for simulation
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

# Initialize vectors
S <- numeric(N + 1)
R10 <- numeric(N + 1)
R11 <- numeric(N + 1)
R20 <- numeric(N + 1)
R21 <- numeric(N + 1)
sigma_t <- numeric(N)

# Initial values
S[1] <- 1
R10[1] <- -0.10392006
R11[1] <- 0.05214729
R20[1] <- 0.00655017
R21[1] <- 0.01580152

# Simulation
for (i in 1:N) {
  R1t <- (1 - true_params$theta1) * R10[i] + true_params$theta1 * R11[i]
  R2t <- (1 - true_params$theta2) * R20[i] + true_params$theta2 * R21[i]
  sigma_t[i] <- true_params$beta0 + true_params$beta1 * R1t + true_params$beta2 * sqrt(max(R2t, 0))
  dW1 <- rnorm(1, mean = 0, sd = sqrt(dt))
  dW2 <- rnorm(1, mean = 0, sd = sqrt(dt))
  dW3 <- rnorm(1, mean = 0, sd = sqrt(dt))
  S[i + 1] <- S[i] * exp(-0.5 * sigma_t[i]^2 * dt + sigma_t[i] * dW1)
  R10[i + 1] <- R10[i] + true_params$lambda10 * (-R10[i]) * dt + true_params$lambda10 * sigma_t[i] * dW2
  R11[i + 1] <- R11[i] + true_params$lambda11 * (-R11[i]) * dt + true_params$lambda11 * sigma_t[i] * dW3
  R20[i + 1] <- R20[i] + true_params$lambda20 * (sigma_t[i]^2 - R20[i]) * dt
  R21[i + 1] <- R21[i] + true_params$lambda21 * (sigma_t[i]^2 - R21[i]) * dt
}
plot(t, S, type = 'l', col = 'blue', main = 'Asset Price Simulation', xlab = 'Time', ylab = 'S(t)')


R1t <- (1 - true_params$theta1) * R10 + true_params$theta1 * R11
R2t <- (1 - true_params$theta2) * R20 + true_params$theta2 * R21

X <- cbind(S, R1t, R10, R2t, R20)

# Negative log-likelihood function
neg_log_likelihood <- function(params) {
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
    Sigmak <- matrix(0, nrow = 5, ncol = 5)
    Sigmak[1, 1] <- beta_term * S_k
    Sigmak[2, 2] <- ((1 - theta1) * lambda10 + theta1 * lambda11) * beta_term
    Sigmak[3, 3] <- lambda10 * beta_term
    diff <- X[k + 1, 1:3] - X[k, 1:3] - dt * Fk
    SigmaSigmaT <- dt * Sigmak %*% t(Sigmak)
    epsilon <- 1e-8 
    SigmaSigmaT_reg <- SigmaSigmaT + epsilon * diag(nrow(SigmaSigmaT))
    inv_Sigma <- solve(SigmaSigmaT_reg)
    det_val <- det(SigmaSigmaT_reg)
    ll <- ll + t(diff) %*% inv_Sigma[1:3,1:3] %*% diff + log(det_val)
  }
  return(as.numeric(ll))
}

init_guess <- c(0.03)
lower_bounds <- c(0.01)
upper_bounds <- c(0.05)

fit <- optim(par = init_guess, fn = neg_log_likelihood, method = "L-BFGS-B",
             lower = lower_bounds, upper = upper_bounds)


fit$par
fit$value

neg_log_likelihood(true_params$beta0)
neg_log_likelihood_vec<-Vectorize(neg_log_likelihood)
beta0_test<-seq(0.03,0.05,by=0.001)
neg_log_likelihood_vec(beta0_test)
plot(beta0_test,neg_log_likelihood_vec(beta0_test))

n_simulations <- 500
beta0_estimates <- numeric(n_simulations)

# Simulation and estimation function
for (sim in 1:n_simulations) {
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
    sigma_t[i] <- true_params$beta0 + true_params$beta1 * R1t + true_params$beta2 * sqrt(pmax(R2t, 0))
    dW1 <- rnorm(1, mean = 0, sd = sqrt(dt))
    dW2 <- rnorm(1, mean = 0, sd = sqrt(dt))
    dW3 <- rnorm(1, mean = 0, sd = sqrt(dt))
    S[i + 1] <- S[i] * exp(-0.5 * sigma_t[i]^2 * dt + sigma_t[i] * dW1)
    R10[i + 1] <- R10[i] + true_params$lambda10 * (-R10[i]) * dt + true_params$lambda10 * sigma_t[i] * dW2
    R11[i + 1] <- R11[i] + true_params$lambda11 * (-R11[i]) * dt + true_params$lambda11 * sigma_t[i] * dW3
    R20[i + 1] <- R20[i] + true_params$lambda20 * (sigma_t[i]^2 - R20[i]) * dt
    R21[i + 1] <- R21[i] + true_params$lambda21 * (sigma_t[i]^2 - R21[i]) * dt
  }
  
  R1t <- (1 - true_params$theta1) * R10 + true_params$theta1 * R11
  R2t <- (1 - true_params$theta2) * R20 + true_params$theta2 * R21
  X <- cbind(S, R1t, R10, R2t, R20)
  
  init_guess <- c(0.03)
  lower_bounds <- c(0.03)
  upper_bounds <- c(0.05)
  
  fit <- optim(par = init_guess, fn = neg_log_likelihood, method = "L-BFGS-B",
               lower = lower_bounds, upper = upper_bounds)
  
  beta0_estimates[sim] <- fit$par
  if (sim %% 50 == 0) cat("Iteration", i, "done\n")
}

# Mean and 95% confidence interval
mean_estimate <- mean(beta0_estimates)
ci_lower <- quantile(beta0_estimates, 0.025)
ci_upper <- quantile(beta0_estimates, 0.975)

cat("Mean Estimate of beta0:", mean_estimate, "\n")
cat("95% Confidence Interval: [", ci_lower, ",", ci_upper, "]\n")

# Plot histogram
hist(beta0_estimates, breaks = 40, col = "lightblue", probability = TRUE,
     main = expression(paste("Histogram of ", hat(beta)[0], " Estimates")),
     xlab = expression(hat(beta)[0]), xlim = c(0.03, 0.05))
lines(density(beta0_estimates), col = "red", lwd = 2)
abline(v = true_params$beta0, col = "red", lwd = 2, lty = 2)
abline(v = ci_lower, col = "darkgreen", lwd = 2, lty = 3)
abline(v = ci_upper, col = "darkgreen", lwd = 2, lty = 3)
legend("topright", legend = c("True Value", "95% CI"), col = c("red", "darkgreen"),
       lty = c(2, 3), lwd = 2)

hist(beta0_estimates, breaks = 50, probability = TRUE,
     main = expression(paste("Histogram of ", hat(beta)[0], " Estimates")),
     xlab = "S(T)", col = "lightblue", border = "white")
lines(density(beta0_estimates), col = "red", lwd = 2)


abline(v = mean_estimate, col = "blue", lwd = 2, lty = 2)
abline(v = ci_lower, col = "darkgreen", lwd = 2, lty = 3)
abline(v = ci_upper, col = "darkgreen", lwd = 2, lty = 3)
abline(v = true_params$beta0, col = "darkred", lwd = 2, lty = 2)


legend("topleft", legend = c("Density", "Mean", "95% CI", "True value"),
       col = c("red", "blue", "darkgreen", "darkred"), lty = c(1, 2, 3, 2), lwd = 2)

