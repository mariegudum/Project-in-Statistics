library(ggplot2)
set.seed(123)

T <- 1       
N <- 50000     
dt <- T/N      
t <- seq(0, T, length.out = N + 1)
beta0 <- 0.04
beta1 <- -0.1
beta2 <- 0.6
theta1 <- 0.43
theta2 <- 0.21
lambda10 <- 80
lambda11 <- 10
lambda20 <- 40
lambda21 <- 3

S <- numeric(N + 1)
R10 <- numeric(N + 1)
R11 <- numeric(N + 1)
R20 <- numeric(N + 1)
R21 <- numeric(N + 1)
sigma_t <- numeric(N + 1)

S[1] <- 1
R10[1] <- -0.10392006
R11[1] <- 0.05214729
R20[1] <- 0.00655017
R21[1] <- 0.01580152

for (i in 1:N) {
  R1t <- (1 - theta1) * R10[i] + theta1 * R11[i]
  R2t <- (1 - theta2) * R20[i] + theta2 * R21[i]
  sigma_t[i] <- beta0 + beta1 * R1t + beta2 * sqrt(max(R2t, 0))
  dW <- rnorm(1, mean = 0, sd = sqrt(dt))
  S[i + 1] <- S[i] * exp(-0.5 * sigma_t[i]^2 * dt + sigma_t[i] * dW)
  dR10 <- lambda10 * (-R10[i]) * dt + lambda10 * sigma_t[i] * dW
  dR11 <- lambda11 * (-R11[i]) * dt + lambda11 * sigma_t[i] * dW
  dR20 <- lambda20 * (sigma_t[i]^2 - R20[i]) * dt
  dR21 <- lambda21 * (sigma_t[i]^2 - R21[i]) * dt
  R10[i + 1] <- R10[i] + dR10
  R11[i + 1] <- R11[i] + dR11
  R20[i + 1] <- R20[i] + dR20
  R21[i + 1] <- R21[i] + dR21
}
sigma_t[N + 1] <- beta0 + beta1 * R1t + beta2 * sqrt(max(R2t, 0))

plot(t, S, type = 'l', col = 'blue', main = 'Asset Price Simulation', xlab = 'Time', ylab = 'S(t)')
plot(t, sigma_t, type = 'l', col = 'darkred', main = 'Stochastic volatility', xlab = 'Time', ylab = expression(sigma[t]))

###### Estimation method 0

X <- cbind(S, R10, R11, R20, R21)

log_likelihood <- function(theta, X, dt) {
  beta0  <- theta[1]
  beta1  <- theta[2]
  beta2  <- theta[3]
  theta1 <- theta[4]
  theta2 <- theta[5]
  lambda10 <- theta[6]
  lambda11 <- theta[7]
  lambda20 <- theta[8]
  lambda21 <- theta[9]
  N <- nrow(X) - 1
  loglik <- 0
  for (k in 1:N) {
    S_k  <- X[k, 1]
    R10k <- X[k, 2]
    R11k <- X[k, 3]
    R20k <- X[k, 4]
    R21k <- X[k, 5]
    R1t <- (1 - theta1) * R10k + theta1 * R11k
    R2t <- (1 - theta2) * R20k + theta2 * R21k
    beta_term <- beta0 + beta1 * R1t + beta2 * sqrt(max(R2t, 0))
    beta_term_sq <- beta_term^2
    
    # F vector
    Fk <- c(
      0,
      -(1 - theta1) * lambda10 * R10k - theta1 * lambda11 * R11k,
      -lambda10 * R10k,
      (1 - theta2) * lambda20 * (beta_term_sq - R20k) + theta2 * lambda21 * (beta_term_sq - R21k),
      lambda20 * (beta_term_sq - R20k)
    )
    
    # Sigma vector
    Sigmak <- matrix(c(
      beta_term * S_k,
      ((1 - theta1) * lambda10 + theta1 * lambda11) * beta_term,
      lambda10 * beta_term,
      0,
      0
    ), ncol = 1)
    
    cov_k <- dt * (Sigmak %*% t(Sigmak))  # 5x5 covariance matrix
    mean_k <- X[k, ] + dt * Fk
    innov <- X[k + 1, ] - mean_k
    
    # Safeguard for numerical stability
    cov_k <- cov_k + 1e-8 * diag(5)
    
    # Mahalanobis term
    mahal <- t(innov) %*% solve(cov_k, innov)
    
    # Log determinant
    log_det <- log(det(cov_k))
    
    loglik <- loglik + mahal + log_det
  }
  return(as.numeric(loglik))
}

# Initial guess
theta_init <- c(0.03, -0.05, 0.5, 0.4, 0.2, 50, 10, 30, 5)

# Use L-BFGS-B to apply bounds (e.g., positivity constraints)
fit <- optim(
  par = theta_init,
  fn = log_likelihood,
  method = "L-BFGS-B",
  lower = c(0, -1, 0, 0, 0, 1, 1, 1, 1),  # some reasonable lower bounds
  upper = c(1, 1, 2, 1, 1, 200, 200, 200, 200),
  X = X,
  dt = dt
)

fit$par      # Estimated parameters
fit$value    # Final log-likelihood



fit <- optim(
  par = theta_init,
  fn = function(theta) -log_likelihood(theta, X, dt),  # negate for MLE
  method = "L-BFGS-B",
  lower = c(0, -1, 0, 0, 0, 1, 1, 1, 1),
  upper = c(1, 1, 2, 1, 1, 200, 200, 200, 200)
)

fit$par      # Estimated parameters
-fit$value    # Final log-likelihood



###### Estimation method 1 

result <- optim(
  par = init_params,
  fn = neg_log_likelihood,
  method = "BFGS",
  S = S,
  R10 = R10,
  R11 = R11,
  R20 = R20,
  R21 = R21,
  dt = dt,
  theta1 = theta1,
  theta2 = theta2
)

print(result$par)

###### Estimation method 2


mle_result <- mle(log_likelihood,
                  start = list(beta0 = 0.04, beta1 = -0.1, beta2 = 0.6,
                               lambda10 = 80, lambda11 = 10,
                               lambda20 = 40, lambda21 = 3),
                  method = "L-BFGS-B",
                  lower = c(0, -1, 0, 1, 1, 1, 1))
summary(mle_result)




### Simulate 1000 paths

library(ggplot2)
set.seed(123)

T <- 1       
N <- 50000     
dt <- T/N      
t <- seq(0, T, length.out = N + 1)

beta0 <- 0.04
beta1 <- -0.1
beta2 <- 0.6
theta1 <- 0.43
theta2 <- 0.21
lambda10 <- 80
lambda11 <- 10
lambda20 <- 40
lambda21 <- 3

n_sim <- 1000

# Storage for simulations (e.g., final value of S)
S_paths <- matrix(0, nrow = N + 1, ncol = n_sim)
sigma_paths <- matrix(0, nrow = N + 1, ncol = n_sim)

for (sim in 1:n_sim) {
  S <- numeric(N + 1)
  sigma_t <- numeric(N + 1)
  R10 <- numeric(N + 1)
  R11 <- numeric(N + 1)
  R20 <- numeric(N + 1)
  R21 <- numeric(N + 1)
  
  S[1] <- 1
  R10[1] <- -0.10392006
  R11[1] <- 0.05214729
  R20[1] <- 0.00655017
  R21[1] <- 0.01580152
  
  for (i in 1:N) {
    R1t <- (1 - theta1) * R10[i] + theta1 * R11[i]
    R2t <- (1 - theta2) * R20[i] + theta2 * R21[i]
    sigma_t[i] <- beta0 + beta1 * R1t + beta2 * sqrt(max(R2t, 0))
    dW <- rnorm(1, mean = 0, sd = sqrt(dt))
    S[i + 1] <- S[i] * exp(-0.5 * sigma_t[i]^2 * dt + sigma_t[i] * dW)
    dR10 <- lambda10 * (-R10[i]) * dt + lambda10 * sigma_t[i] * dW
    dR11 <- lambda11 * (-R11[i]) * dt + lambda11 * sigma_t[i] * dW
    dR20 <- lambda20 * (sigma_t[i]^2 - R20[i]) * dt
    dR21 <- lambda21 * (sigma_t[i]^2 - R21[i]) * dt
    R10[i + 1] <- R10[i] + dR10
    R11[i + 1] <- R11[i] + dR11
    R20[i + 1] <- R20[i] + dR20
    R21[i + 1] <- R21[i] + dR21
  }
  
  sigma_t[N + 1] <- beta0 + beta1 * R1t + beta2 * sqrt(max(R2t, 0))
  S_paths[, sim] <- S
  sigma_paths[, sim] <- sigma_t
}


matplot(t, S_paths[, 1:10], type = 'l', lty = 1, col = rainbow(10),
        main = '10 Sample Simulated Asset Price Paths',
        xlab = 'Time', ylab = 'S(t)')


mean_sigma <- rowMeans(sigma_paths)
plot(t, mean_sigma, type = 'l', col = 'red', lwd = 2,
     main = 'Average Stochastic Volatility Path',
     xlab = 'Time', ylab = expression(sigma[t]))


S_T <- S_paths[N + 1, ]

hist(S_T, breaks = 50, probability = TRUE,
     main = "Distribution of Terminal Asset Price S(T)",
     xlab = "S(T)", col = "lightblue", border = "white")
lines(density(S_T), col = "red", lwd = 2)

mean_ST <- mean(S_T)
sd_ST <- sd(S_T)
error_margin <- qnorm(0.975) * sd_ST / sqrt(length(S_T))
ci_lower <- mean_ST - error_margin
ci_upper <- mean_ST + error_margin

abline(v = mean_ST, col = "blue", lwd = 2, lty = 2)
abline(v = ci_lower, col = "darkgreen", lwd = 2, lty = 3)
abline(v = ci_upper, col = "darkgreen", lwd = 2, lty = 3)

legend("topleft", legend = c("Density", "Mean", "95% CI"),
       col = c("red", "blue", "darkgreen"), lty = c(1, 2, 3), lwd = 2)

# Print statistics
cat("Mean of S(T):", mean_ST, "\n")
cat("Standard deviation of S(T):", sd_ST, "\n")
cat("95% Confidence Interval: [", ci_lower, ",", ci_upper, "]\n")


## F and sigma function

# sigma = beta_term
beta_term <- beta0 + beta1 * R1t + beta2 * sqrt(max(R2t, 0))
beta_term_sq <- beta_term^2

# Vector F
F <- c(
  0,
  -(1 - theta1) * lambda10 * R10[i] - theta1 * lambda11 * R11[i],
  -lambda10 * R10[i],
  (1 - theta2) * lambda20 * (beta_term_sq - R20[i]) + theta2 * lambda21 * (beta_term_sq - R21[i]),
  lambda20 * (beta_term_sq - R20[i])
)


# Vector Sigma
Sigma <- c(
  beta_term * S,
  ((1 - theta1) * lambda10 + theta1 * lambda11) * beta_term,
  lambda10 * beta_term,
  0,
  0
)



