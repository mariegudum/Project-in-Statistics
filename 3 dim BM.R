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
  return(-as.numeric(ll))
}

init_guess <- c(0.04)
lower_bounds <- c(0)
upper_bounds <- c(0.5)

fit <- optim(par = init_guess, fn = neg_log_likelihood, method = "L-BFGS-B",
             lower = lower_bounds, upper = upper_bounds)

neg_log_likelihood(true_params$beta0)
neg_log_likelihood_vec<-Vectorize(neg_log_likelihood)
beta0_test<-seq(0,0.1,by=0.01)
neg_log_likelihood_vec(beta0_test)
plot(beta0_test,neg_log_likelihood_vec(beta0_test))


fit$par
fit$value



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
    dW1 <- rnorm(1, mean = 0, sd = sqrt(dt))
    dW2 <- rnorm(1, mean = 0, sd = sqrt(dt))
    dW3 <- rnorm(1, mean = 0, sd = sqrt(dt))
    S[i + 1] <- S[i] * exp(-0.5 * sigma_t[i]^2 * dt + sigma_t[i] * dW1)
    dR10 <- lambda10 * (-R10[i]) * dt + lambda10 * sigma_t[i] * dW2
    dR11 <- lambda11 * (-R11[i]) * dt + lambda11 * sigma_t[i] * dW3
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

