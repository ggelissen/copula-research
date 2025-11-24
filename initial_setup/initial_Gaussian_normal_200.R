# Initial setup for simulating influence of margin knowledge on copula 
# model performance. We simulate from a 2D copula and estimate the 
# model parameters with known and unknown margins. Choose a sample size 
# (N_sim) and run 100 Monte Carlo simulations (n=100). Based on the 
# resulting parameter values, we compute the bias, variance and MSE for 
# both cases and conclude the effect of the margin knowledge. 

# Configuration
# Copula Family: Gaussian
# Marginal Dist: standard Normal
# Sample Size  : 200

# 1. Load necessary libraries
library(copula)
library(stats)

# 2. Set Simulation Parameters
N_sim <- 100       # Number of Monte Carlo simulations
n_obs <- 200       # Sample size for each simulation
rho_true <- 0.5    # True correlation parameter (Kendall's tau = 0.33)

# 3. Initialize Storage
est_rho_known <- numeric(N_sim)
est_rho_unknown <- numeric(N_sim)

true_copula <- normalCopula(dim = 2, param = rho_true, dispstr = "un")

true_mvdc <- mvdc(
  copula = true_copula,
  margins = c("norm", "norm"),
  paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1))
)

# 4. Run the Simulation Loop
set.seed(1234)
for (i in 1:N_sim) {
  
  X_data <- rMvdc(n_obs, true_mvdc)
  
  # --- Case 1: Known Margins ---
  U_known <- pnorm(X_data)
  fit_known <- fitCopula(normalCopula(dim = 2, dispstr = "un"),
                         data = U_known,
                         method = "ml")
  est_rho_known[i] <- coef(fit_known)
  
  
  # --- Case 2: Unknown Margins ---
  U_unknown <- pobs(X_data)
  fit_unknown <- fitCopula(normalCopula(dim = 2, dispstr = "un"),
                           data = U_unknown,
                           method = "ml")
  est_rho_unknown[i] <- coef(fit_unknown)
  
}

# 5. Analyze Results
bias_known <- mean(est_rho_known) - rho_true
var_known <- var(est_rho_known)
mse_known <- mean((est_rho_known - rho_true)^2) 

bias_unknown <- mean(est_rho_unknown) - rho_true
var_unknown <- var(est_rho_unknown)
mse_unknown <- mean((est_rho_unknown - rho_true)^2)

# 6. Format and Print Summary
output_string <- sprintf(
"--- Simulation Results (N=%d, n=%d) ---

True Parameter (rho): %.3f

Case 1: KNOWN Margins (Parametric ML)
-----------------------------------------
Bias     : %f
Variance : %f
MSE      : %f

Case 2: UNKNOWN Margins (Semi-parametric ML)
-----------------------------------------
Bias     : %f
Variance : %f
MSE      : %f
",
  N_sim, n_obs, rho_true,
  bias_known, var_known, mse_known,
  bias_unknown, var_unknown, mse_unknown
)

cat(output_string)


# 7. Visualize Results
library(ggplot2)

results_df <- data.frame(
  Estimate = c(est_rho_known, est_rho_unknown),
  Case = factor(rep(c("Known Margins", "Unknown Margins"), each = N_sim),
                levels = c("Known Margins", "Unknown Margins"))
)

plot_density <- ggplot(results_df, aes(x = Estimate, fill = Case)) +

  geom_density(alpha = 0.7) +
  
  geom_vline(xintercept = rho_true, linetype = "dashed", color = "black", size = 1) +
  
  labs(
    title = "Distribution of Parameter Estimates (n=200, N=100)",
    x = "Estimated Rho",
    y = "Density",
    fill = "Estimation Case"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(plot_density)
ggsave("figs/initial_setup/density_Gaussian_normal_200.png", plot = plot_density, device = "png", width = 8, height = 5)


plot_boxplot <- ggplot(results_df, aes(x = Case, y = Estimate, fill = Case)) +
  geom_boxplot() +
  
  geom_hline(yintercept = rho_true, linetype = "dashed", color = "red", size = 1) +
  
  labs(
    title = "Comparison of Estimator Spread (n=200, N=100)",
    x = "Estimation Case",
    y = "Estimated Rho"
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(plot_boxplot)
ggsave("figs/initial_setup/boxplot_Gaussian_normal_200.png", plot = plot_boxplot, device = "png", width = 6, height = 5)