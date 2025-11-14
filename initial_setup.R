# Initial setup for simulating influence of margin knowledge on copula 
# model performance. We simulate from a 2D Gaussian copula and estimate
# the model parameters with known and unknown margins. Choose a sample
# size of 200 (N_sim=200) and run 100 Monte Carlo simulations (n=100).
# Based on the resulting parameter values, we compute the bias, variance
# and MSE for both cases and conclude the effect of the margin knowledge. 

# 1. Load necessary libraries
# 'copula' is essential for copula definitions, simulation, and fitting.
# 'stats' is needed for pnorm, var.
library(copula)
library(stats)

# 2. Set Simulation Parameters
N_sim <- 100       # Number of Monte Carlo simulations
n_obs <- 200       # Sample size for each simulation
rho_true <- 0.5    # True correlation parameter (Kendall's tau = 0.33)

# 3. Initialize Storage
# We will store the estimated correlation parameter from each simulation
est_rho_known <- numeric(N_sim)
est_rho_unknown <- numeric(N_sim)

# Define the true copula object
# A 2-dimensional Gaussian (normal) copula
true_copula <- normalCopula(dim = 2, param = rho_true, dispstr = "un")

# Define the true multivariate distribution (MVDC)
# This binds the Gaussian copula to standard normal margins
true_mvdc <- mvdc(
  copula = true_copula,
  margins = c("norm", "norm"),
  paramMargins = list(list(mean = 0, sd = 1), list(mean = 0, sd = 1))
)

# 4. Run the Simulation Loop
set.seed(1234) # for reproducibility
for (i in 1:N_sim) {
  
  # Generate n_obs from the true multivariate distribution
  # X_data is (n_obs x 2) matrix with N(0,1) margins
  X_data <- rMvdc(n_obs, true_mvdc)
  
  # --- Case 1: Known Margins ---
  # We know the margins are N(0,1), so we transform X_data
  # using the true marginal CDF, pnorm().
  U_known <- pnorm(X_data)
  
  # Fit the copula using Maximum Likelihood (ML) on the true uniform data
  # We initialize the copula to be fit (param=0.5 is just a starting value)
  fit_known <- fitCopula(normalCopula(dim = 2, dispstr = "un"),
                         data = U_known,
                         method = "ml")
  est_rho_known[i] <- coef(fit_known)
  
  
  # --- Case 2: Unknown Margins ---
  # We pretend we don't know the margins and must estimate them.
  # We convert X_data to pseudo-observations using the empirical CDF.
  # This is the standard non-parametric transformation.
  U_unknown <- pobs(X_data)
  
  # Fit the copula using ML on the pseudo-observations
  fit_unknown <- fitCopula(normalCopula(dim = 2, dispstr = "un"),
                           data = U_unknown,
                           method = "ml")
  est_rho_unknown[i] <- coef(fit_unknown)
  
} # End of simulation loop

# 5. Analyze Results

# --- Known Margins Metrics ---
bias_known <- mean(est_rho_known) - rho_true
var_known <- var(est_rho_known)
mse_known <- mean((est_rho_known - rho_true)^2) # Bias^2 + Variance

# --- Unknown Margins Metrics ---
bias_unknown <- mean(est_rho_unknown) - rho_true
var_unknown <- var(est_rho_unknown)
mse_unknown <- mean((est_rho_unknown - rho_true)^2) # Bias^2 + Variance

# 6. Format and Print Summary

# Use sprintf() to format the entire output as a single string.
# This provides more control and produces a cleaner block.
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

# Print the single, formatted string to the console.
cat(output_string)


# 7. Visualize Results with ggplot2

# Load the ggplot2 library (you may need to install it: install.packages("ggplot2"))
library(ggplot2)

# --- Create a 'long' format data.frame suitable for ggplot2 ---
# This data frame will have two columns:
# 'Estimate' (the rho value)
# 'Case' (a factor: "Known Margins" or "Unknown Margins")
results_df <- data.frame(
  Estimate = c(est_rho_known, est_rho_unknown),
  Case = factor(rep(c("Known Margins", "Unknown Margins"), each = N_sim),
                levels = c("Known Margins", "Unknown Margins"))
)

# --- Plot 1: Overlaid Density Plots ---
# This plot is excellent for seeing the central tendency and shape
plot_density <- ggplot(results_df, aes(x = Estimate, fill = Case)) +
  # Draw the density curves with transparency
  geom_density(alpha = 0.7) +
  
  # Add a vertical dashed line for the true parameter value
  geom_vline(xintercept = rho_true, linetype = "dashed", color = "black", linewidth = 1) +
  
  # Add labels and a clean theme
  labs(
    title = "Distribution of Parameter Estimates (n=200, N=100)",
    x = "Estimated Rho (ρ)",
    y = "Density",
    fill = "Estimation Case"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# Print the density plot
print(plot_density)


# --- Plot 2: Side-by-Side Boxplots ---
# This plot is excellent for comparing variance (the height of the boxes)
plot_boxplot <- ggplot(results_df, aes(x = Case, y = Estimate, fill = Case)) +
  # Draw the boxplots
  geom_boxplot() +
  
  # Add a horizontal dashed line for the true parameter value
  geom_hline(yintercept = rho_true, linetype = "dashed", color = "red", linewidth = 1) +
  
  # Add labels and a clean theme
  labs(
    title = "Comparison of Estimator Spread (n=200, N=100)",
    x = "Estimation Case",
    y = "Estimated Rho (ρ)"
  ) +
  theme_bw() +
  # Hide the legend as the x-axis is already descriptive
  theme(legend.position = "none")

# Print the boxplot
print(plot_boxplot)