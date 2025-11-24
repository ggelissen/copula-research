# This script performs an analysis of the effect
# of sample size (n) on estimator MSE and MSE Inflation for
# a single, baseline scenario: Normal Copula + Normal Margins.

# Load Libraries
library(copula)
library(stats)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- Global Simulation Parameters ---
N_sim <- 100         # Monte Carlo replications
TRUE_TAU <- 1/3     # True Kendall's tau (constant across families)

# Add fixed parameters (will be used if we add t-copulas later)
T_COPULA_DF <- 4
T_MARGIN_DF <- 3
GAMMA_MARGIN_SHAPE <- 2
GAMMA_MARGIN_RATE <- 1

# --- Parameter Grid ---
# Create a high-resolution grid of n_obs values
n_values <- seq(50, 1000, by = 50)

param_grid <- data.frame(
  n_obs = n_values,
  family = "normal",
  margins = "norm",
  stringsAsFactors = FALSE
)

cat("Starting n-convergence simulation.\n")
cat("Scenarios to run:", nrow(param_grid), "\n")
cat("This will take a significant amount of time (est. 20-30 min).\n")


# --- Main Simulation Function (Identical to V3) ---

run_scenario <- function(N_sim, n_obs, family, margins, true_tau) {
  
  # --- Define Copula ---
  if (family == "normal") {
    true_param <- iTau(ellipCopula("normal"), true_tau)
    copula_obj <- normalCopula(dim = 2, param = true_param, dispstr = "un")
  } else if (family == "t") {
    true_param <- iTau(ellipCopula("t"), true_tau)
    copula_obj <- tCopula(dim = 2, param = true_param, df = T_COPULA_DF, dispstr = "un")
  } else if (family == "gumbel") {
    true_param <- iTau(gumbelCopula(2), true_tau) 
    copula_obj <- gumbelCopula(dim = 2, param = true_param)
  } else if (family == "clayton") {
    true_param <- iTau(claytonCopula(2), true_tau) 
    copula_obj <- claytonCopula(dim = 2, param = true_param)
  }

  # --- Define Margins & MVDC ---
  if (margins == "norm") {
    margin_names <- c("norm", "norm")
    margin_params <- list(list(mean = 0, sd = 1), list(mean = 0, sd = 1))
    p_margin_fn <- function(x) pnorm(x, 0, 1)
  } else if (margins == "t") {
    margin_names <- c("t", "t")
    margin_params <- list(list(df = T_MARGIN_DF), list(df = T_MARGIN_DF))
    p_margin_fn <- function(x) pt(x, T_MARGIN_DF)
  } else if (margins == "gamma") {
    margin_names <- c("gamma", "gamma")
    margin_params <- list(list(shape = GAMMA_MARGIN_SHAPE, rate = GAMMA_MARGIN_RATE),
                          list(shape = GAMMA_MARGIN_SHAPE, rate = GAMMA_MARGIN_RATE))
    p_margin_fn <- function(x) pgamma(x, shape = GAMMA_MARGIN_SHAPE, rate = GAMMA_MARGIN_RATE)
  }
  
  true_mvdc <- mvdc(
    copula = copula_obj,
    margins = margin_names,
    paramMargins = margin_params
  )
  
  est_param_known <- numeric(N_sim)
  est_param_unknown <- numeric(N_sim)

  # --- Run N_sim loops ---
  for (i in 1:N_sim) {
    X_data <- rMvdc(n_obs, true_mvdc)
    U_known <- p_margin_fn(X_data)
    U_unknown <- pobs(X_data)
    
    # --- Fit Copulas ---
    if (family == "normal") { 
      fit_cop_obj <- normalCopula(dim = 2, dispstr = "un")
    } else if (family == "t") {
      fit_cop_obj <- tCopula(dim = 2, dispstr = "un", df = T_COPULA_DF)
    } else if (family == "gumbel") {
      fit_cop_obj <- gumbelCopula(dim = 2, param = 2) 
    } else if (family == "clayton") {
      fit_cop_obj <- claytonCopula(dim = 2, param = 2) 
    }
    
    fit_known <- try(fitCopula(fit_cop_obj, data = U_known, method = "ml"), silent = TRUE)
    if (!inherits(fit_known, "try-error")) {
      est_param_known[i] <- coef(fit_known)[1]
    } else { est_param_known[i] <- NA }
    
    fit_unknown <- try(fitCopula(fit_cop_obj, data = U_unknown, method = "ml"), silent = TRUE)
    if (!inherits(fit_unknown, "try-error")) {
      est_param_unknown[i] <- coef(fit_unknown)[1]
    } else { est_param_unknown[i] <- NA }
  }
  
  # --- Calculate and return metrics ---
  est_param_known <- na.omit(est_param_known)
  est_param_unknown <- na.omit(est_param_unknown)
  
  metrics <- data.frame(
    bias_known     = mean(est_param_known) - true_param,
    var_known      = var(est_param_known),
    mse_known      = mean((est_param_known - true_param)^2),
    bias_unknown   = mean(est_param_unknown) - true_param,
    var_unknown    = var(est_param_unknown),
    mse_unknown    = mean((est_param_unknown - true_param)^2)
  )
  return(metrics)
}


# --- Run the Experiment ---

# Initialize list to store results
results_list <- vector("list", nrow(param_grid))

# Run the main loop
pb <- txtProgressBar(min = 0, max = nrow(param_grid), style = 3)

for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  
  metrics <- run_scenario(
    N_sim = N_sim,
    n_obs = params$n_obs,
    family = as.character(params$family), 
    margins = as.character(params$margins),
    true_tau = TRUE_TAU
  )
  
  results_list[[i]] <- cbind(params, metrics)
  
  setTxtProgressBar(pb, i)
}
close(pb)

cat("\nSimulation complete.\n")

# Combine all results into one master data.frame
results_df <- do.call(rbind, results_list)

# --- Process and Plot Results ---

# Create 'MSE_Inflation' metric
results_df <- results_df %>%
  mutate(
    MSE_Inflation = mse_unknown / mse_known
  )

# Pivot for plotting
results_long_mse <- results_df %>%
  select(n_obs, family, margins, mse_known, mse_unknown) %>%
  pivot_longer(
    cols = c("mse_known", "mse_unknown"),
    names_to = "Case",
    values_to = "MSE",
    names_prefix = "mse_"
  ) %>%
  mutate(Case = factor(Case, levels = c("known", "unknown")))


# --- Absolute MSE Decay (Log Scale) ---
gg_mse_decay <- ggplot(results_long_mse, aes(x = n_obs, y = MSE, color = Case)) +
  geom_line(size = 1) +
  geom_point(alpha = 0.5) +
  # Using a log-scale on y-axis is essential to see the effect
  scale_y_log10(breaks = c(0.1, 0.01, 0.001, 0.0001), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(breaks = seq(0, 1000, by = 200)) +
  labs(
    title = "Absolute MSE Decay vs. Sample Size",
    subtitle = "Scenario: Normal Copula, Normal Margins (True Tau = 1/3)",
    x = "Sample Size (n)",
    y = "Mean Squared Error (MSE, log scale)"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(gg_mse_decay)


# --- Relative Cost (MSE Inflation) Decay ---
gg_cost_decay <- ggplot(results_df, aes(x = n_obs, y = MSE_Inflation)) +
  geom_line(size = 1, color = "blue") +
  geom_point(alpha = 0.5, color = "blue") +
  # Add a dashed red line at 1.0 (no cost)
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "red", size = 1) +
  scale_x_continuous(breaks = seq(0, 1000, by = 200)) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.2)) + # Adjust 'by' as needed
  coord_cartesian(ylim = c(0.8, NA)) + # Zoom in, but allow to see 1.0
  labs(
    title = "Relative Cost of Unknown Margins vs. Sample Size",
    subtitle = "Scenario: Normal Copula, Normal Margins (True Tau = 1/3)",
    x = "Sample Size (n)",
    y = "MSE Inflation (MSE_Unknown / MSE_Known)"
  ) +
  theme_bw()

print(gg_cost_decay)



