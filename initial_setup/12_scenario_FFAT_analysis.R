# This script replaces the two old OFAT bar charts with a
# single, definitive Full Factorial experiment.
#
# It analyzes the interaction between:
# 1. Copula Family (4 types)
# 2. Marginal Distribution (3 types)
#
# It produces grouped bar charts and a heatmap to
# visualize the interaction effects.

# 1. Load Libraries
library(copula)
library(stats)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- 2. Global Simulation Parameters ---
# We fix n=200 as a representative sample size
FIXED_N_OBS <- 200
# We use N_sim=500 for low-noise results
N_sim <- 500
# We fix tau=1/3
TRUE_TAU <- 1/3

# Add fixed parameters
T_COPULA_DF <- 4
T_MARGIN_DF <- 3
GAMMA_MARGIN_SHAPE <- 2
GAMMA_MARGIN_RATE <- 1

# --- 3. Parameter Grid ---
# This is a full factorial (4 x 3 = 12 scenarios)
param_grid <- expand.grid(
  family = factor(c("normal", "t", "gumbel", "clayton")),
  margins = factor(c("norm", "t", "gamma")),
  stringsAsFactors = TRUE
)

# We add n_obs to the grid
param_grid$n_obs <- FIXED_N_OBS

cat("Starting Full Factorial 2D simulation.\n")
cat("N_sim =", N_sim, ", n_obs =", FIXED_N_OBS, "\n")
cat("Scenarios to run:", nrow(param_grid), "\n")
cat("This will take a few minutes.\n")


# --- 4. Main Simulation Function (Identical to V5) ---

run_scenario <- function(N_sim, n_obs, family, margins, true_tau) {
  
  # --- 4.1. Define Copula ---
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

  # --- 4.2. Define Margins & MVDC ---
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

  # --- 4.3. Run N_sim loops ---
  for (i in 1:N_sim) {
    X_data <- rMvdc(n_obs, true_mvdc)
    U_known <- p_margin_fn(X_data)
    U_unknown <- pobs(X_data)
    
    # --- 4.4. Fit Copulas ---
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
  
  # --- 4.5. Calculate and return metrics ---
  est_param_known <- na.omit(est_param_known)
  est_param_unknown <- na.omit(est_param_unknown)
  
  if (length(est_param_known) < 2 || length(est_param_unknown) < 2) {
    return(data.frame(
      bias_known = NA, var_known = NA, mse_known = NA,
      bias_unknown = NA, var_unknown = NA, mse_unknown = NA
    ))
  }
  
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


# --- 5. Run the Full Experiment ---

results_list <- vector("list", nrow(param_grid))
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

results_df <- do.call(rbind, results_list)
results_df <- na.omit(results_df)

# --- 6. Process and Plot Results ---

# 6.1. Create 'MSE_Inflation' metric
results_df <- results_df %>%
  mutate(
    MSE_Inflation = mse_unknown / mse_known,
    # Round for better heatmap labels
    MSE_Inflation_Label = round(MSE_Inflation, 2)
  )

# --- Plot 1: Effect of Copula Family (Grouped Bar Chart) ---
gg_plot1 <- ggplot(results_df, aes(x = family, y = MSE_Inflation, fill = margins)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "red") +
  labs(
    title = "Effect of Copula Family on 'Cost of Ignorance'",
    subtitle = paste("Grouped by Marginal Distribution (n=", FIXED_N_OBS, ", N_sim=", N_sim, ")", sep=""),
    x = "Copula Family",
    y = "MSE Inflation (MSE_Unknown / MSE_Known)",
    fill = "Margins"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(gg_plot1)


# --- Plot 2: Effect of Marginal Distribution (Grouped Bar Chart) ---
gg_plot2 <- ggplot(results_df, aes(x = margins, y = MSE_Inflation, fill = family)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "red") +
  labs(
    title = "Effect of Margins on 'Cost of Ignorance'",
    subtitle = paste("Grouped by Copula Family (n=", FIXED_N_OBS, ", N_sim=", N_sim, ")", sep=""),
    x = "Marginal Distribution",
    y = "MSE Inflation (MSE_Unknown / MSE_Known)",
    fill = "Family"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(gg_plot2)


# --- Plot 3: Interaction Heatmap (The Best Visualization) ---
gg_plot3 <- ggplot(results_df, aes(x = family, y = margins, fill = MSE_Inflation)) +
  # geom_tile() creates the heatmap
  geom_tile(color = "white") +
  
  # Add the text label in the center of the cell
  geom_text(aes(label = MSE_Inflation_Label), color = "black", size = 4) +
  
  # Use a better color scale (e.g., from blue to red)
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 1.0, limit = c(min(results_df$MSE_Inflation), max(results_df$MSE_Inflation))) +
  labs(
    title = "Interaction: Cost of Margin Ignorance",
    subtitle = paste("MSE Inflation (Unknown/Known) at n=", FIXED_N_OBS, ", N_sim=", N_sim, sep=""),
    x = "Copula Family",
    y = "Marginal Distribution",
    fill = "MSE Inflation"
  ) +
  theme_minimal() + # A clean theme for heatmaps
  theme(legend.position = "right")

print(gg_plot3)