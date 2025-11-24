# Investigates the influence of:
# 1. Copula Family (Normal, t, Gumbel, Clayton)
# 2. Marginal Distribution (Normal, t, Gamma)
# 3. Sample Size (n=50, 200, 1000)
# on the cost of unknown margins for copula estimation.
#
# This version uses an 8-scenario "One-Factor-at-a-Time" (OFAT)
# grid for a fast, intermediate analysis.

# 1. Load Libraries
library(copula)
library(stats)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- 2. Global Simulation Parameters ---
N_sim <- 100         # Monte Carlo replications (set >1000 for paper)
TRUE_TAU <- 1/3     # True Kendall's tau (constant across families)

# --- 3. Parameter Grid (OFAT 8-Scenario Version) ---

# Define the baseline scenario
baseline <- data.frame(
  n_obs = 200,
  family = "normal",  
  margins = "norm",
  stringsAsFactors = FALSE
)

# Excursion A: Varying n_obs
grid_n <- data.frame(
  n_obs = c(50, 1000), 
  family = "normal",  
  margins = "norm",
  stringsAsFactors = FALSE
)

# Excursion B: Varying family
grid_family <- data.frame(
  n_obs = 200,
  family = c("t", "gumbel", "clayton"), # "normal" is in baseline
  margins = "norm",
  stringsAsFactors = FALSE
)

# Excursion C: Varying margins
grid_margins <- data.frame(
  n_obs = 200,
  family = "normal",  
  margins = c("t", "gamma"), 
  stringsAsFactors = FALSE
)

# Combine all scenarios
param_grid <- rbind(baseline, grid_n, grid_family, grid_margins)

# Ensure factor levels are set correctly for plotting later
param_grid$family <- factor(param_grid$family, levels = c("normal", "t", "gumbel", "clayton"))
param_grid$margins <- factor(param_grid$margins, levels = c("norm", "t", "gamma"))


# Add fixed parameters for "t" copula and "t" / "gamma" margins
T_COPULA_DF <- 4
T_MARGIN_DF <- 3
GAMMA_MARGIN_SHAPE <- 2
GAMMA_MARGIN_RATE <- 1


# --- 4. Main Simulation Function ---

#' run_scenario
#'
#' Runs a full Monte Carlo simulation (N_sim loops) for a single
#' parameter combination.
#'
#' @param N_sim Number of loops
#' @param n_obs Sample size per loop
#' @param family Copula family name (string)
#' @param margins Margin family name (string)
#' @param true_tau The target Kendall's tau
#'
#' @return A 1-row data.frame with 6 summary metrics
#'         (bias, var, mse for both known/unknown)
#'
run_scenario <- function(N_sim, n_obs, family, margins, true_tau) {
  
  # --- 4.1. Define Copula ---
  # *** CORRECTED SECTION ***
  # iTau() calls are now INSIDE the relevant if-block
  
  if (family == "normal") {
    true_param <- iTau(ellipCopula("normal"), true_tau)
    copula_obj <- normalCopula(dim = 2, param = true_param, dispstr = "un")
  } else if (family == "t") {
    true_param <- iTau(ellipCopula("t"), true_tau)
    copula_obj <- tCopula(dim = 2, param = true_param, df = T_COPULA_DF, dispstr = "un")
  } else if (family == "gumbel") {
    true_param <- iTau(gumbelCopula(2), true_tau) # param > 1
    copula_obj <- gumbelCopula(dim = 2, param = true_param)
  } else if (family == "clayton") {
    true_param <- iTau(claytonCopula(2), true_tau) # param > 0
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
  
  # Storage for this scenario
  est_param_known <- numeric(N_sim)
  est_param_unknown <- numeric(N_sim)

  # --- 4.3. Run N_sim loops ---
  for (i in 1:N_sim) {
    # Generate data from the true MVDC
    X_data <- rMvdc(n_obs, true_mvdc)
    
    # Case 1: Known Margins
    # We apply the true parametric marginal CDF
    U_known <- p_margin_fn(X_data)
    
    # Case 2: Unknown Margins
    # We apply the empirical marginal CDF (pseudo-observations)
    U_unknown <- pobs(X_data)
    
    # --- 4.4. Fit Copulas ---
    # Need to create the "unfitted" copula object for fitCopula
    if (family == "normal") { 
      fit_cop_obj <- normalCopula(dim = 2, dispstr = "un")
    } else if (family == "t") {
      # For t-copula, ML is complex. We fix df and estimate rho.
      fit_cop_obj <- tCopula(dim = 2, dispstr = "un", df = T_COPULA_DF)
    } else if (family == "gumbel") {
      fit_cop_obj <- gumbelCopula(dim = 2, param = 2) # Start value
    } else if (family == "clayton") {
      fit_cop_obj <- claytonCopula(dim = 2, param = 2) # Start value
    }
    
    # Fit Known
    fit_known <- try(fitCopula(fit_cop_obj, data = U_known, method = "ml"), silent = TRUE)
    if (!inherits(fit_known, "try-error")) {
      est_param_known[i] <- coef(fit_known)[1]
    } else {
      est_param_known[i] <- NA # Mark as failed
    }
    
    # Fit Unknown
    fit_unknown <- try(fitCopula(fit_cop_obj, data = U_unknown, method = "ml"), silent = TRUE)
    if (!inherits(fit_unknown, "try-error")) {
      est_param_unknown[i] <- coef(fit_unknown)[1]
    } else {
      est_param_unknown[i] <- NA # Mark as failed
    }
    
  } # End N_sim loop
  
  # --- 4.5. Calculate and return metrics ---
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


# --- 5. Run the Full Experiment ---

cat("Starting OFAT 8-scenario simulation (Corrected Version)... This will be quick.\n")
cat("Scenarios to run: \n")
print(param_grid)

# Initialize list to store results
results_list <- vector("list", nrow(param_grid))

# Run the main loop
pb <- txtProgressBar(min = 0, max = nrow(param_grid), style = 3)

for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  
  metrics <- run_scenario(
    N_sim = N_sim,
    n_obs = params$n_obs,
    family = as.character(params$family), # Use as.character to pass string
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

# --- 6. Process and Plot Results ---

# 6.1. Create 'MSE_Inflation' metric
results_df <- results_df %>%
  mutate(
    MSE_Inflation = mse_unknown / mse_known,
    Var_Inflation = var_unknown / var_known,
    # Ensure factors are set for correct plot ordering
    family = factor(family, levels = c("normal", "t", "gumbel", "clayton")), 
    margins = factor(margins, levels = c("norm", "t", "gamma"))
  )

# 6.2. Pivot for plotting
results_long_mse <- results_df %>%
  select(n_obs, family, margins, mse_known, mse_unknown) %>%
  pivot_longer(
    cols = c("mse_known", "mse_unknown"),
    names_to = "Case",
    values_to = "MSE",
    names_prefix = "mse_"
  ) %>%
  mutate(Case = factor(Case, levels = c("known", "unknown")))


# --- Plot 1: Effect of Sample Size (n) ---
plot_data_n <- results_long_mse %>%
  filter(family == "normal", margins == "norm") 

gg_effect_n <- ggplot(plot_data_n, aes(x = Case, y = MSE, fill = Case)) +
  geom_col(position = "dodge") +
  facet_wrap(~ sprintf("n = %d", n_obs), scales = "free_y") +
  labs(
    title = "Effect of Sample Size on Estimator MSE",
    subtitle = "Fixed: Normal Copula, Normal Margins. True Tau = 1/3.", 
    y = "Mean Squared Error (MSE)"
  ) +
  theme_bw()

print(gg_effect_n)

# --- Plot 1b: MSE Inflation vs Sample Size (n) ---
plot_data_n_inf <- results_df %>%
  filter(family == "normal", margins == "norm") 

gg_effect_n_inf <- ggplot(plot_data_n_inf, aes(x = n_obs, y = MSE_Inflation)) +
  geom_line() +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = c(50, 200, 1000)) + 
  labs(
    title = "Cost of Unknown Margins vs. Sample Size",
    subtitle = "Fixed: Normal Copula, Normal Margins. True Tau = 1/3.", 
    y = "MSE Inflation (MSE_Unknown / MSE_Known)"
  ) +
  theme_bw()

print(gg_effect_n_inf)


# --- Plot 2: Effect of Copula Family ---
plot_data_fam <- results_long_mse %>%
  filter(n_obs == 200, margins == "norm")

gg_effect_fam <- ggplot(plot_data_fam, aes(x = Case, y = MSE, fill = Case)) +
  geom_col(position = "dodge") +
  facet_wrap(~ family, scales = "free_y") +
  labs(
    title = "Effect of Copula Family on Estimator MSE",
    subtitle = "Fixed: n=200, Normal Margins. True Tau = 1/3.",
    y = "Mean Squared Error (MSE)"
  ) +
  theme_bw()

print(gg_effect_fam)

# --- Plot 2b: MSE Inflation vs Copula Family ---
plot_data_fam_inf <- results_df %>%
  filter(n_obs == 200, margins == "norm")

gg_effect_fam_inf <- ggplot(plot_data_fam_inf, aes(x = family, y = MSE_Inflation, fill = family)) +
  geom_col() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = "Cost of Unknown Margins vs. Copula Family",
    subtitle = "Fixed: n=200, Normal Margins. True Tau = 1/3.",
    y = "MSE Inflation (MSE_Unknown / MSE_Known)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(gg_effect_fam_inf)


# --- Plot 3: Effect of Marginal Distribution ---
plot_data_mar <- results_long_mse %>%
  filter(n_obs == 200, family == "normal") 

gg_effect_mar <- ggplot(plot_data_mar, aes(x = Case, y = MSE, fill = Case)) +
  geom_col(position = "dodge") +
  facet_wrap(~ margins, scales = "free_y") +
  labs(
    title = "Effect of Margins on Estimator MSE",
    subtitle = "Fixed: n=200, Normal Copula. True Tau = 1/3.", 
    y = "Mean Squared Error (MSE)"
  ) +
  theme_bw()

print(gg_effect_mar)

# --- Plot 3b: MSE Inflation vs Marginal Distribution ---
plot_data_mar_inf <- results_df %>%
  filter(n_obs == 200, family == "normal") 

gg_effect_mar_inf <- ggplot(plot_data_mar_inf, aes(x = margins, y = MSE_Inflation, fill = margins)) +
  geom_col() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    title = "Cost of Unknown Margins vs. Marginal Distribution",
    subtitle = "Fixed: n=200, Normal Copula. True Tau = 1/3.", 
    y = "MSE Inflation (MSE_Unknown / MSE_Known)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(gg_effect_mar_inf)