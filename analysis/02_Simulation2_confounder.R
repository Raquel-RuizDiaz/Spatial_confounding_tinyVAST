# Simulation 2 (spatial latent variable as confounder)
# The covariate shares the smooth latent field while keeping its own short-range variation:
#   z <- u + lambda * omega,  mu = z + omega

# Reuses run_confounding_simulation() from 03_Simulation3_mediator.R without changing that file.

source_lines <- readLines("analysis/03_Simulation3_mediator.R")
run_start <- grep("^# Run experiment", source_lines) - 2
source_text <- paste(source_lines[1:run_start], collapse = "\n")

replace_code <- function(text, old, new) {
  stopifnot(grepl(old, text, fixed = TRUE))
  gsub(old, new, text, fixed = TRUE)
}

# Covariate shares omega; the induced step is removed
source_text <- replace_code(
  source_text,
  "omega_s <- omega_s + Xreal_to_Z * zW_sj[, 1]",
  "zW_sj[, 1] <- zW_sj[, 1] + Xreal_to_Z * omega_s[, 1]"
)
# Total effect is the in-sample slope of the true mean on z
source_text <- replace_code(
  source_text,
  "covariate_effects[i] + (i == 1) * Xreal_to_Z",
  "covariate_effects[i] + (i == 1) * cov(omega_s[, 1], zW_sj[, 1]) / var(zW_sj[, 1])"
)
source_text <- replace_code(
  source_text,
  "covariate_effects[coef_idx] + (coef_idx == 1) * Xreal_to_Z",
  "covariate_effects[coef_idx] + (coef_idx == 1) * cov(omega_s[, 1], zW_sj[, 1]) / var(zW_sj[, 1])"
)
# SGLMM-selected reports the SGLMM slope, not the RSR slope
source_text <- replace_code(
  source_text,
  "\n            estimates <- as.list(fit$sdrep, what = \"Estimate\", report = TRUE)$alphaprime_j\n            SE        <- as.list(fit$sdrep, what = \"Std. Error\", report = TRUE)$alphaprime_j\n",
  "\n            estimates <- as.list(fit$sdrep, what = \"Estimate\")$alpha_j\n            SE        <- as.list(fit$sdrep, what = \"Std. Error\")$alpha_j\n"
)
# Label SGLMM-selected coefficients by name, not by position
source_text <- replace_code(
  source_text,
  "field_hat      <- as.numeric(pred_with_gmrf - fit$tmb_inputs$tmb_data$X_ij %*% estimates)\n",
  paste0(
    "field_hat      <- as.numeric(pred_with_gmrf - fit$tmb_inputs$tmb_data$X_ij %*% estimates)\n",
    "          if (model_name == \"VAST_Selected\") {\n",
    "            coef_pos  <- match(selected_covs, covariate_names) + 1\n",
    "            full_est  <- c(estimates[1], rep(0, n_covariates))\n",
    "            full_se   <- c(SE[1], rep(NA, n_covariates))\n",
    "            full_est[coef_pos] <- estimates[-1]\n",
    "            full_se[coef_pos]  <- SE[-1]\n",
    "            estimates <- full_est\n",
    "            SE        <- full_se\n",
    "            p_values  <- 1 - pchisq((estimates / SE)^2, df = 1)\n",
    "          }\n"
  )
)
eval(parse(text = source_text))
library(tidyr)

config_shared <- list(
  seed = 101,
  n_simulations = 50,
  n_s = 200,
  n_false_covariates = 3,
  shared_strength = c(0, 0.35, 0.5, 1, 1.5, 2, 4),   # lambda
  decorrelation_dist_omega = 2.0,
  decorrelation_dist_z = 0.5,
  spatial_sd = 0.2,
  measurement_sd = 0.2
)

set.seed(config_shared$seed)

results_shared <- data.frame()

for (cs in config_shared$shared_strength) {
  cat("\n--- Shared strength:", cs, "---\n")
  results <- run_confounding_simulation(
    n_s                  = config_shared$n_s,
    n_false_covariates   = config_shared$n_false_covariates,
    Xreal_to_Z           = cs,
    decorrelation_dist_omega = config_shared$decorrelation_dist_omega,
    decorrelation_dist_z = config_shared$decorrelation_dist_z,
    n_simulations        = config_shared$n_simulations,
    spatial_sd           = config_shared$spatial_sd,
    measurement_sd       = config_shared$measurement_sd,
    include_variable_selection = TRUE
  )
  results_shared <- rbind(results_shared, results)
}

write.csv(results_shared, "results/Simulation2_confounder_results.csv", row.names = FALSE)

###################
# Summaries
###################

true_cov <- results_shared %>% filter(Covariate == "real")

estimand_summary <- true_cov %>%
  group_by(Model, Xreal_to_Z) %>%
  summarise(
    n           = n(),
    mean_corr   = mean(Realized_Corr),
    mean_total  = mean(Total_Effect),
    mean_beta   = mean(Estimate),
    sd_beta     = sd(Estimate),
    bias_direct = mean(Estimate - True_Effect),
    rmse_direct = rmse(Estimate - True_Effect),
    bias_total  = mean(Estimate - Total_Effect),
    rmse_total  = rmse(Estimate - Total_Effect),
    rmse_unconditional = mean(RMSE_latent),
    retained    = mean(!is.na(SE)),
    coverage_direct = mean(abs(Estimate - True_Effect) < 1.96 * SE, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(estimand_summary, "results/Simulation2_confounder_estimand.csv", row.names = FALSE)

field_summary <- true_cov %>%
  filter(Model %in% c("VAST_GMRF", "VAST_RSR")) %>%
  group_by(Model, Xreal_to_Z) %>%
  summarise(
    mean_field_slope_on_X = mean(Field_slope_on_X),
    mean_var_field_hat    = mean(Var_field_hat),
    mean_var_omega        = mean(Var_omega0),
    mean_cor_field_omega  = mean(Cor_field_omega0, na.rm = TRUE),
    n_flat_field          = sum(is.na(Cor_field_omega0)),
    .groups = "drop"
  )
write.csv(field_summary, "results/Simulation2_confounder_field.csv", row.names = FALSE)

# Indirect effect: true covariate, and Euclidean norm across all covariates
delta_all <- results_shared %>%
  filter(Model %in% c("VAST_GMRF", "VAST_RSR")) %>%
  select(Model, Simulation, Xreal_to_Z, Covariate, Estimate) %>%
  pivot_wider(names_from = Model, values_from = Estimate) %>%
  mutate(delta = VAST_RSR - VAST_GMRF)
delta_replicate <- delta_all %>%
  group_by(Simulation, Xreal_to_Z) %>%
  summarise(delta_true = delta[Covariate == "real"], delta_norm = sqrt(sum(delta^2)), .groups = "drop") %>%
  left_join(
    true_cov %>%
      filter(Model %in% c("Covariate_Only", "VAST_GMRF")) %>%
      select(Model, Simulation, Xreal_to_Z, RMSE_latent) %>%
      pivot_wider(names_from = Model, values_from = RMSE_latent) %>%
      mutate(rmse_gap = VAST_GMRF - Covariate_Only),
    by = c("Simulation", "Xreal_to_Z")
  )
delta_summary <- delta_replicate %>%
  group_by(Xreal_to_Z) %>%
  summarise(
    mean_delta    = mean(delta_true),
    sd_delta      = sd(delta_true),
    mean_rmse_gap = mean(rmse_gap),
    cor_gap_norm  = cor(rmse_gap, delta_norm),
    .groups = "drop"
  )
write.csv(delta_summary, "results/Simulation2_confounder_indirect_effect.csv", row.names = FALSE)
cat("Pooled cor(rmse gap, delta norm):", cor(delta_replicate$rmse_gap, delta_replicate$delta_norm), "\n")

# The manuscript figure and tables for Simulation 2 are made by 02b_Simulation2_confounder_figure_tables.R
