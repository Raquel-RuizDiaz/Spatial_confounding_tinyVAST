# Simulation 3 (spatial latent variable as mediator).
# Figure 3 is made from these results by 03b_Simulation3_mediator_figure.R.

library(tinyVAST)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggpubr)
library(viridis)

# Helper functions
rmse <- function(vec) sqrt(mean(vec^2))
# REVISION 2: tinyVAST >= 1.6 keeps rmvnorm_prec() internal with argument `prec`
rmvnorm_prec <- function(Q, n = 1) tinyVAST:::rmvnorm_prec(prec = Q, n = n)

# Configuration
config_fig4a <- list(
  seed = 101,
  n_simulations = 50,
  n_s = 200,                           # Fixed sample size
  n_false_covariates = 3,              # Fixed false covariates
  confounding_strength = c(0, 0.15, 0.25, 0.35, 0.5, 0.75, 1, 2, 4),
  confounding_strength_reported = c(0, 0.15, 0.25, 0.5, 0.75, 1, 4),
  alpha = 0.05,
  decorrelation_dist_omega = 2.0,          # omega: smooth large-scale latent field
  decorrelation_dist_z = 0.5,          # z: local covariate
  measurement_error = "normal",
  spatial_sd = 0.2,
  measurement_sd = 0.2
)

###################
# Simulation function
# DGP: mediator-like setup where correlation between X and Z is induced
###################

run_confounding_simulation <- function(
    n_s = 200,
    n_false_covariates = 3,
    decorrelation_dist_omega = 2.0,
    decorrelation_dist_z = 0.5,
    Xreal_to_Z = 0,
    measurement_error = "normal",
    n_simulations = 50,
    spatial_sd = 0.2,
    measurement_sd = 0.2,
    alpha = 0.05,
    include_variable_selection = TRUE
) {
  
  n_covariates      <- 1 + n_false_covariates
  covariate_effects <- c(1, rep(0, n_false_covariates))
  results           <- data.frame()
  
  # Separate spatial parameters for Z (smooth) and X (local)
  kappa_omega <- sqrt(8) / decorrelation_dist_omega
  tau_omega   <- 1 / kappa_omega / spatial_sd / sqrt(4 * pi)
  kappa_z <- sqrt(8) / decorrelation_dist_z
  tau_z   <- 1 / kappa_z / spatial_sd / sqrt(4 * pi)
  
  cat("Running simulation: n_s =", n_s, ", Xreal_to_Z =", Xreal_to_Z, "\n")
  
  for (sim in 1:n_simulations) {
    
    xy_i <- pracma::poisson2disk(n = n_s)
    mesh <- fmesher::fm_mesh_2d(xy_i)
    spde <- fmesher::fm_fem(mesh)
    A_is <- Matrix::drop0(fmesher::fm_evaluator(mesh, loc = xy_i)$proj$A)
    
    # Separate precision matrices
    Q_omega <- (kappa_omega^4 * spde$c0 + 2 * kappa_omega^2 * spde$g1 + spde$g2) * tau_omega^2
    Q_z <- (kappa_z^4 * spde$c0 + 2 * kappa_z^2 * spde$g1 + spde$g2) * tau_z^2
    
    # omega: smooth large-scale latent field (independent)
    omega_s  <- as.matrix(A_is %*% rmvnorm_prec(Q = Q_omega, n = 1))
    # X: local covariate
    zW_sj <- as.matrix(A_is %*% rmvnorm_prec(Q = Q_z, n = n_covariates))
    
    # Induce correlation between X_true and Z (mediator-like)
    # X_true causes Z, making them correlated
    omega0_s <- omega_s   # REVISION 2: keep the independent part of omega for variance partitioning
    omega_s <- omega_s + Xreal_to_Z * zW_sj[, 1]
    
    # Realized correlation after modification
    realized_corr <- cor(omega_s, zW_sj[, 1])[1, 1]
    
    # Response
    linear_predictor <- omega_s + zW_sj %*% covariate_effects
    
    if (measurement_error == "normal") {
      d_s     <- rnorm(n = length(linear_predictor), mean = linear_predictor, sd = measurement_sd)
      Family  <- gaussian()
      invlink <- identity
    } else {
      d_s     <- rpois(n = length(linear_predictor), lambda = exp(linear_predictor))
      Family  <- poisson()
      invlink <- exp
    }
    
    covariate_names <- c("real", paste0("false_", 1:n_false_covariates))
    dat <- data.frame(X1 = xy_i[,1], X2 = xy_i[,2], d = d_s, dtrue = linear_predictor)
    dat <- cbind(dat, setNames(as.data.frame(zW_sj), covariate_names))
    
    models_to_fit <- c("Covariate_Only", "VAST_GMRF", "VAST_RSR")
    if (include_variable_selection) models_to_fit <- c(models_to_fit, "VAST_Selected")
    
    for (model_name in models_to_fit) {
      tryCatch({
        
        if (model_name == "Covariate_Only") {
          formula_str <- paste("d ~", paste(covariate_names, collapse = " + "))
          fit <- lm(as.formula(formula_str), data = dat)
          coef_summary <- summary(fit)$coefficients
          predictions  <- predict(fit, type = "response")
          
          for (i in 1:n_covariates) {
            coef_name <- covariate_names[i]
            results <- rbind(results, data.frame(
              Model              = model_name,
              Simulation         = sim,
              Sample_Size        = n_s,
              N_False_Covariates = n_false_covariates,
              Xreal_to_Z         = Xreal_to_Z,
              Realized_Corr      = realized_corr,
              Covariate          = coef_name,
              True_Effect        = covariate_effects[i],
              Total_Effect       = covariate_effects[i] + (i == 1) * Xreal_to_Z,  # REVISION 2
              Estimate           = coef_summary[coef_name, "Estimate"],
              SE                 = coef_summary[coef_name, "Std. Error"],
              p_value            = coef_summary[coef_name, 4],
              Significant        = coef_summary[coef_name, 4] < alpha,
              RMSE_latent        = rmse(dat$dtrue - predictions),
              RMSE_with_gmrf     = NA,
              Var_omega0         = var(omega0_s[, 1]),   # REVISION 2
              Var_gammaX         = var(Xreal_to_Z * zW_sj[, 1]),
              Var_field_hat      = NA,
              Field_slope_on_X   = NA,
              Cor_field_omega0   = NA,
              stringsAsFactors   = FALSE
            ))
          }
          
        } else {
          get_rsr <- (model_name == "VAST_RSR")
          
          if (model_name == "VAST_Selected") {
            full_fit <- tinyVAST(
              formula         = as.formula(paste("d ~", paste(covariate_names, collapse = " + "))),
              data            = cbind(dat, "var" = "n"),
              family          = Family, space_term = "",
              variable_column = "var", space_columns = c("X1", "X2"),
              spatial_domain  = mesh,
              control         = tinyVASTcontrol(reml = TRUE, get_rsr = FALSE)
            )
            full_estimates <- as.list(full_fit$sdrep, what = "Estimate")$alpha_j
            full_SE        <- as.list(full_fit$sdrep, what = "Std. Error")$alpha_j
            full_p_values  <- 1 - pchisq((full_estimates / full_SE)^2, df = 1)
            selected_covs  <- covariate_names[full_p_values[-1] < alpha]
            
            formula_selected <- if (length(selected_covs) > 0) {
              as.formula(paste("d ~", paste(selected_covs, collapse = " + ")))
            } else {
              as.formula("d ~ 1")
            }
            
            fit <- tinyVAST(
              formula         = formula_selected,
              data            = cbind(dat, "var" = "n"),
              family          = Family, space_term = "",
              variable_column = "var", space_columns = c("X1", "X2"),
              spatial_domain  = mesh,
              control         = tinyVASTcontrol(reml = TRUE, get_rsr = TRUE)
            )
            estimates <- as.list(fit$sdrep, what = "Estimate", report = TRUE)$alphaprime_j
            SE        <- as.list(fit$sdrep, what = "Std. Error", report = TRUE)$alphaprime_j
            
          } else {
            fit <- tinyVAST(
              formula         = as.formula(paste("d ~", paste(covariate_names, collapse = " + "))),
              data            = cbind(dat, "var" = "n"),
              family          = Family, space_term = "",
              variable_column = "var", space_columns = c("X1", "X2"),
              spatial_domain  = mesh,
              control         = tinyVASTcontrol(reml = TRUE, get_rsr = get_rsr)
            )
            if (get_rsr) {
              estimates <- as.list(fit$sdrep, what = "Estimate", report = TRUE)$alphaprime_j
              SE        <- as.list(fit$sdrep, what = "Std. Error", report = TRUE)$alphaprime_j
            } else {
              estimates <- as.list(fit$sdrep, what = "Estimate")$alpha_j
              SE        <- as.list(fit$sdrep, what = "Std. Error")$alpha_j
            }
          }
          
          p_values       <- 1 - pchisq((estimates / SE)^2, df = 1)
          predictions    <- invlink(fit$tmb_inputs$tmb_data$X_ij %*% estimates)
          pred_with_gmrf <- predict(fit)
          # REVISION 2: estimated spatial-field contribution (conditional minus fixed effects),
          # used to check how variance is partitioned between beta and the field
          field_hat      <- as.numeric(pred_with_gmrf - fit$tmb_inputs$tmb_data$X_ij %*% estimates)

          for (i in 2:length(estimates)) {
            coef_idx <- i - 1
            results <- rbind(results, data.frame(
              Model              = model_name,
              Simulation         = sim,
              Sample_Size        = n_s,
              N_False_Covariates = n_false_covariates,
              Xreal_to_Z         = Xreal_to_Z,
              Realized_Corr      = realized_corr,
              Covariate          = covariate_names[coef_idx],
              True_Effect        = covariate_effects[coef_idx],
              Total_Effect       = covariate_effects[coef_idx] + (coef_idx == 1) * Xreal_to_Z,  # REVISION 2
              Estimate           = estimates[i],
              SE                 = SE[i],
              p_value            = p_values[i],
              Significant        = p_values[i] < alpha,
              RMSE_latent        = rmse(dat$dtrue - predictions),
              RMSE_with_gmrf     = rmse(dat$dtrue - pred_with_gmrf),
              Var_omega0         = var(omega0_s[, 1]),   # REVISION 2
              Var_gammaX         = var(Xreal_to_Z * zW_sj[, 1]),
              Var_field_hat      = var(field_hat),
              Field_slope_on_X   = coef(lm(field_hat ~ zW_sj[, 1]))[2],
              Cor_field_omega0   = cor(field_hat, omega0_s[, 1]),
              stringsAsFactors   = FALSE
            ))
          }
        }
        
      }, error = function(e) {
        cat("Model", model_name, "failed for sim", sim, ":", e$message, "\n")
      })
    }
    
    if (sim %% 10 == 0) cat("Completed", sim, "of", n_simulations, "\n")
  }
  
  return(results)
}

###################
# Run experiment
###################

set.seed(config_fig4a$seed)

results_fig4a <- data.frame()

for (cs in config_fig4a$confounding_strength) {
  cat("\n--- Confounding strength:", cs, "---\n")
  results <- run_confounding_simulation(
    n_s                  = config_fig4a$n_s,
    n_false_covariates   = config_fig4a$n_false_covariates,
    Xreal_to_Z           = cs,
    decorrelation_dist_omega = config_fig4a$decorrelation_dist_omega,
    decorrelation_dist_z = config_fig4a$decorrelation_dist_z,
    n_simulations        = config_fig4a$n_simulations,
    measurement_error    = config_fig4a$measurement_error,
    spatial_sd           = config_fig4a$spatial_sd,
    measurement_sd       = config_fig4a$measurement_sd,
    alpha                = config_fig4a$alpha,
    include_variable_selection = TRUE
  )
  results_fig4a <- rbind(results_fig4a, results)
}

write.csv(results_fig4a, "results/Simulation3_mediator_results.csv", row.names = FALSE)

# summaries and figure use the reported subset of confounding levels
results_fig4a <- results_fig4a %>%
  filter(Xreal_to_Z %in% config_fig4a$confounding_strength_reported)

###################
# bias and RMSE of beta against both estimands
# Direct effect = 1 (beta in the DGP); total effect = 1 + gamma (beta + Xreal_to_Z),
# because omega = omega0 + gamma * X, so a unit change in X shifts E[y] by 1 + gamma.
###################

estimand_summary <- results_fig4a %>%
  filter(Covariate == "real") %>%
  group_by(Model, Xreal_to_Z) %>%
  summarise(
    n            = n(),
    mean_corr    = mean(Realized_Corr),
    mean_beta    = mean(Estimate),
    sd_beta      = sd(Estimate),
    bias_direct  = mean(Estimate - True_Effect),
    rmse_direct  = sqrt(mean((Estimate - True_Effect)^2)),
    bias_total   = mean(Estimate - Total_Effect),
    rmse_total   = sqrt(mean((Estimate - Total_Effect)^2)),
    .groups      = "drop"
  )

print(as.data.frame(estimand_summary), digits = 3)
write.csv(estimand_summary, "results/Simulation3_mediator_estimand.csv", row.names = FALSE)

###################
# variance partitioning diagnostics
# Field_slope_on_X ~ 0  -> the spatial field does not absorb gamma * X (beta carries it)
# Field_slope_on_X ~ gamma -> the field absorbs gamma * X (beta targets the direct effect)
# Cor_field_omega0 -> how well the field recovers the independent part of omega
###################

partition_summary <- results_fig4a %>%
  filter(Covariate == "real", Model != "Covariate_Only") %>%
  group_by(Model, Xreal_to_Z) %>%
  summarise(
    mean_field_slope_on_X = mean(Field_slope_on_X),
    mean_var_field_hat    = mean(Var_field_hat),
    mean_var_omega0       = mean(Var_omega0),
    mean_var_gammaX       = mean(Var_gammaX),
    mean_cor_field_omega0 = mean(Cor_field_omega0),
    .groups               = "drop"
  )

print(as.data.frame(partition_summary), digits = 3)
write.csv(partition_summary, "results/Simulation3_mediator_variance_partition.csv", row.names = FALSE)

# Indirect effect per replicate, delta = beta_RSR - beta_SGLMM (eq. 4-5), and the
# unconditional RMSE gap (SGLMM - GLM) it should explain if partitioning drives the gap
delta_summary <- results_fig4a %>%
  filter(Covariate == "real", Model %in% c("Covariate_Only", "VAST_GMRF", "VAST_RSR")) %>%
  select(Model, Simulation, Xreal_to_Z, Estimate, RMSE_latent) %>%
  tidyr::pivot_wider(names_from = Model, values_from = c(Estimate, RMSE_latent)) %>%
  mutate(
    delta    = Estimate_VAST_RSR - Estimate_VAST_GMRF,
    rmse_gap = RMSE_latent_VAST_GMRF - RMSE_latent_Covariate_Only
  ) %>%
  group_by(Xreal_to_Z) %>%
  summarise(
    mean_delta         = mean(delta),
    sd_delta           = sd(delta),
    mean_rmse_gap      = mean(rmse_gap),
    cor_gap_abs_delta  = cor(rmse_gap, abs(delta)),
    .groups            = "drop"
  )

print(as.data.frame(delta_summary), digits = 3)
write.csv(delta_summary, "results/Simulation3_mediator_indirect_effect.csv", row.names = FALSE)


test <- results_fig4a %>%
  filter(Covariate == "real", Model == "Covariate_Only")

###################
# Realized correlation summary
###################

corr_summary_4a <- results_fig4a %>%
  filter(Covariate == "real", Model == "Covariate_Only") %>%
  group_by(Xreal_to_Z) %>%
  summarise(
    mean_corr = mean(Realized_Corr, na.rm = TRUE),
    sd_corr   = sd(Realized_Corr,   na.rm = TRUE),
    .groups   = "drop"
  )

print(corr_summary_4a)

# The manuscript figure and tables for Simulation 3 are made by 03b_Simulation3_mediator_figure.R and 03c_Simulation3_mediator_tables.R
