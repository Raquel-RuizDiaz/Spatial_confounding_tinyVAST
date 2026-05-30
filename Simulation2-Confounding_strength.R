# Figure 4 -Confounding strength

library(tinyVAST)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggpubr)
library(viridis)

# Helper functions
rmse <- function(vec) sqrt(mean(vec^2))

# Configuration
config_fig4a <- list(
  seed = 101,
  n_simulations = 50,
  n_s = 200,                           # Fixed sample size
  n_false_covariates = 3,              # Fixed false covariates (same as original)
  confounding_strength = c(0, 0.5, 1, 2, 4),  # NEW axis: varying confounding
  alpha = 0.05,
  decorrelation_dist_z = 2.0,          # Z: smooth large-scale latent field
  decorrelation_dist_x = 0.5,          # X: local covariate
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
    decorrelation_dist_z = 2.0,
    decorrelation_dist_x = 0.5,
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
  kappa_z <- sqrt(8) / decorrelation_dist_z
  tau_z   <- 1 / kappa_z / spatial_sd / sqrt(4 * pi)
  kappa_x <- sqrt(8) / decorrelation_dist_x
  tau_x   <- 1 / kappa_x / spatial_sd / sqrt(4 * pi)
  
  cat("Running simulation: n_s =", n_s, ", Xreal_to_Z =", Xreal_to_Z, "\n")
  
  for (sim in 1:n_simulations) {
    
    xy_i <- pracma::poisson2disk(n = n_s)
    mesh <- fmesher::fm_mesh_2d(xy_i)
    spde <- fmesher::fm_fem(mesh)
    A_is <- Matrix::drop0(fmesher::fm_evaluator(mesh, loc = xy_i)$proj$A)
    
    # Separate precision matrices
    Q_z <- (kappa_z^4 * spde$c0 + 2 * kappa_z^2 * spde$g1 + spde$g2) * tau_z^2
    Q_x <- (kappa_x^4 * spde$c0 + 2 * kappa_x^2 * spde$g1 + spde$g2) * tau_x^2
    
    # Z: smooth large-scale latent field (independent)
    z_s  <- as.matrix(A_is %*% rmvnorm_prec(Q = Q_z, n = 1))
    # X: local covariate
    x_sj <- as.matrix(A_is %*% rmvnorm_prec(Q = Q_x, n = n_covariates))
    
    # Induce correlation between X_true and Z (mediator-like)
    # X_true causes Z, making them correlated
    z_s <- z_s + Xreal_to_Z * x_sj[, 1]
    
    # Realized correlation after modification
    realized_corr <- cor(z_s, x_sj[, 1])[1, 1]
    
    # Response
    linear_predictor <- z_s + x_sj %*% covariate_effects
    
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
    dat <- cbind(dat, setNames(as.data.frame(x_sj), covariate_names))
    
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
              Estimate           = coef_summary[coef_name, "Estimate"],
              SE                 = coef_summary[coef_name, "Std. Error"],
              p_value            = coef_summary[coef_name, 4],
              Significant        = coef_summary[coef_name, 4] < alpha,
              RMSE_latent        = rmse(dat$dtrue - predictions),
              RMSE_with_gmrf     = NA,
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
              Estimate           = estimates[i],
              SE                 = SE[i],
              p_value            = p_values[i],
              Significant        = p_values[i] < alpha,
              RMSE_latent        = rmse(dat$dtrue - predictions),
              RMSE_with_gmrf     = rmse(dat$dtrue - pred_with_gmrf),
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
    decorrelation_dist_z = config_fig4a$decorrelation_dist_z,
    decorrelation_dist_x = config_fig4a$decorrelation_dist_x,
    n_simulations        = config_fig4a$n_simulations,
    measurement_error    = config_fig4a$measurement_error,
    spatial_sd           = config_fig4a$spatial_sd,
    measurement_sd       = config_fig4a$measurement_sd,
    alpha                = config_fig4a$alpha,
    include_variable_selection = TRUE
  )
  results_fig4a <- rbind(results_fig4a, results)
}

write.csv(results_fig4a, "results/fig4a_confounding_strength_results.csv", row.names = FALSE)

# Summarize
tmp = subset(results_fig4a, Covariate == "real" )
tapply( tmp$RMSE_latent, INDEX = list(tmp$Model,tmp$Sample_Size, tmp$Confounding_Label), FUN = mean )

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

###################
# Plotting
###################

model_colors <- c(
  "Covariate_Only" = "#440154FF",
  "VAST_GMRF"      = "#31688EFF",
  "VAST_RSR"       = "#35B779FF",
  "VAST_Selected"  = "gold"
)

model_labels <- c(
  "Covariate_Only" = "GLM",
  "VAST_GMRF"      = "SGLMM",
  "VAST_RSR"       = "RSR-SGLMM",
  "VAST_Selected"  = "SGLMM-selected"
)

# Add confounding labels
corr_lookup_4a <- results_fig4a %>%
  filter(Covariate == "real", Model == "Covariate_Only") %>%
  group_by(Xreal_to_Z) %>%
  summarise(mean_corr = mean(Realized_Corr, na.rm = TRUE), .groups = "drop") %>%
  arrange(Xreal_to_Z) %>%
  mutate(Confounding_Label = factor(
    sprintf("r = %.2f", mean_corr),
    levels = sprintf("r = %.2f", mean_corr)
  ))

results_fig4a <- results_fig4a %>%
  select(-any_of("Confounding_Label")) %>%
  left_join(corr_lookup_4a %>% select(Xreal_to_Z, Confounding_Label),
            by = "Xreal_to_Z")

shared_theme <- theme_bw() +
  theme(
    plot.title      = element_text(size = 12, face = "bold", hjust = 0.5),
    strip.text      = element_text(size = 11, face = "bold"),
    axis.title.x    = element_text(size = 11, margin = margin(t = 10)),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 9),
    axis.title.y    = element_text(size = 11, margin = margin(r = 8)),
    axis.text.y     = element_text(size = 9),
    legend.position = "none"
  )

# Panel A: Unconditional RMSE across confounding strength
p_uncond <- results_fig4a %>%
  filter(Covariate == "real") %>%
  ggplot(aes(x = Model, y = RMSE_latent, fill = Model)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.3) +
  facet_wrap(~ Confounding_Label, nrow = 1) +
  scale_fill_manual(values = model_colors, labels = model_labels) +
  scale_x_discrete(labels = model_labels) +
  labs(
    #title = "A) Unconditional prediction RMSE across confounding strength",
    x = "Model", y = "RMSE (unconditional prediction)"
  ) +
  shared_theme

# Panel B: Conditional RMSE across confounding strength
p_cond <- results_fig4a %>%
  filter(Covariate == "real", Model != "Covariate_Only") %>%
  ggplot(aes(x = Model, y = RMSE_with_gmrf, fill = Model)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.3) +
  facet_wrap(~ Confounding_Label, nrow = 1) +
  scale_fill_manual(values = model_colors, labels = model_labels) +
  scale_x_discrete(labels = model_labels) +
  labs(
    title = "B) Conditional prediction RMSE across confounding strength",
    x = "Model", y = "RMSE (conditional prediction)"
  ) +
  shared_theme

# Panel C: Coefficient estimates across confounding strength
p_coef <- results_fig4a %>%
  filter(Covariate == "real") %>%
  ggplot(aes(x = Model, y = Estimate, fill = Model)) +
  geom_violin(alpha = 0.8) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
  facet_wrap(~ Confounding_Label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = model_colors, labels = model_labels) +
  scale_x_discrete(labels = model_labels) +
  labs(
    title = "C) Coefficient estimate across confounding strength",
    x = "Model", y = "Estimate"
  ) +
  shared_theme

# Shared legend
legend_plot <- results_fig4a %>%
  filter(Covariate == "real") %>%
  ggplot(aes(x = Model, y = RMSE_latent, fill = Model)) +
  geom_boxplot() +
  scale_fill_manual(values = model_colors, labels = model_labels) +
  labs(fill = "Model") +
  theme_bw() +
  theme(
    legend.title     = element_text(size = 12, face = "bold"),
    legend.text      = element_text(size = 11),
    legend.position  = "bottom",
    legend.direction = "horizontal"
  )
shared_legend <- ggpubr::get_legend(legend_plot)

fig4a_final <- plot_grid(
  p_uncond, p_cond, shared_legend,
  ncol        = 1,
  rel_heights = c(0.3, 0.3, 0.07)
)

fig4a_final <- plot_grid(
  p_uncond, shared_legend,
  ncol        = 1,
  rel_heights = c(1, 0.08)
)

print(fig4a_final)

ggsave(
  filename = "figures/Figure4a_confounding_strength_A.png",
  plot     = fig4a_final,
  bg       = "white",
  width    = 10, height = 4, units = "in", dpi = 300
)