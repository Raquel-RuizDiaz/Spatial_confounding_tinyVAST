# Sensitivity analysis - factorial design
# Keeps the original sample size x false covariates factorial
# BUT adds induced confounding between X_true and Z

library(tinyVAST)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggpubr)
library(viridis)

# Helper functions
rmse <- function(vec) sqrt(mean(vec^2))

# Configuration — same as original Figure 4 factorial but with fixed confounding level
config_fig4b <- list(
  seed = 101,
  n_simulations = 50,
  sample_sizes       = c(50, 100, 200), 
  n_false_covariates = c(1, 3, 5),      
  Xreal_to_Z         = 1,                 # Fixed moderate confounding (r ~ 0.5-0.7)
  alpha = 0.05,
  decorrelation_dist_omega = 2.0,
  decorrelation_dist_z = 0.5,
  measurement_error = "normal",
  spatial_sd = 0.2,
  measurement_sd = 0.2
)

###################
# Simulation function — original factorial design with confounding
###################

run_factorial_simulation <- function(
    n_s = 100,
    n_false_covariates = 3,
    decorrelation_dist_omega = 2.0,
    decorrelation_dist_z = 0.5,
    Xreal_to_Z = 1,
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
  
  kappa_omega <- sqrt(8) / decorrelation_dist_omega
  tau_omega   <- 1 / kappa_omega / spatial_sd / sqrt(4 * pi)
  kappa_z <- sqrt(8) / decorrelation_dist_z
  tau_z   <- 1 / kappa_z / spatial_sd / sqrt(4 * pi)
  
  cat("Running factorial: n_s =", n_s, ", n_false =", n_false_covariates, "\n")
  
  for (sim in 1:n_simulations) {
    
    xy_i <- pracma::poisson2disk(n = n_s)
    mesh <- fmesher::fm_mesh_2d(xy_i)
    spde <- fmesher::fm_fem(mesh)
    A_is <- Matrix::drop0(fmesher::fm_evaluator(mesh, loc = xy_i)$proj$A)
    
    Q_omega <- (kappa_omega^4 * spde$c0 + 2 * kappa_omega^2 * spde$g1 + spde$g2) * tau_omega^2
    Q_z <- (kappa_z^4 * spde$c0 + 2 * kappa_z^2 * spde$g1 + spde$g2) * tau_z^2
    
    omega_s  <- as.matrix(A_is %*% rmvnorm_prec(Q = Q_omega, n = 1))
    zW_sj <- as.matrix(A_is %*% rmvnorm_prec(Q = Q_z, n = n_covariates))
    
    # Induce confounding between X_true and Z
    omega_s <- omega_s + Xreal_to_Z * zW_sj[, 1]
    
    realized_corr <- cor(omega_s, zW_sj[, 1])[1, 1]
    
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
          formula_str  <- paste("d ~", paste(covariate_names, collapse = " + "))
          fit          <- lm(as.formula(formula_str), data = dat)
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

set.seed(config_fig4b$seed)
cat("Fixed confounding: Xreal_to_Z =", config_fig4b$Xreal_to_Z, "\n")

results_fig4b <- data.frame()

for (n_s in config_fig4b$sample_sizes) {
  for (n_false in config_fig4b$n_false_covariates) {
    results <- run_factorial_simulation(
      n_s                  = n_s,
      n_false_covariates   = n_false,
      Xreal_to_Z           = config_fig4b$Xreal_to_Z,
      decorrelation_dist_omega = config_fig4b$decorrelation_dist_omega,
      decorrelation_dist_z = config_fig4b$decorrelation_dist_z,
      n_simulations        = config_fig4b$n_simulations,
      measurement_error    = config_fig4b$measurement_error,
      spatial_sd           = config_fig4b$spatial_sd,
      measurement_sd       = config_fig4b$measurement_sd,
      alpha                = config_fig4b$alpha,
      include_variable_selection = TRUE
    )
    results_fig4b <- rbind(results_fig4b, results)
  }
}

#write.csv(results_fig4b, "results/fig4b_factorial_with_confounding_results.csv", row.names = FALSE)

###################
# Plotting — original Figure 4 style
# Rows: sample size, Columns: number of false covariates
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

shared_theme <- theme_bw() +
  theme(
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 13, margin = margin(t = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, size =10),
    axis.title.y = element_text(size = 13, margin = margin(r = 15)),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    panel.grid      = element_blank()
  )

# Mean realized correlation for subtitle
mean_corr <- results_fig4b %>%
  filter(Covariate == "real", Model == "Covariate_Only") %>%
  summarise(r = round(mean(Realized_Corr, na.rm = TRUE), 2)) %>%
  pull(r)

# Panel A: Unconditional RMSE — original Figure 4 layout
fig4b_uncond <- results_fig4b %>%
  filter(Covariate == "real") %>%
  ggplot(aes(x = Model, y = RMSE_latent, fill = Model)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.3) +
  facet_grid(
    rows = vars(Sample_Size),
    cols = vars(N_False_Covariates),
    labeller = labeller(
      Sample_Size        = function(x) paste("n =", x),
      N_False_Covariates = function(x) paste("False covariates:", x)
    )
  ) +
  scale_fill_manual(values = model_colors, labels = model_labels) +
  scale_x_discrete(labels = model_labels) +
  labs(
    title    = "Unconditional prediction",
    #subtitle = "Rows: sample size | Columns: number of false covariates",
    x        = "Model",
    y        = "RMSE"
  ) +
  shared_theme

# Panel B: Conditional RMSE
fig4b_cond <- results_fig4b %>%
  filter(Covariate == "real", Model != "Covariate_Only") %>%
  ggplot(aes(x = Model, y = RMSE_with_gmrf, fill = Model)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.3) +
  facet_grid(
    rows = vars(Sample_Size),
    cols = vars(N_False_Covariates),
    labeller = labeller(
      Sample_Size        = function(x) paste("n =", x),
      N_False_Covariates = function(x) paste("False covariates:", x)
    )
  ) +
  scale_fill_manual(values = model_colors, labels = model_labels) +
  scale_x_discrete(labels = model_labels) +
  labs(
    title    = paste0("Conditional prediction RMSE (induced confounding: r \u2248 ", mean_corr, ")"),
    #subtitle = "Rows: sample size | Columns: number of false covariates",
    x        = "Model",
    y        = "RMSE (conditional prediction)"
  ) +
  shared_theme

# Shared legend
legend_plot <- results_fig4b %>%
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

fig4b_final <- plot_grid(
  fig4b_uncond, shared_legend,
  ncol        = 1,
  rel_heights = c(1, 0.08)
)

print(fig4b_final)

ggsave(
  filename = "figures/Figure4b_factorial_with_confounding_final.png",
  plot     = fig4b_uncond,
  bg       = "white",
  width    = 10, height = 8, units = "in", dpi = 300
)
