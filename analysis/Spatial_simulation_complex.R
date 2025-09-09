# Figure 4: Cross-treatment comparison
# Extending the simulation to compare across sample sizes and number of false covariates

library(tinyVAST)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)

# Helper functions
rmse <- function(vec) sqrt(mean(vec^2))

# Configuration for Figure 4 - factorial design
config_fig4 <- list(
  seed = 101,
  n_simulations = 50,  # Reduced for faster computation across treatments
  sample_sizes = c(50, 100, 200),  # Multiple sample sizes
  n_false_covariates = c(1, 3, 5),  # Different numbers of false covariates (excluding the real one)
  alpha = 0.05,
  decorrelation_dist = 0.5,
  measurement_error = "normal"
)

# Modified simulation function for factorial design
run_factorial_simulation <- function(
    n_s = 100,
    n_false_covariates = 3,  # Number of FALSE covariates (real one is always included)
    decorrelation_dist = 0.5,
    measurement_error = "normal",
    n_simulations = 50,
    spatial_sd = 0.2,
    alpha = 0.05,
    include_variable_selection = TRUE
) {
  
  # Set up covariate structure: 1 real + n_false_covariates false
  n_covariates <- 1 + n_false_covariates
  covariate_effects <- c(1, rep(0, n_false_covariates))  # First is real, rest are false
  
  # Storage for results
  results <- data.frame()
  
  # Simulation parameters (exactly as in your working code)
  kappa <- sqrt(8) / decorrelation_dist
  tau <- 1 / kappa / spatial_sd / sqrt(4 * pi)
  
  cat("Running factorial simulation: n_s =", n_s, ", n_false =", n_false_covariates, "\n")
  
  for (sim in 1:n_simulations) {
    
    # Generate spatial structure (exactly as in your working code)
    xy_i <- pracma::poisson2disk(n = n_s)
    mesh <- fmesher::fm_mesh_2d(xy_i)
    spde <- fmesher::fm_fem(mesh)
    
    A_is <- Matrix::drop0(fmesher::fm_evaluator(mesh, loc = xy_i)$proj$A)
    Q_ss <- (kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2) * tau^2
    
    # Generate spatial field and covariates (exactly as in your working code)
    z_s <- as.matrix(A_is %*% rmvnorm_prec(Q = Q_ss, n = 1))
    x_sj <- as.matrix(A_is %*% rmvnorm_prec(Q = Q_ss, n = n_covariates))
    
    # Generate response (exactly as in your working code)
    linear_predictor <- z_s + x_sj %*% covariate_effects
    
    if (measurement_error == "normal") {
      SD <- 0.2
      d_s <- rnorm(n = length(linear_predictor), 
                   mean = linear_predictor, sd = SD)
      Family <- gaussian()
      invlink <- identity
    } else if (measurement_error == "poisson") {
      d_s <- rpois(n = length(linear_predictor), 
                   lambda = exp(linear_predictor))
      Family <- poisson()
      invlink <- exp
    }
    
    # Prepare data (exactly as in your working code)
    covariate_names <- c("real", paste0("false_", 1:n_false_covariates))
    dat <- data.frame(
      X1 = xy_i[,1], X2 = xy_i[,2],
      d = d_s, dtrue = linear_predictor
    )
    dat <- cbind(dat, setNames(as.data.frame(x_sj), covariate_names))
    
    # Fit models with error handling (exactly as in your working code)
    models_to_fit <- c("Covariate_Only", "VAST_GMRF", "VAST_RSR")
    if (include_variable_selection) {
      models_to_fit <- c(models_to_fit, "VAST_Selected")
    }
    
    for (model_name in models_to_fit) {
      
      tryCatch({
        
        if (model_name == "Covariate_Only") {
          # Standard regression (exactly as in your working code)
          formula_str <- paste("d ~", paste(covariate_names, collapse = " + "))
          
          if (measurement_error == "normal") {
            fit <- lm(as.formula(formula_str), data = dat)
            coef_summary <- summary(fit)$coefficients
          } else {
            fit <- glm(as.formula(formula_str), data = dat, family = Family)
            coef_summary <- summary(fit)$coefficients
          }
          
          predictions <- predict(fit, type = "response")
          
          # Store results for each covariate (exactly as in your working code)
          for (i in 1:n_covariates) {
            coef_name <- covariate_names[i]
            results <- rbind(results, data.frame(
              Model = model_name,
              Simulation = sim,
              Sample_Size = n_s,
              N_False_Covariates = n_false_covariates,  # Added this for factorial design
              Covariate = coef_name,
              True_Effect = covariate_effects[i],
              Estimate = coef_summary[coef_name, "Estimate"],
              SE = coef_summary[coef_name, "Std. Error"],
              p_value = coef_summary[coef_name, 4],  # p-value column
              Significant = coef_summary[coef_name, 4] < alpha,
              RMSE_data = rmse(dat$d - predictions),
              RMSE_latent = rmse(dat$dtrue - predictions),
              stringsAsFactors = FALSE
            ))
          }
          
        } else {
          # VAST models (exactly as in your working code)
          get_rsr <- (model_name == "VAST_RSR")
          
          # For variable selection, first fit full model
          if (model_name == "VAST_Selected") {
            # Fit full model first to get p-values
            full_fit <- tinyVAST(
              formula = as.formula(paste("d ~", paste(covariate_names, collapse = " + "))),
              data = cbind(dat, "var" = "n"),
              family = Family,
              space_term = "",
              variable_column = "var",
              space_columns = c("X1", "X2"),
              spatial_domain = mesh,
              control = tinyVASTcontrol(reml = TRUE, get_rsr = FALSE)
            )
            
            # Select significant covariates
            full_estimates <- as.list(full_fit$sdrep, what = "Estimate")$alpha_j
            full_SE <- as.list(full_fit$sdrep, what = "Std. Error")$alpha_j
            full_p_values <- 1 - pchisq((full_estimates / full_SE)^2, df = 1)
            
            # Select covariates (excluding intercept)
            selected_covs <- covariate_names[full_p_values[-1] < alpha]
            
            if (length(selected_covs) > 0) {
              formula_selected <- as.formula(paste("d ~", paste(selected_covs, collapse = " + ")))
            } else {
              formula_selected <- as.formula("d ~ 1")
            }
            
            # Fit selected model
            fit <- tinyVAST(
              formula = formula_selected,
              data = cbind(dat, "var" = "n"),
              family = Family,
              space_term = "",
              variable_column = "var",
              space_columns = c("X1", "X2"),
              spatial_domain = mesh,
              control = tinyVASTcontrol(reml = TRUE, get_rsr = TRUE)
            )
            
            estimates <- as.list(fit$sdrep, what = "Estimate", report = TRUE)$alphaprime_j
            SE <- as.list(fit$sdrep, what = "Std. Error", report = TRUE)$alphaprime_j
            
          } else {
            # Regular VAST fit
            fit <- tinyVAST(
              formula = as.formula(paste("d ~", paste(covariate_names, collapse = " + "))),
              data = cbind(dat, "var" = "n"),
              family = Family,
              space_term = "",
              variable_column = "var",
              space_columns = c("X1", "X2"),
              spatial_domain = mesh,
              control = tinyVASTcontrol(reml = TRUE, get_rsr = get_rsr)
            )
            
            if (get_rsr) {
              estimates <- as.list(fit$sdrep, what = "Estimate", report = TRUE)$alphaprime_j
              SE <- as.list(fit$sdrep, what = "Std. Error", report = TRUE)$alphaprime_j
            } else {
              estimates <- as.list(fit$sdrep, what = "Estimate")$alpha_j
              SE <- as.list(fit$sdrep, what = "Std. Error")$alpha_j
            }
          }
          
          p_values <- 1 - pchisq((estimates / SE)^2, df = 1)
          predictions <- invlink(fit$tmb_inputs$tmb_data$X_ij %*% estimates)
          
          # Store results for each covariate
          coef_names <- c("(Intercept)", covariate_names)
          for (i in 2:length(estimates)) {  # Skip intercept
            coef_idx <- i - 1
            results <- rbind(results, data.frame(
              Model = model_name,
              Simulation = sim,
              Sample_Size = n_s,
              N_False_Covariates = n_false_covariates,  # Added this for factorial design
              Covariate = covariate_names[coef_idx],
              True_Effect = covariate_effects[coef_idx],
              Estimate = estimates[i],
              SE = SE[i],
              p_value = p_values[i],
              Significant = p_values[i] < alpha,
              RMSE_data = rmse(dat$d - predictions),
              RMSE_latent = rmse(dat$dtrue - predictions),
              stringsAsFactors = FALSE
            ))
          }
        }
        
      }, error = function(e) {
        cat("Model", model_name, "failed for simulation", sim, ":", e$message, "\n")
      })
    }
    
    # Progress update
    if (sim %% 20 == 0) {
      cat("Completed", sim, "of", n_simulations, "simulations\n")
    }
  }
  
  return(results)
}

# Enhanced summary function for factorial design
create_factorial_summary <- function(results) {
  
  summary_stats <- results %>%
    group_by(Model, Covariate, Sample_Size, N_False_Covariates) %>%
    summarise(
      n_sims = n(),
      true_effect = first(True_Effect),
      
      # Bias and accuracy
      mean_estimate = mean(Estimate, na.rm = TRUE),
      bias = mean(Estimate - True_Effect, na.rm = TRUE),
      rmse_estimate = sqrt(mean((Estimate - True_Effect)^2, na.rm = TRUE)),
      
      # Standard errors
      mean_SE = mean(SE, na.rm = TRUE),
      empirical_SE = sd(Estimate, na.rm = TRUE),
      
      # Inference
      power_or_type_I = mean(Significant, na.rm = TRUE),
      coverage = mean(
        (Estimate - 1.96 * SE <= True_Effect) & 
          (Estimate + 1.96 * SE >= True_Effect), 
        na.rm = TRUE
      ),
      
      # Prediction
      mean_RMSE_data = mean(RMSE_data, na.rm = TRUE),
      mean_RMSE_latent = mean(RMSE_latent, na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    mutate(
      effect_type = ifelse(true_effect == 0, "No Effect", "Real Effect"),
      performance_metric = ifelse(true_effect == 0, "Type I Error", "Power")
    )
  
  return(summary_stats)
}

# Set seed
set.seed(config_fig4$seed)

# Run factorial simulation
cat("Starting factorial simulation for Figure 4...\n")
factorial_results <- data.frame()

for (n_s in config_fig4$sample_sizes) {
  for (n_false in config_fig4$n_false_covariates) {
    
    cat("\n--- Sample Size:", n_s, "False Covariates:", n_false, "---\n")
    
    results <- run_factorial_simulation(
      n_s = n_s,
      n_false_covariates = n_false,
      n_simulations = config_fig4$n_simulations,
      decorrelation_dist = config_fig4$decorrelation_dist,
      measurement_error = config_fig4$measurement_error,
      alpha = config_fig4$alpha,
      include_variable_selection = TRUE
    )
    
    factorial_results <- rbind(factorial_results, results)
  }
}


# Create summary for factorial results
factorial_summary <- create_factorial_summary(factorial_results)

factorial_results<- read.csv("results/factorial_simulation_results.csv")
factorial_summary<- read.csv("results/factorial_simulation_summary.csv")


model_colors <- c(
  "Covariate_Only" = "#440154FF",
  "VAST_GMRF" = "#31688EFF", 
  "VAST_RSR" = "#35B779FF",
  "VAST_Selected" = "gold"
)

# Custom x-axis labels
model_labels <- c(
  "Covariate_Only" = "GLM",
  "VAST_GMRF" = "SGLMM", 
  "VAST_RSR" = "RSR-SGLMM",
  "VAST_Selected" = "SGLMM-selected"
)

# Figure 4: RMSE_latent metric with sample sizes as rows, false covariates as columns
fig4 <- factorial_results %>%
  filter(Covariate == "real") %>%
  ggplot(aes(x = Model, y = RMSE_latent, fill = Model)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
  #geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 0.8) +
  facet_grid(rows = vars(Sample_Size), 
             cols = vars(N_False_Covariates),
             labeller = labeller(
               Sample_Size = function(x) paste("n =", x),
               N_False_Covariates = function(x) paste("False Covariates:", x)
             )) +
  scale_fill_manual(values = model_colors, labels = model_labels) +
  scale_x_discrete(labels = model_labels) +
  labs(
    #title = "Figure 4: Prediction RMSE Distributions Across Sample Sizes and Covariate Complexity",
    #subtitle = "RMSE for predicting latent density",
    x = "Model",
    y = "RMSE (unconditional prediction)",
    fill = "Model"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, margin = margin(t = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, size =10),
    axis.title.y = element_text(size = 14, margin = margin(r = 15)),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
  )

# Display the plot
print(fig4)


# Save Figure 4
ggsave(
  filename = "figures/Figure4_cross_treatment_RMSElatent_comparison_FINAL.png",
  plot = fig4,
  width = 350, height = 250, units = "mm", dpi = 300,
  bg = "transparent"
)

# Print summary statistics
cat("\n=== Figure 4 Summary Statistics ===\n")
print(head(factorial_summary))

# Save the factorial results
write.csv(factorial_results, "results/factorial_simulation_results.csv", row.names = FALSE)
write.csv(factorial_summary, "results/factorial_simulation_summary.csv", row.names = FALSE)

