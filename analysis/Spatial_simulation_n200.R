library(tinyVAST)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(cv)

# Configuration
config <- list(
  seed = 101,
  n_simulations = 50,
  sample_sizes = 200,
  alpha = 0.05,  # Standard significance level
  decorrelation_dist = 0.75,
  measurement_error = "normal",
  spatial_sd = 0.2,
  measurement_sd = 0.2,
  include_cv_rmse = FALSE
)

# Helper functions
rmse <- function(vec) sqrt(mean(vec^2))

# Main simulation function with enhanced features
run_enhanced_simulation <- function(
    n_s = 200,
    n_covariates = 4,  # 1 real + 3 false (from first code)
    covariate_effects = c(1, 0, 0, 0),  # Effects for each covariate
    decorrelation_dist = 0.75,
    measurement_error = "normal",
    spatial_sd = 0.2,
    measurement_sd = 0.2,
    n_simulations = 50,
    alpha = 0.05,
    include_variable_selection = TRUE,  # Option to include from first code,
    include_cv_rmse = TRUE
) {
  
  # Storage for results
  results <- data.frame()
  
  # Simulation parameters
  kappa <- sqrt(8) / decorrelation_dist
  tau <- 1 / kappa / spatial_sd / sqrt(4 * pi)
  
  cat("Running enhanced simulation with n_s =", n_s, "\n")
  
  for (sim in 1:n_simulations) {
    
    # Generate spatial structure (same as both codes)
    xy_i <- pracma::poisson2disk(n = n_s)
    mesh <- fmesher::fm_mesh_2d(xy_i)
    spde <- fmesher::fm_fem(mesh)
    
    A_is <- Matrix::drop0(fmesher::fm_evaluator(mesh, loc = xy_i)$proj$A)
    Q_ss <- (kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2) * tau^2
    
    # Generate spatial field and covariates
    z_s <- as.matrix(A_is %*% rmvnorm_prec(Q = Q_ss, n = 1))
    x_sj <- as.matrix(A_is %*% rmvnorm_prec(Q = Q_ss, n = n_covariates))
    
    # Generate response
    linear_predictor <- z_s + x_sj %*% covariate_effects
    
    if (measurement_error == "normal") {
      SD <- measurement_sd
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
    
    # Prepare data
    covariate_names <- c("real", paste0("false_", 1:(n_covariates-1)))
    dat <- data.frame(
      X1 = xy_i[,1], X2 = xy_i[,2],
      d = d_s, dtrue = linear_predictor
    )
    dat <- cbind(dat, setNames(as.data.frame(x_sj), covariate_names))
    
    # Fit models with error handling
    models_to_fit <- c("Covariate_Only", "VAST_GMRF", "VAST_RSR")
    if (include_variable_selection) {
      models_to_fit <- c(models_to_fit, "VAST_Selected")
    }
    
    for (model_name in models_to_fit) {
      
      tryCatch({
        
        if (model_name == "Covariate_Only") {
          # Standard regression
          formula_str <- paste("d ~", paste(covariate_names, collapse = " + "))
          
          if (measurement_error == "normal") {
            fit <- lm(as.formula(formula_str), data = dat)
            coef_summary <- summary(fit)$coefficients
          } else {
            fit <- glm(as.formula(formula_str), data = dat, family = Family)
            coef_summary <- summary(fit)$coefficients
          }
          
          predictions <- predict(fit, type = "response")
          if( isTRUE(include_cv_rmse) ){
            CV = cv(fit, k = 5)
          }else{
            CV = list("CV crit" = NA)
          }
          
          # Store results for each covariate
          for (i in 1:n_covariates) {
            coef_name <- covariate_names[i]
            results <- rbind(results, data.frame(
              Model = model_name,
              Simulation = sim,
              Sample_Size = n_s,
              Covariate = coef_name,
              True_Effect = covariate_effects[i],
              Estimate = coef_summary[coef_name, "Estimate"],
              SE = coef_summary[coef_name, "Std. Error"],
              p_value = coef_summary[coef_name, 4],  # p-value column
              Significant = coef_summary[coef_name, 4] < alpha,
              RMSE_data = rmse(dat$d - predictions),
              RMSE_latent = rmse(dat$dtrue - predictions),
              RMSE_cv_with_gmrf = sqrt( CV[['CV crit']] ),
              RMSE_with_gmrf = NA,
              stringsAsFactors = FALSE
            ))
          }
          
        } else {
          # VAST models
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
              control = tinyVASTcontrol(reml = TRUE, get_rsr = TRUE)
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
            # Embedding cv::cv in a function conflicts with the way it expects to pass objects ... assigning to global memory to allow it to work
            if(isTRUE(include_cv_rmse)) formula_selected <<- formula_selected
            
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
              control = tinyVASTcontrol(reml = TRUE, get_rsr = TRUE)
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
          pred_with_gmrf = predict( fit )
          
          # Embedding cv::cv in a function conflicts with the way it expects to pass objects ... assigning to global memory to allow it to work
          if( isTRUE(include_cv_rmse) ){
            mesh <<- mesh
            Family <<- Family
            covariate_names <<- covariate_names
            CV = cv( fit, k = 5 )
          }else{
            CV = list("CV crit" = NA)
          }
          
          # Store results for each covariate
          coef_names <- c("(Intercept)", covariate_names)
          for (i in 2:length(estimates)) {  # Skip intercept
            coef_idx <- i - 1
            results <- rbind(results, data.frame(
              Model = model_name,
              Simulation = sim,
              Sample_Size = n_s,
              Covariate = covariate_names[coef_idx],
              True_Effect = covariate_effects[coef_idx],
              Estimate = estimates[i],
              SE = SE[i],
              p_value = p_values[i],
              Significant = p_values[i] < alpha,
              RMSE_data = rmse(dat$d - predictions),
              RMSE_latent = rmse(dat$dtrue - predictions),
              RMSE_cv_with_gmrf = sqrt( CV[['CV crit']] ),
              RMSE_with_gmrf = rmse(dat$dtrue - pred_with_gmrf),
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

# Enhanced summary function
create_enhanced_summary <- function(results) {
  
  summary_stats <- results %>%
    group_by(Model, Covariate, Sample_Size) %>%
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
set.seed(config$seed)

# FIXED: Loop over all sample sizes from config
all_results <- data.frame()

for (sample_size in config$sample_sizes) {
  cat("\n=== Running simulation for sample size:", sample_size, "===\n")
  
  results <- run_enhanced_simulation(
    n_s = sample_size,  # Use the sample size from the loop
    n_simulations = config$n_simulations,  # Use config values
    decorrelation_dist = config$decorrelation_dist,
    measurement_error = config$measurement_error,
    measurement_sd = config$measurement_sd,
    spatial_sd = config$spatial_sd,
    include_cv_rmse = config$include_cv_rmse,
    alpha = config$alpha,
    include_variable_selection = TRUE
  )
  
  # Combine results
  all_results <- rbind(all_results, results)
}

#write.csv(all_results, "results/spatial_model_comparison_nj4_allresults.csv", row.names = FALSE)

# Create summary for all results
summary_results <- create_enhanced_summary(all_results)
print(head(summary_results))

# Save results
#write.csv(summary_results, "results/spatial_model_comparison_nj4.csv", row.names = FALSE)

# Reload
if( FALSE ){
  results_dir = R'(C:\Users\James.Thorson\Desktop\Git\spatio-temporal-autocorrelation\results)'
  summary_results = read.csv( file.path(results_dir, "spatial_model_comparison_nj4.csv") )
}

summary_results<- read.csv("results/spatial_model_comparison_nj4.csv")
all_results<- read.csv("results/spatial_model_comparison_nj4_allresults.csv")

# Summarize
tmp = subset(all_results, Covariate == "real" )
tapply( tmp$RMSE_latent, INDEX = list(tmp$Model,tmp$Sample_Size), FUN = mean )

create_performance_plots <- function(summary_results) {
  
  # Color palette
  model_colors <- c(
    "Covariate_Only" = "#440154FF",
    "VAST_Normal" = "#31688EFF", 
    "VAST_RSR" = "#35B779FF",
    "VAST_Selected" = "gold"
  )
  
  # Custom x-axis labels
  model_labels <- c(
    "Covariate_Only" = "GLM",
    "VAST_Normal" = "SGLMM", 
    "VAST_RSR" = "RSR-SGLMM",
    "VAST_Selected" = "SGLMM-selected"
  )
  
  # Plot 1: Type I Error Rate (False Discovery)
  covariate_labels <- c("false_1" = "False 1", "false_2" = "False 2", "false_3" = "False 3")
  
  type_I_plot <- summary_results %>%
    filter(true_effect == 0) %>%
    ggplot(aes(x = factor(Sample_Size), y = power_or_type_I, fill = Model)) +
    geom_col(position = "dodge", alpha = 0.8) +
    facet_wrap(~Covariate, scales = "free_y", labeller = labeller(Covariate = covariate_labels)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
    scale_fill_manual(values = model_colors, labels = model_labels) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      #title = "Type I Error Rate (False Discovery)",
      #subtitle = "Proportion of times covariate is detected when no effect exists",
      x = "Sample Size",
      y = "Type I Error Rate",
      fill = "Model"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      strip.text = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 13, margin = margin(t = 15)),
      axis.text.x = element_text(angle = 45, hjust = 1, size =10),
      axis.title.y = element_text(size = 13, margin = margin(r = 15)),
      axis.text.y = element_text(size = 12),
      legend.position = "none",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
    )
  
  # Plot 2: Power (True Discovery)
  power_plot <- summary_results %>%
    filter(true_effect == 1) %>%
    ggplot(aes(x = factor(Sample_Size), y = power_or_type_I, fill = Model)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = model_colors, labels = model_labels) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      title = "Statistical Power (True Discovery)",
      subtitle = "Proportion of times covariate is detected when effect exists",
      x = "Sample Size",
      y = "Power",
      fill = "Model"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom"
    )
  
  # Plot 3: Bias in estimates
  bias_plot <- summary_results %>%
    ggplot(aes(x = factor(Sample_Size), y = bias, fill = Model)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_manual(values = model_colors, labels = model_labels) +
    facet_wrap(~Covariate, scales = "free_y") +
    labs(
      title = "Bias in Covariate Estimates",
      subtitle = "Average difference between estimated and true effect",
      x = "Sample Size",
      y = "Bias",
      fill = "Model"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom"
    )
  
  # Plot 4: Coverage probability
  coverage_plot <- summary_results %>%
    ggplot(aes(x = factor(Sample_Size), y = coverage, fill = Model)) +
    geom_col(position = "dodge", alpha = 0.8) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", size = 1) +
    scale_fill_manual(values = model_colors, labels = model_labels) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    facet_wrap(~effect_type) +
    labs(
      title = "95% Confidence Interval Coverage",
      subtitle = "Proportion of times CI contains true parameter value",
      x = "Sample Size",
      y = "Coverage Probability",
      fill = "Model"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "bottom"
    )
  
  return(list(
    type_I = type_I_plot,
    power = power_plot,
    bias = bias_plot,
    coverage = coverage_plot
  ))
}

# Real effect violin plot
model_colors <- c(
  "Covariate_Only" = "#440154FF",
  "VAST_Normal" = "#31688EFF", 
  "VAST_RSR" = "#35B779FF",
  "VAST_Selected" = "gold"
)

# Custom x-axis labels
model_labels <- c(
  "Covariate_Only" = "GLM",
  "VAST_Normal" = "SGLMM", 
  "VAST_RSR" = "RSR-SGLMM",
  "VAST_Selected" = "SGLMM-selected"
)

rmse_slope_violin_plot <- ggplot(subset(all_results, Covariate=="real")) +
  geom_violin(aes(x = Model, y = Estimate, fill = Model)) +
  # facet_grid removed and title added below
  labs(title = "Coefficient estimate") + 
  scale_fill_manual(values = model_colors, labels = model_labels) +
  scale_x_discrete(labels = model_labels) +
  theme_bw() +
  theme(
    # Title is now centered and sized like the old facet label
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
    plot.subtitle = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 13, margin = margin(t = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, size =10),
    axis.title.y = element_text(size = 13, margin = margin(r = 15)),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

rmse_pred_boxplot <- ggplot(subset(all_results, Covariate=="real")) +
  geom_boxplot(aes(x = Model, y = RMSE_latent, fill = Model)) +
  # facet_grid removed and title added below
  scale_fill_manual(values = model_colors, labels = model_labels) +
  scale_x_discrete(labels = model_labels) +
  labs(
    title = "Unconditional prediction",
    x = "Model",
    y = "RMSE latent",
    fill = "Model"
  ) +
  theme_bw() +
  theme(
    # Title is now centered and sized like the old facet label
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
    plot.subtitle = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 13, margin = margin(t = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, size =10),
    axis.title.y = element_text(size = 13, margin = margin(r = 15)),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

rmse_pred_gmrf_boxplot <- ggplot(subset(all_results, Covariate=="real")) +
  geom_boxplot(aes(x = Model, y = RMSE_with_gmrf, fill = Model)) +
  # facet_grid removed and title added below
  scale_fill_manual(values = model_colors, labels = model_labels) +
  scale_x_discrete(labels = model_labels) +
  labs(
    title = "Conditional prediction",
    x = "Model",
    y = "RMSE with GMRF",
    fill = "Model"
  ) +
  theme_bw() +
  theme(
    # Title is now centered and sized like the old facet label
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 13, margin = margin(t = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, size =10),
    axis.title.y = element_text(size = 13, margin = margin(r = 15)),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )



# Create visualizations
plots <- create_performance_plots(summary_results)

# Display plots
print(plots$type_I)

#print(plots$power)
#print(plots$bias)
#print(plots$coverage)
print(rmse_slope_violin_plot)

fig3 <- grid.arrange(
  plots$type_I,
  rmse_slope_violin_plot,
  rmse_pred_gmrf_boxplot,
  rmse_pred_boxplot,
  ncol = 2 #,
  #top = "Figure 3: Spatial Model Performance Comparison"
)

# Option 1B: Using grid.arrange with shared legend on the side
library(grid)
library(ggpubr)
library(gridExtra)
library(cowplot)

legend_plot <- ggplot(subset(all_results, Covariate=="real")) +
  geom_boxplot(aes(x = Model, y = RMSE_latent, fill = Model)) +
  scale_fill_manual(values = model_colors, labels = model_labels) +
  labs(fill = "Model") +
  theme_bw() +
  theme(
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )

# Extract just the legend
shared_legend <- ggpubr::get_legend(legend_plot)

# Arrange the four plots with labels
fig3_panels <- plot_grid(
  plots$type_I,
  rmse_slope_violin_plot,
  rmse_pred_gmrf_boxplot,
  rmse_pred_boxplot,
  labels = c("A)", "B)", "C)", "D)"),
  label_size = 16,
  ncol = 2
  #label_x = c(0.05, 0.05, 0.05, 0.05)  # move C & D a bit right
  #label_y = c(1, 1, 1, 1)         # keep vertical position same
)


# Add shared legend below
fig3_option1b <- plot_grid(
  fig3_panels,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

# Save
ggsave(
  filename = "figures/Figure2_Spatial_model_performance_comparison2.png",
  plot     = fig3_option1b,
  bg       = "white",
  width    = 350, height = 250, units = "mm", dpi = 300
)
