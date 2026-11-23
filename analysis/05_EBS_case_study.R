library(sf)
library(rnaturalearth)
library(terra)
library(fmesher)
library(ggplot2)
library(tinyVAST)
library(gridExtra)  # For arranging multiple plots
library(mvtnorm)  # sample from covariance
library(abind)  # Manipulate results
cat("tinyVAST version:", as.character(packageVersion("tinyVAST")), "\n")

set.seed(123)
options(warn = 1)
NL <- 2  # two Newton steps after nlminb so every fit converges (max |gradient| < 1e-4)
cat("newton_loops =", NL, "\n")

# setup
region <- "EBS"
root_dir <- "."
data_dir <- file.path(root_dir, "data")
region_data_all <- read.csv(file.path(data_dir, "all_EBS_data_2021.csv"))
bdepth <- rast(file.path(data_dir, "bdepth.tif"))

# Process coordinates
loc <- st_multipoint(as.matrix(region_data_all[, c("lon", "lat")]), dim = "XY")
sf_loc <- st_sfc(loc, crs = st_crs(bdepth))
sf_loc <- st_transform(sf_loc, crs = 4326)
region_data_all[, c("lon", "lat")] = st_coordinates(sf_loc)[,1:2]

# Create mesh
mesh <- fmesher::fm_mesh_2d(region_data_all[, c("lon", "lat")], cutoff = 0.5)

# Define species list
species_list <- c("a_poll", "tanner", "a_pcod")
species_names <- c("Adult Pollock", "Tanner Crab", "Adult Pacific Cod")

# Function to fit models and extract coefficients
fit_species_models <- function(species_col) {
  # Create formula dynamically
  formula_str <- paste(species_col, "~ 1 + btemp + I(btemp^2) + factor(year)")
  
  # Spatial model
  tv_spatial <- tinyVAST(
    formula = as.formula(formula_str),
    data = region_data_all,
    spatial_domain = mesh,
    space_columns = c("lon","lat"),
    space_term = "",
    family = tweedie(),
    control = tinyVASTcontrol(
      trace = 0,
      newton_loops = NL,
      get_rsr = TRUE
    )
  )
  
  # Non-spatial model
  tv_nonspatial <- tinyVAST(
    formula = as.formula(formula_str),
    data = region_data_all,
    family = tweedie(),
    control = tinyVASTcontrol(
      trace = 0,
      newton_loops = NL,
      get_rsr = TRUE
    )
  )
  
  cat(species_col, "| spatial max|grad| =", signif(max(abs(tv_spatial$sdrep$gradient.fixed)),3), "| non-spatial max|grad| =", signif(max(abs(tv_nonspatial$sdrep$gradient.fixed)),3), "| pdHess:", tv_spatial$sdrep$pdHess, tv_nonspatial$sdrep$pdHess, "\n")
  # Extract coefficients and covariances
  alpha0_j <- as.list(tv_nonspatial$sdrep, what = "Estimate")$alpha_j[2:3]
  V_alpha0_j <- tv_nonspatial$sdrep$cov.fixed[2:3,2:3]
  alpha_j <- as.list(tv_spatial$sdrep, what = "Estimate")$alpha_j[2:3]
  V_alpha_j <- tv_spatial$sdrep$cov.fixed[2:3,2:3]
  alphaprime_j <- as.list(tv_spatial$sdrep, what = "Estimate", report = TRUE)$alphaprime_j[2:3]
  # Covariance of the RSR slopes: use the rows of sdrep$cov labelled alphaprime_j (robust to ADREPORT order)
  rep_names <- names(tv_spatial$sdrep$value)
  idx <- which(rep_names == "alphaprime_j")[2:3]
  cat(species_col, ": ADREPORT blocks =", paste(unique(rep_names), collapse = ", "), "| alphaprime_j rows used =", paste(idx, collapse = ","), "| original index [2:3] correct =", identical(as.integer(idx), 2:3), "\n")
  V_alphaprime_j <- tv_spatial$sdrep$cov[idx, idx]
  
  return(list(
    alpha0_j = alpha0_j,
    V_alpha0_j = V_alpha0_j,
    alpha_j = alpha_j,
    V_alpha_j = V_alpha_j,
    alphaprime_j = alphaprime_j,
    V_alphaprime_j = V_alphaprime_j
  ))
}

# Fit models for all species
cat("Fitting models for all species...\n")
species_results <- list()
for(i in 1:length(species_list)) {
  cat(paste("Fitting models for", species_names[i], "...\n"))
  species_results[[i]] <- fit_species_models(species_list[i])
}
names(species_results) <- species_list

# Create temperature sequence
x_z <- seq(quantile(region_data_all$btemp, 0.01), 
           quantile(region_data_all$btemp, 0.99), 
           length = 100)
X_zp <- cbind(x_z, x_z^2)

# Calculate predictions for all species
# z: x-axis value
# m:  model
# i:  species
# p:  parameter
confidence_interval = c(0.025, 0.975)
for(m in 1:3 ){
  pred_i = list()
  for(i in 1:length(species_list)) {
    if(m==1){
      beta_p = species_results[[i]]$alpha0_j
      Cov_pp = species_results[[i]]$V_alpha0_j
    }
    if(m==2){
      beta_p = species_results[[i]]$alpha_j
      Cov_pp = species_results[[i]]$V_alpha_j
    }
    if(m==3){
      beta_p = species_results[[i]]$alphaprime_j
      Cov_pp = species_results[[i]]$V_alphaprime_j
    }
    y_z <- X_zp %*% beta_p
    bhat_pr = t(rmvnorm( n = 2000, mean = beta_p, sigma = Cov_pp )) ## here I replaced Cov with Cov_pp
    y_zr <- X_zp %*% bhat_pr
    # MAYBE STANDARDARIZE FOR THE MEAN ?
    pred_i[[i]] <- cbind( "hat" = y_z, t(apply(y_zr, MARGIN=1, FUN = quantile, prob = confidence_interval)) )
  }
  if(m==1){
    pred_m1 = pred_i
  }
  if(m==2){
    pred_m2 = pred_i
  }
  if(m==3){
    pred_m3 = pred_i
  }
}

# Create the summary figure
run_dir <- file.path(root_dir, "figures")
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)


# Print summary statistics
cat("\nSummary of temperature coefficients:\n")
for(i in 1:length(species_list)) {
  cat(paste("\n", species_names[i], ":\n"))
  cat(paste("  Non-spatial - Linear:", round(species_results[[i]]$alpha0_j[1], 3), 
            "Quadratic:", round(species_results[[i]]$alpha0_j[2], 3), "\n"))
  cat(paste("  Spatial - Linear:", round(species_results[[i]]$alpha_j[1], 3), 
            "Quadratic:", round(species_results[[i]]$alpha_j[2], 3), "\n"))
  cat(paste("  Spatial RSR - Linear:", round(species_results[[i]]$alphaprime_j[1], 3), 
            "Quadratic:", round(species_results[[i]]$alphaprime_j[2], 3), "\n"))
}


##########################################
## Figure
##########################################
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggpubr)

# Build long-format dataframe
model_names <- c("GLM", "SGLMM", "RSR-SGLMM")
model_keys  <- c("m1", "m2", "m3")
pred_list   <- list(pred_m1, pred_m2, pred_m3)

plot_df <- data.frame()

for (i in seq_along(species_list)) {
  for (m in 1:3) {
    
    pred_raw <- exp(pred_list[[m]][[i]])
    max_val  <- max(pred_raw)
    pred_scaled <- pred_raw / max_val
    
    plot_df <- rbind(plot_df, data.frame(
      species   = species_names[i],
      model     = model_names[m],
      temp      = x_z,
      mean      = pred_scaled[, 1],
      lower     = pred_scaled[, 2],
      upper     = pred_scaled[, 3]
    ))
  }
}

# Factor ordering
plot_df$model <- factor(plot_df$model, levels = model_names)

model_colors <- c(
  "GLM"       = "#440154FF",
  "SGLMM"     = "#31688EFF",
  "RSR-SGLMM" = "#35B779FF"
)

p_case <- ggplot(plot_df, aes(x = temp, color = model, fill = model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), 
              alpha = 0.2, color = NA) +
  geom_line(aes(y = mean), linewidth = 1.1) +
  facet_wrap(~ species, nrow = 1) +
  scale_color_manual(values = model_colors) +
  scale_fill_manual(values  = model_colors) +
  labs(
    x     = "Bottom temperature (°C)",
    y     = "Relative density",
    color = "Model",
    fill  = "Model"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    strip.text      = element_text(size = 12, face = "bold"),
    axis.title.x    = element_text(size = 13, margin = margin(t = 15)),
    axis.text.x     = element_text(size = 10),
    axis.title.y    = element_text(size = 13, margin = margin(r = 15)),
    axis.text.y     = element_text(size = 12),
    legend.position = "bottom",
    legend.title    = element_text(size = 14, face = "bold"),
    legend.text     = element_text(size = 12),
    panel.grid      = element_blank()
    
  )

ggsave(
  filename = file.path(run_dir, "Figure4_case_study.png"),
  plot     = p_case,
  bg       = "white",
  width    = 10, height = 4, units = "in", dpi = 300
)


# Save fitted coefficients and plotted curves
out_dir <- file.path(root_dir, "results")
saveRDS(species_results, file.path(out_dir, "EBS_case_study_coefficients.rds"))
coef_df <- do.call(rbind, lapply(seq_along(species_list), function(i) data.frame(
  species = species_names[i], model = rep(c("GLM", "SGLMM", "RSR-SGLMM"), each = 2), term = rep(c("btemp", "btemp^2"), 3),
  estimate = c(species_results[[i]]$alpha0_j, species_results[[i]]$alpha_j, species_results[[i]]$alphaprime_j),
  se = sqrt(c(diag(species_results[[i]]$V_alpha0_j), diag(species_results[[i]]$V_alpha_j), diag(species_results[[i]]$V_alphaprime_j))))))
write.csv(coef_df, file.path(out_dir, "EBS_case_study_coefficients.csv"), row.names = FALSE)
write.csv(plot_df, file.path(out_dir, "EBS_case_study_curves.csv"), row.names = FALSE)
