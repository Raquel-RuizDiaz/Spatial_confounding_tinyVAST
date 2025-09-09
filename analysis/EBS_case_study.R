library(sf)
library(rnaturalearth)
library(terra)
library(fmesher)
library(ggplot2)
library(tinyVAST)
library(gridExtra)  # For arranging multiple plots
library(mvtnorm)  # sample from covariance
library(abind)  # Manipulate results

set.seed(123)

# Your existing setup code remains the same...
region <- "EBS"
root_dir <- here::here()
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
      trace = 1,
      get_rsr = TRUE
    )
  )
  
  # Non-spatial model
  tv_nonspatial <- tinyVAST(
    formula = as.formula(formula_str),
    data = region_data_all,
    family = tweedie(),
    control = tinyVASTcontrol(
      trace = 1,
      get_rsr = TRUE
    )
  )
  
  # Extract coefficients and covariances
  alpha0_j <- as.list(tv_nonspatial$sdrep, what = "Estimate")$alpha_j[2:3]
  V_alpha0_j <- tv_nonspatial$sdrep$cov.fixed[2:3,2:3]
  alpha_j <- as.list(tv_spatial$sdrep, what = "Estimate")$alpha_j[2:3]
  V_alpha_j <- tv_spatial$sdrep$cov.fixed[2:3,2:3]
  alphaprime_j <- as.list(tv_spatial$sdrep, what = "Estimate", report = TRUE)$alphaprime_j[2:3]
  V_alphaprime_j <- tv_spatial$sdrep$cov[2:3,2:3]
  
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

# Option 1: Base R plot with all species in one panel
png(file = file.path(run_dir, "EBS_3_species_comparison_renamed2.png"), 
    width = 15, height = 5, res = 300, units = "in")

par(mfrow = c(1, 3), mar = c(3,4,2,1), mgp = c(2,0.5,0), tck = -0.02)

# Individual plots for each species
model_labels <- c("GLM", "SGLMM", "RSR-SGLMM")

for(i in 1:length(species_list)) {
  Y_zzm_log = abind( pred_m1[[i]], pred_m2[[i]], pred_m3[[i]], along = 3 )
  Y_zzm_density <- exp(Y_zzm_log)
  
  # RESCALE EACH CURVE TO HAVE MAXIMUM VALUE OF 1.0
  Y_zzm_rescaled <- Y_zzm_density
  for(m in 1:3) {
    # Find the maximum value across all columns for this model
    max_val <- max(Y_zzm_density[,,m])
    # Rescale all values for this model
    Y_zzm_rescaled[,,m] <- Y_zzm_density[,,m] / max_val
  }
  
  plot( x=1, y = 1, type = "n", xlim = range(x_z), ylim = range(Y_zzm_rescaled),
        ylab = "Relative Density",
        xlab = "Bottom temperature (°C)",
        main = species_names[i],
        cex.lab = 1.5,    # Axis labels size
        cex.axis = 1.3,   # Axis tick labels size
        cex.main = 1.6    # Main title size
  )
  for(m in 1:3){
    lines(
      x = x_z,
      y = Y_zzm_rescaled[,1,m],
      type = "l",
      lty = "solid",
      col = viridisLite::viridis(3)[m],
      lwd = 2
    )
    polygon(
      x = c( x_z, rev(x_z) ),
      y = c( Y_zzm_rescaled[,2,m], rev(Y_zzm_rescaled[,3,m]) ),
      border = NA,
      col = viridisLite::viridis(3, alpha = 0.2)[m]
    )
  }
  
  if(i == 1) {
    legend("bottomleft",
           fill = viridisLite::viridis(3),
           legend = model_labels,
           bty = "n",
           cex = 1.5)
  }
}

dev.off()

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