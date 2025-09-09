# Figure 1: Conceptual figure showing GMRF vs RSR basis functions
# and covariate correlation patterns

library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)
library(Matrix)

# Set up 1D spatial domain
n_locations <- 50
locations <- seq(0, 10, length.out = n_locations)

# Create adjacency matrix for 1D case (neighbors are adjacent points)
create_1d_adjacency <- function(n) {
  A <- Matrix(0, nrow = n, ncol = n)
  for(i in 1:(n-1)) {
    A[i, i+1] <- 1
    A[i+1, i] <- 1
  }
  return(A)
}

# Create GMRF precision matrix
create_gmrf_precision <- function(A, rho = 0.8) {
  I <- Diagonal(nrow(A))
  Q <- t(I - rho * A) %*% (I - rho * A)
  return(Q)
}

# Function to extract basis functions from precision matrix
extract_basis_functions <- function(Q, n_basis = 5) {
  # Eigendecomposition of precision matrix
  eigen_result <- eigen(as.matrix(Q))
  
  # Select a few representative eigenvectors (basis functions)
  # We'll take every 10th eigenvector to show different spatial scales
  indices <- seq(1, min(ncol(eigen_result$vectors), 50), by = 10)[1:n_basis]
  basis_functions <- eigen_result$vectors[, indices]
  eigenvalues <- eigen_result$values[indices]
  
  return(list(
    basis = basis_functions,
    values = eigenvalues,
    indices = indices
  ))
}

# Create data for plotting
A_1d <- create_1d_adjacency(n_locations)
Q_gmrf <- create_gmrf_precision(A_1d)

# Extract GMRF basis functions
gmrf_basis <- extract_basis_functions(Q_gmrf, n_basis = 4)

# Create a smooth covariate that correlates with some basis functions
set.seed(123)
covariate <- sin(2 * pi * locations / 5) + 0.3 * cos(4 * pi * locations / 5) + 
  rnorm(n_locations, 0, 0.1)
covariate <- scale(covariate)[,1]  # Standardize

# RSR approach: Same precision matrix, but project the spatial effects
X <- cbind(1, covariate)  # Design matrix with intercept and covariate
P <- X %*% solve(t(X) %*% X) %*% t(X)  # Projection matrix
I_n <- diag(n_locations)

# RSR uses the same basis functions as GMRF, but projects them
# to be orthogonal to the covariate space
rsr_basis_functions <- matrix(0, nrow = n_locations, ncol = 4)
for(i in 1:4) {
  # Apply orthogonal projection: (I - P) * basis_function
  # This removes any component that's in the covariate space
  rsr_basis_functions[,i] <- (I_n - P) %*% gmrf_basis$basis[,i]
}

# Prepare data for plotting
plot_data_gmrf <- data.frame(
  location = rep(locations, 4),
  basis_value = as.vector(gmrf_basis$basis),
  basis_function = rep(paste("Basis", 1:4), each = n_locations),
  type = "GMRF"
)

plot_data_rsr <- data.frame(
  location = rep(locations, 4),
  basis_value = as.vector(rsr_basis_functions),
  basis_function = rep(paste("Basis", 1:4), each = n_locations),
  type = "RSR"
)

plot_data_combined <- rbind(plot_data_gmrf, plot_data_rsr)

covariate_data <- data.frame(
  location = locations,
  covariate = covariate
)

# Calculate correlations for annotation
correlations_gmrf <- sapply(1:4, function(i) cor(covariate, gmrf_basis$basis[,i]))
correlations_rsr <- sapply(1:4, function(i) cor(covariate, rsr_basis_functions[,i]))

# Create plots
# Panel A: Covariate
plot_covariate <- ggplot(covariate_data, aes(x = location, y = covariate)) +
  geom_line(size = 1.2, color = "#440154FF") +
  geom_point(size = 2, color = "#440154FF", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "A) Simulated Covariate",
    x = "Spatial Location",
    y = "Covariate Value"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

# Panel B: GMRF Basis Functions
plot_gmrf <- plot_data_gmrf %>%
  ggplot(aes(x = location, y = basis_value, color = basis_function)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~basis_function, scales = "free_y", nrow = 2) +
  scale_color_viridis_d(option = "D") +
  theme_minimal() +
  labs(
    title = "B) GMRF Basis Functions",
    subtitle = "Some basis functions correlate with covariate",
    x = "Spatial Location",
    y = "Basis Function Value"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10)
  )

# create a dataframe with correlation annotations per facet
cor_df <- data.frame(
  basis_function = unique(plot_data_gmrf$basis_function),
  correlation = paste0("r = ", round(correlations_gmrf, 3)),
  x = 8, # or change to a position that suits your x range
  y = Inf # or pick a reasonable y
)

plot_gmrf <- plot_gmrf + 
  geom_text(
    data = cor_df,
    aes(x = x, y = y, label = correlation),
    vjust = 1.2, hjust = 1, fontface = "bold", size = 3.5, col="black"
  )

# Panel C: RSR Basis Functions  
plot_rsr <- plot_data_rsr %>%
  ggplot(aes(x = location, y = basis_value, color = basis_function)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(~basis_function, scales = "free_y", nrow = 2) +
  scale_color_viridis_d(option = "D") +
  theme_minimal() +
  labs(
    title = "C) RSR Basis Functions",
    subtitle = "Basis functions orthogonal to covariate (r ≈ 0)",
    x = "Spatial Location",
    y = "Basis Function Value"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10)
  )

# create a dataframe with correlation annotations per facet
cor_df_rsr <- data.frame(
  basis_function = unique(plot_data_rsr$basis_function),
  correlation = paste0("r = ", round(correlations_rsr, 3)),
  x = 8, # or change to a position that suits your x range
  y = Inf # or pick a reasonable y
)

plot_rsr <- plot_rsr + 
  geom_text(
    data = cor_df_rsr,
    aes(x = x, y = y, label = correlation),
    vjust = 1.2, hjust = 1, fontface = "bold", size = 3.5, col="black"
  )

# Combine plots
figure_1 <- grid.arrange(
  plot_covariate,
  plot_gmrf,
  plot_rsr,
  layout_matrix = rbind(c(1, 1),
                        c(2, 3)),
  heights = c(1, 2)
)

ggsave(
  filename = "figures/Conceptual_figure1.png",
  plot     = figure_1,
  width    = 270, height = 200, units = "mm", dpi = 300,
  bg       = "transparent"
)

# Print correlation summary
cat("\n=== CORRELATION SUMMARY ===\n")
cat("GMRF Basis Function Correlations with Covariate:\n")
for(i in 1:4) {
  cat(sprintf("  Basis %d: r = %.3f\n", i, correlations_gmrf[i]))
}

cat("\nRSR Basis Function Correlations with Covariate:\n")
for(i in 1:4) {
  cat(sprintf("  Basis %d: r = %.3f\n", i, correlations_rsr[i]))
}

cat("\nKey Insight:")
cat("\n- GMRF basis functions can be correlated with covariates")
cat("\n- RSR basis functions are orthogonal to covariates (correlation ≈ 0)")
cat("\n- This orthogonality prevents confounding between spatial and covariate effects\n")

# Optional: Create a summary comparison plot
correlation_comparison <- data.frame(
  Basis = rep(paste("Basis", 1:4), 2),
  Method = rep(c("GMRF", "RSR"), each = 4),
  Correlation = c(correlations_gmrf, correlations_rsr)
)

plot_correlation_comparison <- ggplot(correlation_comparison, 
                                      aes(x = Basis, y = Correlation, fill = Method)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_viridis_d(option = "D", begin = 0.3, end = 0.8) +
  theme_minimal() +
  labs(
    title = "D) Correlation Summary: Basis Functions vs Covariate",
    x = "Basis Function",
    y = "Correlation with Covariate",
    fill = "Method"
  ) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )

print(plot_correlation_comparison)
