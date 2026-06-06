library(RTMB)
library(Matrix)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(cowplot)
library(ggpubr)

# ============================================================
# Helper: simulate from GMRF precision
# ============================================================
rmvnorm_prec <- function(mu, prec, n.sims) {
  z0 <- matrix(rnorm(length(mu) * n.sims), ncol = n.sims)
  L  <- Cholesky(prec, super = TRUE)
  z  <- solve(L, z0, system = "Lt")
  z  <- solve(L, z,  system = "Pt")
  return(mu + as.matrix(z))
}

# ============================================================
# Spatial setup (1D transect for clarity)
# ============================================================
set.seed(123)
n_loc     <- 100
locations <- seq(0, 10, length.out = n_loc)

# 1D adjacency
A_ss <- Matrix(0, nrow = n_loc, ncol = n_loc)
for (i in 1:(n_loc - 1)) {
  A_ss[i, i + 1] <- 1
  A_ss[i + 1, i] <- 1
}

# CAR precision matrix (for illustration)
rho  <- 0.95
I_ss <- Diagonal(n = n_loc)
Q_ss <- t(I_ss - rho * A_ss) %*% (I_ss - rho * A_ss)

# ============================================================
# Simulate covariate z_s and latent field omega_s
# Both from same precision matrix so they share spatial scale
# => meaningful confounding arises by chance
# ============================================================
omega_s <- rmvnorm_prec(mu = rep(0, n_loc), prec = Q_ss, n.sims = 1)
z_s     <- rmvnorm_prec(mu = rep(0, n_loc), prec = Q_ss, n.sims = 1)

# Scale both for comparability
omega_s <- scale(omega_s)[, 1]
z_s     <- scale(z_s)[, 1]

# ============================================================
# Fit SGLMM via RTMB
# ============================================================
beta_true <- 1
y_s <- z_s * beta_true + omega_s + rnorm(n_loc, 0, 0.2)

X   <- cbind(1, z_s)   # design matrix: intercept + covariate

jnll_sglmm <- function(parlist) {
  getAll(parlist)
  Q   <- exp(ln_tau) * t(I_ss - exp(ln_rho) * A_ss) %*%
    (I_ss - exp(ln_rho) * A_ss)
  nll <- -1 * dgmrf(omega_s_est, Q = Q, log = TRUE)
  nll <- nll - sum(dnorm(y_s,
                         mean = X %*% betas + omega_s_est,
                         sd   = exp(ln_sd),
                         log  = TRUE))
  return(nll)
}

parlist <- list(
  omega_s_est = rep(0, n_loc),
  ln_tau      = 0,
  ln_rho      = 0,
  ln_sd       = 0,
  betas       = rep(0, ncol(X))
)

obj <- RTMB::MakeADFun(jnll_sglmm, parlist, random = "omega_s_est")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdrep <- sdreport(obj)

# Extract SGLMM estimates
beta_hat    <- summary(sdrep, "fixed")[grep("betas", rownames(summary(sdrep, "fixed"))), 1]
omega_hat   <- summary(sdrep, "random")[, 1]   # empirical Bayes omega

# ============================================================
# Post-hoc RSR adjustment (Equations 4-6 of manuscript)
# ============================================================

# Projection matrix onto covariate space
P_ss <- X %*% solve(t(X) %*% X) %*% t(X)

# Indirect effect: component of A*omega_hat correlated with X (Eq. 4)
# Here A = I (1D, no interpolation needed)
delta_hat <- solve(t(X) %*% X) %*% t(X) %*% omega_hat   # (X'X)^-1 X' A*omega

# RSR total effect (Eq. 5)
beta_star <- beta_hat + delta_hat

# RSR-adjusted random effect (Eq. 6)
omega_star <- omega_hat - X %*% solve(t(X) %*% X) %*% t(X) %*% omega_hat

# Verify: fitted values are identical
fitted_sglmm <- X %*% beta_hat  + omega_hat
fitted_rsr   <- X %*% beta_star + omega_star
cat("Max difference in fitted values (should be ~0):",
    max(abs(fitted_sglmm - fitted_rsr)), "\n")

# Correlations
r_omega_z      <- cor(omega_hat,   z_s)
r_omega_star_z <- cor(omega_star,  z_s)
r_delta_x      <- cor(as.vector(X %*% delta_hat), z_s)

cat(sprintf("Cor(omega_hat,   z_s) = %.3f\n", r_omega_z))
cat(sprintf("Cor(omega_star,  z_s) = %.3f\n", r_omega_star_z))

# ============================================================
# Build plotting data
# ============================================================

plot_df <- data.frame(
  location   = locations,
  z_s        = as.vector(z_s),
  omega_hat  = as.vector(omega_hat),
  delta_comp = as.vector(X %*% delta_hat),   # correlated component
  omega_star = as.vector(omega_star)          # orthogonal residual
)

# ============================================================
# Panel A: Simulated covariate
# ============================================================
p_A <- ggplot(plot_df, aes(x = location, y = z_s)) +
  geom_line(color = "#440154FF", size = 1.1) +
  geom_point(color = "#440154FF", size = 1.5, alpha = 0.6) +
  theme_bw() +
  labs(
    title = "A) Simulated Covariate",
    x     = "Spatial Location",
    y     = "Covariate Value"
  ) +
  theme(plot.title = element_text(size = 13, face = "bold"),
        panel.grid      = element_blank())

# ============================================================
# Panel B: Realized SGLMM random effect A*omega_hat
# annotated with its correlation with z_s
# ============================================================
p_B <- ggplot(plot_df, aes(x = location, y = omega_hat)) +
  geom_line(color = "#31688EFF", size = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
  annotate("text",
           x     = max(locations) * 0.75,
           y     = max(plot_df$omega_hat) * 0.9,
           label = sprintf("r = %.3f", r_omega_z),
           fontface = "bold", size = 4) +
  theme_bw() +
  labs(
    title    = "B) SGLMM Random Effect",
    subtitle = expression(paste("Realized ", hat(omega),
                                " — partially correlated with z")),
    x = "Spatial Location",
    y = expression(hat(omega))
  ) +
  theme(plot.title = element_text(size = 13, face = "bold"),
        panel.grid      = element_blank())

# ============================================================
# Panel C: Decomposition
# omega_hat = delta_component (correlated) + omega_star (orthogonal)
# ============================================================
# Use plain strings as factor levels
component_levels <- c("Unadjusted", "indirect", "orthogonal")

component_labels <- c(
  "Unadjusted"      = "A*hat(omega)~(Unadjusted~random~effect)",
  "indirect"   = "X*hat(delta)~(indirect~-~attributed~to~hat(beta)^'*')",
  "orthogonal" = "A*hat(omega)^'*'~(residual~-~orthogonal~to~z)"
)

decomp_df <- data.frame(
  location  = rep(locations, 3),
  value     = c(plot_df$omega_hat,
                plot_df$delta_comp,
                plot_df$omega_star),
  component = factor(
    rep(component_levels, each = n_loc),
    levels = component_levels
  )
)

# Annotation data frame using the same plain string keys
annot_df <- data.frame(
  component = factor(component_levels, levels = component_levels),
  label     = c(
    sprintf("r = %.3f", r_omega_z),
    sprintf("r = %.3f", r_delta_x),
    sprintf("r = %.3f", r_omega_star_z)
  ),
  x = max(locations) * 0.75,
  y = Inf
)

p_C <- ggplot(decomp_df, aes(x = location, y = value, color = component)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
  facet_wrap(
    ~ component,
    nrow     = 3,
    scales   = "free_y",
    labeller = as_labeller(component_labels, default = label_parsed)
  ) +
  geom_text(
    data        = annot_df,
    aes(x = x, y = y, label = label),
    vjust       = 1.4,
    fontface    = "bold",
    size        = 3.8,
    color       = "black",
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = c("#440154FF", "#FDE725FF", "#35B779FF")) +
  theme_bw() +
  labs(
    title    = "C) RSR Decomposition of Random Effect",
    #subtitle = expression(paste(
    #  "A", hat(omega)^"*", " = A", hat(omega), " \u2212 X", hat(delta))),
    x = "Spatial Location",
    y = "Value"
  ) +
  theme(
    plot.title      = element_text(size = 13, face = "bold"),
    legend.position = "none",
    strip.text      = element_text(size = 10),
    panel.grid      = element_blank()
  )

# ============================================================
# Assemble figure
# ============================================================
fig2 <- plot_grid(
  plot_grid(p_A, p_B, ncol = 2, labels = c("", "")),
  p_C,
  ncol        = 1,
  rel_heights = c(1, 1.6)
)

ggsave(
  filename = "figures/Figure2_RSR_conceptual_final.png",
  plot     = fig2,
  bg       = "white",
  width    = 270, height = 220, units = "mm", dpi = 300
)
