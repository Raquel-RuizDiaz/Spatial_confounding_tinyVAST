# Simulation 2 (spatial latent variable as confounder): Figure 2 and its two tables
# Reads the saved results from 02_Simulation2_confounder.R (no simulation is run).

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

rmse <- function(vec) sqrt(mean(vec^2))

results_shared <- read.csv("results/Simulation2_confounder_results.csv") %>%
  filter(Model != "VAST_Selected")
true_cov <- results_shared %>% filter(Covariate == "real")

# Tables
slope_summary <- true_cov %>%
  group_by(Model, Xreal_to_Z) %>%
  summarise(
    mean_corr   = mean(Realized_Corr),
    mean_total  = mean(Total_Effect),
    mean_beta   = mean(Estimate),
    sd_beta     = sd(Estimate),
    bias_direct = mean(Estimate - True_Effect),
    rmse_direct = rmse(Estimate - True_Effect),
    bias_total  = mean(Estimate - Total_Effect),
    rmse_total  = rmse(Estimate - Total_Effect),
    rmse_uncond = mean(RMSE_latent),
    slope_field = mean(Field_slope_on_X),
    .groups = "drop"
  )

delta_rep <- results_shared %>%
  filter(Model %in% c("VAST_GMRF", "VAST_RSR")) %>%
  select(Model, Simulation, Xreal_to_Z, Covariate, Estimate) %>%
  pivot_wider(names_from = Model, values_from = Estimate) %>%
  mutate(delta = VAST_RSR - VAST_GMRF) %>%
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
cat("Pooled cor(rmse gap, delta norm):", cor(delta_rep$rmse_gap, delta_rep$delta_norm), "\n")

delta_summary <- delta_rep %>%
  group_by(Xreal_to_Z) %>%
  summarise(mean_delta = mean(delta_true), sd_delta = sd(delta_true), mean_gap = mean(rmse_gap), .groups = "drop")

num <- function(x, k = 2) {
  s <- sprintf(paste0("%.", k, "f"), x)
  paste0("$", sub("^-(0\\.0+)$", "\\1", s), "$")
}
lab <- c(Covariate_Only = "GLM", VAST_GMRF = "SGLMM", VAST_RSR = "RSR-SGLMM")

tex <- c(
  "% Requires \\usepackage{booktabs}", "",
  "\\begin{table}[ht]", "\\centering", "\\small",
  "\\caption{Slope estimates for the true covariate in the shared-field sensitivity analysis (50 replicates per level). $\\lambda$ is the strength with which the covariate shares the latent field, and $r$ is the mean realized correlation between $z_s$ and $\\omega_s$. Total is the mean total effect, calculated in each replicate as $1 + \\mathrm{cov}(\\omega, z)/\\mathrm{var}(z)$. Mean $\\hat{\\beta}$ and SD are the mean and standard deviation of the estimated slope across replicates. Bias is the mean of $\\hat{\\beta}$ minus the target, and RMSE is the root mean squared difference between $\\hat{\\beta}$ and the target. GLM and RSR-SGLMM slope estimates are identical.}",
  "\\label{tab:si_shared_slope}", "\\begin{tabular}{llrrrrrrrr}", "\\toprule",
  " & & & & & & \\multicolumn{2}{c}{Direct effect} & \\multicolumn{2}{c}{Total effect} \\\\",
  "\\cmidrule(lr){7-8}\\cmidrule(lr){9-10}",
  "Estimator & $\\lambda$ & $r$ & Total & Mean $\\hat{\\beta}$ & SD & Bias & RMSE & Bias & RMSE \\\\", "\\midrule"
)
for (m in names(lab)) {
  x <- slope_summary[slope_summary$Model == m, ]
  for (j in seq_len(nrow(x))) {
    tex <- c(tex, paste0(if (j == 1) lab[[m]] else "", " & ", x$Xreal_to_Z[j], " & ", num(x$mean_corr[j]), " & ",
      num(x$mean_total[j]), " & ", num(x$mean_beta[j]), " & ", num(x$sd_beta[j]), " & ", num(x$bias_direct[j]), " & ",
      num(x$rmse_direct[j]), " & ", num(x$bias_total[j]), " & ", num(x$rmse_total[j]), " \\\\"))
  }
  if (m != "VAST_RSR") tex <- c(tex, "\\addlinespace")
}
glm <- slope_summary[slope_summary$Model == "Covariate_Only", ]
sg  <- slope_summary[slope_summary$Model == "VAST_GMRF", ]
tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\end{table}", "",
  "\\begin{table}[ht]", "\\centering", "\\small",
  "\\caption{Partitioning and unconditional prediction in the shared-field sensitivity analysis (50 replicates per level). Slope on $z$ is the mean slope from regressing the SGLMM's estimated field on $z$. Mean $\\hat{\\delta}$ and SD $\\hat{\\delta}$ are the mean and standard deviation across replicates of the indirect effect on the true covariate (RSR-SGLMM slope minus SGLMM slope). The final columns give mean unconditional RMSE.}",
  "\\label{tab:si_shared_partition}", "\\begin{tabular}{lrrrrrr}", "\\toprule",
  " & & & & & \\multicolumn{2}{c}{Unconditional RMSE} \\\\", "\\cmidrule(lr){6-7}",
  "$\\lambda$ & $r$ & Slope on $z$ & Mean $\\hat{\\delta}$ & SD $\\hat{\\delta}$ & SGLMM & GLM/RSR-SGLMM \\\\", "\\midrule")
for (j in seq_len(nrow(delta_summary))) {
  tex <- c(tex, paste0(delta_summary$Xreal_to_Z[j], " & ", num(glm$mean_corr[j]), " & ", num(sg$slope_field[j]), " & ",
    num(delta_summary$mean_delta[j]), " & ", num(delta_summary$sd_delta[j]), " & ", num(sg$rmse_uncond[j], 3), " & ",
    num(glm$rmse_uncond[j], 3), " \\\\"))
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex, "results/Tables_Simulation2_confounder.tex")

# Figure
model_colors <- c("Covariate_Only" = "#440154FF", "VAST_GMRF" = "#31688EFF", "VAST_RSR" = "#35B779FF")

plot_dat <- true_cov %>%
  left_join(
    glm %>%
      transmute(Xreal_to_Z, Confounding_Label = factor(sprintf("r = %.2f", mean_corr), levels = sprintf("r = %.2f", mean_corr))),
    by = "Xreal_to_Z"
  )

shared_theme <- theme_bw() +
  theme(
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 13, margin = margin(t = 15)),
    axis.text.x  = element_text(size = 11),
    axis.title.y = element_text(size = 13, margin = margin(r = 15)),
    axis.text.y  = element_text(size = 12),
    legend.position = "none",
    panel.grid   = element_blank()
  )

panel <- function(y, title, ylab, zero_line = TRUE) {
  p <- ggplot(plot_dat, aes(x = Confounding_Label, y = {{ y }}, fill = Model))
  if (zero_line) p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.7)
  p +
    geom_boxplot(position = position_dodge(0.8), width = 0.7, alpha = 0.8, outlier.size = 0.3) +
    stat_summary(aes(group = Model), fun = mean, geom = "point", shape = 23, size = 2,
                 fill = "white", position = position_dodge(0.8)) +
    scale_fill_manual(values = model_colors, labels = lab) +
    labs(title = title, x = "Confounding strength", y = ylab) +
    shared_theme
}

p_direct <- panel(Estimate - True_Effect, "Coefficient estimate error relative to the direct effect", expression(hat(beta) - 1))
p_total  <- panel(Estimate - Total_Effect, "Coefficient estimate error relative to the total effect", expression(hat(beta) - beta[total]))
p_uncond <- panel(RMSE_latent, "Unconditional prediction", "RMSE", zero_line = FALSE)

legend_plot <- ggplot(plot_dat, aes(x = Model, y = RMSE_latent, fill = Model)) +
  geom_boxplot() +
  scale_fill_manual(values = model_colors, labels = lab) +
  labs(fill = "Model") +
  theme_bw() +
  theme(legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
        legend.position = "bottom", legend.direction = "horizontal")

fig <- plot_grid(
  p_direct, p_total, p_uncond, ggpubr::get_legend(legend_plot),
  labels = c("A)", "B)", "C)"), label_size = 16, ncol = 1,
  rel_heights = c(0.3, 0.3, 0.3, 0.07)
)

ggsave("figures/Figure2_Simulation2.png", fig, bg = "white", width = 14, height = 11, units = "in", dpi = 300)
