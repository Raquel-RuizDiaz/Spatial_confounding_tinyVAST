# Simulation 3 (spatial latent variable as mediator)
# Figure 3. Reads results/Simulation3_mediator_results.csv from 03_Simulation3_mediator.R

suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(patchwork)})
res <- read.csv("results/Simulation3_mediator_results.csv")
rep_levels <- c(0, 0.15, 0.25, 0.5, 0.75, 1, 4)
res <- res %>% filter(Covariate == "real", Xreal_to_Z %in% rep_levels)

model_colors <- c(Covariate_Only = "#440154FF", VAST_GMRF = "#31688EFF", VAST_RSR = "#35B779FF", VAST_Selected = "gold")
model_labels <- c(Covariate_Only = "GLM", VAST_GMRF = "SGLMM", VAST_RSR = "RSR-SGLMM", VAST_Selected = "SGLMM-selected")

lab <- res %>% filter(Model == "Covariate_Only") %>% group_by(Xreal_to_Z) %>%
  summarise(r = mean(Realized_Corr), .groups = "drop") %>% arrange(Xreal_to_Z) %>%
  mutate(Label = factor(sprintf("r = %.2f", r), levels = sprintf("r = %.2f", r)),
         block = factor(ifelse(Xreal_to_Z < 4, "low", "high"), levels = c("low", "high")))
res <- res %>% left_join(lab, by = "Xreal_to_Z") %>% mutate(Model = factor(Model, levels = names(model_colors)))

shared_theme <- theme_bw() + theme(
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
  axis.title.x = element_text(size = 13, margin = margin(t = 10)),
  axis.title.y = element_text(size = 13, margin = margin(r = 10)),
  axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 12),
  panel.grid = element_blank(), legend.position = "none")

# Reference marks: direct effect (beta = 1) and total effect (1 + gamma) for each level
refs <- lab %>% select(Label, block, Xreal_to_Z) %>%
  tidyr::crossing(Effect = c("Direct effect (\u03b2 = 1)", "Total effect (1 + \u03b3)")) %>%
  mutate(y = ifelse(grepl("Direct", Effect), 1, 1 + Xreal_to_Z), Effect = factor(Effect))

mkA <- function(d, rf, ylab = TRUE) {
  ggplot(d, aes(x = Label, y = Estimate)) +
    geom_errorbar(data = rf, aes(x = Label, ymin = y, ymax = y, linetype = Effect), inherit.aes = FALSE,
                  width = 0.95, linewidth = 0.7, colour = "black") +
    geom_boxplot(aes(fill = Model), position = position_dodge(0.8), width = 0.7, alpha = 0.8, outlier.size = 0.3) +
    stat_summary(aes(group = Model), fun = mean, geom = "point", shape = 23, size = 2, fill = "white",
                 position = position_dodge(0.8)) +
    scale_fill_manual(values = model_colors, guide = "none") +
    scale_linetype_manual(values = c("solid", "dashed"), name = NULL, drop = FALSE) +
    labs(x = NULL, y = if (ylab) expression(hat(beta)) else NULL) +
    shared_theme + theme(legend.text = element_text(size = 12), legend.key.width = unit(1.5, "cm"))
}
# Split axis: gamma <= 1 on the left, gamma = 4 on the right (separate y scales)
lo <- res %>% filter(block == "low") %>% droplevels()
hi <- res %>% filter(block == "high") %>% droplevels()
pA <- (mkA(lo, refs %>% filter(block == "low") %>% mutate(Label = droplevels(Label))) +
       mkA(hi, refs %>% filter(block == "high") %>% mutate(Label = droplevels(Label)), ylab = FALSE) +
       plot_layout(widths = c(6, 1.15), guides = "collect")) & theme(legend.position = "top")
pA <- wrap_elements(full = pA + plot_annotation(
  title = "Coefficient estimates relative to the direct and total effects",
  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))))

pB <- ggplot(res, aes(x = Label, y = Estimate - Total_Effect, fill = Model)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.7) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7, alpha = 0.8, outlier.size = 0.3) +
  stat_summary(aes(group = Model), fun = mean, geom = "point", shape = 23, size = 2, fill = "white", position = position_dodge(0.8)) +
  scale_fill_manual(values = model_colors, labels = model_labels) +
  labs(title = "Coefficient estimate error relative to the total effect", x = NULL, y = expression(hat(beta) - (1 + gamma))) +
  shared_theme

pC <- ggplot(res, aes(x = Label, y = RMSE_latent, fill = Model)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7, alpha = 0.8, outlier.size = 0.3) +
  stat_summary(aes(group = Model), fun = mean, geom = "point", shape = 23, size = 2, fill = "white", position = position_dodge(0.8)) +
  scale_fill_manual(values = model_colors, labels = model_labels, name = "Model") +
  labs(title = "Unconditional prediction", x = "Confounding strength", y = "RMSE") +
  shared_theme + theme(legend.position = "bottom", legend.title = element_text(size = 12, face = "bold"),
                       legend.text = element_text(size = 11))

fig <- wrap_plots(pA, pB, pC, ncol = 1, heights = c(1.1, 1, 1)) +
  plot_annotation(tag_levels = list(c("A)", "B)", "C)"))) & theme(plot.tag = element_text(size = 16, face = "bold"))
ggsave("figures/Figure3_Simulation3.png", fig, width = 14, height = 11.5, units = "in", dpi = 300, bg = "white")
