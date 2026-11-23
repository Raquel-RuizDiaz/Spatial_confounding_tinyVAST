# Simulation 3 (spatial latent variable as mediator)
# the two LaTeX tables. Reads the summaries saved by 03_Simulation3_mediator.R
# Writes results/Tables_Simulation3_mediator.tex (requires \usepackage{booktabs}).

estimand  <- read.csv("results/Simulation3_mediator_estimand.csv")
partition <- read.csv("results/Simulation3_mediator_variance_partition.csv")

num <- function(x, k = 2) {
  s <- sprintf(paste0("%.", k, "f"), x)
  paste0("$", sub("^-(0\\.0+)$", "\\1", s), "$")
}
lab <- c(Covariate_Only = "GLM", VAST_GMRF = "SGLMM", VAST_RSR = "RSR-SGLMM", VAST_Selected = "SGLMM-selected")

tex <- c(
  "\\begin{table}[ht]", "\\centering", "\\small",
  "\\caption{Slope estimates for the true covariate in \\textit{Simulation 3} (spatial latent variable as mediator; 50 replicates per level). $\\gamma$ is the strength with which the latent variable includes the covariate ($\\omega^*_s = \\omega_s + \\gamma z_s$), and $r$ is the mean realized correlation between $z_s$ and $\\omega^*_s$. Mean $\\hat{\\beta}$ and SD are the mean and standard deviation of the estimated slope across replicates. Bias is the mean of $\\hat{\\beta}$ minus the target, and RMSE is the root mean squared difference between $\\hat{\\beta}$ and the target, where the target is either the direct effect ($\\beta = 1$) or the total effect ($1 + \\gamma$).}",
  "\\label{tab:sim3_estimand}", "\\begin{tabular}{llrrrrrrr}", "\\toprule",
  " & & & & & \\multicolumn{2}{c}{Direct effect} & \\multicolumn{2}{c}{Total effect} \\\\",
  "\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}",
  "Estimator & $\\gamma$ & $r$ & Mean $\\hat{\\beta}$ & SD & Bias & RMSE & Bias & RMSE \\\\", "\\midrule"
)
for (m in names(lab)) {
  x <- estimand[estimand$Model == m, ]
  for (j in seq_len(nrow(x))) {
    tex <- c(tex, paste0(if (j == 1) lab[[m]] else "", " & ", x$Xreal_to_Z[j], " & ", num(x$mean_corr[j]), " & ",
      num(x$mean_beta[j]), " & ", num(x$sd_beta[j]), " & ", num(x$bias_direct[j]), " & ", num(x$rmse_direct[j]), " & ",
      num(x$bias_total[j]), " & ", num(x$rmse_total[j]), " \\\\"))
  }
  if (m != "VAST_Selected") tex <- c(tex, "\\addlinespace")
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\end{table}", "",
  "\\begin{table}[ht]", "\\centering", "\\small",
  "\\caption{Diagnostics of the estimated spatial latent variable in \\textit{Simulation 3} (means over 50 replicates per level). $\\mathbf{A}\\hat{\\boldsymbol{\\omega}}$ is the estimated spatial latent variable at the sampled locations, $\\omega$ is the innovation of the latent variable (independent of $z$), and $\\gamma z$ is the component added to the latent variable by the true covariate. Slope on $z$ is the slope from regressing $\\mathbf{A}\\hat{\\boldsymbol{\\omega}}$ on $z$, which indicates how much of $\\gamma z$ the estimated latent variable absorbed. $\\mathrm{var}(\\cdot)$ is the variance across sampled locations, and $\\mathrm{cor}(\\mathbf{A}\\hat{\\boldsymbol{\\omega}}, \\omega)$ is the Pearson correlation across sampled locations, which indicates how closely the estimated latent variable matches $\\omega$.}",
  "\\label{tab:sim3_field}", "\\begin{tabular}{llrrrrr}", "\\toprule",
  "Estimator & $\\gamma$ & Slope on $z$ & $\\mathrm{var}(\\mathbf{A}\\hat{\\boldsymbol{\\omega}})$ & $\\mathrm{var}(\\omega)$ & $\\mathrm{var}(\\gamma z)$ & $\\mathrm{cor}(\\mathbf{A}\\hat{\\boldsymbol{\\omega}}, \\omega)$ \\\\", "\\midrule")
for (m in c("VAST_GMRF", "VAST_RSR")) {
  x <- partition[partition$Model == m, ]
  for (j in seq_len(nrow(x))) {
    tex <- c(tex, paste0(if (j == 1) lab[[m]] else "", " & ", x$Xreal_to_Z[j], " & ", num(x$mean_field_slope_on_X[j]), " & ",
      num(x$mean_var_field_hat[j], 3), " & ", num(x$mean_var_omega0[j], 3), " & ", num(x$mean_var_gammaX[j], 3), " & ",
      num(x$mean_cor_field_omega0[j]), " \\\\"))
  }
  if (m == "VAST_GMRF") tex <- c(tex, "\\addlinespace")
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex, "results/Tables_Simulation3_mediator.tex")
