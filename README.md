# Analysis scripts

Code for the simulations and case study in *When and how to use restricted spatial regression to
separate environmental effects from spatial confounding*. Run every script from the project root
(the folder that contains `analysis/`, `data/`, `figures/` and `results/`).

## Scripts

All scripts are in `analysis/`.

| Script | What it does | Manuscript output |
|---|---|---|
| `00_Conceptual_figure.R` | Conceptual illustration of the RSR decomposition | Figure S1 |
| `01_Simulation1_independent.R` | Simulation 1: spatial latent variable independent of the covariate | Figure 1 |
| `02_Simulation2_confounder.R` | Simulation 2: spatial latent variable as confounder. Runs the simulation and saves `results/Simulation2_confounder_*.csv` | — |
| `02b_Simulation2_confounder_figure_tables.R` | Figure and tables for Simulation 2, from the saved results | Figure 2, Tables 1–2 |
| `03_Simulation3_mediator.R` | Simulation 3: spatial latent variable as mediator. Runs the simulation and saves `results/Simulation3_mediator_results.csv` and summaries | — |
| `03b_Simulation3_mediator_figure.R` | Figure for Simulation 3, from the saved results | Figure 3 |
| `03c_Simulation3_mediator_tables.R` | Tables for Simulation 3, from the saved summaries | Tables 3–4 |
| `04_Sensitivity_factorial.R` | Sensitivity analysis across sample sizes and numbers of false covariates | Figure S2 |
| `05_EBS_case_study.R` | Eastern Bering Sea case study (uses `data/`) | Figure 4 |

`02_Simulation2_confounder.R` reuses the simulation function defined in `03_Simulation3_mediator.R`,
so keep both files in `analysis/`. Scripts ending in `b` or `c` only read saved results and run in
seconds; the others refit the models.

## Outputs

Figures are written to `figures/` under the names used in the manuscript:

| Manuscript | Script | File |
|---|---|---|
| Figure 1 | `01_Simulation1_independent.R` | `figures/Figure1_Simulation1.png` |
| Figure 2 | `02b_Simulation2_confounder_figure_tables.R` | `figures/Figure2_Simulation2.png` |
| Figure 3 | `03b_Simulation3_mediator_figure.R` | `figures/Figure3_Simulation3.png` |
| Figure 4 | `05_EBS_case_study.R` | `figures/Figure4_case_study.png` |
| Figure S1 | `00_Conceptual_figure.R` | `figures/FigureS1_concept.png` |
| Figure S2 | `04_Sensitivity_factorial.R` | `figures/FigureS2_factorial.png` |

Results are written to `results/`, named after the script that writes them:

| Script | Files in `results/` |
|---|---|
| `01_Simulation1_independent.R` | `Simulation1_independent_results.csv`, `Simulation1_independent_summary.csv` |
| `02_Simulation2_confounder.R` | `Simulation2_confounder_results.csv`, `Simulation2_confounder_estimand.csv`, `Simulation2_confounder_field.csv`, `Simulation2_confounder_indirect_effect.csv` |
| `02b_Simulation2_confounder_figure_tables.R` | `Tables_Simulation2_confounder.tex` (Tables 1–2) |
| `03_Simulation3_mediator.R` | `Simulation3_mediator_results.csv`, `Simulation3_mediator_estimand.csv`, `Simulation3_mediator_variance_partition.csv`, `Simulation3_mediator_indirect_effect.csv` |
| `03c_Simulation3_mediator_tables.R` | `Tables_Simulation3_mediator.tex` (Tables 3–4) |
| `04_Sensitivity_factorial.R` | `Sensitivity_factorial_results.csv` |
| `05_EBS_case_study.R` | `EBS_case_study_coefficients.csv`, `EBS_case_study_curves.csv`, `EBS_case_study_coefficients.rds` |

## Software

Results were produced with tinyVAST 1.6.3 (CRAN). In this version `rmvnorm_prec()` is internal, so
the simulation scripts define a small wrapper that calls `tinyVAST:::rmvnorm_prec(prec = Q, n = n)`.
Other packages: RTMB, fmesher, pracma, Matrix, mvtnorm, abind, cv, dplyr, tidyr, ggplot2, cowplot,
ggpubr, gridExtra, patchwork, viridis, sf, terra, rnaturalearth.
