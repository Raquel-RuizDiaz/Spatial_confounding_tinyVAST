# When and how to use restricted spatial regression to separate environmental effects from spatial confounding

R code to reproduce the simulations, case study, figures and tables in:

> Ruiz Diaz, R. & Thorson, J. T. *When and how to use restricted spatial regression to separate environmental effects from spatial confounding.* Preprint: EcoEvoRxiv, https://doi.org/10.32942/X28351

The paper compares a spatial generalized linear mixed model (SGLMM) with the restricted spatial regression (RSR) estimator, as implemented in the [tinyVAST](https://cran.r-project.org/package=tinyVAST) R package. The SGLMM slope targets the direct effect of a covariate, and the RSR slope targets its total effect. Three simulations, in which an unmeasured spatial process is independent of the covariate, confounds its effect or mediates it, and a case study of three Bering Sea species show which estimator targets which effect.

## Repository structure

```
.
├── analysis/   R scripts, numbered in the order they are run
├── data/       Bering Sea survey data used in the case study
├── figures/    Figures as they appear in the paper
└── results/    Simulation and case-study outputs used in the paper
```

## Requirements

The results were produced with R 4.5.3, tinyVAST 1.6.3, fmesher 0.8.0, RTMB 2.0 and TMB 1.9.23. tinyVAST compiles C++ code, so a working compiler is needed (Rtools on Windows, the Xcode command line tools on macOS).

```r
install.packages(c(
  "tinyVAST", "fmesher", "RTMB", "Matrix", "pracma", "mvtnorm", "abind", "cv",
  "dplyr", "tidyr", "ggplot2", "cowplot", "ggpubr", "gridExtra", "patchwork",
  "viridis", "scales", "sf", "terra", "rnaturalearth"
))

# To use the same tinyVAST release as the paper:
# remotes::install_version("tinyVAST", version = "1.6.3")
```

### Outputs

Each script writes its outputs to `results/` and `figures/`, named after the script:

| Script | `results/` | `figures/` |
|---|---|---|
| `00` | | `FigureS1_concept.png` |
| `01` | `Simulation1_independent_results.csv`, `Simulation1_independent_summary.csv` | `Figure1_Simulation1.png` |
| `02` | `Simulation2_confounder_results.csv`, `Simulation2_confounder_estimand.csv`, `Simulation2_confounder_field.csv`, `Simulation2_confounder_indirect_effect.csv` | |
| `02b` | `Tables_Simulation2_confounder.tex` | `Figure2_Simulation2.png` |
| `03` | `Simulation3_mediator_results.csv`, `Simulation3_mediator_estimand.csv`, `Simulation3_mediator_variance_partition.csv`, `Simulation3_mediator_indirect_effect.csv` | |
| `03b` | | `Figure3_Simulation3.png` |
| `03c` | `Tables_Simulation3_mediator.tex` | |
| `04` | `Sensitivity_factorial_results.csv` | `FigureS2_factorial.png` |
| `05` | `EBS_case_study_coefficients.csv`, `EBS_case_study_curves.csv`, `EBS_case_study_coefficients.rds` | `Figure4_case_study.png` |

The `*_results.csv` files contain one row per model, covariate and simulation replicate. The other CSV files are summaries across replicates.

## Using the RSR estimator in your own analysis

The RSR estimator is computed after fitting, without changing the SGLMM. In tinyVAST, request it with `get_rsr = TRUE`:

```r
library(tinyVAST)

fit <- tinyVAST(
  formula        = y ~ x1 + x2,
  data           = dat,             # data frame with columns y, x1, x2, X, Y
  space_term     = "",              # spatial random effect with no covariance among variables
  space_columns  = c("X", "Y"),
  spatial_domain = fmesher::fm_mesh_2d(dat[, c("X", "Y")]),
  control        = tinyVASTcontrol(get_rsr = TRUE)
)

# SGLMM slopes (direct effect)
as.list(fit$sdrep, what = "Estimate")$alpha_j
# RSR slopes (total effect) and their standard errors
as.list(fit$sdrep, what = "Estimate",   report = TRUE)$alphaprime_j
as.list(fit$sdrep, what = "Std. Error", report = TRUE)$alphaprime_j
```

## Questions

Please open an issue in this repository.
