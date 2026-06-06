
library(RTMB)
library(Matrix)

# Simulate GMRF from precision ... from textbook
rmvnorm_prec <-
function( mu, # estimated fixed and random effects
          prec, # estimated joint precision
          n.sims) {

  require(Matrix)
  # Simulate values
  z0 = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  # Q = t(P) * L * t(L) * P
  L = Cholesky(prec, super=TRUE)
  # Calcualte t(P) * solve(t(L)) * z0 in two steps
  z = solve(L, z0, system = "Lt") # z = Lt^-1 * z
  z = solve(L, z, system = "Pt") # z = Pt    * z
  return(mu + as.matrix(z))
}

###################
# Simulate latent variable z_s and covariates x_s
# data y_s = z_s + error
###################

set.seed(101)

# Create grid of locations
xy_i = expand.grid( x=seq_len(10), y=seq_len(10) )

# Create adjacency matrix
D_ss = dist( xy_i, diag=TRUE, upper=TRUE )
A_ss = Matrix(ifelse( as.matrix(D_ss)==1, 1, 0 ))

# Make precision matrix
rho = 0.95
I_ss = Diagonal( n=nrow(A_ss) )
Q_ss = t(I_ss - rho*A_ss) %*% (I_ss - rho*A_ss)

# Simulate latent variable
z_s = rmvnorm_prec( mu=rep(0,nrow(A_ss)), prec = Q_ss, n = 1 )

# Irrelevant covariates
x_sj = rmvnorm_prec( mu=rep(0,nrow(A_ss)), prec = Q_ss, n = 3 )
x_sj = cbind( 1, x_sj )

# Simulate data
yhat_s = z_s
y_s = rnorm( n=length(yhat_s), mean=yhat_s, sd=1 )

###################
# Fit wrong linear model, y_s ~ x_s
# Tends to think one or more covariates is significant, anti-conservative tests 
###################

Lm = lm( y_s ~ 0 + x_sj )
summary(Lm)

#################
# Fit standard SAR
# Tends to correctly identify insignificant covariates, lower slopes, higher SEs
# BUT: adjusts slopes too
#################

# SAR
jnll_sar = function(parlist){
  getAll(parlist)                   
  Q = exp(ln_tau) * t(I_ss - exp(ln_rho)*A_ss) %*% (I_ss - exp(ln_rho)*A_ss)
  nll = -1 * dgmrf( omega_s, Q = Q, log=TRUE )
  nll = nll - sum(dnorm(y_s, mean=x_sj %*% betas + omega_s, sd=exp(ln_sd), log=TRUE))
  return(nll)
}

# 
parlist = list(
  omega_s = rep(0, nrow(A_ss)),
  ln_tau = 0,
  ln_rho = 0,
  ln_sd = 0,
  betas = rep(0, ncol(x_sj))
)

# Estimate SAR
obj_sar = RTMB::MakeADFun( jnll_sar, parlist, random = "omega_s" )
opt_sar = nlminb( obj_sar$par, obj_sar$fn, obj_sar$gr )
sdrep_sar = sdreport( obj_sar )

#################
# Fit Restricted Spatial Regression (RSR)
#################

# RSR
jnll_rsr = function(parlist){
  getAll(parlist)                   
  P_ss = x_sj %*% solve(t(x_sj) %*% x_sj) %*% t(x_sj)
  #P_ss = Matrix(0, nrow=nrow(A_ss), ncol=nrow(A_ss))
  Q = exp(ln_tau) * t(I_ss - exp(ln_rho)*A_ss) %*% (I_ss - exp(ln_rho)*A_ss)
  nll = -1 * dgmrf( omega_s, Q = Q, log=TRUE )
  omegaprime_s = (diag(nrow(A_ss)) - P_ss) %*% omega_s
  nll = nll - sum(dnorm(y_s, mean=x_sj %*% betas + omegaprime_s, sd=exp(ln_sd), log=TRUE))
  return(nll)
}

# Estimate SAR
obj_rsr = RTMB::MakeADFun( jnll_rsr, parlist, random = "omega_s" )
opt_rsr = nlminb( obj_rsr$par, obj_rsr$fn, obj_rsr$gr )
sdrep_rsr = sdreport( obj_rsr )

##################
# Compare betas with slopes from linear model
##################

summary(Lm)
summary(sdrep_sar, "fixed")
summary(sdrep_rsr, "fixed")

###################
# visualization
###################

# 1. Create spatial visualization of the latent field
plot_latent <- ggplot(data.frame(x = xy_i$x, 
                                 y = xy_i$y, 
                                 z = as.vector(z_s)), 
                      aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_viridis(option = "D", name = "Value") +
  theme_minimal() +
  labs(title = "A) Latent Spatial Field",
       x = "X coordinate",
       y = "Y coordinate") +
  theme(text = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))

# 2. Create spatial visualization of the observed data
plot_observed <- ggplot(data.frame(x = xy_i$x, 
                                   y = xy_i$y, 
                                   z = as.vector(y_s)), 
                        aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_viridis(option = "D", name = "Value") +
  theme_minimal() +
  labs(title = "B) Observed Data",
       x = "X coordinate",
       y = "Y coordinate") +
  theme(text = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"))

# 3. Create coefficient comparison plot
# Extract coefficients and standard errors correctly
lm_coef <- data.frame(
  Model = "LM",
  Coefficient = paste0("β", 0:3),
  Estimate = coef(Lm),
  SE = summary(Lm)$coefficients[,2]
)

sar_coef <- data.frame(
  Model = "SAR",
  Coefficient = paste0("β", 0:3),
  Estimate = summary(sdrep_sar, "fixed")[4:7,1],
  SE = summary(sdrep_sar, "fixed")[4:7,2]
)

rsr_coef <- data.frame(
  Model = "RSR",
  Coefficient = paste0("β", 0:3),
  Estimate = summary(sdrep_rsr, "fixed")[4:7,1],
  SE = summary(sdrep_rsr, "fixed")[4:7,2]
)

# Combine all coefficients
coef_data <- rbind(lm_coef, sar_coef, rsr_coef)

# Add confidence intervals
coef_data <- coef_data %>%
  mutate(
    lower = Estimate - 1.96 * SE,
    upper = Estimate + 1.96 * SE
  )

plot_coef <- ggplot(coef_data, aes(x = Coefficient, y = Estimate, color = Model)) +
  geom_point(position = position_dodge(width = 0.4), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), 
                position = position_dodge(width = 0.4),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  scale_color_viridis_d(option = "D") +
  labs(title = "C) Coefficient Estimates with 95% CI",
       x = "Coefficient",
       y = "Estimate") +
  theme(text = element_text(size = 12),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom")

# Arrange plots in a grid
grid.arrange(plot_latent, plot_observed, plot_coef,
             layout_matrix = rbind(c(1,2),
                                   c(3,3)),
             widths = c(1,1),
             heights = c(1,1))

# Print numerical results for verification
print("Linear Model Coefficients:")
print(summary(Lm)$coefficients)
print("\nSAR Model Coefficients:")
print(summary(sdrep_sar, "fixed")[4:7,])
print("\nRSR Model Coefficients:")
print(summary(sdrep_rsr, "fixed")[4:7,])

