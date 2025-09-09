
library(RTMB)
library(Matrix)
# library(tweedie)  # loading conflicts with RTMB

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

rescale <-
function( IminusP ){
  invIminusP = solve(IminusP)     # solve(adsparse) returns dense-matrix

  # Hadamard squared LU-decomposition
  # See: https://eigen.tuxfamily.org/dox/group__QuickRefPage.html
  squared_invIminusP = AD(invIminusP)
  squared_invIminusP@x = squared_invIminusP@x^2
  targetvar_k = rep(1,nrow(IminusP))

  # Rescale Gamma
  margvar_k = solve(squared_invIminusP, targetvar_k)
  #D_kk = AD(Matrix(nrow=length(margvar_k), ncol=length(margvar_k)))
  D_kk = AD(Matrix(diag(rep(1,length(margvar_k)))))
  D_kk@x = 1/sqrt(margvar_k)
  out = D_kk %*% IminusP
  return(out)
}

###################
# Simulate latent variable z_s and covariates x_s
# data d_s = z_s + error
###################

#set.seed(101)
n_j = 1 + 3 # Intercept + slopes
names_j = paste0("x", seq_len(n_j)-1 )
measurement_error = c("normal", "lognormal", "tweedie", "poisson")[3]
method = c("LS", "WLS", "IWLS")[3]  # LS seems to work
random = c("omega_s","beta_j")[1:2]   # opt_sar$obj = opt_rsr$obj  ~only if~ using REML (i.e., slopes are random)
gmrf_type = c( "sar", "spde" )[2]

if( gmrf_type == "sar" ){
  # derived
  n_x = n_y = 10
  n_s = n_x * n_y

  # Create grid of locations
  xy_i = expand.grid( x=seq_len(n_x), y=seq_len(n_y) )

  # Create adjacency matrix
  D_ss = dist( xy_i, diag=TRUE, upper=TRUE )
  A_ss = Matrix(ifelse( as.matrix(D_ss)==1, 1, 0 ))

  # Make precision matrix
  rho_s = 0.20
  I_ss = Diagonal( n=nrow(A_ss) )
  IminusP_ss = I_ss - rho_s*A_ss
  IminusP_ss = rescale( IminusP_ss)
  Q_ss = t(IminusP_ss) %*% IminusP_ss
  A_is = Diagonal( n = nrow(Q_ss) )
}
if( gmrf_type == "spde" ){
  n_s = 50
  xy_i = pracma::poisson2disk( n=n_s )

  # Create mesh
  decorrelation_dist = 0.5
  SD = 1
  kappa = sqrt(8) / decorrelation_dist   # Range = sqrt(8) / kappa
  tau = 1 / kappa / SD / sqrt(4 * pi) # SD = 1 / sqrt(4 * pi * kappa^2 * tau^2)
  mesh = fmesher::fm_mesh_2d( xy_i )
  spde = fmesher::fm_fem( mesh )

  # make precision
  A_is = drop0(fmesher::fm_evaluator( mesh, loc=xy_i )$proj$A)
  Q_ss = (kappa^4*spde$c0 + 2*kappa^2*spde$g1 + spde$g2) * tau^2
}

# Simulate latent variable
z_s = as.matrix(A_is %*% rmvnorm_prec( mu=rep(0,nrow(Q_ss)), prec = Q_ss, n = 1 ))

# Irrelevant covariates
x_sj = as.matrix(A_is %*% rmvnorm_prec( mu=rep(0,nrow(Q_ss)), prec = Q_ss, n = n_j - 1 ))
x_sj = cbind( 1, x_sj )
colnames(x_sj) = names_j

# Simulate data
# f = inverse-link
# g = weighting term (variance as function of mean)
if( measurement_error=="normal" ){
  SD = 0.5
  f = identity
  g = \(x) x^0
  dhat_s = f(z_s)
  d_s = rnorm( n=length(dhat_s), mean=dhat_s, sd=SD )
}
if( measurement_error=="tweedie" ){
  phi = 1
  power = 1.5
  f = exp
  g = function(v) v^power     # ACTUALLY
  dhat_s = f(z_s)
  d_s = tweedie::rtweedie( n=length(dhat_s), mu=dhat_s[,1], phi = phi, power = power )
}
if( measurement_error=="lognormal" ){
  sdlog = 1
  f = exp
  g = identity
  dhat_s = f(z_s)
  d_s = rlnorm( n=length(dhat_s), meanlog=log(dhat_s[,1]), sdlog = 1 )
}
if( measurement_error=="poisson" ){
  f = exp
  g = identity
  dhat_s = f(z_s)
  d_s = rpois( n=length(dhat_s), lambda=dhat_s[,1] )
}

# Compile
dat = data.frame( xy_i,
             dhat = dhat_s,
             d = d_s,
             x_sj )

###################
# Fit wrong linear model, d_s ~ x_s
# Tends to think one or more covariates is significant, anti-conservative tests 
###################

devtools::install_github("vast-lib/tinyVAST@dev", force=TRUE)
library(tinyVAST)

#
if( measurement_error=="normal" ){
  form = as.formula(paste0("d ~ 0 + ", paste0( names_j, collapse=" + " )))
  #Lm1 = tinyVAST( form, data = dat, family = tweedie(link="log") )
  Lm1 = lm( form, data = dat )
  Coef = summary(Lm1)$coef
  family = gaussian()
}
if( measurement_error=="tweedie" ){
  form = as.formula(paste0("d ~ 0 + ", paste0( names_j, collapse=" + " )))
  #Lm1 = tinyVAST( form, data = dat, family = tweedie(link="log") )
  Lm1 = mgcv::gam( form, data = dat, family = mgcv::Tweedie(p = power) )
  Coef = summary(Lm1)$p.table
  family = tweedie()
}
if( measurement_error=="lognormal" ){
  form = as.formula(paste0("log(d) ~ 0 + ", paste0( names_j, collapse=" + " )))
  Lm1 = tinyVAST( form, data = dat, family = gaussian() )
  Coef = summary(Lm1)$coef
  family = gaussian()
}
if( measurement_error=="poisson" ){
  form = as.formula(paste0("d ~ 0 + ", paste0( names_j, collapse=" + " )))
  #Lm1 = tinyVAST( form, data = dat, family = poisson(link="log") )
  Lm1 = glm( form, data = dat, family = poisson(link="log") )
  Coef = summary(Lm1)$coef
  family = poisson()
}

#################
# Fit standard SAR
# Tends to correctly identify insignificant covariates, lower slopes, higher SEs
# BUT: adjusts slopes too
#
# betaprimeLS_j seems to work fine, but only gives identical results when using a linear mixed model (identity link and normal distribution)
# when using a identity link and normal distribution, betaprimeLS_j in obj_sar and beta_j in obj_rsr are identical only when using REML
#################

# SAR
jnll_sar = function(parlist){
  getAll(parlist)                   
  if( gmrf_type == "sar" ){
    IminusP_ss = Diagonal(n=n_s) - exp(ln_rho_s)*A_ss
    IminusP_ss = rescale( IminusP_ss)
    Q_ss = (t(IminusP_ss) %*% IminusP_ss) * exp(2*ln_tau)
  }
  if( gmrf_type == "spde" ){
    Q_ss = (exp(ln_kappa*4)*spde$c0 + 2*exp(ln_kappa*2)*spde$g1 + spde$g2) * exp(2*ln_tau)
  }
  nll = -1 * dgmrf( omega_s, Q = Q_ss, log=TRUE )
  omega_i = (A_is %*% omega_s)[,1]
  dhat_s = as.vector(f(omega_i + (x_sj %*% beta_j)))
  if( measurement_error=="normal" ){
    nll = nll - sum(dnorm(d_s, mean=dhat_s, sd=exp(ln_sd), log=TRUE))
  }
  if( measurement_error=="tweedie" ){
    nll = nll - sum(dtweedie(d_s, mu=dhat_s, phi=exp(ln_phi), p=1 + plogis(finv_power), log=TRUE))
  }
  if( measurement_error=="lognormal" ){
    nll = nll - sum(dlnorm( d_s, meanlog=log(dhat_s), sdlog=exp(ln_sdlog), log=TRUE))
  }
  if( measurement_error=="poisson" ){
    nll = nll - sum(dpois( d_s, lambda=dhat_s, log=TRUE))
  }
  # Backtransform betas
  W = AD(Diagonal( n = length(dhat_s), x = seq_along(dhat_s) ))
  W@x = g(dhat_s) #  * constant # When  W = constant * mu, constant drops out
  betaprimeLS_j = beta_j + (solve(t(x_sj) %*% x_sj) %*% (t(x_sj) %*% omega_i))
  betaprimeWLS_j = beta_j + ( solve(t(x_sj) %*% W %*% x_sj) %*% ( t(x_sj) %*% W %*% omega_i ) )[,1]
  REPORT( beta_j )
  REPORT( betaprimeLS_j )
  ADREPORT( betaprimeLS_j )
  REPORT( betaprimeWLS_j )
  ADREPORT( betaprimeWLS_j )
  return(nll)
}

# 
map = NULL
parlist = list(
  omega_s = rep(0, nrow(Q_ss)),
  ln_tau = 0,
  beta_j = rep(0, ncol(x_sj))
)
if( gmrf_type == "sar" ){
  parlist$ln_rho_s = log(0.2)
}
if( gmrf_type == "spde" ){
  parlist$ln_kappa = log( 1 / 0.2)
}
if( measurement_error=="normal" ){
  parlist$ln_sd = log(1)
}
if( measurement_error=="tweedie" ){
  parlist$ln_phi = 0
  parlist$finv_power = qlogis( power - 1 )
  #map$finv_power = factor(NA)
}
if( measurement_error=="lognormal" ){
  parlist$ln_sdlog = log(1)
}
# Estimate SAR
obj_sar = MakeADFun( jnll_sar,
                     parlist,
                     random = random,
                     #random = "omega_s",
                     #profile = "beta_j",
                     silent = TRUE,
                     map = map )
opt_sar = nlminb( obj_sar$par,
                  obj_sar$fn,
                  obj_sar$gr )
sdrep_sar = sdreport( obj_sar, bias.correct=TRUE )

# Compare
summary(sdrep_sar, "random", p.value = TRUE)
summary(sdrep_sar, "report", p.value = TRUE)
Coef

#################
# tinyVAST
#################

fit = tinyVAST(
  formula = form,
  data = cbind(dat,"var"="n"),
  family = family,
  space_term = "",
  variable_column = "var",
  space_columns = c("X1","X2"),
  spatial_domain = mesh,
  control = tinyVASTcontrol( reml = TRUE )
)

# fit$rep$alphaprime_j  ~should match~  betaprimeLS_j
cbind(
  "tinyVAST" = fit$rep$alphaprime_j,
  "bespoke_RTMB" = as.list(sdrep_sar, report=TRUE, what = "Estimate")$betaprimeLS_j
)

#################
# Fit Restricted Spatial Regression (RSR)
# Breaks sparsity in inner Hessian matrix
#################

# RSR
jnll_rsr = function(parlist){
  getAll(parlist)                   
  # Densities
  IminusP_ss = Diagonal(n=n_s) - exp(ln_rho_s)*A_ss
  IminusP_ss = rescale( IminusP_ss)
  Q_ss = t(IminusP_ss) %*% IminusP_ss
  nll = -1 * dgmrf( omega_s, Q = Q_ss, scale = exp(ln_tau), log=TRUE )
  # Rotation
  #P_ss = x_sj %*% solve(t(x_sj) %*% x_sj) %*% t(x_sj)
  #omegaprime_s = (diag(nrow(A_ss)) - P_ss) %*% omega_s
  if( method == "LS" ){
    omegaprime_s = omega_s - ( x_sj %*% ( solve(t(x_sj) %*% x_sj) %*% ( t(x_sj) %*% omega_s ) ) )
  }
  if( method == "WLS" ){   # Weighting based only on covariates
    mu_s = f( (x_sj %*% beta_j)[,1] )
    W = AD(Diagonal( n = length(mu_s), x = seq_along(mu_s) ))
    W@x = g(mu_s) #  * constant # When  W = constant * mu, constant drops out
    omegaprime_s = omega_s - ( x_sj %*% ( solve(t(x_sj) %*% W %*% x_sj) %*% ( t(x_sj) %*% W %*% omega_s ) ) )[,1]
  }
  if( method == "IWLS" ){    # Weighting based on covariates and rotated GMRF
    omegaprime_s = omega_s
    for(i in seq_len(3)){
      mu_s = f( omegaprime_s + (x_sj %*% beta_j)[,1] )
      W = AD(Diagonal( n = length(mu_s), x = seq_along(mu_s) ))
      W@x = g(mu_s) #  * constant # When  W = constant * mu, constant drops out
      omegaprime_s = omega_s - ( x_sj %*% ( solve(t(x_sj) %*% W %*% x_sj) %*% ( t(x_sj) %*% W %*% omega_s ) ) )[,1]
    }
  }
  # Data likelihood
  dhat_s = f(omegaprime_s + (x_sj %*% beta_j)[,1] )
  if( measurement_error=="tweedie" ){
    nll = nll - sum(dtweedie(d_s, mu=dhat_s, phi=exp(ln_phi), p=1 + plogis(finv_power), log=TRUE))
  }
  if( measurement_error=="lognormal" ){
    nll = nll - sum(dlnorm( d_s, meanlog=log(dhat_s), sdlog=exp(ln_sdlog), log=TRUE))
  }
  if( measurement_error=="poisson" ){
    nll = nll - sum(dpois( d_s, lambda=dhat_s, log=TRUE))
  }
  REPORT( beta_j )
  REPORT( omegaprime_s )
  return(nll)
}
jnll_rsr( parlist )

# Estimate SAR
obj_rsr = MakeADFun( jnll_rsr,
                     parlist,
                     random = random,
                     silent = TRUE,
                     map = map )
opt_rsr = nlminb( obj_rsr$par,
                  obj_rsr$fn,
                  obj_rsr$gr,
                  control = list(trace=1) )
sdrep_rsr = sdreport( obj_rsr )

#
H = obj_rsr$env$spHess(random=TRUE)
image(H)

# Check covariance
omegaprime_s = obj_rsr$report()$omegaprime_s
#cov( cbind(omegaprime_s,x_sj) )
t(omegaprime_s) %*% x_sj

##################
# Compare betas with slopes from linear model
##################

Coef
summary(sdrep_sar, "fixed")
summary(sdrep_rsr, "fixed")
