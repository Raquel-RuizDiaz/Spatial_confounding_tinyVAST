

library(RTMB)
library(Matrix)
# library(tweedie)  # loading conflicts with RTMB
#library(mgcv)

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
# data y_s = z_s + error
###################

#set.seed(101)
rho_s = 0.2 #  Adjacency = 4 cells
rho_t = 0.45  # Adjacency = 2 times
n_fit = 3
n_forecast = 2
n_x = n_y = 10
n_j = 1 + 1  # First is the intercept
beta_j = c( 2, rep(0,n_j-1) ) # Response for intercept and n_j covariates
phi = 1
power = 1.5
random = c("omega_st","beta_j")[1:2]   # opt_sar$obj = opt_rsr$obj  ~only if~ using REML (i.e., slopes are random)

# Derived
n_t = n_fit + n_forecast
n_s = n_x * n_y
names_j = paste0( "x", seq_len(n_j)-1 )

# Create spatial precision
xy_i = expand.grid( x=seq_len(n_x), y=seq_len(n_y) )
D_ss = dist( xy_i, diag=TRUE, upper=TRUE )
A_ss = Matrix(ifelse( as.matrix(D_ss)==1, 1, 0 ))
I_ss = Diagonal( n=nrow(A_ss) )
IminusP_ss = I_ss - rho_s*A_ss
IminusP_ss = rescale( IminusP_ss)
Q_ss = t(IminusP_ss) %*% IminusP_ss

#
D_tt = dist( seq_len(n_t), diag=TRUE, upper=TRUE )
A_tt = Matrix(ifelse( as.matrix(D_tt)==1, 1, 0 ))
I_tt = Diagonal( n=nrow(A_tt) )
Q_tt = t(I_tt - rho_t*A_tt) %*% (I_tt - rho_t*A_tt)
IminusP_tt = I_tt - rho_t*A_tt
IminusP_tt = rescale( IminusP_tt )
Q_tt = t(IminusP_tt) %*% IminusP_tt

# Simulate latent variable
Q_kk = kronecker( Q_tt, Q_ss )
z_st = matrix( rmvnorm_prec( mu=rep(0,nrow(Q_kk)), prec = Q_kk, n = 1 ), ncol=n_t )

# Irrelevant covariates
x_stj = array( rmvnorm_prec( mu=rep(0,nrow(Q_kk)), prec = Q_kk, n = n_j-1 ), dim=c(n_s, n_t, n_j-1) )
x_stj = sweep( x_stj, MARGIN=2, FUN = "+", STAT = seq(-1,1,length=n_t) )
x_stj = abind::abind( array(1,dim=c(n_s,n_t,1)), x_stj )
apply(x_stj, MARGIN=2, FUN=mean )

#
x_kj = matrix(x_stj, ncol=n_j, dimnames=list(NULL,names_j))

# Simulate data
dhat_k = exp(as.vector(z_st) + (x_kj %*% beta_j ))
d_k = tweedie::rtweedie( n=length(dhat_k), mu=dhat_k[,1], phi = phi, power = power )

# Compile
dat = expand.grid( x = seq_len(n_x),
                   y = seq_len(n_y),
                   t = seq_len(n_t) )
dat = cbind( dat,
             dhat = dhat_k,
             d = d_k,
             x_kj )

# Fitted
dat1 = subset( dat, t <= n_fit )
#dat1 = dat
d1_k = dat1$d
x1_kj = as.matrix(dat1[,colnames(x_kj)])

###################
# Fit wrong linear model, y_s ~ x_s
# Tends to think one or more covariates is significant, anti-conservative tests
###################

library(tinyVAST)

#
form = as.formula(paste0("d ~ 0 + ", paste0( names_j, collapse=" + " )))
Lm1 = tinyVAST( form, data = dat1, family = tweedie() )

# Without covariates
Lm2 = tinyVAST( d ~ 1, data = dat1, family = tweedie() )

#
dat$lm1 = predict( Lm1, newdata = dat)
dat$lm2 = predict( Lm2, newdata = dat)

summary(Lm1, what="fixed")
summary(Lm2, what="fixed")

#################
# Fit standard SAR
# Tends to correctly identify insignificant covariates, lower slopes, higher SEs
# BUT: adjusts slopes too
#################

# SAR ... predict to all years
jnll_sar = function(parlist){
  getAll(parlist)
  # Precisions
  #Q_ss = t(I_ss - exp(ln_rho_s)*A_ss) %*% (I_ss - exp(ln_rho_s)*A_ss)
  #Q_tt = t(I_tmp - exp(ln_rho_t)*A_tt[seq_len(n_fit),seq_len(n_fit)]) %*% (I_tmp - exp(ln_rho_t)*A_tt[seq_len(n_fit),seq_len(n_fit)])
  IminusP_ss = Diagonal(n=n_s) - exp(ln_rho_s)*A_ss
  IminusP_ss = rescale( IminusP_ss)
  Q_ss = t(IminusP_ss) %*% IminusP_ss
  IminusP_tt = Diagonal(n=n_t) - exp(ln_rho_t)*A_tt
  IminusP_tt = rescale( IminusP_tt)
  Q_tt = t(IminusP_tt) %*% IminusP_tt
  # Densities
  dens_s <- function(x) dgmrf(x, Q = Q_ss, scale = exp(ln_tau), log=TRUE)
  dens_t <- function(x) dgmrf(x, Q = Q_tt, scale = 1, log=TRUE)
  nll = -1 * dseparable(dens_s, dens_t)(omega_st)
  # Data likelihood
  omega_k = as.vector( omega_st )
  dhat_k = exp(omega_k + (x_kj %*% beta_j))
  which_fit = which( dat$t <= n_fit )
  nll = nll - sum(dtweedie(d_k[which_fit], mu=dhat_k[which_fit], phi=exp(ln_phi), p=1 + plogis(finv_power), log=TRUE))
  REPORT( dhat_k )
  # Backtransform betas
  which_fit = which( dat$t <= n_fit )
  betaprime_j = beta_j + (solve(t(x_kj[which_fit,]) %*% x_kj[which_fit,]) %*% (t(x_kj[which_fit,]) %*% omega_k[which_fit]))
  REPORT( beta_j )
  REPORT( betaprime_j )
  ADREPORT( betaprime_j )
  return(nll)
}

#
parlist = list(
  omega_st = array(0, dim=c(n_s,n_t)),
  ln_tau = 0,
  ln_rho_s = log(0.2),
  ln_rho_t = log(0.2),
  ln_phi = 0,
  finv_power = 0,
  beta_j = rep(0, n_j)
)

# Estimate SAR
jnll_sar( parlist )
start_time = Sys.time()
obj_sar = MakeADFun( jnll_sar,
                     parlist,
                     random = random,
                     silent = TRUE )    # profile = "beta_j"
opt_sar = nlminb( obj_sar$par,
                  obj_sar$fn,
                  obj_sar$gr,
                  control = list(iter.max = 1e4, eval.max = 1e4, trace = 1) )
sdrep_sar = sdreport( obj_sar,
                      getReportCovariance = FALSE )
time_sar = Sys.time() - start_time
summary(sdrep_sar, "fixed", p.value = TRUE)

dat$sar = obj_sar$report()$dhat_k

#################
# Compare them
#################

#
m_t = tapply( dat$d,
                 INDEX = dat$t,
                 FUN = \(x)median((x)) )
mhat_t = tapply( dat$dhat,
                 INDEX = dat$t,
                 FUN = \(x)median((x)) )
m1_t = tapply( dat$lm1,
                 INDEX = dat$t,
                 FUN = \(x)median((x)) )
m2_t = tapply( dat$lm2,
                 INDEX = dat$t,
                 FUN = \(x)median((x)) )
m3_t = tapply( dat$sar,
                 INDEX = dat$t,
                 FUN = \(x)median((x)) )

matplot( cbind(m_t, mhat_t, m1_t, m2_t,m3_t),
         col = c("black","black","blue","blue","red"),
         type = "l",
         lty = c("solid","dotted","solid","dotted","solid"),
         log = "y" )



#################
# Fit Restricted Spatial Regression (RSR)
#################

# RSR
jnll_rsr = function(parlist){
  getAll(parlist)
  # Precisions
  #Q_ss = t(I_ss - exp(ln_rho_s)*A_ss) %*% (I_ss - exp(ln_rho_s)*A_ss)
  #Q_tt = t(I_tt - exp(ln_rho_t)*A_tt) %*% (I_tt - exp(ln_rho_t)*A_tt)
  IminusP_ss = Diagonal(n=n_s) - exp(ln_rho_s)*A_ss
  IminusP_ss = rescale( IminusP_ss)
  Q_ss = t(IminusP_ss) %*% IminusP_ss
  IminusP_tt = Diagonal(n=n_t) - exp(ln_rho_t)*A_tt
  IminusP_tt = rescale( IminusP_tt)
  Q_tt = t(IminusP_tt) %*% IminusP_tt
  # Densities
  dens_s <- function(x) dgmrf(x, Q = Q_ss, scale = exp(ln_tau), log=TRUE)
  dens_t <- function(x) dgmrf(x, Q = Q_tt, scale = 1, log=TRUE)
  nll = -1 * dseparable( dens_s, dens_t)(omega_st)
  # Rotation
  which_fit = which( dat$t <= n_fit )
  omega_k = as.vector( omega_st )
  #P_kk = x_kj %*% solve(t(x_kj) %*% x_kj) %*% t(x_kj)
  #omegaprime_k = (diag(nrow(P_kk)) - P_kk) %*% omega_s
  omegaprime_k = omega_k - ( x_kj %*% ( solve(t(x_kj[which_fit,]) %*% x_kj[which_fit,]) %*% ( t(x_kj[which_fit,]) %*% omega_k[which_fit] ) ) )
  # Data likelihood
  dhat_k = exp(omegaprime_k + (x_kj %*% beta_j))
  nll = nll - sum(dtweedie(d_k[which_fit], mu=dhat_k[which_fit], phi=exp(ln_phi), p=1 + plogis(finv_power), log=TRUE))
  REPORT( dhat_k )
  REPORT( omegaprime_k )
  REPORT( beta_j )
  return(nll)
}
jnll_rsr( parlist )

# Estimate SAR
start_time = Sys.time()
obj_rsr = MakeADFun( jnll_rsr,
                     parlist,
                     random = random,
                     silent = TRUE )
opt_rsr = nlminb( obj_rsr$par,
                  obj_rsr$fn,
                  obj_rsr$gr,
                  control = list(iter.max = 1e4, eval.max = 1e4, trace = 1) )
sdrep_rsr = sdreport( obj_rsr,
                      getReportCovariance = FALSE )
time_rsr = Sys.time() - start_time
#summary(opt_rsr, "fixed", p.value = TRUE)

# Check covariance
omegaprime_k = obj_rsr$report()$omegaprime_k
#cov( cbind(omegaprime_k,x_kj) )
t(omegaprime_k) %*% x_kj[which( dat$t <= n_fit ),]

##################
# Compare betas with slopes from linear model
##################

summary(Lm1, "fixed")
summary(sdrep_rsr, "random")
summary(sdrep_sar, "report")
