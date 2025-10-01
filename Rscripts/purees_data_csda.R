# Set wd ------------------------------------------------------------------

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Load libraries ----------------------------------------------------------

library("tidyverse")
library("ACutils") # devtools::install_github("https://github.com/alessandrocolombi/ACutils")
library("mvtnorm")
library("salso")
library("FGM") #  devtools::install_github("alessandrocolombi/FGMpackage")
library("gmp") # che fa?
library("mcclust")
library("mcclust.ext")
library("logr")
library("tidygraph")
library("ggraph")
library("igraph")
library("Rcpp")
library("RcppArmadillo")
library("RcppEigen")
library('RcppGSL')
library("fda")
library("coda")
library("lattice")

# Load custom functions ---------------------------------------------------

sourceCpp("../src/wade.cpp")
Rcpp::sourceCpp('../src/UpdateParamsGSL.cpp')
source("../R/utility_functions.R")
source("../R/bulky_functions.R")
source("../R/data_generation.R")



# Load data and plot ------------------------------------------------------
load("../data/purees.Rdat")
curves = purees$data
wavelengths = purees$wavelengths
strawberry <- curves[which(curves$Group == "Strawberry"), ]
data = strawberry[,-1]*100
data = as.matrix(data)    # n x r

#Generate basis
p =  40
n = dim(data)[1]
r = dim(data)[2]
range_x = range(wavelengths)

data_W <- t(as.matrix(data))
basis <- create.bspline.basis(rangeval=range(wavelengths), nbasis=p, norder=3) 
data.fd <- Data2fd(y = data_W, argvals = wavelengths, basisobj = basis)
plot.fd(data.fd, main="B-splines")

BaseMat <- eval.basis(wavelengths,basis) # matrix r x p 



# Initialization ----------------------------------------------------------
n = dim(data)[1]
r = dim(data)[2]
# Compute quantities for function UpdateParamGSL
tbase_base = t(BaseMat)%*%BaseMat  # p x p (phi_t * phi)
tbase_data = t(BaseMat)%*%t(data)  # p x n (phi_t * Y_t) 
Sdata = sum(diag(data%*%t(data))) 

set_UpdateParamsGSL_list = set_UpdateParamsGSL(
  tbase_base  = tbase_base,
  tbase_data  = tbase_data,
  Sdata       = Sdata,
  a_tau_eps   = 10,
  b_tau_eps   = 0.001,
  sigma_mu    = 100,
  r           = r,
  Update_Beta = TRUE,
  Update_Mu   = TRUE,
  Update_Tau  = TRUE
)


# Exp_theta = 50
# Var_theta =  1
# a_theta = Exp_theta*Exp_theta/Var_theta
# b_theta = Exp_theta/Var_theta
# rho0 = c(3,5,3,5,4)

# Set the initialization values of the chains
initialization_values = set_initialization(
  Beta          = matrix(rnorm(n=p*n), nrow = p, ncol = n),
  mu            = rnorm(n=p),
  tau_eps       = 100,
  K             = matrix(0,p,p),
  G             = matrix(0,p,p),
  z             = rep(1,p), 
  rho           = p,#rho0,
  a_sigma       = 1,
  b_sigma       = 1,
  c_sigma       = 1,
  d_sigma       = 1,
  c_theta       = 1,#a_theta,
  d_theta       = 1,#b_theta,
  sigma         = 0.0001,
  theta         = 2,
  weights_a0    = rep(1,p-1),
  weights_d0    = rep(1,p-1),
  total_weights = 0,
  total_K       = matrix(0,p,p),
  total_graphs  = matrix(0,p,p),
  graph_start   = NULL,
  graph_density = 0.8,
  beta_sig2     = 0.01,
  d             = 3
)



# Run ---------------------------------------------------------------------

# Set the number of iterations and burn-in
niter   <- 50000
burn_in <- 10000

chains = Gibbs_sampler_update(
  set_UpdateParamsGSL_list,
  niter, 
  initialization_values,
  alpha_target           = 0.234,
  alpha_add              = 0.5,
  adaptation_step        = 1/(1000*p),
  seed                   = 22111996,
  update_sigma_prior     = FALSE,
  update_theta_prior     = FALSE,
  update_weights         = FALSE,
  update_partition       = FALSE,
  update_graph           = TRUE,
  perform_shuffle        = FALSE
)

beepr::beep()

## Graph analysis ----------------------------------------------------------
# Extract last plinks
last_plinks = tail(chains$G, n=1)[[1]]


#Criterion 2 to select the threshold
bfdr_select = BFDR_selection(last_plinks, tol = seq(0.1, 1, by = 0.001))

# Inspect the threshold and assign final graph
bfdr_select$best_treshold
G_est = bfdr_select$best_truncated_graph


### Plot estimated matrices
ACutils::ACheatmap(
  last_plinks+diag(p),
  use_x11_device = F,
  horizontal = F,
  main = "Estimated plinks matrix",
  center_value = 0.84,
  col.upper = "darkred",
  col.center = "white",
  col.lower = "skyblue"
)

pal = hcl.colors(n = 3, palette = "Reds")
ACutils::ACheatmap(
  last_plinks+diag(p),
  use_x11_device = F,
  horizontal = F,
  main = "Estimated plinks matrix",
  center_value = 0.84,
  col.upper = "darkred",
  col.n_breaks = 21,
  col.center = "skyblue",
  col.lower = "white"
)


ACutils::ACheatmap(
  G_est+diag(p),
  use_x11_device = F,
  horizontal = F,
  main = "Estimated Graph",
  center_value = NULL,
  col.upper = "darkred",
  col.center = "grey50",
  col.lower = "white"
)


ACutils::ACheatmap(
  tail(chains$K,n=1)[[1]]-diag(diag(tail(chains$K,n=1)[[1]])),
  use_x11_device = F,
  horizontal = T,
  main = "Estimated Precision matrix",
  center_value = 0,
  col.upper = "darkred",
  col.center = "white",
  col.lower = "skyblue"
)



