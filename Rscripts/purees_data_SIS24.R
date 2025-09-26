# Set wd ------------------------------------------------------------------

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Load libraries ----------------------------------------------------------

library("tidyverse")
library("ACutils") # devtools::install_github("https://github.com/alessandrocolombi/ACutils")
library("mvtnorm")
library("salso")
library("FGM") # devtools::install_github("SmearyTundra/FGM")
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
  theta         = 1,
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
niter   <- 10000#100000
burn_in <-  1000#10000

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
  update_weights         = TRUE,
  update_partition       = TRUE,
  update_graph           = TRUE,
  perform_shuffle        = TRUE
)


## 1. Plot smoothed curves -------------------------------------------------
# Compute the mean of Beta in order to have data_post
sum_Beta <- matrix(0, p, n)
for(i in (burn_in+1):niter){
  sum_Beta <- sum_Beta + chains$Beta[[i]]
}
mean_Beta <- sum_Beta/(niter-burn_in)
data_post <- BaseMat %*% mean_Beta

# Compute the x value, create the basis and the functional object
x <- seq(0, 1, length.out=r)
basis <- create.bspline.basis(rangeval=range(x), nbasis=p, norder=3) 
data.fd <- Data2fd(y = data_post, argvals = x, basisobj = basis)

# Plot smoothed curves (1st one only)
plot.fd(data.fd[1,], main="smoothed curves", ylim=c(-1,3))
# plot(x, data_post[,1], type='l', ylim=c(-2,4))
lines(x,data[1,], main="smoothed curves", col='red')



## 2. Plot of the final Beta matrix ----------------------------------------

ACutils::ACheatmap(
  chains$Beta[[niter]],  
  use_x11_device = F,
  horizontal = F,
  main = "Estimated Beta matrix",
  center_value = NULL,
  col.upper = "darkred",
  col.center = "orange",
  col.lower = "lightyellow"
)



## 3. Traceplots (tau_eps, mu, beta) ---------------------------------------

# tau_eps
tau_plot <- as.vector(chains$tau_eps)
tau_plot <- tau_plot[(burn_in+1):niter]
plot(ts(mcmc(tau_plot)))


# mu
mu_plot <- matrix(0, niter, p)
for(i in 1:niter){
  mu_plot[i, ] <- chains$mu[[i]]
}
par(mfrow=c(1,3))
plot(ts(mcmc(mu_plot[, 1])), ylim=c(0.2,0.5))
plot(ts(mcmc(mu_plot[, 2])), ylim=c(0.2,0.5))
plot(ts(mcmc(mu_plot[, 3])), ylim=c(0.2,0.5))


# first element of first beta
beta1_plot <- rep(0, niter)
for(i in 1:niter){
  beta1_plot[i] <- chains$Beta[[i]][1,1]
}
# and 10th element of first beta
beta10_plot <- rep(0, niter)
for(i in 1:niter){
  beta10_plot[i] <- chains$Beta[[i]][10,1]
}

par(mfrow=c(1,2))
plot(ts(mcmc(beta1_plot)))
plot(ts(mcmc(beta10_plot)))


# 4. PYP parameters (sigma/theta) -----------------------------------------
# dsig = density(chains$sigma[9000:10000])
# plot(dsig$x, dsig$y)

plot(chains$theta, type = "l")
abline(h = 1, lty = 2, col = "red")

dtheta = density(chains$theta)
plot(dtheta$x, dtheta$y, type = "l")

# Posterior analysis ------------------------------------------------------
## Recomputing the partition in other forms and the number of groups
rho <- chains$rho
r = do.call(rbind, lapply(chains$rho, rho_to_r))
z = do.call(rbind, lapply(chains$rho, rho_to_z))
num_clusters = do.call(rbind, lapply(chains$rho, length))
num_clusters = as.vector(num_clusters)



## Barplot of changepoints -------------------------------------------------
bar_heights = colSums(r)
barplot(
  bar_heights,
  names = seq_along(bar_heights),
  border = "NA",
  space = 0,
  yaxt = "n",
  main="Changepoint frequency distribution",
  #col = color,
  cex.names=.6,
  las=2
)



## Evolution of the number of clusters -------------------------------------
par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
plot(
  x = seq_along(num_clusters),
  y = num_clusters,
  type = "n",
  xlab = "Iterations",
  ylab = "Number of groups",
  main = "Number of groups - Traceplot"
)
lines(x = seq_along(num_clusters), y = num_clusters)

par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
barplot(prop.table(table(num_clusters)), main = "Number of groups")


## Best partition ----------------------------------------------------------

sim_matrix <- salso::psm(z)
heatmap(sim_matrix)
image(sim_matrix)
est_part = minbinder(sim_matrix)$cl #binder
#est_part = minVI(sim_matrix)$cl #VI
est_part
rho_est_part = z_to_rho(est_part)
rho_est_part
prima = (data.fd$basis$rangeval[1]+data.fd$basis$params[1])/2
ultima = (data.fd$basis$rangeval[2]+tail(data.fd$basis$params,1))/2
nodi = c(prima,data.fd$basis$params,ultima)

# Plot smoothed curves
par(mfrow = c(1,1), mar = c(2,2,2,1), bty = "l")
plot.fd(data.fd[1:n,], main="smoothed curves", ylim=c(-1,3))
# abline(v = nodi, lwd = 2, lty = 2)
abline(v = nodi[cumsum(rho_est_part)], lwd = 2, lty = 2, col = "red")
abline(h=0,col = "white")


### Retrieving best partition using Binder only on visited ones (order is guaranteed here)  
# dists <- VI_LB(z, psm_mat = sim_matrix)
# best_partition_index = which.min(dists)
# rho_est = rho[[best_partition_index]]
# z_est = z[best_partition_index,]
# 
# # VI loss
# dists[best_partition_index]  
# 
# # select best partition
# unname(z_est)



# Graph analysis ----------------------------------------------------------
# Extract last plinks
last_plinks = tail(chains$G, n=1)[[1]]

# Criterion 1 to select the threshold (should not work very well) and assign final graph
threshold = 0.5
G_est <- matrix(0,p,p)
G_est[which(last_plinks>threshold)] = 1

#Criterion 2 to select the threshold
bfdr_select = BFDR_selection(last_plinks, tol = seq(0.1, 1, by = 0.001))

# Inspect the threshold and assign final graph
bfdr_select$best_treshold
G_est = bfdr_select$best_truncated_graph


### Plot estimated matrices
ACutils::ACheatmap(
  last_plinks,
  use_x11_device = F,
  horizontal = F,
  main = "Estimated plinks matrix",
  center_value = NULL,
  col.upper = "black",
  col.center = "grey50",
  col.lower = "white"
)


ACutils::ACheatmap(
  G_est+diag(p),
  use_x11_device = F,
  horizontal = F,
  main = "Estimated Graph",
  center_value = NULL,
  col.upper = "black",
  col.center = "grey50",
  col.lower = "white"
)


ACutils::ACheatmap(
  tail(chains$K,n=1)[[1]]-diag(diag(tail(chains$K,n=1)[[1]])),
  use_x11_device = F,
  horizontal = F,
  main = "Estimated Precision matrix",
  center_value = 0,
  col.upper = "black",
  col.center = "grey50",
  col.lower = "white"
)

