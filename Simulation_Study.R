# wd ----------------------------------------------------------------------
wd_pc_ale = "C:/Users/colom/FGMRandomBlocks/"
wd_pc_luciano = "C://Users//lucia//Desktop//PhD//my collaborations//Alessandro Colombi work//FGMRandomBlocks"
wd_bocconi = "/home/colombi/FGMRandomBlocks/"
wd_vec = c(wd_pc_ale,wd_pc_luciano,wd_bocconi)
choose_wd = wd_vec[3] # <--- modify here to select the wd according to the user
wd = paste0(choose_wd,"./")
setwd(wd)

# Load libraries ----------------------------------------------------------

library("tidyverse")
#library("ACutils") # devtools::install_github("https://github.com/alessandrocolombi/ACutils")
library("mvtnorm")
library("FGM") #  devtools::install_github("alessandrocolombi/FGMpackage")
library("gmp") 
library("logr")
library("tidygraph")
library("ggraph")
library("igraph")
library("Rcpp")
library("RcppArmadillo")
library("RcppEigen")
library('RcppGSL')
# library("fda")
library("coda")
library("lattice")
library("parallel")

cat("R version: ", R.version.string, "\n", sep = "")
cat("FGM version: ", as.character(utils::packageVersion("FGM")), "\n", sep = "")
cat("FGM path: ", find.package("FGM"), "\n", sep = "")
cat(
  "Thread env: ",
  paste(
    c(
      "OMP_NUM_THREADS",
      "OMP_THREAD_LIMIT",
      "OPENBLAS_NUM_THREADS",
      "MKL_NUM_THREADS",
      "MKL_DOMAIN_NUM_THREADS",
      "RCPP_PARALLEL_NUM_THREADS"
    ),
    Sys.getenv(
      c(
        "OMP_NUM_THREADS",
        "OMP_THREAD_LIMIT",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "MKL_DOMAIN_NUM_THREADS",
        "RCPP_PARALLEL_NUM_THREADS"
      ),
      unset = ""
    ),
    sep = "=",
    collapse = ", "
  ),
  "\n",
  sep = ""
)

# Load custom functions ---------------------------------------------------

source("./utility_functions.R")
source("./bulky_functions.R")
source("./get_things.R")
source("./bdgraph.R")
source("./Gibbs_memory_optimized.R")

# Simulation Set up --------------------------------------------------------
## Partition --------------------------------------------------------------
rho_true = c(5,6,5,7,7,4,6)
Ktrue = length(rho_true)
p = sum(rho_true)

rho0_1 = rho_true
rho0_shift = c(3,7,7,5,7,3,8)
rho0_SM = c(11,5,3,4,3,4,9,1) 
sum(rho0_shift)
sum(rho0_SM)
## Graph ------------------------------------------------------------------
Gsmall = matrix(0,nrow = Ktrue, ncol = Ktrue)
Gsmall[1,2] <- Gsmall[2,5] <- Gsmall[2,6] <- 1
Gsmall[3,6] <- Gsmall[3,7] <- Gsmall[5,7] <- 1
diag(Gsmall) = rep(1,Ktrue)
rho_true_starts = c(1,rho_true[1:(Ktrue-1)])
rho_true_ends = rho_true
cum_starts = cumsum(rho_true_starts)
cum_ends = cumsum(rho_true_ends)
Gtrue = matrix(0,p,p)
for(i in 1:(Ktrue-1)){
  for(j in (i+1):Ktrue){
    if(Gsmall[i,j] > 0){
      Gtrue[cum_starts[i]:cum_ends[i],cum_starts[j]:cum_ends[j]] <- 1
    }
  }
}
Gtrue = Gtrue + t(Gtrue)
for(i in 1:Ktrue){
  if(Gsmall[i,i] > 0)
    Gtrue[cum_starts[i]:cum_ends[i],cum_starts[i]:cum_ends[i]] <- 1
}

# # visualize graph and three different partitions
# z0_1 = rho_to_z(rho0_1) 
# z0_shift = rho_to_z(rho0_shift) 
# z0_SM = rho_to_z(rho0_SM) 
# cp_true <- which(diff(z0_1) != 0)
# cp_shift <- which(diff(z0_shift) != 0)
# cp_SM <- which(diff(z0_SM) != 0)
# ACheatmap_nolegend(
#   Gtrue,
#   use_x11_device = FALSE,
#   center_value = NULL,
#   col.lower = "white"
# )
# for(k in cp_true){
#   abline(v = (k - 0.5)/(p-1), col = "red", lwd = 2)
#   abline(h = (k - 0.5)/(p-1), col = "red", lwd = 2)
# }
# for(k in cp_shift){
#   abline(v = (k - 0.5)/(p-1), col = "red", lwd = 2)
#   abline(h = (k - 0.5)/(p-1), col = "red", lwd = 2)
# }
# for(k in cp_SM){
#   abline(v = (k - 0.5)/(p-1), col = "red", lwd = 2)
#   abline(h = (k - 0.5)/(p-1), col = "red", lwd = 2)
# }

## Beta simulation -----------------------------------------------------------
seed = 123131
set.seed(seed)
n = 500
d = 3
U = diag(p)
Nrep = 2 # <---

Omega_true_arr = BDgraph::rgwish(n = Nrep, adj = Gtrue, b = d, D = U)
dim(Omega_true_arr)

# Nrep datasets, each one n x p
data_list <- lapply(seq_len(Nrep), function(b) {
  Omega_b <- Omega_true_arr[, , b]
  Sigma_b <- solve(Omega_b)
  
  MASS::mvrnorm(
    n = n,
    mu = rep(0, p),
    Sigma = Sigma_b
  )
})

# Initialization ----------------------------------------------------------
mat_ones <- matrix(1, nrow = 40, ncol = 40)  # start with all ones
diag(mat_ones) <- 0                          # set diagonal to 0

a_sigma = 1
b_sigma = 1
initialization_values_h = set_initialization_h(
  Beta          = matrix(rnorm(n=p*n), nrow = p, ncol = n),
  mu            = rep(0,p),
  tau_eps       = 0,
  K             = matrix(0,p,p),
  G             = matrix(0,p,p),
  z             = rep(1,p), 
  rho           = p,
  a_sigma       = a_sigma,
  b_sigma       = b_sigma,
  c_sigma       = 0.87908, 
  d_sigma       = 0.93759, 
  c_theta       = 0.87908, 
  d_theta       = 0.93759, 
  sigma         = 0.5,
  theta         = 3,
  weights_a0    = rep(1,p-1),
  weights_d0    = rep(1,p-1),
  total_weights = 0,
  total_K       = matrix(0,p,p),
  total_graphs  = matrix(0,p,p),
  graph_start   = mat_ones,
  graph_density = 0.5, 
  beta_sig2     = 0.1, 
  d             = 3,
  gamma = c(rep(0,p-1),1)
) 

algorithm_graph <- Sys.getenv("FGM_GRAPH_ALGORITHM", unset = "rjmcmc")
if(!algorithm_graph %in% c("rjmcmc", "bdmcmc", "rjmcmc.mpl"))
  stop("FGM_GRAPH_ALGORITHM must be one of: rjmcmc, bdmcmc, rjmcmc.mpl")
etas <- c(0, 0.5, 0.75, 0.9)

## Run options ------------------------------------------------------------
# Set each option to TRUE/FALSE to decide which initial partitions are run
# and saved for every pair (data_idx, eta).
run_rho0_1     <- TRUE
run_rho0_shift <- FALSE
run_rho0_SM    <- FALSE
debug_sampler <- identical(Sys.getenv("FGM_DEBUG_SAMPLER"), "1")

is_eta_zero = function(eta){
  isTRUE(all.equal(eta, 0))
}

should_run_eta_for_config = function(config_name, eta){
  config_name == "rho0_1" || !is_eta_zero(eta)
}

niter   <- 20#0000 # <---
burn_in <- 0
thin = 1#0
algorithm <- "rjmcmc"
(niter-burn_in)/thin

# Output directory
out_dir <- "chains"
dir.create(out_dir, showWarnings = FALSE)


# Random partition --------------------------------------------------------
cat("\nConfigurations to run:\n")
cat("  rho0_1     = ", run_rho0_1, "\n", sep = "")
cat("  rho0_shift = ", run_rho0_shift, "\n", sep = "")
cat("  rho0_SM    = ", run_rho0_SM, "\n", sep = "")
cat("  debug_sampler = ", debug_sampler, "\n", sep = "")
cat("  algorithm_graph = ", algorithm_graph, "\n", sep = "")

for(data_idx in 1:Nrep){
  cat("\n","--- Repetitions ",data_idx,"/",Nrep," --- ","\n")
  initialization_values_h$Beta = t(data_list[[data_idx]])
  for (eta in etas) {
    
    eta_chr <- as.character(eta)
    
    if(run_rho0_1 && should_run_eta_for_config("rho0_1", eta)){
      message("  eta = ", eta, " | rho0_1")
      
      ## First run: rho_0 is the true partition
      chain_rho0_1 <- Gibbs_sampler_update_h_optimized(
        set_UpdateParamsGSL_list = NULL,
        niter,
        initialization_values_h,
        alpha_target       = 0.234,
        alpha_add          = 0.5,
        adaptation_step    = 1 / (10 * p),
        seed               = 22111996,
        update_sigma_prior = TRUE,
        update_theta_prior = TRUE,
        update_weights     = TRUE,
        update_partition   = TRUE,
        update_graph       = TRUE,
        perform_shuffle    = TRUE,
        update_gamma       = TRUE,
        rho_0              = rho0_1,
        eta                = eta,
        compute_partition_update_info = FALSE,
        sample_eta         = FALSE,
        algorithm_graph    = algorithm_graph,
        rj_iters           = 1,
        thin_save = thin,
        keep_beta = F,
        show_progress = !debug_sampler,
        debug_sampler = debug_sampler
      )
      message("  sampler returned | eta = ", eta, " | rho0_1")
      
      fit_file <- file.path(
        out_dir,
        paste0(
          "SS_rho0_1_eta_", eta_chr,"_Nsim_",data_idx,
          "_niter", niter,
          "_thin", thin,
          ".rds"
        )
      )
      saveRDS(
        chain_rho0_1,
        file = fit_file
      )
      message("  saved ", fit_file)
      
      rm(chain_rho0_1)
    }
    
    if(run_rho0_shift && should_run_eta_for_config("rho0_shift", eta)){
      message("  eta = ", eta, " | rho0_shift")
      
      ## Second run: rho_0 is the shifted wrt the true partition
      chain_rho0_shift <- Gibbs_sampler_update_h_optimized(
        set_UpdateParamsGSL_list = NULL,
        niter,
        initialization_values_h,
        alpha_target       = 0.234,
        alpha_add          = 0.5,
        adaptation_step    = 1 / (10 * p),
        seed               = 22111996,
        update_sigma_prior = TRUE,
        update_theta_prior = TRUE,
        update_weights     = TRUE,
        update_partition   = TRUE,
        update_graph       = TRUE,
        perform_shuffle    = TRUE,
        update_gamma       = TRUE,
        rho_0              = rho0_shift,
        eta                = eta,
        compute_partition_update_info = FALSE,
        sample_eta         = FALSE,
        algorithm_graph    = algorithm_graph,
        rj_iters           = 1,
        thin_save = thin,
        keep_beta = F,
        show_progress = !debug_sampler,
        debug_sampler = debug_sampler
      )
      message("  sampler returned | eta = ", eta, " | rho0_shift")
      
      fit_file <- file.path(
        out_dir,
        paste0(
          "SS_rho0_shift_eta_", eta_chr,"_Nsim_",data_idx,
          "_niter", niter,
          "_thin", thin,
          ".rds"
        )
      )
      saveRDS(
        chain_rho0_shift,
        file = fit_file
      )
      message("  saved ", fit_file)
      
      rm(chain_rho0_shift)
    }
    
    if(run_rho0_SM && should_run_eta_for_config("rho0_SM", eta)){
      message("  eta = ", eta, " | rho0_SM")
      
      ## Third run: rho_0 is a split-merge modification of the true partition
      chain_rho0_SM <- Gibbs_sampler_update_h_optimized(
        set_UpdateParamsGSL_list = NULL,
        niter,
        initialization_values_h,
        alpha_target       = 0.234,
        alpha_add          = 0.5,
        adaptation_step    = 1 / (10 * p),
        seed               = 22111996,
        update_sigma_prior = TRUE,
        update_theta_prior = TRUE,
        update_weights     = TRUE,
        update_partition   = TRUE,
        update_graph       = TRUE,
        perform_shuffle    = TRUE,
        update_gamma       = TRUE,
        rho_0              = rho0_SM,
        eta                = eta,
        compute_partition_update_info = FALSE,
        sample_eta         = FALSE,
        algorithm_graph    = algorithm_graph,
        rj_iters           = 1,
        thin_save = thin,
        keep_beta = F,
        show_progress = !debug_sampler,
        debug_sampler = debug_sampler
      )
      message("  sampler returned | eta = ", eta, " | rho0_SM")
      
      fit_file <- file.path(
        out_dir,
        paste0(
          "SS_rho0_SM_eta_", eta_chr,"_Nsim_",data_idx,
          "_niter", niter,
          "_thin", thin,
          ".rds"
        )
      )
      saveRDS(
        chain_rho0_SM,
        file = fit_file
      )
      message("  saved ", fit_file)
      
      rm(chain_rho0_SM)
    }
  }
}



