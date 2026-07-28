# wd ----------------------------------------------------------------------
wd_pc_ale = "C:/Users/colom/FGMRandomBlocks/"
wd_pc_luciano = "C://Users//lucia//Desktop//PhD//my collaborations//Alessandro Colombi work//FGMRandomBlocks"
wd_bocconi = "/home/colombi/FGMRandomBlocks/"
wd_vec = c(wd_pc_ale,wd_pc_luciano,wd_bocconi)
choose_wd = wd_vec[1] # <--- modify here to select the wd according to the user
wd = paste0(choose_wd,"./")
setwd(wd)
# Load libraries ----------------------------------------------------------

library("tidyverse")
#library("ACutils") # devtools::install_github("https://github.com/alessandrocolombi/ACutils")
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
# library("fda")
library("coda")
library("lattice")
# library(pheatmap)
library('mclust')
# Load custom functions ---------------------------------------------------

compute_ARI <- function(est, z_expert, z_true){
  c(
    ARI_true = adjustedRandIndex(z_true, est),
    ARI_expert = adjustedRandIndex(z_expert, est)
  )
}

matrix_sqrt <- function(A) {
  eig <- eigen(A, symmetric = TRUE)
  eig$vectors %*% diag(sqrt(pmax(eig$values, 0))) %*% t(eig$vectors)
}

wasserstein_gaussian <- function(Sigma1, Sigma2) {
  Sigma1_half <- matrix_sqrt(Sigma1)
  middle <- Sigma1_half %*% Sigma2 %*% Sigma1_half
  middle_half <- matrix_sqrt(middle)
  W2_sq <- sum(diag(Sigma1 + Sigma2 - 2 * middle_half))
  sqrt(max(W2_sq, 0))
}


source("utility_functions.R")
source("bulky_functions.R")

# this one because it does not load from FGM the get_things script, I do not know why!
source("get_things.R")
# this one for modified graph sampling function
source("bdgraph.R")
# this one is to load optimized Gibbs
source("Gibbs_memory_optimized.R")

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

idx <- upper.tri(Gtrue)
g_true <- Gtrue[idx]


# # visualize graph and three different partitions
z0_1 = rho_to_z(rho0_1)
z0_shift = rho_to_z(rho0_shift)
z0_SM = rho_to_z(rho0_SM)
cp_true <- which(diff(z0_1) != 0)
cp_shift <- which(diff(z0_shift) != 0)
cp_SM <- which(diff(z0_SM) != 0)
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

true_cpn_nodes = c(cp_true,40)
z_true = rho_to_z(rho0_1)

run_one_chain = function(rho0, eta, update_bool_list){
  if(eta=='CSDA'){eta_ = 0}else{eta_ = as.numeric(eta)}
  chain = Gibbs_sampler_update_h_optimized(
    set_UpdateParamsGSL_list = NULL,
    niter,
    initialization_values_h,
    alpha_target       = 0.234,
    alpha_add          = 0.5,
    adaptation_step    = 1 / (10 * p),
    seed               = 22111996,
    update_sigma_prior = update_bool_list$update_sigma_prior,
    update_theta_prior = update_bool_list$update_theta_prior,
    update_weights     = update_bool_list$update_weights,
    update_partition   = update_bool_list$update_partition,
    update_graph       = update_bool_list$update_graph,
    perform_shuffle    = update_bool_list$perform_shuffle,
    update_gamma       = update_bool_list$update_gamma,
    rho_0              = rho0,
    eta                = eta_,
    compute_partition_update_info = FALSE,
    sample_eta         = FALSE,
    algorithm_graph    = algorithm_graph,
    rj_iters           = 1,
    thin_save = thin,
    keep_beta = F
  ) 
  return(chain)
}


evaluate_chain <- function(chains,rho_type,eta,nsimul,niter,burn_in,thin,
                           g_true, save_img_graph = F, graphplot_dir = NULL,
                           z_true){
  chains = keep_last(chains, n = (niter-burn_in)/thin)
  
  ### Precision & Covariance Evaluation
  Omega_est  <- chains$K_mean[, , 1]
  Omega_true <- Omega_true_arr[, , nsimul]
  Sigma_est  <- solve(Omega_est)
  Sigma_true <- solve(Omega_true)
  # Wasserstein 
  Wasserstein <- wasserstein_gaussian(Sigma_est, Sigma_true)
  # Frobenius 
  Frobenius <- norm(Omega_est - Omega_true, type = "F")
  Relative_Frobenius <- Frobenius / norm(Omega_true, type = "F")
  
  ### Graph estimation
  bfdr_select <- BFDR_selection(chains$plinks[,,1],tol = seq(0.1, 1, by = 0.001))
  G_est <- bfdr_select$best_truncated_graph
  g_est <- G_est[idx]
  
  if(save_img_graph){
      png(filename = 
            file.path(graphplot_dir,
                      paste0("graph_eta_", eta,"_", rho_type,"_sim_", nsimul,".png")
                      ),
          width = 1200,
          height = 1200,
          res = 200
      )
      ACheatmap_nolegend(
        G_est,
        use_x11_device = FALSE,
        center_value = NULL,
        col.lower = "white",
        main = paste("eta =", eta,"|", rho_type,"| simulation", nsimul)
      )
      ## add true changepoints
      for(k in cp_true){
        abline(v = (k-0.5)/(p-1), col="red", lwd=2)
        abline(h = (k-0.5)/(p-1), col="red", lwd=2)
      }
      dev.off()
    } # END IF plot graph

  # Confusion matrix
  TP <- sum(g_true == 1 & g_est == 1)
  FP <- sum(g_true == 0 & g_est == 1)
  FN <- sum(g_true == 1 & g_est == 0)
  TN <- sum(g_true == 0 & g_est == 0)
  # Metrics
  SHD <- FP + FN
  # between / within edges
  TP_within <- sum(g_true[within_block] == 1 & g_est [within_block] == 1)
  FP_within <- sum(g_true[within_block] == 0 & g_est [within_block] == 1)
  FN_within <- sum(g_true[within_block] == 1 & g_est [within_block] == 0)
  den <- 2*TP_within + FP_within + FN_within
  F1_within <- if(den == 0) NA else (2*TP_within)/den
  
  TP_between <- sum(g_true[between_block] == 1 & g_est [between_block] == 1)
  FP_between <- sum(g_true[between_block] == 0 & g_est [between_block] == 1)
  FN_between <- sum(g_true[between_block] == 1 & g_est [between_block] == 0)
  den <- 2*TP_between + FP_between + FN_between
  F1_between <- if(den == 0) NA else (2*TP_between)/den
  
  den <- 2*TP + FP + FN
  F1 <- if (den == 0) NA else (2*TP)/den
  
  ### partition estimation
  if(rho_type == 'None'){
    results_single_chain <- data.frame(
      nsimul = nsimul,
      eta = eta,
      rho_type = rho_type,
      
      ARI_true_bind   = 0,
      ARI_true_vi     = 0,
      ARI_true_p05    = 0,
      ARI_true_argmax = 0,
      
      ARI_exp_bind   = 0,
      ARI_exp_vi     = 0,
      ARI_exp_p05    = 0,
      ARI_exp_argmax = 0,
      
      threshold = bfdr_select$best_treshold,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      SHD = SHD,
      STD_SHD = SHD / choose(p,2),
      F1 = F1,
      Sensitivity = if (TP+FN==0) NA else TP/(TP+FN),
      Specificity = if (TN+FP==0) NA else TN/(TN+FP),
      
      F1_within = F1_within,
      F1_between = F1_between,
      
      Wasserstein = Wasserstein,
      Frobenius = Frobenius,
      Relative_Frobenius = Relative_Frobenius
    )
  }else{
    if(rho_type == 'rho0_1'){
      # expert_cpn_nodes = c(cp_true,40)
      z_expert <- rho_to_z(rho0_1)
    }else if(rho_type == 'rho0_shift'){
      # expert_cpn_nodes = c(cp_shift,40)
      z_expert <- rho_to_z(rho0_shift)
    }else if(rho_type == 'rho0_SM'){
      ## Expert partition = rho0_SM
      # expert_cpn_nodes = c(cp_SM,40)
      z_expert <- rho_to_z(rho0_SM)
    }else{
      stop("Unknown rho_type")
    }
    
    rho <- chains$rho
    r_cpn = do.call(rbind, lapply(rho, rho_to_r))
    r_cpn = cbind(r_cpn, rep(1,nrow(r_cpn)) )
    z = do.call(rbind, lapply(rho, rho_to_z))
    # num_clusters = do.call(rbind, lapply(rho, length))
    # num_clusters = as.vector(num_clusters)
    # maxK = max(num_clusters)
    # TabK = tabulate(num_clusters, nbins = maxK)/(niter-burn_in)
    r_no_final = r_cpn[,1:(p-1)]
    bar_heights_norm = colSums(r_no_final)/nrow(r_no_final)
    # bars_names = as.character(1:(p-1))
    sim_matrix <- salso::psm(z)
    rownames(sim_matrix) <- 1:p
    colnames(sim_matrix) <- 1:p
    # pheatmap(sim_matrix,cluster_rows = FALSE,cluster_cols = FALSE)
    
    est_part_bind = minbinder(sim_matrix)$cl #binder
    est_part_vi = minVI(sim_matrix)$cl #VI
    est_p05 <- r_to_z( c(as.integer(bar_heights_norm > 0.5),1) )
    
    K <- apply(z,1,max)
    argmax_K = which.max(table(K)) %>% names() %>% as.integer()
    idx_K <- which(K==argmax_K)
    B_argmaxK <- (z[idx_K,-40] != z[idx_K,-1])*1
    p_argmaxK <- colMeans(B_argmaxK)
    risk_argmaxK <- B_argmaxK %*% (1-p_argmaxK) + (1-B_argmaxK) %*% p_argmaxK
    est_argmaxK <- z[idx_K[which.min(risk_argmaxK)],]
    
    ari_bind   <- compute_ARI(est_part_bind, z_expert, z_true)
    ari_vi     <- compute_ARI(est_part_vi, z_expert, z_true)
    ari_p05    <- compute_ARI(est_p05, z_expert, z_true)
    ari_argmax <- compute_ARI(est_argmaxK, z_expert, z_true)
    
    results_single_chain <- data.frame(
      nsimul = nsimul,
      eta = eta,
      rho_type = rho_type,
      
      ARI_true_bind   = ari_bind[1],
      ARI_true_vi     = ari_vi[1],
      ARI_true_p05    = ari_p05[1],
      ARI_true_argmax = ari_argmax[1],
      
      ARI_exp_bind   = ari_bind[2],
      ARI_exp_vi     = ari_vi[2],
      ARI_exp_p05    = ari_p05[2],
      ARI_exp_argmax = ari_argmax[2],
      
      threshold = bfdr_select$best_treshold,
      TP = TP,
      FP = FP,
      FN = FN,
      TN = TN,
      SHD = SHD,
      STD_SHD = SHD / choose(p,2),
      F1 = F1,
      Sensitivity = if (TP+FN==0) NA else TP/(TP+FN),
      Specificity = if (TN+FP==0) NA else TN/(TN+FP),
      
      F1_within = F1_within,
      F1_between = F1_between,
      
      Wasserstein = Wasserstein,
      Frobenius = Frobenius,
      Relative_Frobenius = Relative_Frobenius
    )
  } # END IF rho type
  # results_single_chain <- do.call(rbind, results_single_chain)
  rownames(results_single_chain) <- NULL
  return(results_single_chain)
}

## Beta simulation -----------------------------------------------------------
seed = 123131

n = 500
d = 10
U = diag(p)
Nrep = 50

set.seed(seed)
Omega_true_arr = BDgraph::rgwish(n = Nrep, adj = Gtrue, b = d, D = U)

set.seed(seed)
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

algorithm_graph <- "rjmcmc"
etas <- c('CSDA', 0, 0.5, 0.75, 0.9)

# B_idx <- 1:2   # 1 = (0.5, 0.1), 2 = (0.6, 0.03)
# 
# EBeta = c(0.5,0.6)
# SigBeta = c(0.1,0.03)

niter   <- 500
burn_in <- 20
thin = 2
(niter-burn_in)/thin

# Output directory
out_dir <- "chains"
dir.create(out_dir, showWarnings = FALSE)

## TRUE if the pair belongs to the same partition block
within_block <- outer(z_true, z_true, "==")[idx]
## TRUE if the pair belongs to different blocks
between_block <- !within_block

rho_list <- list(
  rho0_1     = rho0_1,
  rho0_shift = rho0_shift,
  rho0_SM    = rho0_SM
)

# Gibbs -------------------------------------------------------------------

# graph_path = 'C:\\Users\\lucia\\Desktop\\PhD\\my collaborations\\Alessandro Colombi work\\Simulation Study\\Graph_plots'
# graphplot_dir <- file.path(graph_path)

# files_n = expand.grid(
#   data_idx = 1:Nrep,
#   eta = etas,
#   rho_type = names(rho_list),
#   # B_idx = B_idx,
#   KEEP.OUT.ATTRS = FALSE#,
#   # stringsAsFactors = FALSE
# )

nchains <- Nrep * (2 + (length(etas)-2)*length(rho_list))
# nchains <- 2 * (2 + (length(etas)-2)*length(rho_list))

results_graph_part_dbig <- vector("list", nchains)

counter = 1

# # Output directory
# setwd("C:\\Users\\lucia\\Desktop\\PhD\\my collaborations\\Alessandro Colombi work\\Simulation Study")
# out_dir <- "chains"
# dir.create(out_dir, showWarnings = FALSE)


for(eta in etas){
  message('eta: ', eta)
  for(data_idx in 1:Nrep){
    message('N simul: ', data_idx)
    initialization_values_h$Beta = t(data_list[[data_idx]])
    if(eta == 'CSDA'){
      rho_name_csda = "None"
      update_bool_list = list( update_sigma_prior = F,
                               update_theta_prior = F,
                               update_weights     = F,
                               update_partition   = F,
                               update_graph       = TRUE,
                               perform_shuffle    = F,
                               update_gamma       = F)
      chain = run_one_chain(rho0=p, eta, update_bool_list)
      # results for CSDA
      results_graph_part_dbig[[counter]] <-
        evaluate_chain(
          chains=chain, rho_type=rho_name_csda, eta=eta, nsimul=data_idx, 
          niter=niter, burn_in=burn_in, thin=thin, g_true=g_true, z_true=z_true
        )
      counter <- counter + 1
      rm(chain)
    }else if(eta == 0){
      rho_name_freepart = "rho0_1"
      update_bool_list = list( update_sigma_prior = T,
                               update_theta_prior = T,
                               update_weights     = T,
                               update_partition   = T,
                               update_graph       = TRUE,
                               perform_shuffle    = T,
                               update_gamma       = T)
      chain = run_one_chain(rho0=p, eta, update_bool_list)
      # results for CSDA
      results_graph_part_dbig[[counter]] <-
        evaluate_chain(
          chains=chain, rho_type=rho_name_freepart, eta=eta, nsimul=data_idx, 
          niter=niter, burn_in=burn_in, thin=thin, g_true=g_true, z_true=z_true
        )
      counter <- counter + 1
      rm(chain)
    }else{
      for(rho_name in names(rho_list)){
        message('rho type: ', rho_name)
        rho0 = rho_list[[rho_name]]
        update_bool_list = list( update_sigma_prior = T,
                                 update_theta_prior = T,
                                 update_weights     = T,
                                 update_partition   = T,
                                 update_graph       = TRUE,
                                 perform_shuffle    = T,
                                 update_gamma       = T)
        chain = run_one_chain(rho0, eta, update_bool_list)
        # results for eta in 0, 0.5, 0.75, 0.9
        results_graph_part_dbig[[counter]] <-
          evaluate_chain(
            chains=chain, rho_type=rho_name, eta=eta, nsimul=data_idx, 
            niter=niter, burn_in=burn_in, thin=thin, g_true=g_true, z_true=z_true
          )
        counter <- counter + 1
        rm(chain)
        } # END FOR rho types
    } # END IF etas is csda or others
  } # END FOR data
} # END FOR etas


# class(results_graph_part_dbig)
# length(results_graph_part_dbig)

results_graph_part_dbig <- do.call(rbind, results_graph_part_dbig)
results_graph_part_dbig


saveRDS(results_graph_part_dbig, file = 'results_graph_part_dbig')
# results_graph_part_dbig = readRDS('C:\\Users\\lucia\\Desktop\\PhD\\my collaborations\\Alessandro Colombi work\\Simulation Study\\results_graph_part_dbig.rds')


