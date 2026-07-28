# wd ----------------------------------------------------------------------
wd_pc_ale = "C:/Users/colom/FGMRandomBlocks/"
wd_pc_luciano = "C://Users//lucia//Desktop//PhD//my collaborations//Alessandro Colombi work//FGMRandomBlocks"
wd_bocconi = "/home/colombi/FGMRandomBlocks/"
wd_vec = c(wd_pc_ale, wd_pc_luciano, wd_bocconi)
choose_wd = wd_vec[3] # <--- modify here to select the wd according to the user
wd = paste0(choose_wd, "./")
setwd(wd)

# Thread limits -----------------------------------------------------------
# Keep each parallel R worker single-threaded. This avoids using
# n_cores * BLAS/OpenMP threads on shared servers.
limit_threaded_libraries = function(n_threads = 1L){
  n_threads = as.character(n_threads)
  Sys.setenv(
    OMP_NUM_THREADS = n_threads,
    OMP_THREAD_LIMIT = n_threads,
    OPENBLAS_NUM_THREADS = n_threads,
    MKL_NUM_THREADS = n_threads,
    MKL_DOMAIN_NUM_THREADS = n_threads,
    BLIS_NUM_THREADS = n_threads,
    VECLIB_MAXIMUM_THREADS = n_threads,
    RCPP_PARALLEL_NUM_THREADS = n_threads,
    NUMEXPR_NUM_THREADS = n_threads,
    GOTO_NUM_THREADS = n_threads,
    OMP_DYNAMIC = "FALSE",
    MKL_DYNAMIC = "FALSE"
  )

  if(requireNamespace("RhpcBLASctl", quietly = TRUE)){
    RhpcBLASctl::blas_set_num_threads(as.integer(n_threads))
    RhpcBLASctl::omp_set_num_threads(as.integer(n_threads))
  }

  invisible(NULL)
}

limit_threaded_libraries(1L)

# Load libraries ----------------------------------------------------------

library("tidyverse")
#library("ACutils") # devtools::install_github("https://github.com/alessandrocolombi/ACutils")
library("mvtnorm")
# library("salso")
library("FGM") #  devtools::install_github("alessandrocolombi/FGMpackage")
library("gmp")
library("mcclust")
library("mcclust.ext")
library("logr")
library("tidygraph")
library("ggraph")
library("igraph")
library("Rcpp")
library("RcppArmadillo")
library("RcppEigen")
library("RcppGSL")
# library("fda")
library("coda")
library("lattice")
# library(pheatmap)
library("mclust")
library("parallel")

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

source("./utility_functions.R")
source("./bulky_functions.R")
source("./get_things.R")
source("./bdgraph.R")
source("./Gibbs_memory_optimized.R")

# Simulation Set up -------------------------------------------------------
## Partition --------------------------------------------------------------
rho_true = c(5, 6, 5, 7, 7, 4, 6)
Ktrue = length(rho_true)
p = sum(rho_true)

rho0_1 = rho_true
rho0_shift = c(3, 7, 7, 5, 7, 3, 8)
rho0_SM = c(11, 5, 3, 4, 3, 4, 9, 1)
sum(rho0_shift)
sum(rho0_SM)

## Graph ------------------------------------------------------------------
Gsmall = matrix(0, nrow = Ktrue, ncol = Ktrue)
Gsmall[1, 2] <- Gsmall[2, 5] <- Gsmall[2, 6] <- 1
Gsmall[3, 6] <- Gsmall[3, 7] <- Gsmall[5, 7] <- 1
diag(Gsmall) = rep(1, Ktrue)
rho_true_starts = c(1, rho_true[1:(Ktrue - 1)])
rho_true_ends = rho_true
cum_starts = cumsum(rho_true_starts)
cum_ends = cumsum(rho_true_ends)
Gtrue = matrix(0, p, p)
for(i in 1:(Ktrue - 1)){
  for(j in (i + 1):Ktrue){
    if(Gsmall[i, j] > 0){
      Gtrue[cum_starts[i]:cum_ends[i], cum_starts[j]:cum_ends[j]] <- 1
    }
  }
}
Gtrue = Gtrue + t(Gtrue)
for(i in 1:Ktrue){
  if(Gsmall[i, i] > 0)
    Gtrue[cum_starts[i]:cum_ends[i], cum_starts[i]:cum_ends[i]] <- 1
}

idx <- upper.tri(Gtrue)
g_true <- Gtrue[idx]

z0_1 = rho_to_z(rho0_1)
z0_shift = rho_to_z(rho0_shift)
z0_SM = rho_to_z(rho0_SM)
cp_true <- which(diff(z0_1) != 0)
cp_shift <- which(diff(z0_shift) != 0)
cp_SM <- which(diff(z0_SM) != 0)

true_cpn_nodes = c(cp_true, 40)
z_true = rho_to_z(rho0_1)

# Beta simulation ---------------------------------------------------------
seed = 123131

n = 500
d = 10
U = diag(p)
Nrep = 15 # <---

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
mat_ones <- matrix(1, nrow = 40, ncol = 40)
diag(mat_ones) <- 0

a_sigma = 1
b_sigma = 1
initialization_values_h = set_initialization_h(
  Beta          = matrix(rnorm(n = p * n), nrow = p, ncol = n),
  mu            = rep(0, p),
  tau_eps       = 0,
  K             = matrix(0, p, p),
  G             = matrix(0, p, p),
  z             = rep(1, p),
  rho           = p,
  a_sigma       = a_sigma,
  b_sigma       = b_sigma,
  c_sigma       = 0.87908,
  d_sigma       = 0.93759,
  c_theta       = 0.87908,
  d_theta       = 0.93759,
  sigma         = 0.5,
  theta         = 3,
  weights_a0    = rep(1, p - 1),
  weights_d0    = rep(1, p - 1),
  total_weights = 0,
  total_K       = matrix(0, p, p),
  total_graphs  = matrix(0, p, p),
  graph_start   = mat_ones,
  graph_density = 0.5,
  beta_sig2     = 0.1,
  d             = 3,
  gamma = c(rep(0, p - 1), 1)
)

algorithm_graph <- Sys.getenv("FGM_GRAPH_ALGORITHM", unset = "rjmcmc")
if(!algorithm_graph %in% c("rjmcmc", "bdmcmc", "rjmcmc.mpl"))
  stop("FGM_GRAPH_ALGORITHM must be one of: rjmcmc, bdmcmc, rjmcmc.mpl")

etas <- c("CSDA", 0, 0.5, 0.75, 0.9)

## MCMC options -----------------------------------------------------------
niter   <- 500
burn_in <- 20
thin = 2
(niter - burn_in) / thin

sampler_seed <- 22111996
alpha_target <- 0.234
alpha_add <- 0.5
adaptation_step <- 1 / (10 * p)
rj_iters <- 1
keep_beta <- FALSE

## Parallel options -------------------------------------------------------
avail_cores = parallel::detectCores(logical = TRUE)
if(is.na(avail_cores))
  avail_cores = 1L
requested_cores = 15L # <---
n_cores = min(requested_cores, avail_cores, Nrep)

# Output directories ------------------------------------------------------
out_dir <- "chains"
result_dir <- file.path(out_dir, "results_parallel_2")
log_dir <- file.path(out_dir, "logs_parallel_2")
status_dir <- file.path(out_dir, "status_parallel_2")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)

## TRUE if the pair belongs to the same partition block
within_block <- outer(z_true, z_true, "==")[idx]
## TRUE if the pair belongs to different blocks
between_block <- !within_block

rho_list <- list(
  rho0_1     = rho0_1,
  rho0_shift = rho0_shift,
  rho0_SM    = rho0_SM
)

if(identical(Sys.getenv("SIMULATION_STUDY_2_PARALLEL_SMOKE"), "1")){
  Nrep = min(Nrep, 1L)
  data_list = data_list[seq_len(Nrep)]
  Omega_true_arr = Omega_true_arr[, , seq_len(Nrep), drop = FALSE]
  etas = c("CSDA", 0.5)
  n_cores = 1L
  niter = 4
  burn_in = 0
  thin = 2
  out_dir = file.path(tempdir(), "Simulation_Study_2_parallel_smoke")
  result_dir = file.path(out_dir, "results_parallel_2")
  log_dir = file.path(out_dir, "logs_parallel_2")
  status_dir = file.path(out_dir, "status_parallel_2")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)
}

# Parallel runner helpers -------------------------------------------------

project_dir = normalizePath(wd, winslash = "/", mustWork = TRUE)
source_paths = file.path(
  project_dir,
  c(
    "utility_functions.R",
    "bulky_functions.R",
    "get_things.R",
    "bdgraph.R",
    "Gibbs_memory_optimized.R"
  )
)

is_eta_zero = function(eta){
  eta_num = suppressWarnings(as.numeric(eta))
  !is.na(eta_num) && isTRUE(all.equal(eta_num, 0))
}

sanitize_tag_value = function(x){
  gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
}

make_update_bool_list = function(update_partition){
  list(
    update_sigma_prior = update_partition,
    update_theta_prior = update_partition,
    update_weights     = update_partition,
    update_partition   = update_partition,
    update_graph       = TRUE,
    perform_shuffle    = update_partition,
    update_gamma       = update_partition
  )
}

make_tasks_for_rep = function(data_idx, etas, rho_list, p){
  tasks = list()

  for(eta in etas){
    if(eta == "CSDA"){
      tasks[[length(tasks) + 1]] = list(
        data_idx = data_idx,
        eta = eta,
        rho_type = "None",
        rho0 = p,
        update_bool_list = make_update_bool_list(update_partition = FALSE)
      )
    }else if(is_eta_zero(eta)){
      tasks[[length(tasks) + 1]] = list(
        data_idx = data_idx,
        eta = eta,
        rho_type = "rho0_1",
        rho0 = p,
        update_bool_list = make_update_bool_list(update_partition = TRUE)
      )
    }else{
      for(rho_name in names(rho_list)){
        tasks[[length(tasks) + 1]] = list(
          data_idx = data_idx,
          eta = eta,
          rho_type = rho_name,
          rho0 = rho_list[[rho_name]],
          update_bool_list = make_update_bool_list(update_partition = TRUE)
        )
      }
    }
  }

  tasks
}

make_task_tag = function(rho_type, eta, data_idx, niter, n_obs){
  paste0(
    "SS2_", sanitize_tag_value(rho_type),
    "_eta_", sanitize_tag_value(eta),
    "_Nsim_", data_idx,
    "_niter", niter,
    "_n", n_obs
  )
}

make_result_file = function(result_dir, rho_type, eta, data_idx, niter, n_obs){
  file.path(result_dir, paste0(make_task_tag(rho_type, eta, data_idx, niter, n_obs), ".rds"))
}

make_status_file = function(status_dir, rho_type, eta, data_idx, niter, n_obs, status){
  file.path(status_dir, paste0(make_task_tag(rho_type, eta, data_idx, niter, n_obs), ".", status))
}

fmt_path = function(x){
  normalizePath(x, winslash = "/", mustWork = FALSE)
}

append_log = function(log_file, ...){
  cat(
    sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0(..., collapse = ""),
    "\n",
    file = log_file,
    append = TRUE,
    sep = ""
  )
}

combine_data_frames = function(rows){
  rows = rows[!vapply(rows, is.null, logical(1))]
  if(length(rows) == 0)
    return(data.frame())
  do.call(rbind, rows)
}

run_one_chain = function(init_values, rho0, eta, update_bool_list,
                          niter, thin, p,
                          algorithm_graph, sampler_seed,
                          alpha_target, alpha_add, adaptation_step,
                          rj_iters, keep_beta){
  eta_ = if(eta == "CSDA") 0 else as.numeric(eta)

  Gibbs_sampler_update_h_optimized(
    set_UpdateParamsGSL_list = NULL,
    niter,
    init_values,
    alpha_target       = alpha_target,
    alpha_add          = alpha_add,
    adaptation_step    = adaptation_step,
    seed               = sampler_seed,
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
    rj_iters           = rj_iters,
    thin_save          = thin,
    keep_beta          = keep_beta,
    show_progress      = FALSE,
    debug_sampler      = FALSE
  )
}

evaluate_chain <- function(chains, rho_type, eta, nsimul, niter, burn_in, thin,
                           g_true, save_img_graph = FALSE, graphplot_dir = NULL,
                           z_true){
  chains = keep_last(chains, n = (niter - burn_in) / thin)

  ### Precision & Covariance Evaluation
  Omega_est  <- chains$K_mean[, , 1]
  Omega_true <- Omega_true_arr[, , nsimul]
  Sigma_est  <- solve(Omega_est)
  Sigma_true <- solve(Omega_true)
  Wasserstein <- wasserstein_gaussian(Sigma_est, Sigma_true)
  Frobenius <- norm(Omega_est - Omega_true, type = "F")
  Relative_Frobenius <- Frobenius / norm(Omega_true, type = "F")

  ### Graph estimation
  bfdr_select <- BFDR_selection(chains$plinks[, , 1], tol = seq(0.1, 1, by = 0.001))
  G_est <- bfdr_select$best_truncated_graph
  g_est <- G_est[idx]

  if(save_img_graph){
    png(filename =
          file.path(
            graphplot_dir,
            paste0("graph_eta_", eta, "_", rho_type, "_sim_", nsimul, ".png")
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
      main = paste("eta =", eta, "|", rho_type, "| simulation", nsimul)
    )
    for(k in cp_true){
      abline(v = (k - 0.5) / (p - 1), col = "red", lwd = 2)
      abline(h = (k - 0.5) / (p - 1), col = "red", lwd = 2)
    }
    dev.off()
  }

  # Confusion matrix
  TP <- sum(g_true == 1 & g_est == 1)
  FP <- sum(g_true == 0 & g_est == 1)
  FN <- sum(g_true == 1 & g_est == 0)
  TN <- sum(g_true == 0 & g_est == 0)
  SHD <- FP + FN

  TP_within <- sum(g_true[within_block] == 1 & g_est[within_block] == 1)
  FP_within <- sum(g_true[within_block] == 0 & g_est[within_block] == 1)
  FN_within <- sum(g_true[within_block] == 1 & g_est[within_block] == 0)
  den <- 2 * TP_within + FP_within + FN_within
  F1_within <- if(den == 0) NA else (2 * TP_within) / den

  TP_between <- sum(g_true[between_block] == 1 & g_est[between_block] == 1)
  FP_between <- sum(g_true[between_block] == 0 & g_est[between_block] == 1)
  FN_between <- sum(g_true[between_block] == 1 & g_est[between_block] == 0)
  den <- 2 * TP_between + FP_between + FN_between
  F1_between <- if(den == 0) NA else (2 * TP_between) / den

  den <- 2 * TP + FP + FN
  F1 <- if (den == 0) NA else (2 * TP) / den

  if(rho_type == "None"){
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
      STD_SHD = SHD / choose(p, 2),
      F1 = F1,
      Sensitivity = if (TP + FN == 0) NA else TP / (TP + FN),
      Specificity = if (TN + FP == 0) NA else TN / (TN + FP),

      F1_within = F1_within,
      F1_between = F1_between,

      Wasserstein = Wasserstein,
      Frobenius = Frobenius,
      Relative_Frobenius = Relative_Frobenius
    )
  }else{
    if(rho_type == "rho0_1"){
      z_expert <- rho_to_z(rho0_1)
    }else if(rho_type == "rho0_shift"){
      z_expert <- rho_to_z(rho0_shift)
    }else if(rho_type == "rho0_SM"){
      z_expert <- rho_to_z(rho0_SM)
    }else{
      stop("Unknown rho_type")
    }

    rho <- chains$rho
    r_cpn = do.call(rbind, lapply(rho, rho_to_r))
    r_cpn = cbind(r_cpn, rep(1, nrow(r_cpn)))
    z = do.call(rbind, lapply(rho, rho_to_z))
    r_no_final = r_cpn[, 1:(p - 1)]
    bar_heights_norm = colSums(r_no_final) / nrow(r_no_final)

    sim_matrix <- mcclust::comp.psm(z)  #salso::psm(z)
    rownames(sim_matrix) <- 1:p
    colnames(sim_matrix) <- 1:p

    est_part_bind = minbinder(sim_matrix)$cl
    est_part_vi = minVI(sim_matrix)$cl
    est_p05 <- r_to_z(c(as.integer(bar_heights_norm > 0.5), 1))

    K <- apply(z, 1, max)
    argmax_K = which.max(table(K)) %>% names() %>% as.integer()
    idx_K <- which(K == argmax_K)
    B_argmaxK <- (z[idx_K, -40] != z[idx_K, -1]) * 1
    p_argmaxK <- colMeans(B_argmaxK)
    risk_argmaxK <- B_argmaxK %*% (1 - p_argmaxK) + (1 - B_argmaxK) %*% p_argmaxK
    est_argmaxK <- z[idx_K[which.min(risk_argmaxK)], ]

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
      STD_SHD = SHD / choose(p, 2),
      F1 = F1,
      Sensitivity = if (TP + FN == 0) NA else TP / (TP + FN),
      Specificity = if (TN + FP == 0) NA else TN / (TN + FP),

      F1_within = F1_within,
      F1_between = F1_between,

      Wasserstein = Wasserstein,
      Frobenius = Frobenius,
      Relative_Frobenius = Relative_Frobenius
    )
  }

  rownames(results_single_chain) <- NULL
  results_single_chain
}

run_single_task = function(task, init_values,
                           niter, burn_in, thin, p, n_obs,
                           result_dir, status_dir, log_file,
                           algorithm_graph, sampler_seed,
                           alpha_target, alpha_add, adaptation_step,
                           rj_iters, keep_beta,
                           g_true, z_true){
  data_idx = task$data_idx
  eta = task$eta
  rho_type = task$rho_type
  task_tag = make_task_tag(rho_type, eta, data_idx, niter, n_obs)
  result_file = make_result_file(result_dir, rho_type, eta, data_idx, niter, n_obs)
  started_file = make_status_file(status_dir, rho_type, eta, data_idx, niter, n_obs, "started")
  finished_file = make_status_file(status_dir, rho_type, eta, data_idx, niter, n_obs, "finished")
  failed_file = make_status_file(status_dir, rho_type, eta, data_idx, niter, n_obs, "failed")

  out = tryCatch({
    if(file.exists(finished_file))
      unlink(finished_file)
    if(file.exists(failed_file))
      unlink(failed_file)

    cat(
      "started_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
      "tag: ", task_tag, "\n",
      "result_file: ", fmt_path(result_file), "\n",
      file = started_file,
      sep = ""
    )

    append_log(log_file, "START task ", task_tag, " | result_file=", fmt_path(result_file))
    task_start = Sys.time()

    chain = run_one_chain(
      init_values = init_values,
      rho0 = task$rho0,
      eta = eta,
      update_bool_list = task$update_bool_list,
      niter = niter,
      thin = thin,
      p = p,
      algorithm_graph = algorithm_graph,
      sampler_seed = sampler_seed,
      alpha_target = alpha_target,
      alpha_add = alpha_add,
      adaptation_step = adaptation_step,
      rj_iters = rj_iters,
      keep_beta = keep_beta
    )

    result_row = evaluate_chain(
      chains = chain,
      rho_type = rho_type,
      eta = eta,
      nsimul = data_idx,
      niter = niter,
      burn_in = burn_in,
      thin = thin,
      g_true = g_true,
      z_true = z_true
    )

    saveRDS(result_row, file = result_file)
    rm(chain)
    gc(verbose = FALSE)

    cat(
      "finished_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
      "tag: ", task_tag, "\n",
      "result_file: ", fmt_path(result_file), "\n",
      file = finished_file,
      sep = ""
    )

    append_log(
      log_file,
      "SAVED result ", task_tag,
      " | elapsed_min=", round(as.numeric(difftime(Sys.time(), task_start, units = "mins")), 3),
      " | result_file=", fmt_path(result_file)
    )

    manifest = data.frame(
      data_idx = data_idx,
      eta = as.character(eta),
      rho_type = rho_type,
      result_file = result_file,
      status = "success",
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )

    list(manifest = manifest, result = result_row)
  }, error = function(e){
    cat(
      "failed_at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
      "tag: ", task_tag, "\n",
      "error: ", conditionMessage(e), "\n",
      file = failed_file,
      sep = ""
    )

    append_log(log_file, "FAILED task ", task_tag, " | error=", conditionMessage(e))

    manifest = data.frame(
      data_idx = data_idx,
      eta = as.character(eta),
      rho_type = rho_type,
      result_file = result_file,
      status = "failed",
      error_message = conditionMessage(e),
      stringsAsFactors = FALSE
    )

    list(manifest = manifest, result = NULL)
  })

  out
}

run_single_rep = function(data_idx, data_list, initialization_values_h,
                          etas, rho_list,
                          niter, burn_in, thin, p, n_obs,
                          result_dir, log_dir, status_dir,
                          algorithm_graph, sampler_seed,
                          alpha_target, alpha_add, adaptation_step,
                          rj_iters, keep_beta,
                          g_true, z_true){
  log_file = file.path(log_dir, paste0("Simulation_Study_2_Nsim_", data_idx, ".log"))
  append_log(log_file, "START repetition ", data_idx)

  init_values = initialization_values_h
  init_values$Beta = t(data_list[[data_idx]])

  tasks = make_tasks_for_rep(data_idx, etas, rho_list, p)
  append_log(log_file, "TASKS repetition ", data_idx, " | n_tasks=", length(tasks))

  task_outputs = lapply(tasks, function(task){
    run_single_task(
      task = task,
      init_values = init_values,
      niter = niter,
      burn_in = burn_in,
      thin = thin,
      p = p,
      n_obs = n_obs,
      result_dir = result_dir,
      status_dir = status_dir,
      log_file = log_file,
      algorithm_graph = algorithm_graph,
      sampler_seed = sampler_seed,
      alpha_target = alpha_target,
      alpha_add = alpha_add,
      adaptation_step = adaptation_step,
      rj_iters = rj_iters,
      keep_beta = keep_beta,
      g_true = g_true,
      z_true = z_true
    )
  })

  append_log(log_file, "END repetition ", data_idx)

  list(
    manifest = combine_data_frames(lapply(task_outputs, function(x) x$manifest)),
    results = combine_data_frames(lapply(task_outputs, function(x) x$result))
  )
}

all_tasks = unlist(
  lapply(seq_len(Nrep), make_tasks_for_rep, etas = etas, rho_list = rho_list, p = p),
  recursive = FALSE
)
nchains = length(all_tasks)

cat("\nSimulation Study 2 parallel settings:\n")
cat("  Nrep            = ", Nrep, "\n", sep = "")
cat("  n_cores         = ", n_cores, "\n", sep = "")
cat("  niter           = ", niter, "\n", sep = "")
cat("  burn_in         = ", burn_in, "\n", sep = "")
cat("  thin            = ", thin, "\n", sep = "")
cat("  etas            = ", paste(etas, collapse = ", "), "\n", sep = "")
cat("  expected tasks  = ", nchains, "\n", sep = "")
cat("  algorithm_graph = ", algorithm_graph, "\n", sep = "")

cat(sprintf("[%s] START whole script: %s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Simulation_Study_2_parallel.R"))
flush.console()

cat(sprintf("[%s] Creating PSOCK cluster with %s workers\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), n_cores))
flush.console()

cl = parallel::makeCluster(n_cores)
on.exit(parallel::stopCluster(cl), add = TRUE)

cat(sprintf("[%s] Cluster created\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
flush.console()

cat(sprintf("[%s] Exporting objects to workers\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
flush.console()

parallel::clusterExport(
  cl,
  varlist = c(
    "wd", "source_paths",
    "limit_threaded_libraries",
    "data_list", "Omega_true_arr", "initialization_values_h",
    "etas", "rho_list",
    "p", "n", "idx", "g_true", "within_block", "between_block",
    "z_true", "rho0_1", "rho0_shift", "rho0_SM", "cp_true",
    "niter", "burn_in", "thin",
    "result_dir", "log_dir", "status_dir",
    "algorithm_graph", "sampler_seed",
    "alpha_target", "alpha_add", "adaptation_step",
    "rj_iters", "keep_beta",
    "compute_ARI", "matrix_sqrt", "wasserstein_gaussian",
    "is_eta_zero", "sanitize_tag_value", "make_update_bool_list",
    "make_tasks_for_rep", "make_task_tag", "make_result_file", "make_status_file",
    "fmt_path", "append_log", "combine_data_frames",
    "run_one_chain", "evaluate_chain", "run_single_task", "run_single_rep"
  ),
  envir = environment()
)

cat(sprintf("[%s] Worker object export completed\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
flush.console()

for(worker_id in seq_along(cl)){
  cat(sprintf("[%s] Loading packages/sources on worker %s/%s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), worker_id, length(cl)))
  flush.console()
  parallel::clusterCall(cl[worker_id], function(wd, source_paths, limit_threaded_libraries){
    setwd(wd)
    limit_threaded_libraries(1L)
    suppressPackageStartupMessages({
      library("tidyverse")
      #library("ACutils") # devtools::install_github("https://github.com/alessandrocolombi/ACutils")
      library("mvtnorm")
      # library("salso")
      library("FGM") #  devtools::install_github("alessandrocolombi/FGMpackage")
      library("gmp")
      library("mcclust")
      library("mcclust.ext")
      library("logr")
      library("tidygraph")
      library("ggraph")
      library("igraph")
      library("Rcpp")
      library("RcppArmadillo")
      library("RcppEigen")
      library("RcppGSL")
      # library("fda")
      library("coda")
      library("lattice")
      # library(pheatmap)
      library("mclust")
      library("parallel")
    })

    for(path in source_paths)
      source(path)

    NULL
  }, wd = wd, source_paths = source_paths, limit_threaded_libraries = limit_threaded_libraries)
  cat(sprintf("[%s] Worker %s/%s ready\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), worker_id, length(cl)))
  flush.console()
}

cat(sprintf("[%s] Starting parallel repetitions\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
flush.console()

parallel_outputs = parallel::parLapplyLB(cl, seq_len(Nrep), function(data_idx){
  run_single_rep(
    data_idx = data_idx,
    data_list = data_list,
    initialization_values_h = initialization_values_h,
    etas = etas,
    rho_list = rho_list,
    niter = niter,
    burn_in = burn_in,
    thin = thin,
    p = p,
    n_obs = n,
    result_dir = result_dir,
    log_dir = log_dir,
    status_dir = status_dir,
    algorithm_graph = algorithm_graph,
    sampler_seed = sampler_seed,
    alpha_target = alpha_target,
    alpha_add = alpha_add,
    adaptation_step = adaptation_step,
    rj_iters = rj_iters,
    keep_beta = keep_beta,
    g_true = g_true,
    z_true = z_true
  )
})

manifest_df = combine_data_frames(lapply(parallel_outputs, function(x) x$manifest))
results_graph_part_dbig = combine_data_frames(lapply(parallel_outputs, function(x) x$results))

results_file = file.path(out_dir, paste0(
  "results_graph_part_dbig_parallel",
  "_Nrep", Nrep,
  "_niter", niter,
  "_n", n,
  ".rds"
))

results_csv = file.path(out_dir, paste0(
  "results_graph_part_dbig_parallel",
  "_Nrep", Nrep,
  "_niter", niter,
  "_n", n,
  ".csv"
))

manifest_file = file.path(out_dir, paste0(
  "Simulation_Study_2_parallel_manifest",
  "_Nrep", Nrep,
  "_niter", niter,
  "_n", n,
  ".csv"
))

saveRDS(results_graph_part_dbig, file = results_file)
write.csv(results_graph_part_dbig, file = results_csv, row.names = FALSE)
write.csv(manifest_df, file = manifest_file, row.names = FALSE)

cat("\nSaved objects / output manifest\n")
cat("results_rds:   ", fmt_path(results_file), "\n", sep = "")
cat("results_csv:   ", fmt_path(results_csv), "\n", sep = "")
cat("manifest_csv:  ", fmt_path(manifest_file), "\n", sep = "")
cat("result_dir:    ", fmt_path(result_dir), "\n", sep = "")
cat("log_dir:       ", fmt_path(log_dir), "\n", sep = "")
cat("status_dir:    ", fmt_path(status_dir), "\n", sep = "")

cat("\nStatus counts:\n")
print(table(manifest_df$status))

if(any(manifest_df$status != "success")){
  cat("\nFailed tasks:\n")
  print(manifest_df[manifest_df$status != "success", , drop = FALSE])
}

cat(sprintf("[%s] END whole script: %s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Simulation_Study_2_parallel.R"))
flush.console()

results_graph_part_dbig
