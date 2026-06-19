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

## Beta simulation -----------------------------------------------------------
seed = 123131
set.seed(seed)
n = 500
d = 3
U = diag(p)
Nrep = 50 # <---

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

algorithm_graph <- "rjmcmc"
etas <- c(0, 0.5, 0.75, 0.9)

## Run options ------------------------------------------------------------
# Set each option to TRUE/FALSE to decide which initial partitions are run
# and saved for every pair (data_idx, eta).
run_rho0_1     <- TRUE
run_rho0_shift <- FALSE
run_rho0_SM    <- FALSE

rho0_configs = list(
  list(name = "rho0_1",     rho0 = rho0_1,     run = run_rho0_1),
  list(name = "rho0_shift", rho0 = rho0_shift, run = run_rho0_shift),
  list(name = "rho0_SM",    rho0 = rho0_SM,    run = run_rho0_SM)
)
rho0_configs = rho0_configs[vapply(rho0_configs, function(x) isTRUE(x$run), logical(1))]

if(length(rho0_configs) == 0)
  stop("At least one run_rho0_* option must be TRUE")

## Parallel options -------------------------------------------------------
avail_cores = parallel::detectCores(logical = TRUE)
if(is.na(avail_cores))
  avail_cores = 1L
requested_cores = 30L # <---
n_cores = min(requested_cores, avail_cores, Nrep)

## MCMC options -----------------------------------------------------------
niter   <- 200000 # <---
burn_in <- 0
thin = 10
algorithm <- "rjmcmc"
(niter-burn_in)/thin

sampler_seed <- 22111996
alpha_target <- 0.234
alpha_add <- 0.5
adaptation_step <- 1 / (10 * p)
rj_iters <- 1
keep_beta <- FALSE

# Output directory --------------------------------------------------------
out_dir <- "chains"
log_dir <- file.path(out_dir, "logs_parallel")

if(!dir.exists(out_dir))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if(!dir.exists(log_dir))
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

if(identical(Sys.getenv("SIMULATION_STUDY_PARALLEL_SMOKE"), "1")){
  Nrep = min(Nrep, 1L)
  data_list = data_list[seq_len(Nrep)]
  etas = etas[1]
  n_cores = 1L
  niter = 3
  burn_in = 0
  thin = 1
  out_dir = file.path(tempdir(), "Simulation_Study_parallel_smoke")
  log_dir = file.path(out_dir, "logs_parallel")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
}

cat("\nConfigurations to run:\n")
cat("  rho0_1     = ", run_rho0_1, "\n", sep = "")
cat("  rho0_shift = ", run_rho0_shift, "\n", sep = "")
cat("  rho0_SM    = ", run_rho0_SM, "\n", sep = "")
cat("\nParallel settings:\n")
cat("  Nrep       = ", Nrep, "\n", sep = "")
cat("  n_cores    = ", n_cores, "\n", sep = "")
cat("  etas       = ", paste(etas, collapse = ", "), "\n", sep = "")

# Parallel MCMC runner ----------------------------------------------------

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

make_chain_file = function(out_dir, config_name, eta, data_idx, niter, thin){
  eta_chr = as.character(eta)
  file.path(
    out_dir,
    paste0(
      "SS_", config_name, "_eta_", eta_chr, "_Nsim_", data_idx,
      "_niter", niter,
      "_thin", thin,
      ".rds"
    )
  )
}

fmt_path = function(x){
  normalizePath(x, winslash = "/", mustWork = FALSE)
}

run_single_chain = function(config, eta, data_idx, init_values,
                            p, niter, thin, out_dir,
                            algorithm_graph, sampler_seed,
                            alpha_target, alpha_add, adaptation_step,
                            rj_iters, keep_beta){
  fit_file = make_chain_file(out_dir, config$name, eta, data_idx, niter, thin)

  result = tryCatch({
    cat("Start chain\n")
    cat("  data_idx: ", data_idx, "\n", sep = "")
    cat("  eta:      ", eta, "\n", sep = "")
    cat("  config:   ", config$name, "\n", sep = "")
    cat("  fit_file: ", fmt_path(fit_file), "\n", sep = "")

    chain = Gibbs_sampler_update_h_optimized(
      set_UpdateParamsGSL_list = NULL,
      niter,
      init_values,
      alpha_target       = alpha_target,
      alpha_add          = alpha_add,
      adaptation_step    = adaptation_step,
      seed               = sampler_seed,
      update_sigma_prior = TRUE,
      update_theta_prior = TRUE,
      update_weights     = TRUE,
      update_partition   = TRUE,
      update_graph       = TRUE,
      perform_shuffle    = TRUE,
      update_gamma       = TRUE,
      rho_0              = config$rho0,
      eta                = eta,
      compute_partition_update_info = FALSE,
      sample_eta         = FALSE,
      algorithm_graph    = algorithm_graph,
      rj_iters           = rj_iters,
      thin_save          = thin,
      keep_beta          = keep_beta
    )

    saveRDS(chain, file = fit_file)
    rm(chain)
    gc(verbose = FALSE)

    cat("\nCompleted chain\n")

    list(
      data_idx = data_idx,
      eta = eta,
      config = config$name,
      fit_file = fit_file,
      status = "success",
      error_message = NA_character_
    )
  }, error = function(e){
    cat("\nChain failed\n")
    cat("  data_idx: ", data_idx, "\n", sep = "")
    cat("  eta:      ", eta, "\n", sep = "")
    cat("  config:   ", config$name, "\n", sep = "")
    cat("  error:    ", conditionMessage(e), "\n", sep = "")

    list(
      data_idx = data_idx,
      eta = eta,
      config = config$name,
      fit_file = fit_file,
      status = "failed",
      error_message = conditionMessage(e)
    )
  })

  result
}

run_single_rep = function(data_idx, data_list, initialization_values_h,
                          rho0_configs, etas,
                          p, niter, thin, out_dir, log_dir,
                          algorithm_graph, sampler_seed,
                          alpha_target, alpha_add, adaptation_step,
                          rj_iters, keep_beta){
  log_file = file.path(log_dir, paste0("Simulation_Study_Nsim_", data_idx, ".log"))
  log_open = FALSE
  on.exit({
    if(log_open)
      try(sink(), silent = TRUE)
  }, add = TRUE)

  sink(log_file, split = FALSE)
  log_open = TRUE

  cat(sprintf("[%s] START repetition %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), data_idx))

  init_values = initialization_values_h
  init_values$Beta = t(data_list[[data_idx]])

  rep_results = list()
  for(eta in etas){
    for(config in rho0_configs){
      rep_results[[length(rep_results) + 1]] = run_single_chain(
        config = config,
        eta = eta,
        data_idx = data_idx,
        init_values = init_values,
        p = p,
        niter = niter,
        thin = thin,
        out_dir = out_dir,
        algorithm_graph = algorithm_graph,
        sampler_seed = sampler_seed,
        alpha_target = alpha_target,
        alpha_add = alpha_add,
        adaptation_step = adaptation_step,
        rj_iters = rj_iters,
        keep_beta = keep_beta
      )
    }
  }

  cat(sprintf("[%s] END repetition %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), data_idx))

  if(log_open){
    try(sink(), silent = TRUE)
    log_open = FALSE
  }

  out = do.call(rbind, lapply(rep_results, as.data.frame))
  out$log_file = log_file
  out
}

cl = parallel::makeCluster(n_cores)
on.exit(parallel::stopCluster(cl), add = TRUE)

cat(sprintf("[%s] START whole script: %s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Simulation_Study_parallel.R"))
flush.console()

parallel::clusterExport(
  cl,
  varlist = c(
    "wd", "source_paths",
    "data_list", "initialization_values_h",
    "rho0_configs", "etas",
    "p", "niter", "thin",
    "out_dir", "log_dir",
    "algorithm_graph", "sampler_seed",
    "alpha_target", "alpha_add", "adaptation_step",
    "rj_iters", "keep_beta",
    "make_chain_file", "fmt_path",
    "run_single_chain", "run_single_rep"
  ),
  envir = environment()
)

for(worker_id in seq_along(cl)){
  parallel::clusterCall(cl[worker_id], function(wd, source_paths){
    setwd(wd)
    suppressPackageStartupMessages({
      library("tidyverse")
      library("mvtnorm")
      library("salso")
      library("FGM")
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
      library("coda")
      library("lattice")
    })

    for(path in source_paths)
      source(path)

    NULL
  }, wd = wd, source_paths = source_paths)
}

results = parallel::parLapplyLB(cl, seq_len(Nrep), function(data_idx){
  run_single_rep(
    data_idx = data_idx,
    data_list = data_list,
    initialization_values_h = initialization_values_h,
    rho0_configs = rho0_configs,
    etas = etas,
    p = p,
    niter = niter,
    thin = thin,
    out_dir = out_dir,
    log_dir = log_dir,
    algorithm_graph = algorithm_graph,
    sampler_seed = sampler_seed,
    alpha_target = alpha_target,
    alpha_add = alpha_add,
    adaptation_step = adaptation_step,
    rj_iters = rj_iters,
    keep_beta = keep_beta
  )
})

results_df = do.call(rbind, results)

summary_file = file.path(out_dir, paste0(
  "Simulation_Study_parallel_summary",
  "_Nrep", Nrep,
  "_niter", niter,
  "_thin", thin,
  ".csv"
))

write.csv(
  results_df,
  file = summary_file,
  row.names = FALSE
)

cat("\nSaved objects / output manifest\n")
cat("summary_csv: ", fmt_path(summary_file), "\n", sep = "")
cat("output_dir:  ", fmt_path(out_dir), "\n", sep = "")
cat("log_dir:     ", fmt_path(log_dir), "\n", sep = "")

for(i in seq_len(nrow(results_df))){
  cat("\nRun ", i, "/", nrow(results_df), "\n", sep = "")
  cat("status:   ", results_df$status[i], "\n", sep = "")
  cat("Nsim:     ", results_df$data_idx[i], "\n", sep = "")
  cat("eta:      ", results_df$eta[i], "\n", sep = "")
  cat("config:   ", results_df$config[i], "\n", sep = "")
  cat("fit_rds:  ", fmt_path(results_df$fit_file[i]), "\n", sep = "")
  cat("log_file: ", fmt_path(results_df$log_file[i]), "\n", sep = "")
  if(!is.na(results_df$error_message[i]))
    cat("error:    ", results_df$error_message[i], "\n", sep = "")
}

cat(sprintf("[%s] END whole script: %s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "Simulation_Study_parallel.R"))
flush.console()

results_df
