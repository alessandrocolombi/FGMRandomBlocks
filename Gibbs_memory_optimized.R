# data = t(Beta_current - mu_current)
# options = options
# s = s
# rho_0 = rho_0
# algorithm = algorithm_graph

Gibbs_sampler_h_fast = function(
    data,
    options,
    s,
    rho_0,
    algorithm,
    debug_sampler = FALSE){
  
  if(debug_sampler){
    message(sprintf("[%s] inner iter %s: enter Gibbs_sampler_h_fast",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s))
    flush.console()
  }
  
  n = nrow(data)
  p = ncol(data)
  
  # --- STATE VARIABLES (no history)
  sigma_prior = options$sigma_prior_0
  theta_prior = options$theta_prior_0
  rho         = options$rho
  weights_a   = options$weights_a0
  weights_d   = options$weights_d0
  total_weights = options$total_weights0
  total_K       = options$total_K0
  total_graphs  = options$total_graphs0
  graph         = options$graph
  
  adaptation_step = options$adaptation_step
  
  gamma = options$gamma
  eta   = options$eta
  
  beta_params = estimate_Beta_params(options$beta_mu, options$beta_sig2)
  
  # constant parameters
  alpha_add    = options$alpha_add
  alpha_target = options$alpha_target
  d = options$d

  
  t_over_p = s / p
  
  # ===================== GRAPH UPDATE =====================

  if (options$update_graph) {
    if(debug_sampler){
      message(sprintf("[%s] inner iter %s: start post_graph_sampling",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s))
      flush.console()
    }

    graph_print = suppressWarnings(
      as.integer(Sys.getenv("FGM_GRAPH_PRINT", unset = "1"))
    )
    if(is.na(graph_print) || graph_print <= 0)
      graph_print = 1L

    graph_cores = suppressWarnings(
      as.integer(Sys.getenv("FGM_GRAPH_CORES", unset = "1"))
    )
    if(is.na(graph_cores) || graph_cores <= 0)
      graph_cores = 1L

    if(identical(Sys.getenv("FGM_DUMP_GRAPH_INPUT"), "1")){
      dump_dir = Sys.getenv("FGM_DUMP_GRAPH_INPUT_DIR", unset = "")
      if(!nzchar(dump_dir))
        dump_dir = file.path(getwd(), "debug_graph_inputs")
      dir.create(dump_dir, recursive = TRUE, showWarnings = FALSE)

      dump_file = file.path(
        dump_dir,
        sprintf(
          "graph_input_pid%s_iter%s_%s.rds",
          Sys.getpid(),
          s,
          format(Sys.time(), "%Y%m%d_%H%M%S")
        )
      )

      saveRDS(
        list(
          data = data,
          rho = rho,
          n = n,
          algorithm = algorithm,
          rj_steps = options$rj_steps,
          graph_print = graph_print,
          graph_cores = graph_cores,
          graph = graph,
          beta_params = beta_params,
          threshold = 1e-8,
          timestamp = Sys.time(),
          session_info = sessionInfo()
        ),
        file = dump_file
      )
      message(sprintf("[%s] inner iter %s: dumped graph input to %s",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s, dump_file))
      flush.console()
    }
    
    output = post_graph_sampling(
      data, rho, n,
      method = "ggm",
      algorithm = algorithm,
      iter = options$rj_steps,
      burnin = 0,
      g.start = graph,
      save = TRUE,
      print = graph_print,
      cores = graph_cores,
      threshold = 1e-8,
      beta_params = beta_params
    )
    
    if(debug_sampler){
      message(sprintf("[%s] inner iter %s: end post_graph_sampling",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s))
      flush.console()
    }
    
    graph  = output$last_graph
    last_K = output$last_K
    
    total_graphs  = total_graphs + graph
    total_K       = total_K + last_K
    total_weights = total_weights + 1
  }else{
    stop("This Gibbs currently supports only graph updating")
  }
  
  # ===================== GAMMA UPDATE =====================
  if(debug_sampler){
    message(sprintf("[%s] inner iter %s: start gamma update",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s))
    flush.console()
  }
  
  if(options$update_gamma){
    eta_sampled = if(options$sample_eta) runif(1,0,0.99) else eta
    
    rvec_current = c(rho_to_r(rho),1)
    
    gamma_updated = update_gamma_h_optimized(
      rvec = c(rho_to_r(rho_0),1),
      rvec_current = rvec_current,
      endpoints = cumsum(rho_0)[-length(rho_0)],
      modify_endpoint = FALSE,
      gamma = gamma,
      eta = eta_sampled,
      theta = theta_prior,
      sigma = sigma_prior
    )$gamma_updated
  } else {
    gamma_updated <- gamma
    eta_sampled <- eta
  }
  
  # ===================== PARTITION UPDATE =====================
  if(debug_sampler){
    message(sprintf("[%s] inner iter %s: start partition update",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s))
    flush.console()
  }
  
  if(options$update_partition){
    
    # ===== KEEP ORIGINAL SAFE LOGIC =====
    weights_a_ext = c(weights_a, 1)
    weights_d_ext = c(weights_d, 1)
    
    rho_current_list = split_by_gamma(rho, gamma_updated)
    wa_list  = split_weights_by_gamma(weights_a_ext, gamma_updated)
    wd_list  = split_weights_by_gamma(weights_d_ext, gamma_updated)
    
    H = length(rho_current_list)
    
    # ===== CRUCIAL: REMOVE LAST ARTIFICIAL NODE =====
    wa_list[[H]] = wa_list[[H]][-length(wa_list[[H]])]
    wd_list[[H]] = wd_list[[H]][-length(wd_list[[H]])]
    
    # ===== PREALLOC LISTS (FAST & SAFE) =====
    rho_list_updated <- vector("list", H)
    wa_h_updated <- vector("list", H)
    wd_h_updated <- vector("list", H)
    
    for(h in seq_len(H)){
      
      rho_h = rho_current_list[[h]]
      wa_h  = wa_list[[h]]
      wd_h  = wd_list[[h]]
      
      uph = update_partition_h(
        rho_current_list,
        rho_h,
        options$alpha_add,
        wa_h,
        wd_h,
        theta_prior,
        sigma_prior,
        graph,
        beta_params,
        h = h, H = H, p = p
      )
      
      # ===== UPDATE WEIGHTS =====
      if(options$update_weights){
        if(uph$accepted){ 
          if(uph$choose_add){
            wa_h = update_weight(
              weights = wa_h,
              uph$candidate,
              adaptation_step,
              t_over_p,
              uph$alpha_accept,
              alpha_target
            )
          } else {
            wd_h = update_weight(
              weights = wd_h,
              uph$candidate,
              adaptation_step,
              t_over_p,
              uph$alpha_accept,
              alpha_target
            )
          }
        }
      }
      
      rho_list_updated[[h]] = uph$rho_updated
      wa_h_updated[[h]] = wa_h
      wd_h_updated[[h]] = wd_h
    }
    
    # ===== SAFE RECONSTRUCTION =====
    rho <- unlist(rho_list_updated, use.names = FALSE)
    weights_a <- unlist(wa_h_updated, use.names = FALSE)
    weights_d <- unlist(wd_h_updated, use.names = FALSE)
  }
  
  # ===================== PYP =====================
  if(debug_sampler){
    message(sprintf("[%s] inner iter %s: start PYP updates",
                    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s))
    flush.console()
  }
  
  if(options$update_sigma_prior){
    candidate <- runif(1, max(0,-theta_prior), 1)
    alpha_MH <- full_conditional_sigma(candidate, theta_prior, rho,
                                       options$sigma_prior_parameters$a,
                                       options$sigma_prior_parameters$b,
                                       options$sigma_prior_parameters$c,
                                       options$sigma_prior_parameters$d) -
      full_conditional_sigma(sigma_prior, theta_prior, rho,
                             options$sigma_prior_parameters$a,
                             options$sigma_prior_parameters$b,
                             options$sigma_prior_parameters$c,
                             options$sigma_prior_parameters$d)
    
    if(log(runif(1)) <= min(alpha_MH,0)){
      sigma_prior <- candidate
    }
  }
  
  if(options$update_theta_prior){
    theta_prior = full_conditional_theta(
      options$theta_prior_parameters$c,
      options$theta_prior_parameters$d,
      theta_prior,
      length(rho),
      p,
      sigma_prior
    )
  }
  
  last_S = if(options$update_graph) get_S_from_G_rho(graph,rho) else NULL
  
  # ===================== RETURN ONLY CURRENT STATE =====================
  return(list(
    K              = last_K,
    #plinks         = total_graphs / total_weights,
    graph_samples  = graph,
    total_K        = total_K,
    total_graphs   = total_graphs,
    total_weights  = total_weights,
    rho            = rho,
    weights_a      = weights_a,
    weights_d      = weights_d,
    gamma          = gamma_updated,
    eta            = eta_sampled,
    sigma          = sigma_prior,
    theta          = theta_prior,
    last_S         = last_S
  ))
}



# set_UpdateParamsGSL_list
# niter = 500
# init_vals_h = initialization_values_h
# alpha_target       = 0.234
# alpha_add          = 0.5
# adaptation_step    = 1 / (10 * p)
# seed               = 22111996
# update_sigma_prior = TRUE
# update_theta_prior = TRUE
# update_weights     = TRUE
# update_partition   = TRUE
# update_graph       = TRUE
# perform_shuffle    = TRUE
# update_gamma       = TRUE
# rho_0              = rho_0
# eta                = 0.5
# compute_partition_update_info = FALSE
# sample_eta         = FALSE
# algorithm_graph    = "rjmcmc"
# rj_iters           = 1
# thin_save = 100
# keep_beta = FALSE

# OLD one with no simulation study option
# Gibbs_sampler_update_h_optimized = function(
#     set_UpdateParamsGSL_list,
#     niter,
#     initialization_values_h,
#     alpha_target,
#     alpha_add,
#     adaptation_step,
#     seed,
#     update_sigma_prior,
#     update_theta_prior,
#     update_weights,
#     update_partition,
#     update_graph,
#     perform_shuffle,
#     update_gamma,
#     rho_0,
#     eta,
#     compute_partition_update_info,
#     sample_eta,
#     algorithm_graph,
#     rj_iters,
#     thin_save = 100,
#     keep_beta = FALSE){
#   set.seed(seed)
#   
#   # ===== DIMENSIONS =====
#   p <- sum(initialization_values_h$rho)    
#   n <- ncol(initialization_values_h$Beta)
#   
#   # ===== THINNING =====
#   n_save <- max(1, niter %/% thin_save)    
#   
#   # ===== PREALLOCATION =====
#   chains <- list(
#     K = array(NA_real_, c(p,p,n_save)),
#     K_mean = array(NA_real_, c(p,p,1)),
#     plinks = array(NA_real_, c(p,p,1)),
#     graph_samples = array(FALSE, c(p,p,n_save)),
#     
#     total_weights = numeric(n_save),
#     
#     sigma = numeric(n_save),
#     theta = numeric(n_save),
#     
#     mu = matrix(NA_real_, n_save, p),
#     tau_eps = numeric(n_save),
#     
#     weights_a = numeric(p),
#     weights_d = numeric(p),
#     
#     rho = vector("list", n_save),
#     gamma = vector("list", n_save)
#   )
#   
#   if(keep_beta){
#     chains$Beta <- array(NA_real_, c(p, n, n_save))
#   }
#   
#   # ===== CURRENT STATE =====
#   Beta_current    <- initialization_values_h$Beta
#   mu_current      <- initialization_values_h$mu
#   tau_eps_current <- initialization_values_h$tau_eps
#   
#   rho            <- initialization_values_h$rho
#   gamma          <- initialization_values_h$gamma
#   weights_a      <- initialization_values_h$weights_a
#   weights_d      <- initialization_values_h$weights_d
#   total_weights  <- initialization_values_h$total_weights
#   total_K        <- initialization_values_h$total_K
#   total_graphs   <- initialization_values_h$total_graphs
#   graph_start    <- initialization_values_h$graph_start
#   
#   K_current <- initialization_values_h$K   
#   
#   sigma_current <- initialization_values_h$sigma
#   theta_current <- initialization_values_h$theta
#   
#   seeds = sample.int(.Machine$integer.max, niter)
#   
#   save_idx <- 0
#   
#   pb = txtProgressBar(min = 2, max = niter, style = 3)
#   
#   for(s in 2:niter){
#     
#     # ===================== PARAMETER UPDATE =====================
#     fit = UpdateParamsGSL(
#       Beta_current,
#       mu_current,
#       tau_eps_current,
#       K_current,   
#       set_UpdateParamsGSL_list$tbase_base,
#       set_UpdateParamsGSL_list$tbase_data,
#       set_UpdateParamsGSL_list$Sdata,
#       set_UpdateParamsGSL_list$a_tau_eps,
#       set_UpdateParamsGSL_list$b_tau_eps,
#       set_UpdateParamsGSL_list$sigma_mu,
#       set_UpdateParamsGSL_list$r,
#       set_UpdateParamsGSL_list$Update_Beta,
#       set_UpdateParamsGSL_list$Update_Mu,
#       set_UpdateParamsGSL_list$Update_Tau,
#       seeds[s]
#     )
#     
#     # ===== UPDATE CURRENT STATE =====
#     Beta_current    <- fit$Beta
#     mu_current      <- fit$mu
#     tau_eps_current <- fit$tau_eps
#     
#     # ===================== OPTIONS =====================
#     options = set_options_h(
#       sigma_prior_0 = sigma_current,
#       sigma_prior_parameters = list(
#         a=initialization_values_h$a_sigma,
#         b=initialization_values_h$b_sigma,
#         c=initialization_values_h$c_sigma,
#         d=initialization_values_h$d_sigma
#       ),
#       theta_prior_0 = theta_current,
#       theta_prior_parameters=list(
#         c=initialization_values_h$c_theta,
#         d=initialization_values_h$d_theta
#       ),
#       rho = rho,
#       weights_a0 = weights_a,
#       weights_d0 = weights_d,
#       total_weights0 = total_weights,
#       total_K0 = total_K,
#       total_graphs0 = total_graphs,
#       graph = graph_start,
#       alpha_target = alpha_target,
#       beta_mu = initialization_values_h$graph_density,
#       beta_sig2 = initialization_values_h$beta_sig2,
#       d = initialization_values_h$d,
#       alpha_add = alpha_add,
#       adaptation_step = adaptation_step,
#       update_sigma_prior = update_sigma_prior,
#       update_theta_prior = update_theta_prior,
#       update_weights = update_weights,
#       update_partition = update_partition,
#       update_graph = update_graph,
#       perform_shuffle = perform_shuffle,
#       gamma = gamma,
#       update_gamma = update_gamma,
#       eta = eta,
#       sample_eta = sample_eta
#     )
#     
#     options$rj_steps = rj_iters
#     
#     # ===================== INNER GIBBS =====================
#     res <- Gibbs_sampler_h_fast(
#       data = t(Beta_current - mu_current),
#       options = options,
#       s = s,
#       rho_0 = rho_0,
#       algorithm = algorithm_graph
#     )
#     
#     # ===================== UPDATE STATE =====================
#     rho <- res$rho
#     gamma <- res$gamma
#     
#     weights_a <- res$weights_a
#     weights_d <- res$weights_d
#     
#     total_weights <- res$total_weights
#     total_K <- res$total_K
#     total_graphs <- res$total_graphs
#     
#     graph_start <- res$graph_samples
#     
#     K_current <- res$K  
#     
#     sigma_current <- res$sigma
#     theta_current <- res$theta
#     
#     # ===================== THINNED STORAGE =====================
#     if(s %% thin_save == 0){
#       save_idx <- save_idx + 1
#       
#       chains$K[,,save_idx] <- res$K
#       chains$graph_samples[,,save_idx] <- res$graph_samples
#       
#       chains$total_weights[save_idx] = res$total_weights
#       
#       chains$mu[save_idx, ] <- mu_current
#       chains$tau_eps[save_idx] <- tau_eps_current
#       
#       chains$sigma[save_idx] <- sigma_current
#       chains$theta[save_idx] <- theta_current
#       
#       if(keep_beta){
#         chains$Beta[,,save_idx] <- Beta_current
#       }
#       
#       if(s == niter){
#         chains$weights_a = res$weights_a
#         chains$weights_d = res$weights_d
#         chains$plinks[,,1] <- total_graphs / total_weights
#         chains$K_mean[,,1] <- total_K / total_weights
#       }
#       
#       chains$rho[[save_idx]] <- rho
#       chains$gamma[[save_idx]] <- gamma
#     }
#     
#     setTxtProgressBar(pb, s)
#   }
#   
#   close(pb)
#   return(chains)
# }


# NEW one with simulation study option (by setting set_UpdateParamsGSL_list = NULL)
# when calling
Gibbs_sampler_update_h_optimized = function(
    set_UpdateParamsGSL_list,
    niter,
    initialization_values_h,
    alpha_target,
    alpha_add,
    adaptation_step,
    seed,
    update_sigma_prior,
    update_theta_prior,
    update_weights,
    update_partition,
    update_graph,
    perform_shuffle,
    update_gamma,
    rho_0,
    eta,
    compute_partition_update_info,
    sample_eta,
    algorithm_graph,
    rj_iters,
    thin_save = 100,
    keep_beta = FALSE,
    show_progress = TRUE,
    debug_sampler = FALSE){
  set.seed(seed)
  
  # ===== DIMENSIONS =====
  p <- sum(initialization_values_h$rho)    
  n <- ncol(initialization_values_h$Beta)
  
  # ===== THINNING =====
  n_save <- max(1, niter %/% thin_save)    
  
  # ===== PREALLOCATION =====
  if(is.null(set_UpdateParamsGSL_list)){
    chains <- list(
      K = array(NA_real_, c(p,p,n_save)),
      K_mean = array(NA_real_, c(p,p,1)),
      plinks = array(NA_real_, c(p,p,1)),
      graph_samples = array(FALSE, c(p,p,n_save)),
      
      total_weights = numeric(n_save),
      
      sigma = numeric(n_save),
      theta = numeric(n_save),
      
      weights_a = numeric(p),
      weights_d = numeric(p),
      
      rho = vector("list", n_save),
      gamma = vector("list", n_save)
    )
  }else{
    chains <- list(
      K = array(NA_real_, c(p,p,n_save)),
      K_mean = array(NA_real_, c(p,p,1)),
      plinks = array(NA_real_, c(p,p,1)),
      graph_samples = array(FALSE, c(p,p,n_save)),
      
      total_weights = numeric(n_save),
      
      sigma = numeric(n_save),
      theta = numeric(n_save),
      
      mu = matrix(NA_real_, n_save, p),
      tau_eps = numeric(n_save),
      
      weights_a = numeric(p),
      weights_d = numeric(p),
      
      rho = vector("list", n_save),
      gamma = vector("list", n_save)
    )
  }
  
  if(keep_beta){
    chains$Beta <- array(NA_real_, c(p, n, n_save))
  }
  
  # ===== CURRENT STATE =====
  # current state of "functional update" 
  # NB: in simulation or if don't want functional update 
  # just let set_UpdateParamsGSL_list = NULL so these remains as initialization values
  Beta_current    <- initialization_values_h$Beta
  mu_current      <- initialization_values_h$mu
  tau_eps_current <- initialization_values_h$tau_eps
  
  # current state of "non functional update"
  rho            <- initialization_values_h$rho
  gamma          <- initialization_values_h$gamma
  weights_a      <- initialization_values_h$weights_a
  weights_d      <- initialization_values_h$weights_d
  total_weights  <- initialization_values_h$total_weights
  total_K        <- initialization_values_h$total_K
  total_graphs   <- initialization_values_h$total_graphs
  graph_start    <- initialization_values_h$graph_start
  
  K_current <- initialization_values_h$K   
  
  sigma_current <- initialization_values_h$sigma
  theta_current <- initialization_values_h$theta
  
  seeds = sample.int(.Machine$integer.max, niter)
  
  save_idx <- 0
  
  if(show_progress)
    pb = txtProgressBar(min = 2, max = niter, style = 3)
  
  for(s in 2:niter){
    if(debug_sampler){
      message(sprintf("[%s] outer iter %s/%s: start",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s, niter))
      flush.console()
    }
    
    if(!is.null(set_UpdateParamsGSL_list)){
      # ===================== PARAMETER UPDATE =====================
      fit = UpdateParamsGSL(
        Beta_current,
        mu_current,
        tau_eps_current,
        K_current,   
        set_UpdateParamsGSL_list$tbase_base,
        set_UpdateParamsGSL_list$tbase_data,
        set_UpdateParamsGSL_list$Sdata,
        set_UpdateParamsGSL_list$a_tau_eps,
        set_UpdateParamsGSL_list$b_tau_eps,
        set_UpdateParamsGSL_list$sigma_mu,
        set_UpdateParamsGSL_list$r,
        set_UpdateParamsGSL_list$Update_Beta,
        set_UpdateParamsGSL_list$Update_Mu,
        set_UpdateParamsGSL_list$Update_Tau,
        seeds[s]
      )
      
      # ===== UPDATE CURRENT STATE =====
      Beta_current    <- fit$Beta
      mu_current      <- fit$mu
      tau_eps_current <- fit$tau_eps 
    }
    
    # ===================== OPTIONS =====================
    options = set_options_h(
      sigma_prior_0 = sigma_current,
      sigma_prior_parameters = list(
        a=initialization_values_h$a_sigma,
        b=initialization_values_h$b_sigma,
        c=initialization_values_h$c_sigma,
        d=initialization_values_h$d_sigma
      ),
      theta_prior_0 = theta_current,
      theta_prior_parameters=list(
        c=initialization_values_h$c_theta,
        d=initialization_values_h$d_theta
      ),
      rho = rho,
      weights_a0 = weights_a,
      weights_d0 = weights_d,
      total_weights0 = total_weights,
      total_K0 = total_K,
      total_graphs0 = total_graphs,
      graph = graph_start,
      alpha_target = alpha_target,
      beta_mu = initialization_values_h$graph_density,
      beta_sig2 = initialization_values_h$beta_sig2,
      d = initialization_values_h$d,
      alpha_add = alpha_add,
      adaptation_step = adaptation_step,
      update_sigma_prior = update_sigma_prior,
      update_theta_prior = update_theta_prior,
      update_weights = update_weights,
      update_partition = update_partition,
      update_graph = update_graph,
      perform_shuffle = perform_shuffle,
      gamma = gamma,
      update_gamma = update_gamma,
      eta = eta,
      sample_eta = sample_eta
    )
    
    options$rj_steps = rj_iters
    
    # ===================== INNER GIBBS =====================
    res <- Gibbs_sampler_h_fast(
      data = t(Beta_current - mu_current),
      options = options,
      s = s,
      rho_0 = rho_0,
      algorithm = algorithm_graph,
      debug_sampler = debug_sampler
    )
    
    if(debug_sampler){
      message(sprintf("[%s] outer iter %s/%s: end inner Gibbs",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"), s, niter))
      flush.console()
    }
    
    # ===================== UPDATE STATE =====================
    rho <- res$rho
    gamma <- res$gamma
    
    weights_a <- res$weights_a
    weights_d <- res$weights_d
    
    total_weights <- res$total_weights
    total_K <- res$total_K
    total_graphs <- res$total_graphs
    
    graph_start <- res$graph_samples
    
    K_current <- res$K  
    
    sigma_current <- res$sigma
    theta_current <- res$theta
    
    # ===================== THINNED STORAGE =====================
    if(s %% thin_save == 0){
      save_idx <- save_idx + 1
      
      chains$K[,,save_idx] <- res$K
      chains$graph_samples[,,save_idx] <- res$graph_samples
      
      chains$total_weights[save_idx] = res$total_weights
      
      if(!is.null(set_UpdateParamsGSL_list)){
        chains$mu[save_idx, ] <- mu_current
        chains$tau_eps[save_idx] <- tau_eps_current 
      }
      
      chains$sigma[save_idx] <- sigma_current
      chains$theta[save_idx] <- theta_current
      
      if(keep_beta){
        chains$Beta[,,save_idx] <- Beta_current
      }
      
      if(s == niter){
        chains$weights_a = res$weights_a
        chains$weights_d = res$weights_d
        chains$plinks[,,1] <- total_graphs / total_weights
        chains$K_mean[,,1] <- total_K / total_weights
      }
      
      chains$rho[[save_idx]] <- rho
      chains$gamma[[save_idx]] <- gamma
    }
    
    if(show_progress)
      setTxtProgressBar(pb, s)
  }
  
  if(show_progress)
    close(pb)
  return(chains)
}
