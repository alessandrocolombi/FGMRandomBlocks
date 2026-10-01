# Simulazione marginale di K, con p=40.
# Esecuzione: Rscript PriorK_simulation/PriorK_simulation.R
# Usa i coefficienti corretti di Lijoi et al. (2007), tramite log_C_PY.cpp.
# VM: ./run_job.sh PriorK_simulation/PriorK_simulation.R choose_wd=2

# Directory selezionabile da command line, come nello script di esempio.
# Senza choose_wd si mantiene la directory corrente (anche con Ctrl+Invio).
set_prior_working_directory <- function(args = character()) {
  wd_vec <- c("C:/Users/colom/FGMRandomBlocks",  # choose_wd=1: PC
              "/home/colombi/FGMRandomBlocks") # choose_wd=2: VM
  if (!length(args)) return(invisible(getwd()))
  if (length(args) != 1L || !grepl("^choose_wd=[12]$", args))
    stop("Opzione ammessa: choose_wd=1 (PC) oppure choose_wd=2 (VM).")
  choose_wd <- as.integer(sub("^choose_wd=", "", args))
  wd <- wd_vec[choose_wd]
  if (!dir.exists(wd)) stop("Directory non trovata: ", wd)
  setwd(wd)
  message("Working directory: ", getwd())
  invisible(getwd())
}
set_prior_working_directory(
  if (sys.nframe() == 0L) commandArgs(trailingOnly = TRUE) else character())

# ==================== CONFIGURAZIONE UTENTE ==========================
p <- 40L
B <- 1000L
seed <- 20260930L
theta_grid <- c(5,10)#c(1,5,10,20)
sigma_grid <- c(1e-3,0.1)#c(1e-3, 0.1, 0.5, 0.8)  # 0 <= sigma < 1; theta > -sigma

# Sequenza richiesta dall'utente, modificabile per altre configurazioni.
# Eta_j>0 identifica una posizione candidata a changepoint esperto.
# Per rispettare il paper, assegnare 0 alle posizioni non esperte.
eta <- numeric(p)
eta[c(4, 6, 9, 13, 18, 22, 28, 33, 40)] <- 1
# eta[p]=1 garantisce gamma_p=1. Con questi eta, gamma=eta quasi certamente.
# H=9 e p_tilde=c(4, 2, 3, 4, 5, 4, 6, 5, 7).
# ====================================================================

# Eseguibile anche con Ctrl+Invio dalla cartella progetto o da questa cartella.
.simulation_dir <- if (file.exists("log_C_PY.cpp")) "." else "PriorK_simulation"
# Carica PMF e media teorica; load_logC() compila il nuovo log_C_PY.cpp.
source(file.path(.simulation_dir, "PriorKh_test.R"), local = TRUE)

# Righe/colonne nell'ordine delle griglie; sd e' la deviazione standard di K,
# non l'errore Monte Carlo della media. Celle arrotondate a due decimali.
prior_K_mean_sd_table <- function(result) {
  theta_values <- unique(result$theta)
  sigma_values <- unique(result$sigma)
  cells <- matrix(NA_character_, nrow = length(theta_values),
                  ncol = length(sigma_values))
  cells[cbind(match(result$theta, theta_values), match(result$sigma, sigma_values))] <-
    sprintf("%.2f (%.2f)", result$mean, sqrt(result$variance))
  colnames(cells) <- paste0("sigma=", as.character(sigma_values))
  data.frame(theta = theta_values, cells, check.names = FALSE)
}

simulate_prior_K <- function(theta_grid, sigma_grid, eta, B = 1000L,
                             seed = 20260930L,
                             output_dir = file.path(.simulation_dir, "results")) {
  p <- 40L
  if (length(eta) != p || any(!is.finite(eta)) || any(eta < 0 | eta > 1) || eta[p] != 1)
    stop("eta deve avere 40 probabilita' in [0,1], con eta[40]=1.")
  if (length(B) != 1L || !is.finite(B) || B < 2 || B != floor(B))
    stop("B deve essere un intero >=2.")
  if (!length(theta_grid) || !length(sigma_grid)) stop("Le griglie non possono essere vuote.")
  grid <- unique(expand.grid(theta = theta_grid, sigma = sigma_grid))
  for (i in seq_len(nrow(grid))) check_py_parameters(p, grid$theta[i], grid$sigma[i])

  # Nuova implementazione: log(C_mathcal/sigma^k), senza fattori di segno.
  # Compilazione una sola volta; include il limite Dirichlet sigma=0.
  load_logC()
  # Verifica massa, intervallo e media (11) delle PMF non normalizzate.
  # La cache evita il ricalcolo dei coefficienti ad ogni replica Monte Carlo.
  pmfs <- lapply(seq_len(nrow(grid)), function(i) {
    lapply(seq_len(p), function(n) prior_Kh_pmf(n, grid$theta[i], grid$sigma[i]))
  })
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)
  summaries <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    theta <- grid$theta[i]
    sigma <- grid$sigma[i]
    gamma_draws <- matrix(0L, nrow = B, ncol = p)
    sizes_draws <- Kh_draws <- vector("list", B)
    K <- H <- integer(B)
    conditional_mean <- numeric(B)
    for (b in seq_len(B)) {
      gamma <- rbinom(p, size = 1L, prob = eta)
      sizes <- diff(c(0L, which(gamma == 1L)))
      Kh <- vapply(sizes, function(n) {
        sample.int(n, size = 1L, prob = pmfs[[i]][[n]]$probability)
      }, integer(1))
      gamma_draws[b, ] <- gamma
      sizes_draws[[b]] <- sizes
      Kh_draws[[b]] <- Kh
      K[b] <- sum(Kh)
      H[b] <- length(sizes)
      conditional_mean[b] <- sum(vapply(sizes, prior_Kh_mean, numeric(1),
                                        theta = theta, sigma = sigma))
    }
    stopifnot(all(K >= H), all(K <= p))
    q <- quantile(K, probs = c(0.025, 0.975), type = 1, names = FALSE)
    summaries[[i]] <- data.frame(
      p = p, B = B, seed = seed, theta = theta, sigma = sigma,
      mean = mean(K), sd = sd(K), variance = var(K), median = median(K),
      q025 = q[1], q975 = q[2], min = min(K), max = max(K),
      mcse_mean = sd(K) / sqrt(B),
      mean_conditional_expectation = mean(conditional_mean))
    # Indice univoco anche per parametri con rappresentazioni testuali simili.
    tag <- sprintf("PriorK_%03d_theta_%s_sigma_%s", i,
                   format(theta, digits = 16, trim = TRUE),
                   format(sigma, digits = 16, trim = TRUE))
    write.csv(summaries[[i]], file.path(output_dir, paste0(tag, "_summary.csv")),
              row.names = FALSE)
    saveRDS(list(config = list(p = p, B = B, seed = seed, theta = theta,
                              sigma = sigma, eta = eta, grid_index = i,
                              theta_grid = theta_grid, sigma_grid = sigma_grid,
                              quantile_type = 1L, session = sessionInfo()),
                 summary = summaries[[i]], K = K, H = H,
                 gamma = gamma_draws, p_tilde = sizes_draws, Kh = Kh_draws,
                 conditional_mean = conditional_mean),
            file.path(output_dir, paste0(tag, "_draws.rds")))
  }
  result <- do.call(rbind, summaries)
  write.csv(result, file.path(output_dir, "PriorK_grid_summary.csv"), row.names = FALSE)
  write.csv(data.frame(j = seq_len(p), eta = eta),
            file.path(output_dir, "PriorK_eta.csv"), row.names = FALSE)
  mean_sd_table <- prior_K_mean_sd_table(result)
  write.csv(mean_sd_table, file.path(output_dir, "PriorK_mean_sd_table.csv"),
            row.names = FALSE)
  cat("\nDistribuzione di K: media (sd)\n")
  print(mean_sd_table, row.names = FALSE, quote = FALSE, right = TRUE)
  invisible(result)
}

if (sys.nframe() == 0L) simulate_prior_K(theta_grid, sigma_grid, eta, B, seed)
