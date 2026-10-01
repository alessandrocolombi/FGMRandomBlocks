# Prior marginale di K con sigma e theta random, p=40.
# Ctrl+Invio dalla cartella progetto o PriorK_simulation.
# VM: ./run_job.sh PriorK_simulation/PriorK_sim_ab_theta.R choose_wd=2
# Unico file prodotto: tabella CSV con celle media (sd) di K.

set_ab_prior_working_directory <- function(args = character()) {
  wd_vec <- c("C:/Users/colom/FGMRandomBlocks", "/home/colombi/FGMRandomBlocks")
  if (!length(args)) return(invisible(getwd()))
  if (length(args) != 1L || !grepl("^choose_wd=[12]$", args))
    stop("Opzione ammessa: choose_wd=1 (PC) oppure choose_wd=2 (VM).")
  wd <- wd_vec[as.integer(sub("^choose_wd=", "", args))]
  if (!dir.exists(wd)) stop("Directory non trovata: ", wd)
  setwd(wd)
  message("Working directory: ", getwd())
}
set_ab_prior_working_directory(
  if (sys.nframe() == 0L) commandArgs(trailingOnly = TRUE) else character())

# ==================== CONFIGURAZIONE UTENTE ==========================
p <- 40L
B <- 10000L
seed <- 20261001L
a_sigma <- b_sigma <- 1           # FISSI: sigma ~ Uniforme(0,1)
a_theta_grid <- c(0.5,0.88,1.5)   # Esempio modificabile: shape Gamma
b_theta_grid <- c(0.5,0.94,1.5)   # Esempio modificabile: RATE Gamma

eta <- numeric(p)
eta[c(4, 6, 9, 13, 18, 22, 28, 33, 40)] <- 1
# Con questi eta: H=9, p_tilde=(4,2,3,4,5,4,6,5,7).
.ab_prior_dir <- if (file.exists("log_C_PY.cpp")) "." else "PriorK_simulation"
output_file <- file.path(.ab_prior_dir, "results", "PriorK_ab_theta_mean_sd_table.csv")
# ====================================================================

# Riutilizza la PMF verificata e il nuovo C++; non esegue i test al source().
source(file.path(.ab_prior_dir, "PriorKh_test.R"), local = TRUE)

simulate_prior_K_ab_theta <- function(a_theta_grid, b_theta_grid, eta,
                                      B = 1000L, seed = 20261001L,
                                      output_file = file.path(.ab_prior_dir, "results",
                                                             "PriorK_ab_theta_mean_sd_table.csv")) {
  p <- 40L
  a_sigma <- b_sigma <- 1
  if (!is.numeric(a_theta_grid) || !is.numeric(b_theta_grid) ||
      !length(a_theta_grid) || !length(b_theta_grid))
    stop("Le griglie devono essere numeriche e non vuote.")
  if (length(B) != 1L || !is.finite(B) || B < 2 || B != floor(B))
    stop("B deve essere un intero >=2.")
  if (length(eta) != p || any(!is.finite(eta)) || any(eta < 0 | eta > 1) || eta[p] != 1)
    stop("eta deve contenere 40 probabilita' in [0,1], con eta[40]=1.")
  a_values <- unique(a_theta_grid)
  b_values <- unique(b_theta_grid)
  cells <- matrix(NA_character_, nrow = length(a_values), ncol = length(b_values),
                  dimnames = list(NULL, paste0("b_theta=", as.character(b_values))))
  if (any(is.finite(a_values) & a_values > 0) &&
      any(is.finite(b_values) & b_values > 0)) load_logC()
  set.seed(seed)

  for (j in seq_along(b_values)) for (i in seq_along(a_values)) {
    a_theta <- a_values[i]
    b_theta <- b_values[j]
    if (!is.finite(a_theta) || !is.finite(b_theta) || a_theta <= 0 || b_theta <= 0) {
      message(sprintf("Salto a_theta=%g, b_theta=%g: richiesti parametri Gamma positivi e finiti.",
                      a_theta, b_theta))
      next
    }
    K <- integer(B)
    for (b in seq_len(B)) {
      # UNA coppia (theta,sigma) per replica, condivisa da tutti i sottogruppi.
      sigma <- rbeta(1L, shape1 = a_sigma, shape2 = b_sigma)
      U <- rgamma(1L, shape = a_theta, rate = b_theta)
      theta <- U - sigma
      # Nessun clipping o ricampionamento: non alterare la prior in caso di
      # underflow/arrotondamento con iperparametri estremi.
      if (!is.finite(theta) || sigma <= 0 || sigma >= 1 || U <= 0 || theta <= -sigma)
        stop(sprintf(paste0("Precisione numerica insufficiente: a_theta=%g, b_theta=%g, ",
                            "replica=%d. Impossibile rappresentare theta>-sigma."),
                     a_theta, b_theta, b))
      gamma <- rbinom(p, size = 1L, prob = eta)
      sizes <- diff(c(0L, which(gamma == 1L)))

      # Cache solo per questa replica: theta e sigma cambiano alla successiva.
      # Stessa PMF e stessi controlli massa/media di PriorK_simulation.R.
      distinct_sizes <- unique(sizes)
      pmfs <- lapply(distinct_sizes, function(n) prior_Kh_pmf(n, theta, sigma)$probability)
      Kh <- vapply(sizes, function(n) {
        sample.int(n, size = 1L, prob = pmfs[[match(n, distinct_sizes)]])
      }, integer(1))
      K[b] <- sum(Kh)
      stopifnot(K[b] >= length(sizes), K[b] <= p)
    }
    # sd della distribuzione marginale di K, non errore Monte Carlo della media.
    cells[i, j] <- sprintf("%.2f (%.2f)", mean(K), sd(K))
    message(sprintf("Completato a_theta=%g, b_theta=%g (%d repliche).", a_theta, b_theta, B))
  }
  result <- data.frame(a_theta = a_values, cells, check.names = FALSE)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  write.csv(result, output_file, row.names = FALSE, na = "NA")
  cat("\nPrior marginale di K: media (sd); sigma ~ Beta(1,1); b_theta = rate\n")
  print(result, row.names = FALSE, quote = FALSE, right = TRUE)
  invisible(result)
}

# Eseguire anche questa riga con Ctrl+Invio per avviare la simulazione.
if (sys.nframe() == 0L)
  simulate_prior_K_ab_theta(a_theta_grid, b_theta_grid, eta, B, seed, output_file)
