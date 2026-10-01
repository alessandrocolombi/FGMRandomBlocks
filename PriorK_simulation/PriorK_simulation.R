# Prior di K con theta e sigma fissati. Unico output: tabella media (sd).
# VM: ./run_job.sh PriorK_simulation/PriorK_simulation.R choose_wd=2 eta_usr=0.5
# Ctrl+Invio: modificare eta_usr nella configurazione e proseguire dall'alto.

# ==================== CONFIGURAZIONE UTENTE ==========================
p <- 40L
B <- 1000L
seed <- 20260930L
eta_usr <- 1
# Griglie utente: coppie non ammissibili -> NA, senza interrompere.
theta_grid <- c(-0.5,0,1,2,5)
sigma_grid <- c(1e-5,0.1,0.51,0.75,0.99)
# ====================================================================

.simulation_dir <- if (file.exists("PriorK_helpers.R")) "." else "PriorK_simulation"
source(file.path(.simulation_dir, "PriorK_helpers.R"), local = TRUE)
.cli_args <- if (sys.nframe() == 0L) commandArgs(trailingOnly = TRUE) else character()
.options <- read_prior_options(.cli_args, eta_usr)
eta_usr <- .options$eta_usr
.simulation_dir <- if (file.exists("PriorK_helpers.R")) "." else "PriorK_simulation"
source(file.path(.simulation_dir, "PriorKh_test.R"), local = TRUE)
eta <- make_prior_eta(eta_usr) # Eta delle posizioni esperte; eta[40]=1 sempre.

simulate_prior_K <- function(theta_grid, sigma_grid, eta, B = 1000L,
                             seed = 20260930L,
                             output_dir = file.path(.simulation_dir, "results")) {
  check_prior_simulation(eta, B)
  if (!is.numeric(theta_grid) || !is.numeric(sigma_grid) ||
      !length(theta_grid) || !length(sigma_grid))
    stop("Le griglie devono essere numeriche e non vuote.")
  theta_values <- unique(theta_grid)
  sigma_values <- unique(sigma_grid)
  cells <- matrix(NA_character_, length(theta_values), length(sigma_values),
                  dimnames = list(NULL, paste0("sigma=", as.character(sigma_values))))
  set.seed(seed)
  for (j in seq_along(sigma_values)) for (i in seq_along(theta_values)) {
    theta <- theta_values[i]; sigma <- sigma_values[j]
    if (!is.finite(theta) || !is.finite(sigma) || sigma < 0 || sigma >= 1 || theta <= -sigma)
      next
    load_logC()
    pmfs <- lapply(seq_len(40L), function(n) prior_Kh_pmf(n, theta, sigma)$probability)
    K <- vapply(seq_len(B), function(b) draw_K_given_parameters(theta, sigma, eta, pmfs), integer(1))
    cells[i, j] <- sprintf("%.2f (%.2f)", mean(K), sd(K))
  }
  result <- data.frame(theta = theta_values, cells, check.names = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(result, file.path(output_dir, "PriorK_mean_sd_table.csv"), row.names = FALSE, na = "NA")
  print(result, row.names = FALSE, quote = FALSE, right = TRUE)
  invisible(result)
}

if (sys.nframe() == 0L) simulate_prior_K(theta_grid, sigma_grid, eta, B, seed)
