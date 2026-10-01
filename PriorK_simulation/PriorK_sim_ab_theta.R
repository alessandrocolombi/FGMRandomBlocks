# Prior di K con sigma~Beta(1,1), theta+sigma~Gamma(shape=a_theta, rate=b_theta).
# Unico output: tabella media (sd). Eseguibile con Ctrl+Invio.
# VM: ./run_job.sh PriorK_simulation/PriorK_sim_ab_theta.R choose_wd=2 eta_usr=0.5

# ==================== CONFIGURAZIONE UTENTE ==========================
p <- 40L
B <- 10000L
seed <- 20261001L
eta_usr <- 1
a_sigma <- b_sigma <- 1  # Fissi anche nel campionatore condiviso.
a_theta_grid <- c(0.5,0.88,1.5)
b_theta_grid <- c(0.5,0.94,1.5) # RATE
# ====================================================================

.ab_prior_dir <- if (file.exists("PriorK_helpers.R")) "." else "PriorK_simulation"
source(file.path(.ab_prior_dir, "PriorK_helpers.R"), local = TRUE)
.cli_args <- if (sys.nframe() == 0L) commandArgs(trailingOnly = TRUE) else character()
.options <- read_prior_options(.cli_args, eta_usr)
eta_usr <- .options$eta_usr
.ab_prior_dir <- if (file.exists("PriorK_helpers.R")) "." else "PriorK_simulation"
source(file.path(.ab_prior_dir, "PriorKh_test.R"), local = TRUE)
eta <- make_prior_eta(eta_usr) # Eta delle posizioni esperte; eta[40]=1 sempre.
output_file <- file.path(.ab_prior_dir, "results", "PriorK_ab_theta_mean_sd_table.csv")

simulate_prior_K_ab_theta <- function(a_theta_grid, b_theta_grid, eta,
                                      B = 1000L, seed = 20261001L,
                                      output_file = file.path(.ab_prior_dir, "results",
                                                             "PriorK_ab_theta_mean_sd_table.csv")) {
  check_prior_simulation(eta, B)
  if (!is.numeric(a_theta_grid) || !is.numeric(b_theta_grid) ||
      !length(a_theta_grid) || !length(b_theta_grid))
    stop("Le griglie devono essere numeriche e non vuote.")
  a_values <- unique(a_theta_grid); b_values <- unique(b_theta_grid)
  cells <- matrix(NA_character_, length(a_values), length(b_values),
                  dimnames = list(NULL, paste0("b_theta=", as.character(b_values))))
  set.seed(seed)
  for (j in seq_along(b_values)) for (i in seq_along(a_values)) {
    a_theta <- a_values[i]; b_theta <- b_values[j]
    if (!is.finite(a_theta) || !is.finite(b_theta) || a_theta <= 0 || b_theta <= 0) next
    K <- draw_K_hyperprior(a_theta, b_theta, eta, B)
    cells[i, j] <- sprintf("%.2f (%.2f)", mean(K), sd(K))
  }
  result <- data.frame(a_theta = a_values, cells, check.names = FALSE)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  write.csv(result, output_file, row.names = FALSE, na = "NA")
  print(result, row.names = FALSE, quote = FALSE, right = TRUE)
  invisible(result)
}

if (sys.nframe() == 0L)
  simulate_prior_K_ab_theta(a_theta_grid, b_theta_grid, eta, B, seed, output_file)
