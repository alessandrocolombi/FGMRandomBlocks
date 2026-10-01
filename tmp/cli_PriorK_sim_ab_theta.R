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

stopifnot(eta_usr == 0.5, eta[4] == 0.5, eta[40] == 1)
