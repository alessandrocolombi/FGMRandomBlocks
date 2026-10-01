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

stopifnot(eta_usr == 0.5, eta[4] == 0.5, eta[40] == 1)
