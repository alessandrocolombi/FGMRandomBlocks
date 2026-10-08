# Distribuzione prior marginale di K su una griglia di eta.
# Eseguibile dall'alto con Ctrl+Invio.
# Working directory: FGMRandomBlocks oppure PriorK_simulation.
# Nessuna tabella o estrazione salvata: solo figure PDF e PNG, se richiesto.

# ==================== CONFIGURAZIONE UTENTE ==========================
a_theta <- 0.88
b_theta <- 0.94                 # RATE Gamma, non scale
a_sigma <- 1                   # shape1 Beta
b_sigma <- 1                   # shape2 Beta
eta_grid <- c(0, 0.5, 0.75, 0.9)
B <- 100000L
seed <- 20261001L
# sigma ~ Beta(a_sigma,b_sigma), p=40; eta[40]=1 sempre.
save_figures <- TRUE
# ====================================================================

.plot_dir <- if (file.exists("PriorK_helpers.R")) "." else "PriorK_simulation"
source(file.path(.plot_dir, "PriorK_helpers.R"), local = TRUE)
source(file.path(.plot_dir, "PriorKh_test.R"), local = TRUE)
if (!is.numeric(eta_grid) || !length(eta_grid) || any(!is.finite(eta_grid)) ||
    any(eta_grid < 0 | eta_grid > 1)) stop("eta_grid deve contenere valori in [0,1].")
eta_grid <- unique(eta_grid)
set.seed(seed)

# Un campione e una PMF per ciascun eta; nessuna KDE continua.
prior_K_distributions <- lapply(eta_grid, function(eta_usr) {
  eta <- make_prior_eta(eta_usr)
  K <- draw_K_hyperprior(a_theta, b_theta, eta, B, a_sigma = a_sigma, b_sigma = b_sigma)
  list(eta_usr = eta_usr, probability = tabulate(K, nbins = 40L) / B,
       K_mean = mean(K))
})
# Stessi limiti e tick in tutte le figure per facilitare il confronto.
y_max <- max(vapply(prior_K_distributions, function(x) max(x$probability), numeric(1))) * 1.12
y_ticks <- pretty(c(0, y_max), n = 5)
y_ticks <- y_ticks[y_ticks >= 0 & y_ticks <= y_max]
x_ticks <- c(1, seq(5, 40, by = 5))

plot_prior_K_distribution <- function(distribution) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(4.8, 6.8, 3.2, 1), mgp = c(2.5, 0.8, 0),
      cex.axis = 2, cex.lab = 2, cex.main = 2, las = 1, bty = "l")
  k <- seq_len(40L)
  plot(k, distribution$probability, type = "n", xlim = c(0.5, 40.5),
       ylim = c(0, y_max), xaxs = "i", yaxs = "i",
       axes = FALSE, ann = FALSE, bty = "l")
  # Griglia verticale e orizzontale disegnata PRIMA delle barre.
  abline(v = x_ticks, h = y_ticks, col = "grey88", lwd = 0.8)
  rect(k - 0.38, 0, k + 0.38, distribution$probability, col = "grey30", border = NA)
  axis(1, at = x_ticks)
  axis(2, at = y_ticks, las = 1)
  box(bty = "l")
  abline(v = distribution$K_mean, col = "black", lwd = 2, lty = 2)
  title(main = bquote(eta == .(distribution$eta_usr)), line = 1)
  title(xlab = "K", line = 2.6)
  title(ylab = "P(K)", line = 4.6)
  invisible(NULL)
}

save_prior_K_figure <- function(path, format, distribution) {
  if (format == "pdf") {
    # Cairo incorpora i font e preserva eta greca nei lettori PDF.
    if (capabilities("cairo")) {
      cairo_pdf(path, width = 9, height = 6)
    } else {
      pdf(path, width = 9, height = 6, useDingbats = FALSE)
    }
  } else {
    png(path, width = 1800, height = 1200, res = 200)
  }
  on.exit(dev.off())
  plot_prior_K_distribution(distribution)
  invisible(NULL)
}

# Una figura per eta (le precedenti restano nella cronologia del pannello Plots).
# Per rivederne una: plot_prior_K_distribution(prior_K_distributions[[1]])
for (distribution in prior_K_distributions) {
  plot_prior_K_distribution(distribution)
  if (save_figures) {
    figure_dir <- file.path(.plot_dir, "results", "figures")
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
    figure_name <- paste0("PriorK_a_theta_", prior_number_label(a_theta),
                          "_b_theta_", prior_number_label(b_theta),
                          "_eta_usr_", prior_number_label(distribution$eta_usr),
                          "_a_sigma_", prior_number_label(a_sigma),
                          "_b_sigma_", prior_number_label(b_sigma))
    save_prior_K_figure(file.path(figure_dir, paste0(figure_name, ".pdf")), "pdf", distribution)
    save_prior_K_figure(file.path(figure_dir, paste0(figure_name, ".png")), "png", distribution)
  }
}
