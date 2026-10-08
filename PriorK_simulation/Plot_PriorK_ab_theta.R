# Distribuzione prior marginale di K: eseguire dall'alto con Ctrl+Invio.
# Working directory: FGMRandomBlocks oppure PriorK_simulation.
# Nessuna tabella o estrazione salvata: solo figure PDF e PNG.

# ==================== CONFIGURAZIONE UTENTE ==========================
a_theta <- 0.88
b_theta <- 0.94                 # RATE Gamma, non scale
a_sigma <- 1                   # shape1 Beta
b_sigma <- 1                   # shape2 Beta
eta_usr <- 0.5                  # Probabilita' dei changepoint esperti interni
B <- 10000L
seed <- 20261001L
# sigma ~ Beta(a_sigma,b_sigma), p=40; l'ultimo nodo ha sempre eta[40]=1.
save_figures <- FALSE
# ====================================================================

.plot_dir <- if (file.exists("PriorK_helpers.R")) "." else "PriorK_simulation"
source(file.path(.plot_dir, "PriorK_helpers.R"), local = TRUE)
source(file.path(.plot_dir, "PriorKh_test.R"), local = TRUE)
eta <- make_prior_eta(eta_usr)
set.seed(seed)
K <- draw_K_hyperprior(a_theta, b_theta, eta, B, a_sigma = a_sigma, b_sigma = b_sigma)

# K e' discreto: si rappresentano frequenze relative, senza KDE continua.
k <- seq_len(40L)
probability <- tabulate(K, nbins = 40L) / B
K_mean <- mean(K)
K_sd <- sd(K)
K_quantiles <- quantile(K, c(0.025, 0.975), type = 1, names = FALSE)
minimum_K <- sum(eta == 1)

plot_prior_K_distribution <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(5, 5, 6, 1), las = 1)
  plot(k, probability, type = "n", xlim = c(minimum_K - 0.5, 40.5),
       ylim = c(0, max(probability) * 1.16), xaxs = "i", yaxs = "i",
       xlab = "Numero totale di cluster K", ylab = "Probabilita' stimata P(K = k)",
       xaxt = "n", bty = "l")
  abline(h = axTicks(2), col = "grey90", lwd = 0.6)
  rect(k - 0.38, 0, k + 0.38, probability, col = "#327FAD", border = NA)
  ticks <- sort(unique(c(minimum_K, seq(5, 40, by = 5))))
  axis(1, at = ticks[ticks >= minimum_K])
  abline(v = K_mean, col = "#C57422", lwd = 2, lty = 2)
  title(main = "Distribuzione prior marginale di K", line = 4)
  # Testo ASCII per PDF portabili anche senza il font matematico Symbol.
  mtext(sprintf("a_theta = %g; b_theta = %g (rate); eta_usr = %g; sigma ~ Beta(%g,%g)",
                a_theta, b_theta, eta_usr, a_sigma, b_sigma), side = 3, line = 2.6, cex = 0.9)
  mtext(sprintf("B = %s | media (sd) = %.2f (%.2f) | quantili 2.5%% / 97.5%% = %d / %d",
                format(B, big.mark = " ", scientific = FALSE), K_mean, K_sd,
                K_quantiles[1], K_quantiles[2]), side = 3, line = 1.3, cex = 0.85)
  mtext("p = 40; eta[40] = 1. Linea tratteggiata: media di K.", side = 1, line = 3.6, cex = 0.8)
}

# Eseguire questa riga con Ctrl+Invio per mostrare il grafico.
plot_prior_K_distribution()

if (save_figures) {
  figure_dir <- file.path(.plot_dir, "results", "figures")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  figure_name <- paste0("PriorK_a_theta_", prior_number_label(a_theta),
                        "_b_theta_", prior_number_label(b_theta),
                        "_eta_usr_", prior_number_label(eta_usr),
                        "_a_sigma_", prior_number_label(a_sigma),
                        "_b_sigma_", prior_number_label(b_sigma))
  save_prior_K_figure <- function(path, format) {
    if (format == "pdf") {
      pdf(path, width = 9, height = 6, useDingbats = FALSE)
    } else {
      png(path, width = 1800, height = 1200, res = 200)
    }
    on.exit(dev.off())
    plot_prior_K_distribution()
    invisible(NULL)
  }
  save_prior_K_figure(file.path(figure_dir, paste0(figure_name, ".pdf")), "pdf")
  save_prior_K_figure(file.path(figure_dir, paste0(figure_name, ".png")), "png")
}
