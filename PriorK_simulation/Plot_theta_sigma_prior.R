# Prior: sigma ~ Beta(a_sigma,b_sigma), theta+sigma | sigma ~ Gamma(a_theta,b_theta).
# Eseguibile con Ctrl+Invio, senza pacchetti aggiuntivi.
# Working directory: FGMRandomBlocks oppure PriorK_simulation.
# b_theta e' RATE (come rgamma(n, alpha, beta) in utility_functions.R).
# Se b_theta rappresenta invece SCALE, impostare gamma_parameterization="scale".

# ==================== CONFIGURAZIONE UTENTE ==========================
# Valori di esempio: sostituire con gli iperparametri desiderati.
a_theta <- 0.88
b_theta <- 0.94
a_sigma <- 1
b_sigma <- 1
gamma_parameterization <- "rate"  # "rate" oppure "scale"
N <- 1000L                     # numero di estrazioni Monte Carlo
seed <- 20261001L
sigma_conditional <- c(0.1, 0.5, 0.9)
upper_plot_quantile <- 0.995      # limite destro dei grafici di theta
save_plot <- FALSE
script_dir <- if (dir.exists("PriorK_simulation")) "PriorK_simulation" else "."
output_dir <- file.path(script_dir, "results", "theta_sigma_prior")
# ====================================================================

stopifnot(all(is.finite(c(a_theta, b_theta, a_sigma, b_sigma))),
          all(c(a_theta, b_theta, a_sigma, b_sigma) > 0),
          length(N) == 1L, is.finite(N), N >= 1000, N == floor(N),
          gamma_parameterization %in% c("rate", "scale"),
          length(sigma_conditional) > 0L, all(is.finite(sigma_conditional)),
          all(sigma_conditional > 0 & sigma_conditional < 1),
          upper_plot_quantile > 0.5, upper_plot_quantile < 1)
gamma_rate <- if (gamma_parameterization == "rate") b_theta else 1 / b_theta

# Campionamento gerarchico: U e sigma indipendenti, theta=U-sigma.
set.seed(seed)
sigma_draws <- rbeta(N, shape1 = a_sigma, shape2 = b_sigma)
U_draws <- rgamma(N, shape = a_theta, rate = gamma_rate)
theta_draws <- U_draws - sigma_draws
prior_draws <- data.frame(theta = theta_draws, sigma = sigma_draws)

# Il supporto e' 0<sigma<1, theta>-sigma; theta puo' essere negativa.
# Nessuna troncatura delle estrazioni. Solo i grafici limitano la coda destra.
theta_limits <- c(-1, max(0.1, unname(quantile(theta_draws, upper_plot_quantile))))
theta_density <- density(theta_draws, from = theta_limits[1],
                         to = theta_limits[2], n = 1024)
# La marginale di theta e' una KDE Monte Carlo, non una densita' esatta.
# La KDE puo' avere bias vicino al bordo theta=-1.
sigma_axis <- seq(0.001, 0.999, length.out = 300)
theta_axis <- seq(theta_limits[1], theta_limits[2], length.out = 400)

# Densita' congiunta ESATTA: f(theta,sigma)=f_Gamma(theta+sigma)*f_Beta(sigma).
# Il cambio di variabili (U,sigma)->(theta,sigma) ha Jacobiano unitario.
joint_density <- outer(theta_axis, sigma_axis, function(theta, sigma) {
  u <- theta + sigma
  value <- numeric(length(u))
  inside <- u > 0
  value[inside] <- exp(dgamma(u[inside], shape = a_theta, rate = gamma_rate, log = TRUE) +
                        dbeta(sigma[inside], a_sigma, b_sigma, log = TRUE))
  value
})

plot_theta_sigma_prior <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(2, 2), mar = c(4.2, 4.4, 3, 1), oma = c(0, 0, 3, 0))

  plot(sigma_axis, dbeta(sigma_axis, a_sigma, b_sigma), type = "l", lwd = 2,
       col = "#246A9B", xlab = expression(sigma), ylab = "Densita'",
       main = "Marginale di sigma (esatta)", xlim = c(0, 1))

  plot(theta_density$x, theta_density$y, type = "l", lwd = 2,
       col = "#246A9B", xlab = expression(theta), ylab = "Densita'",
       main = "Marginale di theta (Monte Carlo)")
  abline(v = 0, lty = 3, col = "grey50")

  colors <- hcl.colors(length(sigma_conditional), "Dark 3")
  conditional_density <- vapply(sigma_conditional, function(sigma) {
    d <- dgamma(theta_axis + sigma, shape = a_theta, rate = gamma_rate)
    d[!is.finite(d)] <- NA_real_  # singolarita' al bordo se shape<1
    d
  }, numeric(length(theta_axis)))
  matplot(theta_axis, conditional_density, type = "l", lty = 1, lwd = 2,
          col = colors, xlab = expression(theta), ylab = "Densita'",
          main = "Theta | sigma (densita' esatte)")
  legend("topright", legend = paste0("sigma = ", sigma_conditional),
         col = colors, lty = 1, lwd = 2, bty = "n", cex = 0.85)

  # Contorni etichettati con i valori della densita', non log-densita'.
  image(theta_axis, sigma_axis, joint_density,
        col = hcl.colors(80, "YlOrRd", rev = TRUE),
        xlab = expression(theta), ylab = expression(sigma),
        main = "Densita' congiunta (esatta)")
  contour(theta_axis, sigma_axis, joint_density, add = TRUE,
          nlevels = 7, col = "grey25", labcex = 0.7)
  lines(-sigma_axis, sigma_axis, lty = 2, lwd = 2)

  mtext(sprintf("Gamma(shape=%g, %s=%g); Beta(%g, %g); N=%s",
                a_theta, gamma_parameterization, b_theta, a_sigma, b_sigma,
                format(N, big.mark = " ", scientific = FALSE)),
        outer = TRUE, line = 1, cex = 1)
}

# Grafico nel pannello Plots: eseguire questa riga con Ctrl+Invio.
plot_theta_sigma_prior()

# Controllo dei momenti: E(theta)=a_theta/rate - E(sigma).
mean_sigma <- a_sigma / (a_sigma + b_sigma)
var_sigma <- a_sigma*b_sigma / ((a_sigma+b_sigma)^2*(a_sigma+b_sigma+1))
moment_check <- data.frame(
  parameter = c("sigma", "theta"),
  mean_MC = c(mean(sigma_draws), mean(theta_draws)),
  mean_theoretical = c(mean_sigma, a_theta/gamma_rate - mean_sigma),
  sd_MC = c(sd(sigma_draws), sd(theta_draws)),
  sd_theoretical = c(sqrt(var_sigma), sqrt(a_theta/gamma_rate^2 + var_sigma)))
print(moment_check, row.names = FALSE)
cat(sprintf("Grafici di theta limitati a %.1f%% della coda cumulata Monte Carlo.\n",
            100*upper_plot_quantile))

if (save_plot) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  save_prior_png <- function(path) {
    png(path, width = 1800, height = 1400, res = 150)
    on.exit(dev.off())
    plot_theta_sigma_prior()
  }
  save_prior_png(file.path(output_dir, "theta_sigma_prior.png"))
  write.csv(moment_check, file.path(output_dir, "moment_check.csv"), row.names = FALSE)
  saveRDS(list(a_theta=a_theta, b_theta=b_theta, a_sigma=a_sigma, b_sigma=b_sigma,
               gamma_parameterization=gamma_parameterization, N=N, seed=seed,
               sigma_conditional=sigma_conditional,
               upper_plot_quantile=upper_plot_quantile),
          file.path(output_dir, "config.rds"))
}
