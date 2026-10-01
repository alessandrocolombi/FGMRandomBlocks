# PMF di K_h (Lijoi et al. 2007, p.774) e media (11) del Draft.
# Ctrl+Invio/source/Rscript dalla cartella del progetto o PriorK_simulation.
.prior_dir <- if (file.exists("log_C_PY.cpp")) "." else "PriorK_simulation"

check_py_parameters <- function(n, theta, sigma) {
  if (length(n) != 1L || !is.finite(n) || n < 1 || n > 40 || n != floor(n))
    stop("n deve essere intero tra 1 e 40.")
  if (length(sigma) != 1L || !is.finite(sigma) || sigma < 0 || sigma >= 1)
    stop("Richiesto 0 <= sigma < 1.")
  if (length(theta) != 1L || !is.finite(theta) || theta <= -sigma)
    stop("Richiesto theta > -sigma (theta > 0 se sigma=0).")
}

load_logC <- function() {
  if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Installare Rcpp.")
  if (!exists("log_C_PY", envir = environment(load_logC), inherits = FALSE))
    Rcpp::sourceCpp(file.path(.prior_dir, "log_C_PY.cpp"), env = environment(load_logC))
}

prior_Kh_mean <- function(n, theta, sigma) {
  check_py_parameters(n, theta, sigma)
  if (sigma == 0) return(sum(theta / (theta + seq_len(n) - 1)))
  # Formula (11), senza sottrazione di due termini grandi; valida anche theta=0.
  1 + (theta + sigma) *
    (expm1(sum(log1p(sigma / (theta + seq_len(n - 1L))))) / sigma)
}

prior_Kh_logweights <- function(n, theta, sigma) {
  check_py_parameters(n, theta, sigma)
  load_logC()
  # S=C_mathcal/sigma^k: evita cancellazione per sigma vicino a zero.
  log_C_scaled_PY(n, sigma)[seq_len(n) + 1L] +
    c(0, cumsum(log(theta + seq_len(n - 1L) * sigma))) -
    sum(log(theta + seq_len(n - 1L)))
}

prior_Kh_pmf <- function(n, theta, sigma, tolerance = 1e-9) {
  prob <- exp(prior_Kh_logweights(n, theta, sigma))
  target <- prior_Kh_mean(n, theta, sigma)
  # Controlli sulla formula grezza: NON dividere per sum(prob).
  if (any(!is.finite(prob)) || any(prob < 0 | prob > 1 + tolerance) ||
      abs(sum(prob) - 1) > tolerance ||
      abs(sum(seq_len(n) * prob) - target) > tolerance * max(1, target))
    stop("PMF non valida: controllo massa, intervallo o media fallito.")
  data.frame(k = seq_len(n), probability = prob)
}

# Verifica indipendente tramite le probabilita' predittive del Pitman-Yor.
# Non usata per i coefficienti o per la simulazione.
prior_Kh_predictive_check <- function(n, theta, sigma) {
  prob <- 1
  if (n > 1) for (m in seq_len(n - 1L)) {
    k <- seq_len(m)
    next_prob <- numeric(m + 1L)
    next_prob[k] <- prob * (m - sigma * k) / (theta + m)
    next_prob[k + 1L] <- next_prob[k + 1L] + prob * (theta + sigma * k) / (theta + m)
    prob <- next_prob
  }
  prob
}

run_prior_Kh_tests <- function(output_dir = file.path(.prior_dir, "results")) {
  p <- 40L
  test_grid <- expand.grid(theta = c(0, 1, 5), sigma = c(0.2, 0.5, 0.8))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  load_logC()
  distributions <- list()
  checks <- lapply(seq_len(nrow(test_grid)), function(i) {
    theta <- test_grid$theta[i]; sigma <- test_grid$sigma[i]
    prob <- prior_Kh_pmf(p, theta, sigma)$probability
    target <- prior_Kh_mean(p, theta, sigma)
    ref <- prior_Kh_predictive_check(p, theta, sigma)
    distributions[[i]] <<- data.frame(p=p, theta=theta, sigma=sigma,
                                      k=seq_len(p), probability=prob)
    data.frame(p=p, theta=theta, sigma=sigma, raw_mass=sum(prob),
               min_probability=min(prob), max_probability=max(prob),
               numerical_mean=sum(seq_len(p)*prob), theoretical_mean=target,
               mean_error=sum(seq_len(p)*prob)-target,
               max_predictive_error=max(abs(prob-ref)),
               passed=max(abs(prob-ref)) < 1e-9)
  })
  checks <- do.call(rbind, checks)
  write.csv(checks, file.path(output_dir, "PriorKh_validation.csv"), row.names=FALSE)
  write.csv(do.call(rbind, distributions), file.path(output_dir, "PriorKh_probabilities.csv"), row.names=FALSE)

  # Formula esplicita dell'appendice per sigma=1/2, tutti gli n<=40 e k<=n.
  half_error <- 0
  for (n in seq_len(p)) {
    k <- seq_len(n)
    explicit <- (k-2*n)*log(2) + lchoose(2*n-k-1, n-1) + lgamma(n)-lgamma(k)
    half_error <- max(half_error, abs(log_C_PY(n, 0.5)[k+1L] - explicit))
  }
  stopifnot(half_error < 1e-10, all(checks$passed))

  # Bordi del dominio, sigma=0 e theta negativi ammessi, n=1,...,40.
  edge_grid <- data.frame(theta=c(1, 0.1, 1, -0.49, -0.499999, 0, 5),
                          sigma=c(0, 0, 1e-10, 0.5, 0.5, 0.999999, 0.8))
  edge_error <- 0
  for (i in seq_len(nrow(edge_grid))) for (n in seq_len(p)) {
    th <- edge_grid$theta[i]; si <- edge_grid$sigma[i]
    pr <- prior_Kh_pmf(n, th, si)$probability
    edge_error <- max(edge_error, abs(pr-prior_Kh_predictive_check(n, th, si)))
  }
  stopifnot(edge_error < 1e-9)
  writeLines(c(sprintf("Massimo errore log C vs formula esplicita sigma=1/2: %.16g", half_error),
               sprintf("Massimo errore PMF vs predittiva sui casi limite: %.16g", edge_error),
               "Nessuna PMF e' stata normalizzata artificialmente."),
             file.path(output_dir, "PriorKh_additional_checks.txt"))
  print(checks, digits=15, row.names=FALSE)
  invisible(checks)
}

if (sys.nframe() == 0L) run_prior_Kh_tests()
