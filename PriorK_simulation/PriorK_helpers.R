# Funzioni condivise: nessuna simulazione, stampa o scrittura al source().
read_prior_options <- function(args = character(), eta_default = 1) {
  opt <- list(choose_wd = NULL, eta_usr = eta_default)
  keys <- sub("=.*$", "", args)
  if (any(!grepl("=", args, fixed = TRUE)) || anyDuplicated(keys) ||
      any(!keys %in% names(opt)))
    stop("Opzioni ammesse, una volta ciascuna: choose_wd=1|2 eta_usr=VALORE_IN_[0,1].")
  for (i in seq_along(args)) opt[[keys[i]]] <- sub("^[^=]*=", "", args[i])
  opt$eta_usr <- suppressWarnings(as.numeric(opt$eta_usr))
  if (length(opt$eta_usr) != 1L || !is.finite(opt$eta_usr) ||
      opt$eta_usr < 0 || opt$eta_usr > 1)
    stop("eta_usr deve essere un numero in [0,1].")
  if (!is.null(opt$choose_wd)) {
    if (!opt$choose_wd %in% c("1", "2")) stop("choose_wd deve essere 1 oppure 2.")
    wd <- c("C:/Users/colom/FGMRandomBlocks",
            "/home/colombi/FGMRandomBlocks")[as.integer(opt$choose_wd)]
    if (!dir.exists(wd)) stop("Directory non trovata: ", wd)
    setwd(wd)
  }
  opt
}

make_prior_eta <- function(eta_usr) {
  if (length(eta_usr) != 1L || !is.finite(eta_usr) || eta_usr < 0 || eta_usr > 1)
    stop("eta_usr deve essere in [0,1].")
  eta <- numeric(40L)
  eta[c(4, 6, 9, 13, 18, 22, 28, 33)] <- eta_usr
  eta[40] <- 1  # Vincolo strutturale del paper: chiusura dell'ultimo sottogruppo.
  eta
}

check_prior_simulation <- function(eta, B) {
  if (length(B) != 1L || !is.finite(B) || B < 2 || B != floor(B))
    stop("B deve essere un intero >=2.")
  if (length(eta) != 40L || any(!is.finite(eta)) || any(eta < 0 | eta > 1) || eta[40] != 1)
    stop("eta deve contenere 40 probabilita' in [0,1], con eta[40]=1.")
}

# Una estrazione di K; una coppia theta/sigma condivisa da tutti i sottogruppi.
draw_K_given_parameters <- function(theta, sigma, eta, pmfs = NULL) {
  gamma <- rbinom(40L, size = 1L, prob = eta)
  sizes <- diff(c(0L, which(gamma == 1L)))
  if (is.null(pmfs)) {
    pmfs <- vector("list", 40L)
    for (n in unique(sizes)) pmfs[[n]] <- prior_Kh_pmf(n, theta, sigma)$probability
  }
  Kh <- vapply(sizes, function(n) sample.int(n, 1L, prob = pmfs[[n]]), integer(1))
  K <- sum(Kh)
  stopifnot(sum(sizes) == 40L, K >= length(sizes), K <= 40L)
  K
}

# Riutilizzata sia dalla tabella gerarchica sia dal grafico.
# Il chiamante imposta il seme; questa funzione non stampa e non salva.
draw_K_hyperprior <- function(a_theta, b_theta, eta, B, a_sigma = 1, b_sigma = 1) {
  check_prior_simulation(eta, B)
  if (length(a_theta) != 1L || length(b_theta) != 1L ||
      !is.finite(a_theta) || !is.finite(b_theta) || a_theta <= 0 || b_theta <= 0)
    stop("a_theta e b_theta devono essere positivi e finiti (b_theta = rate).")
  if (length(a_sigma) != 1L || length(b_sigma) != 1L ||
      !is.finite(a_sigma) || !is.finite(b_sigma) || a_sigma <= 0 || b_sigma <= 0)
    stop("a_sigma e b_sigma devono essere positivi e finiti.")
  load_logC()
  K <- integer(B)
  for (b in seq_len(B)) {
    sigma <- rbeta(1L, shape1 = a_sigma, shape2 = b_sigma)
    U <- rgamma(1L, shape = a_theta, rate = b_theta)
    theta <- U - sigma
    if (!is.finite(theta) || sigma <= 0 || sigma >= 1 || U <= 0 || theta <= -sigma)
      stop(sprintf("Precisione numerica insufficiente per theta>-sigma alla replica %d.", b))
    K[b] <- draw_K_given_parameters(theta, sigma, eta)
  }
  K
}

prior_number_label <- function(x) {
  # Etichette leggibili: evita le code binarie come 0.9399999999999999.
  vapply(x, function(z) format(z, digits = 15, trim = TRUE, decimal.mark = "."), character(1))
}
