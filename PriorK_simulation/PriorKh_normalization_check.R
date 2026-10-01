# Controllo riga per riga con Ctrl+Invio, senza normalizzare le probabilita'.
# Working directory: cartella del progetto oppure PriorK_simulation.
script_dir <- if (file.exists("log_C_PY.cpp")) "." else "PriorK_simulation"
Rcpp::sourceCpp(file.path(script_dir, "log_C_PY.cpp"))
Rcpp::sourceCpp(file.path(script_dir, "compute_logC_bridge.cpp"))

p <- 40L
theta <- 1
sigma <- 0.5
stopifnot(p == 40L, sigma > 0, sigma < 1, theta > -sigma)
k <- seq_len(p)

# Prefattore esattamente come nel riferimento, p.774.
log_prefactor <- c(0, cumsum(log(theta + seq_len(p - 1L) * sigma))) -
  k * log(sigma) - sum(log(theta + seq_len(p - 1L)))
# C_mathcal del riferimento: nessun (-1)^k o (-1)^p da aggiungere.
log_C <- log_C_PY(p, sigma)[k + 1L]
probability <- exp(log_prefactor + log_C)
# Confronto con il vecchio C++, mantenuto intatto.
old_values <- exp(log_prefactor + original_logC(p, -sigma, 0)[k + 1L])

# Formula (11), in forma numericamente stabile.
mean_theoretical <- 1 + (theta + sigma) / sigma *
  expm1(sum(log1p(sigma / (theta + seq_len(p - 1L)))))
values <- data.frame(p=p, theta=theta, sigma=sigma, k=k, log_C=log_C,
                     probability=probability, old_value_invalid=old_values)
check <- data.frame(p=p, theta=theta, sigma=sigma,
                    total_mass=sum(probability),
                    min_probability=min(probability), max_probability=max(probability),
                    mean_numerical=sum(k*probability), mean_theoretical=mean_theoretical,
                    mean_error=sum(k*probability)-mean_theoretical,
                    old_total_mass=sum(old_values))
out <- file.path(script_dir, "results")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
write.csv(values, file.path(out, "PriorKh_normalization_values.csv"), row.names=FALSE)
write.csv(check, file.path(out, "PriorKh_normalization_summary.csv"), row.names=FALSE)
print(values, digits=12, row.names=FALSE)
print(check, digits=15, row.names=FALSE)
stopifnot(all(probability >= 0 & probability <= 1),
          abs(sum(probability)-1) < 1e-10,
          abs(sum(k*probability)-mean_theoretical) < 1e-10)
