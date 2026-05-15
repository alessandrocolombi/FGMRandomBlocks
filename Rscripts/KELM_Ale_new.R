# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100,wd_bocconi)
choose_wd = wd_vec[4] # <--- modify here
wd = paste0(choose_wd,"FGMRandomBlocks/Rscripts/")
setwd(wd)

library(fda)
suppressWarnings(suppressPackageStartupMessages(library(parallel)))
suppressWarnings(suppressPackageStartupMessages(library(doSNOW)))
# Functions ---------------------------------------------------------------
rbf_kernel_fast <- function(X1, X2, gamma = 1) {
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  p  <- ncol(X1)
  
  # Squared norms of each row
  x1_sq <- rowSums(X1^2)      # length n1
  x2_sq <- rowSums(X2^2)      # length n2
  
  # Compute squared Euclidean distances using (x - y)^2 = x^2 - 2x'y + y^2
  D2 <- matrix(x1_sq, n1, n2, byrow = FALSE) -
    2 * tcrossprod(X1, X2) +
    matrix(x2_sq, n1, n2, byrow = TRUE)
  
  # Avoid numerical issues
  D2[D2 < 0] <- 0
  
  # Apply RBF kernel
  K <- exp(-gamma * D2)
  return(K)
}

mahalanobis_pairwise <- function(X, Y, K) {
  # K is the precision (inverse covariance) matrix
  p <- ncol(X)
  nX <- nrow(X)
  nY <- nrow(Y)
  
  xK <- X %*% K
  yK <- Y %*% K
  
  # (x_i - y_j)' K (x_i - y_j) = 
  # x_i' K x_i + y_j' K y_j - 2 x_i' K y_j
  xKx <- rowSums(X * xK)      # x_i' K x_i
  yKy <- rowSums(Y * yK)      # y_j' K y_j
  xKy <- xK %*% t(Y)          # x_i' K y_j
  
  D2 <- matrix(xKx, nX, nY, byrow = FALSE) +
    matrix(yKy, nX, nY, byrow = TRUE) -
    2 * xKy
  
  # numerics
  D2[D2 < 0] <- 0
  return(D2)
}

custom_kernel_fast <- function(X1, X2,
                               y1, y2,
                               K_s, K_n,
                               gamma = 1) {
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  p  <- ncol(X1)
  
  # indices
  s1 <- which(y1 ==  1)    # strawberry in X1
  n1_ <- which(y1 == -1)   # non in X1
  s2 <- which(y2 ==  1)    # strawberry in X2
  n2_ <- which(y2 == -1)   # non in X2
  
  # preallocate
  K <- matrix(0, n1, n2)
  
  # 1. sâ€“s block: X1[s1,] vs X2[s2,] with Mahalanobis K_s
  if (length(s1) > 0 && length(s2) > 0) {
    d1_s <- X1[s1, , drop = FALSE]
    d2_s <- X2[s2, , drop = FALSE]
    # compute pairwise squared Mahalanobis distances
    D2_s <- mahalanobis_pairwise(d1_s, d2_s, K_s)
    K_ss <- exp(-gamma * D2_s)
    K[s1, s2] <- K_ss
  }
  
  # 2. nâ€“n block: X1[n1_,] vs X2[n2_,] with Mahalanobis K_n
  if (length(n1_) > 0 && length(n2_) > 0) {
    d1_n <- X1[n1_, , drop = FALSE]
    d2_n <- X2[n2_, , drop = FALSE]
    D2_n <- mahalanobis_pairwise(d1_n, d2_n, K_n)
    K_nn <- exp(-gamma * D2_n)
    K[n1_, n2_] <- K_nn
  }
  
  # 3. mixed blocks: sâ€“n and nâ€“s â†’ ordinary RBF (Euclidean)
  if (length(s1) > 0 && length(n2_) > 0) {
    K_sn <- rbf_kernel_fast(X1[s1, , drop = FALSE],
                            X2[n2_, , drop = FALSE],
                            gamma = gamma)
    K[s1, n2_] <- K_sn
  }
  
  if (length(n1_) > 0 && length(s2) > 0) {
    K_ns <- rbf_kernel_fast(X1[n1_, , drop = FALSE],
                            X2[s2, , drop = FALSE],
                            gamma = gamma)
    K[n1_, s2] <- K_ns
  }
  
  return(K)
}


custom_kernel_test_fast <- function(
    X_test,
    X_train,
    y_train,
    K_s,
    K_n,
    gamma = 1
){
  
  n_test  <- nrow(X_test)
  n_train <- nrow(X_train)
  
  Kmat <- matrix(0, n_test, n_train)
  
  idx_s <- which(y_train == 1)
  idx_n <- which(y_train == -1)
  
  # strawberry columns
  if(length(idx_s) > 0){
    
    D2_s <- mahalanobis_pairwise(
      X_test,
      X_train[idx_s, , drop = FALSE],
      K_s
    )
    
    Kmat[, idx_s] <- exp(-gamma * D2_s)
  }
  
  # non-strawberry columns
  if(length(idx_n) > 0){
    
    D2_n <- mahalanobis_pairwise(
      X_test,
      X_train[idx_n, , drop = FALSE],
      K_n
    )
    
    Kmat[, idx_n] <- exp(-gamma * D2_n)
  }
  
  return(Kmat)
}

fit_kelm <- function(
    X_train,
    y_train,
    T_train,
    X_test,
    y_test,
    K_s,
    K_n,
    gamma,
    C,
    zero_crossterms = F
) {
  
  # TRAIN KERNEL
  Omega <- custom_kernel_fast(
    X_train,
    X_train,
    y_train,
    y_train,
    K_s,
    K_n,
    gamma = gamma
  )
  
  if(zero_crossterms){
    
    idx_s <- which(y_train == 1)
    idx_n <- which(y_train == -1)
    
    Omega[idx_s, idx_n] <- 0
    Omega[idx_n, idx_s] <- 0
  }
  
  A <- Omega + diag(1 / C, nrow(Omega))
  B <- solve(A, T_train)
  
  # TEST KERNEL
  K_test <- custom_kernel_test_fast(
    X_test,
    X_train,
    y_train,
    K_s,
    K_n,
    gamma = gamma
  )
  
  # PREDICTION
  f_test <- K_test %*% B
  
  pred_idx <- max.col(f_test)
  
  y_pred <- ifelse(pred_idx == 1, 1, -1)
  
  acc <- mean(y_pred == y_test)
  
  conf <- table(
    Predicted = y_pred,
    True = y_test
  )
  
  return(list(
    accuracy = acc,
    confusion = conf,
    y_pred = y_pred
  ))
}

grid_search_kelm <- function(
    X_train,
    X_test,
    y_train,
    y_test,
    T_train,
    K_s,
    K_n,
    gamma_grid,
    C_grid,
    label,
    zero_crossterms = F
) {
  
  results <- expand.grid(
    gamma = gamma_grid,
    C = C_grid
  )
  
  results$test_accuracy <- NA_real_
  results$setting <- label
  results$confusion <- vector(
    "list",
    nrow(results)
  )
  
  for(i in 1:nrow(results)) {
    
    gamma_i <- results$gamma[i]
    C_i <- results$C[i]
    
    fit <- fit_kelm(
      X_train = X_train,
      y_train = y_train,
      T_train = T_train,
      X_test = X_test,
      y_test = y_test,
      K_s = K_s,
      K_n = K_n,
      gamma = gamma_i,
      C = C_i,
      zero_crossterms = zero_crossterms
    )
    
    results$test_accuracy[i] <- fit$accuracy
    
    results$confusion[[i]] <- fit$confusion
    
    cat(
      "[",
      label,
      "]",
      i,
      "/",
      nrow(results),
      "\n"
    )
  }
  
  results <- results[
    order(-results$test_accuracy),
    ,
    drop = FALSE
  ]
  
  return(results)
}


# Load data  ------------------------------------------------------
load("../data/purees.Rdat")

wavelengths = purees$wavelengths
curves = purees$data

strawberry <- curves[which(curves$Group == "Strawberry"), ]
NONstrawberry <- curves[which(curves$Group != "Strawberry"), ]

data_s = strawberry[,-1]*100
data_n = NONstrawberry[,-1]*100

data_s = as.matrix(data_s)    # n x r
data_n = as.matrix(data_n)    # n x r

n_s = dim(data_s)[1]
n_n = dim(data_n)[1]

r = dim(data_n)[2]

range_x = range(wavelengths)

data_W_s <- t(as.matrix(data_s))
data_W_n <- t(as.matrix(data_n))

# Settings ----------------------------------------------------------------
test_obs_perclass <- 100
Ns_train <- 100  # n_s - test_obs_perclass
Nns_train <- 100 # n_n - test_obs_perclass

Nrep <- 100
seed0 <- 13212
set.seed(seed0)
seeds <- sample(1:999999, size = Nrep)

n_components_seq <- seq(4,40,3)
p_seq <- seq(4,40,3)

C_fixed <- 1e4
B <- 10
gamma_one <- TRUE
gamma_fixed <- 1
num_cores = 33 # <---

# Utilities ----------------------------------------------------------------
make_stratified_folds <- function(y, B) {
  folds <- integer(length(y))
  classes <- sort(unique(y))
  
  for (class_i in classes) {
    idx_class <- which(y == class_i)
    idx_class <- sample(idx_class)
    folds[idx_class] <- rep(seq_len(B), length.out = length(idx_class))
  }
  
  folds
}

confusion_2x2 <- function(y_pred, y_true) {
  table(
    Predicted = factor(y_pred, levels = c(-1, 1)),
    True = factor(y_true, levels = c(-1, 1))
  )
}

fit_rbf_kelm <- function(
    X_train,
    y_train,
    T_train,
    X_eval,
    y_eval,
    gamma,
    C
) {
  
  Omega <- rbf_kernel_fast(
    X_train,
    X_train,
    gamma = gamma
  )
  
  idx_s_train <- which(y_train == 1)
  idx_n_train <- which(y_train == -1)
  
  Omega[idx_s_train, idx_n_train] <- 0
  Omega[idx_n_train, idx_s_train] <- 0
  
  A <- Omega + diag(1 / C, nrow(Omega))
  B_mat <- solve(A, T_train)
  
  K_eval <- rbf_kernel_fast(
    X_eval,
    X_train,
    gamma = gamma
  )
  
  f_eval <- K_eval %*% B_mat
  pred_idx <- max.col(f_eval)
  y_pred <- ifelse(pred_idx == 1, 1, -1)
  acc <- mean(y_pred == y_eval)
  conf <- confusion_2x2(y_pred, y_eval)
  
  list(
    accuracy = acc,
    confusion = conf,
    y_pred = y_pred
  )
}

run_replication <- function(ii) {
  
  seed_ii <- seeds[ii]
  
  set.seed(seed_ii)
  idx_s_train <- sample(1:nrow(data_s), size = Ns_train, replace = FALSE)
  idx_n_train <- sample(1:nrow(data_n), size = Nns_train, replace = FALSE)
  
  idx_s_test <- setdiff(seq_len(nrow(data_s)), idx_s_train)
  idx_n_test <- setdiff(seq_len(nrow(data_n)), idx_n_train)
  
  data_s_train <- data_s[idx_s_train,]
  data_n_train <- data_n[idx_n_train,]
  data_s_test <- data_s[idx_s_test,]
  data_n_test <- data_n[idx_n_test,]
  
  y_train <- c(rep(1, nrow(data_s_train)),
               rep(-1, nrow(data_n_train)))
  
  y_test <- c(rep(1, nrow(data_s_test)),
              rep(-1, nrow(data_n_test)))
  
  T_train <- matrix(c(y_train, -y_train), nrow = length(y_train), ncol = 2)
  
  set.seed(seed_ii + 1000000)
  fold_id <- make_stratified_folds(y_train, B)
  
  rbf_results <- vector("list", length(p_seq))
  custom_results <- vector("list", length(p_seq))
  cv_results <- vector("list", length(p_seq) * B)
  counter_cv <- 1
  
  for(pp in seq_along(p_seq)){
    
    p <- p_seq[pp]
    
    basis <- create.bspline.basis(rangeval=range(wavelengths), nbasis=p, norder=3)
    
    fd_s_train <- Data2fd(y = t(data_s_train), argvals = wavelengths, basisobj = basis)
    fd_n_train <- Data2fd(y = t(data_n_train), argvals = wavelengths, basisobj = basis)
    
    fd_s_test <- Data2fd(y = t(data_s_test), argvals = wavelengths, basisobj = basis)
    fd_n_test <- Data2fd(y = t(data_n_test), argvals = wavelengths, basisobj = basis)
    
    X_s_train <- t(fd_s_train$coefs)
    X_n_train <- t(fd_n_train$coefs)
    
    X_s_test <- t(fd_s_test$coefs)
    X_n_test <- t(fd_n_test$coefs)
    
    X_train <- rbind(X_s_train,X_n_train)
    X_test <- rbind(X_s_test,X_n_test)
    
    train_center <- colMeans(X_train)
    train_scale <- apply(X_train, 2, sd)
    
    X_train_sc <- scale(X_train, center = train_center, scale = train_scale)
    X_test_sc <- scale(X_test, center = train_center, scale = train_scale)
    
    # Splines RBF Kernel with CV gamma selection ---------------------------
    D2_vals_sc <- as.vector(dist(X_train_sc)^2)
    d_med_sc <- median(D2_vals_sc)
    gamma_ref_rbf <- 1 / d_med_sc
    gamma_grid <- gamma_ref_rbf * 10^seq(-1.5, 1.8, length.out = 12)
    
    selected_gamma <- numeric(B)
    selected_accuracy <- numeric(B)
    
    for(b in seq_len(B)){
      
      idx_valid <- which(fold_id == b)
      idx_train_cv <- which(fold_id != b)
      
      X_train_cv <- X_train_sc[idx_train_cv, , drop = FALSE]
      X_valid_cv <- X_train_sc[idx_valid, , drop = FALSE]
      
      y_train_cv <- y_train[idx_train_cv]
      y_valid_cv <- y_train[idx_valid]
      
      T_train_cv <- matrix(
        c(y_train_cv, -y_train_cv),
        nrow = length(y_train_cv),
        ncol = 2
      )
      
      acc_gamma <- numeric(length(gamma_grid))
      
      for(g in seq_along(gamma_grid)){
        
        fit_cv <- fit_rbf_kelm(
          X_train = X_train_cv,
          y_train = y_train_cv,
          T_train = T_train_cv,
          X_eval = X_valid_cv,
          y_eval = y_valid_cv,
          gamma = gamma_grid[g],
          C = C_fixed
        )
        
        acc_gamma[g] <- fit_cv$accuracy
      }
      
      best_g_idx <- which(acc_gamma == max(acc_gamma))[1]
      
      selected_gamma[b] <- gamma_grid[best_g_idx]
      selected_accuracy[b] <- acc_gamma[best_g_idx]
      
      cv_results[[counter_cv]] <- data.frame(
        rep = ii,
        seed = seed_ii,
        p_basis = p,
        fold = b,
        selected_gamma = selected_gamma[b],
        selected_accuracy = selected_accuracy[b]
      )
      
      counter_cv <- counter_cv + 1
    }
    
    gamma_hat <- mean(selected_gamma)
    
    final_fit_rbf <- fit_rbf_kelm(
      X_train = X_train_sc,
      y_train = y_train,
      T_train = T_train,
      X_eval = X_test_sc,
      y_eval = y_test,
      gamma = gamma_hat,
      C = C_fixed
    )
    
    conf_rbf <- final_fit_rbf$confusion
    
    rbf_results[[pp]] <- data.frame(
      rep = ii,
      seed = seed_ii,
      p_basis = p,
      gamma = gamma_hat,
      C = C_fixed,
      test_accuracy = final_fit_rbf$accuracy,
      cv_mean_gamma = mean(selected_gamma),
      cv_sd_gamma = sd(selected_gamma),
      cv_mean_accuracy = mean(selected_accuracy),
      cv_sd_accuracy = sd(selected_accuracy),
      TN = as.integer(conf_rbf["-1", "-1"]),
      FP = as.integer(conf_rbf["1",  "-1"]),
      FN = as.integer(conf_rbf["-1", "1"]),
      TP = as.integer(conf_rbf["1",  "1"])
    )
    
    # Splines Custom Kernel ------------------------------------------------
    Sigma_s <- cov(t(fd_s_train$coefs)) + diag(1e-4,p)
    Sigma_n <- cov(t(fd_n_train$coefs)) + diag(1e-4,p)
    
    K_s <- solve(Sigma_s)
    K_n <- solve(Sigma_n)
    
    Dmat <- diag(train_scale)
    
    K_s_sc <- Dmat %*% K_s %*% Dmat
    K_n_sc <- Dmat %*% K_n %*% Dmat
    
    final_fit_custom <- fit_kelm(
      X_train = X_train_sc,
      y_train = y_train,
      T_train = T_train,
      X_test = X_test_sc,
      y_test = y_test,
      K_s = K_s_sc,
      K_n = K_n_sc,
      gamma = gamma_fixed,
      C = C_fixed,
      zero_crossterms = TRUE
    )
    
    conf_custom <- confusion_2x2(final_fit_custom$y_pred, y_test)
    
    custom_results[[pp]] <- data.frame(
      rep = ii,
      seed = seed_ii,
      p_basis = p,
      gamma = gamma_fixed,
      gamma_one = gamma_one,
      C = C_fixed,
      test_accuracy = final_fit_custom$accuracy,
      TN = as.integer(conf_custom["-1", "-1"]),
      FP = as.integer(conf_custom["1",  "-1"]),
      FN = as.integer(conf_custom["-1", "1"]),
      TP = as.integer(conf_custom["1",  "1"])
    )
  }
  
  list(
    rbf = do.call(rbind, rbf_results),
    custom = do.call(rbind, custom_results),
    cv = do.call(rbind, cv_results)
  )
}

cl <- parallel::makeCluster(num_cores)
registerDoSNOW(cl)

parallel::clusterExport(
  cl,
  varlist = c(
    "seeds",
    "data_s",
    "data_n",
    "wavelengths",
    "Ns_train",
    "Nns_train",
    "p_seq",
    "B",
    "C_fixed",
    "gamma_one",
    "gamma_fixed",
    "make_stratified_folds",
    "confusion_2x2",
    "rbf_kernel_fast",
    "fit_rbf_kelm",
    "mahalanobis_pairwise",
    "custom_kernel_fast",
    "custom_kernel_test_fast",
    "fit_kelm",
    "run_replication"
  )
)

pb <- txtProgressBar(max = Nrep, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

replication_results <- tryCatch(
  {
    foreach(
      ii = seq_len(Nrep),
      .packages = c("fda"),
      .options.snow = opts,
      .inorder = TRUE
    ) %dopar% {
      run_replication(ii)
    }
  },
  finally = {
    close(pb)
    parallel::stopCluster(cl)
  }
)

# FINAL RESULTS -----------------------------------------------------------
results_all <- do.call(
  rbind,
  lapply(replication_results, function(x) x$rbf)
)
rownames(results_all) <- NULL

cv_gamma_by_fold <- do.call(
  rbind,
  lapply(replication_results, function(x) x$cv)
)
rownames(cv_gamma_by_fold) <- NULL

results_custom_kernel <- do.call(
  rbind,
  lapply(replication_results, function(x) x$custom)
)
rownames(results_custom_kernel) <- NULL

accuracy_curve_rbf_cv <- aggregate(
  test_accuracy ~ p_basis,
  data = results_all,
  FUN = function(x) c(mean = mean(x), sd = sd(x))
)
accuracy_curve_rbf_cv <- data.frame(
  p_basis = accuracy_curve_rbf_cv$p_basis,
  mean_accuracy = accuracy_curve_rbf_cv$test_accuracy[, "mean"],
  sd_accuracy = accuracy_curve_rbf_cv$test_accuracy[, "sd"]
)
rownames(accuracy_curve_rbf_cv) <- NULL

accuracy_curve_custom <- aggregate(
  test_accuracy ~ p_basis,
  data = results_custom_kernel,
  FUN = function(x) c(mean = mean(x), sd = sd(x))
)
accuracy_curve_custom <- data.frame(
  p_basis = accuracy_curve_custom$p_basis,
  mean_accuracy = accuracy_curve_custom$test_accuracy[, "mean"],
  sd_accuracy = accuracy_curve_custom$test_accuracy[, "sd"]
)
rownames(accuracy_curve_custom) <- NULL

# Save results -------------------------------------------------------------
output_dir <- file.path(wd, "save")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_file <- file.path(output_dir, "KELM_Ale_new_results.RData")
output_rds <- file.path(output_dir, "KELM_Ale_new_results.rds")

save(
  results_all,
  cv_gamma_by_fold,
  results_custom_kernel,
  accuracy_curve_rbf_cv,
  accuracy_curve_custom,
  replication_results,
  seeds,
  seed0,
  Nrep,
  test_obs_perclass,
  Ns_train,
  Nns_train,
  p_seq,
  C_fixed,
  B,
  gamma_one,
  gamma_fixed,
  file = output_file
)

saveRDS(
  list(
    results_all = results_all,
    cv_gamma_by_fold = cv_gamma_by_fold,
    results_custom_kernel = results_custom_kernel,
    accuracy_curve_rbf_cv = accuracy_curve_rbf_cv,
    accuracy_curve_custom = accuracy_curve_custom,
    replication_results = replication_results,
    seeds = seeds,
    seed0 = seed0,
    Nrep = Nrep,
    test_obs_perclass = test_obs_perclass,
    Ns_train = Ns_train,
    Nns_train = Nns_train,
    p_seq = p_seq,
    C_fixed = C_fixed,
    B = B,
    gamma_one = gamma_one,
    gamma_fixed = gamma_fixed
  ),
  file = output_rds
)

cat("\nSaved results to:\n")
cat(normalizePath(output_file, winslash = "/", mustWork = FALSE), "\n")
cat(normalizePath(output_rds, winslash = "/", mustWork = FALSE), "\n")
