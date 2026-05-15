# wd ----------------------------------------------------------------------
wd_pc = "C:/Users/colom/"
wd_unicatt = "C:/Users/alessandro.colombi/"
wd_g100 = "/g100/home/userexternal/acolombi/"
wd_bocconi = "/home/colombi/"
wd_vec = c(wd_pc,wd_unicatt,wd_g100,wd_bocconi)
choose_wd = wd_vec[1] # <--- modify here
wd = paste0(choose_wd,"FGMRandomBlocks/Rscripts/")
setwd(wd)

# Load results ------------------------------------------------------------
results_file <- file.path(wd, "save", "KELM_Ale_new_results.RData")
load(results_file)

# Plot options ------------------------------------------------------------
save_img <- TRUE
img_dir <- file.path(wd, "img")
if (!dir.exists(img_dir)) {
  dir.create(img_dir, recursive = TRUE)
}

plot_accuracy_curves <- function() {
  ylim_acc <- range(
    c(
      accuracy_curve_rbf_cv$mean_accuracy - accuracy_curve_rbf_cv$sd_accuracy,
      accuracy_curve_rbf_cv$mean_accuracy + accuracy_curve_rbf_cv$sd_accuracy,
      accuracy_curve_custom$mean_accuracy - accuracy_curve_custom$sd_accuracy,
      accuracy_curve_custom$mean_accuracy + accuracy_curve_custom$sd_accuracy
    ),
    na.rm = TRUE
  )
  
  plot(
    x = accuracy_curve_rbf_cv$p_basis,
    y = accuracy_curve_rbf_cv$mean_accuracy,
    type = "b",
    pch = 16,
    lwd = 3,
    ylim = ylim_acc,
    xlab = "Dimension",
    ylab = "Mean accuracy",
    main = "Mean Test Accuracy over Repetitions"
  )
  
  arrows(
    x0 = accuracy_curve_rbf_cv$p_basis,
    y0 = accuracy_curve_rbf_cv$mean_accuracy - accuracy_curve_rbf_cv$sd_accuracy,
    x1 = accuracy_curve_rbf_cv$p_basis,
    y1 = accuracy_curve_rbf_cv$mean_accuracy + accuracy_curve_rbf_cv$sd_accuracy,
    angle = 90,
    code = 3,
    length = 0.04,
    col = "black"
  )
  
  lines(
    x = accuracy_curve_custom$p_basis,
    y = accuracy_curve_custom$mean_accuracy,
    type = "b",
    pch = 17,
    lwd = 3,
    col = "darkred"
  )
  
  arrows(
    x0 = accuracy_curve_custom$p_basis,
    y0 = accuracy_curve_custom$mean_accuracy - accuracy_curve_custom$sd_accuracy,
    x1 = accuracy_curve_custom$p_basis,
    y1 = accuracy_curve_custom$mean_accuracy + accuracy_curve_custom$sd_accuracy,
    angle = 90,
    code = 3,
    length = 0.04,
    col = "darkred"
  )
  
  legend(
    "bottomright",
    legend = c("RBF CV", "Custom gamma=1"),
    pch = c(16, 17),
    lwd = 3,
    col = c("black", "darkred"),
    bty = "n"
  )
}

plot_accuracy_boxplots <- function() {
  acc_all <- rbind(
    data.frame(model = "RBF CV", p_basis = results_all$p_basis,
               test_accuracy = results_all$test_accuracy),
    data.frame(model = "Custom gamma=1", p_basis = results_custom_kernel$p_basis,
               test_accuracy = results_custom_kernel$test_accuracy)
  )
  
  split_acc <- split(
    acc_all$test_accuracy,
    paste(acc_all$model, acc_all$p_basis, sep = " | p=")
  )
  
  boxplot(
    split_acc,
    las = 2,
    ylab = "Test accuracy",
    main = "Accuracy Distribution by Model and Dimension",
    outline = FALSE,
    col = rep(c("grey85", "mistyrose"), each = length(unique(acc_all$p_basis)))
  )
}

plot_gamma_curve <- function() {
  gamma_curve <- aggregate(
    gamma ~ p_basis,
    data = results_all,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  )
  gamma_curve <- data.frame(
    p_basis = gamma_curve$p_basis,
    mean_gamma = gamma_curve$gamma[, "mean"],
    sd_gamma = gamma_curve$gamma[, "sd"]
  )
  rownames(gamma_curve) <- NULL
  
  plot(
    x = gamma_curve$p_basis,
    y = gamma_curve$mean_gamma,
    type = "b",
    pch = 16,
    lwd = 3,
    log = "y",
    xlab = "Dimension",
    ylab = "Mean selected gamma",
    main = "CV-selected Gamma over Repetitions"
  )
}

# Render plots ------------------------------------------------------------
if (save_img) {
  png(file.path(img_dir, "KELM_accuracy_curves.png"), width = 1400, height = 900, res = 140)
  plot_accuracy_curves()
  dev.off()
  
  png(file.path(img_dir, "KELM_accuracy_boxplots.png"), width = 1800, height = 1000, res = 140)
  plot_accuracy_boxplots()
  dev.off()
  
  png(file.path(img_dir, "KELM_gamma_curve.png"), width = 1400, height = 900, res = 140)
  plot_gamma_curve()
  dev.off()
}

if (interactive()) {
  plot_accuracy_curves()
  plot_accuracy_boxplots()
  plot_gamma_curve()
}

cat("Plots generated from:\n")
cat(normalizePath(results_file, winslash = "/", mustWork = FALSE), "\n")
if (save_img) {
  cat("Images saved in:\n")
  cat(normalizePath(img_dir, winslash = "/", mustWork = FALSE), "\n")
}
