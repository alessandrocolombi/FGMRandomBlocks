# wd ----------------------------------------------------------------------
wd_pc_ale = "C:/Users/colom/FGMRandomBlocks/"
wd_pc_luciano = "C://Users//lucia//Desktop//PhD//my collaborations//Alessandro Colombi work//FGMRandomBlocks"
wd_bocconi = "/home/colombi/FGMRandomBlocks/"
wd_vec = c(wd_pc_ale, wd_pc_luciano, wd_bocconi)
choose_wd = wd_vec[1] # <--- modify here to select the wd according to the user
wd = paste0(choose_wd, "./")
setwd(wd)

# Read --------------------------------------------------------------------

result_file <- file.path(
  wd,
  "chains",
  "results_graph_part_dbig_parallel_Nrep50_niter80000_n500.csv"
)

res <- read.csv(result_file)
dim(res)

# Plot options ------------------------------------------------------------

save_plot <- TRUE
show_title <- FALSE
img_dir <- file.path(wd, "img")
dir.create(img_dir, recursive = TRUE, showWarnings = FALSE)

pdf_width <- 6
pdf_height <- 8

n_label <- sub("^.*_n([0-9]+)\\.[^.]+$", "n\\1", basename(result_file))
if(identical(n_label, basename(result_file)))
  n_label <- "n_unknown"

cex_axis <- 2
cex_lab <- 2
cex_main <- 2
cex_names <- 2

plot_mgp <- c(5, 1.4, 0)
plot_mar <- c(9, 9, if(show_title) 4 else 1, 2) + 0.1

# Dark colors used for model competitors.
model_colors <- c(
  "CSDA" = "gray25",
  "eta = 0" = "darkgreen",
  "eta = 0.5" = "darkblue",
  "eta = 0.75" = "darkred",
  "eta = 0.9" = "darkorange"
)

eta_levels <- c("0", "0.5", "0.75", "0.9")

rho_scenarios <- c(
  "rho0_1",
  "rho0_shift",
  "rho0_SM"
)

# Helpers -----------------------------------------------------------------

normalize_eta <- function(x){
  x_chr <- as.character(x)
  x_num <- suppressWarnings(as.numeric(x_chr))
  idx_num <- !is.na(x_num)
  x_chr[idx_num] <- sub(
    "\\.$",
    "",
    sub(
      "0+$",
      "",
      sprintf("%.2f", x_num[idx_num])
    )
  )
  x_chr
}

res$eta_chr <- normalize_eta(res$eta)

metric_specs <- list(
  ARI_true_vi = list(
    column = "ARI_true_vi",
    label = "ARI true VI",
    include_csda = FALSE
  ),
  ARI_true_binder = list(
    column = "ARI_true_bind",
    label = "ARI true Binder",
    include_csda = FALSE
  ),
  ARI_true_p05 = list(
    column = "ARI_true_p05",
    label = "ARI true p > 0.5",
    include_csda = FALSE
  ),
  threshold = list(
    column = "threshold",
    label = "BFDR threshold",
    include_csda = TRUE
  ),
  STD_SHD = list(
    column = "STD_SHD",
    label = "Standardized SHD",
    include_csda = TRUE
  ),
  F1 = list(
    column = "F1",
    label = "F1 score",
    include_csda = TRUE
  ),
  Wasserstein = list(
    column = "Wasserstein",
    label = "Wasserstein distance",
    include_csda = TRUE
  )
)

get_boxplot_values <- function(data, rho_type, metric_col, include_csda){
  values <- list()
  rho_eta_levels <- eta_levels[
    eta_levels %in% data$eta_chr[data$rho_type == rho_type]
  ]

  if(include_csda){
    values[["CSDA"]] <- data[
      data$rho_type == "None" & data$eta_chr == "CSDA",
      metric_col
    ]
  }

  for(eta_value in rho_eta_levels){
    values[[paste("eta =", eta_value)]] <- data[
      data$rho_type == rho_type & data$eta_chr == eta_value,
      metric_col
    ]
  }

  lapply(values, function(x) x[!is.na(x)])
}

make_plot_file <- function(rho_type, metric_name){
  file.path(
    img_dir,
    paste0("boxplot_", rho_type, "_", metric_name, "_", n_label, ".pdf")
  )
}

get_axis_labels <- function(group_names){
  labels <- sub("^eta = ", "eta == ", group_names)
  labels[labels == "CSDA"] <- "'CSDA'"
  parse(text = labels)
}

plot_metric_boxplot <- function(data, rho_type, metric_name, spec){
  metric_col <- spec$column

  if(!metric_col %in% names(data))
    stop("Column not found in results: ", metric_col)

  values <- get_boxplot_values(
    data = data,
    rho_type = rho_type,
    metric_col = metric_col,
    include_csda = spec$include_csda
  )

  empty_groups <- names(values)[vapply(values, length, integer(1)) == 0]
  if(length(empty_groups) > 0){
    warning(
      "Empty groups for ",
      metric_name,
      " / ",
      rho_type,
      ": ",
      paste(empty_groups, collapse = ", ")
    )
  }

  values <- values[vapply(values, length, integer(1)) > 0]
  plot_colors <- model_colors[names(values)]

  if(save_plot){
    pdf(
      file = make_plot_file(rho_type, metric_name),
      width = pdf_width,
      height = pdf_height
    )
    on.exit(dev.off(), add = TRUE)
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(
    mar = plot_mar,
    mgp = plot_mgp,
    las = 2,
    cex.axis = cex_axis,
    cex.lab = cex_lab,
    cex.main = cex_main
  )

  boxplot(
    values,
    col = plot_colors,
    border = "gray15",
    ylab = spec$label,
    main = if(show_title) paste(spec$label, "|", rho_type) else "",
    xaxt = "n",
    cex.axis = cex_axis,
    cex.lab = cex_lab,
    cex.main = cex_main,
    names = names(values),
    outline = TRUE
  )

  axis(
    side = 1,
    at = seq_along(values),
    labels = FALSE
  )
  text(
    x = seq_along(values),
    y = par("usr")[3] - 0.055 * diff(par("usr")[3:4]),
    labels = get_axis_labels(names(values)),
    srt = 45,
    adj = 1,
    xpd = TRUE,
    cex = cex_names
  )

  invisible(values)
}

# Boxplots ----------------------------------------------------------------

rho_types_to_plot <- c(
  "rho0_1",
  "rho0_shift",
  "rho0_SM"
)

metrics_to_plot <- c(
  "ARI_true_vi",
  "ARI_true_binder",
  "ARI_true_p05",
  "threshold",
  "STD_SHD",
  "F1",
  "Wasserstein"
)

for(rho_type_to_plot in rho_types_to_plot){
  for(metric_name in metrics_to_plot){
    plot_metric_boxplot(
      data = res,
      rho_type = rho_type_to_plot,
      metric_name = metric_name,
      spec = metric_specs[[metric_name]]
    )
  }
}
