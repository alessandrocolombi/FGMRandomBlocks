# Replay one dumped graph-sampling input.
#
# Usage on the VM:
#   timeout 60s env FGM_GRAPH_ALGORITHM=rjmcmc Rscript Debug_Replay_Graph_Input.R debug_graph_inputs/<file>.rds

args = commandArgs(trailingOnly = TRUE)
if(length(args) < 1)
  stop("Usage: Rscript Debug_Replay_Graph_Input.R <graph_input.rds> [algorithm]")

dump_file = args[[1]]
if(!file.exists(dump_file))
  stop("Dump file does not exist: ", dump_file)

library("FGM")

source("./utility_functions.R")
source("./bulky_functions.R")
source("./get_things.R")
source("./bdgraph.R")

dump = readRDS(dump_file)

algorithm = if(length(args) >= 2) {
  args[[2]]
} else {
  Sys.getenv("FGM_GRAPH_ALGORITHM", unset = dump$algorithm)
}

if(!algorithm %in% c("rjmcmc", "bdmcmc", "rjmcmc.mpl"))
  stop("algorithm must be one of: rjmcmc, bdmcmc, rjmcmc.mpl")

print_default = if(!is.null(dump$graph_print)) dump$graph_print else 1L
print_period = as.integer(Sys.getenv("FGM_GRAPH_PRINT", unset = as.character(print_default)))
if(is.na(print_period) || print_period <= 0)
  print_period = 1L

cores_default = if(!is.null(dump$graph_cores)) dump$graph_cores else 1L
graph_cores = as.integer(Sys.getenv("FGM_GRAPH_CORES", unset = as.character(cores_default)))
if(is.na(graph_cores) || graph_cores <= 0)
  graph_cores = 1L

cat("Replay file: ", normalizePath(dump_file), "\n", sep = "")
cat("R version: ", R.version.string, "\n", sep = "")
cat("FGM version: ", as.character(utils::packageVersion("FGM")), "\n", sep = "")
cat("FGM path: ", find.package("FGM"), "\n", sep = "")
cat("algorithm: ", algorithm, "\n", sep = "")
cat("rj_steps: ", dump$rj_steps, "\n", sep = "")
cat("n: ", dump$n, "\n", sep = "")
cat("p: ", ncol(dump$data), "\n", sep = "")
cat("rho: ", paste(dump$rho, collapse = ","), "\n", sep = "")
cat("graph edges: ", sum(dump$graph) / 2, "\n", sep = "")
cat("beta_params: ", paste(dump$beta_params, collapse = ","), "\n", sep = "")
cat("print: ", print_period, "\n", sep = "")
cat("cores: ", graph_cores, "\n", sep = "")
cat("Starting post_graph_sampling at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
flush.console()

out = post_graph_sampling(
  dump$data,
  dump$rho,
  dump$n,
  method = "ggm",
  algorithm = algorithm,
  iter = dump$rj_steps,
  burnin = 0,
  g.start = dump$graph,
  save = TRUE,
  print = print_period,
  cores = graph_cores,
  threshold = dump$threshold,
  beta_params = dump$beta_params
)

cat("Finished post_graph_sampling at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
cat("last graph edges: ", sum(out$last_graph) / 2, "\n", sep = "")
cat("has last_K: ", !is.null(out$last_K), "\n", sep = "")

save_file = Sys.getenv("FGM_REPLAY_SAVE", unset = "")
if(nzchar(save_file)){
  saveRDS(out, save_file)
  cat("Saved replay output: ", save_file, "\n", sep = "")
}
