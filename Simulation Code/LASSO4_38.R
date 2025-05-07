source("GC.R")
load("data.RData")
numCores <- as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE")) - 1
options(mc.cores = numCores)
RNGkind("L'Ecuyer-CMRG")
set.seed(8448)
LASSO4_38 <-  mclapply(data[1333:1368], function(zz) lasso(dat = zz$data, type = 2, hdim = 100, omega = 1, activation = 1, y_select = 1, z_select = 2, x_select = 3, max_lag = 2, k = 10, r = 25, penalty = 4, lambdas = 2^seq(-1, 7, 0.1), chooselambda = TRUE), mc.preschedule = TRUE, mc.cores = getOption("mc.cores", numCores))
save(LASSO4_38, file = "LASSO4_38.RData")
print("COMPLETE")
