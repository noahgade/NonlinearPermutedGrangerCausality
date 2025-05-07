source("GC.R")
load("data.RData")
numCores <- as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE")) - 1
options(mc.cores = numCores)
RNGkind("L'Ecuyer-CMRG")
set.seed(8210)
LASSO5_16 <-  mclapply(data[541:576], function(zz) lasso(dat = zz$data, type = 2, hdim = 100, omega = 1, activation = 1, y_select = 1, z_select = 2, x_select = 3, max_lag = 2, k = 10, r = 25, penalty = 5, lambdas = 2^seq(-1, 7, 0.1), chooselambda = TRUE), mc.preschedule = TRUE, mc.cores = getOption("mc.cores", numCores))
save(LASSO5_16, file = "LASSO5_16.RData")
print("COMPLETE")
