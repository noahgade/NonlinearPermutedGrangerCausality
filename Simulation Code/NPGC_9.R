source("GC.R")
load("data.RData")
numCores <- as.integer(Sys.getenv("SLURM_NTASKS_PER_NODE")) - 1
options(mc.cores = numCores)
RNGkind("L'Ecuyer-CMRG")
set.seed(6042)
NPGC_9 <-  mclapply(data, function(zz) npgc(dat = zz$data, type = 0, hdim = 25, omega = 1, activation = 2, y_select = 1, z_select = 2, x_select = 3, max_lag = 8, m = 199, k = 10, r = 25), mc.preschedule = TRUE, mc.cores = getOption("mc.cores", numCores))
save(NPGC_9, file = "NPGC_9.RData")
print("COMPLETE")
