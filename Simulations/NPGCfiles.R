library(tidyverse)
library(data.table)
library(stringr)

# ELM Simulations for NPGC
# changing N, lags, activation, and W0

Ns <- c(25, 50, 100)
lags <- c(2, 4, 8)
gs <- c("tanh" = 0, "relu" = 1, "lrelu" = 2)
NPGCkey <- tibble(names = paste0("NPGC", rep(Ns, each = length(lags) * length(gs)), "_",
                                 rep(rep(lags, each = length(gs)), length(Ns)), "_",
                                 rep(names(gs), length(lags) * length(Ns))),
                  N = rep(Ns, each = length(lags) * length(gs)),
                  lag = rep(rep(lags, each = length(gs)), length(Ns)),
                  g = rep(gs, length(lags) * length(Ns)))
save(NPGCkey, file = "NPGCkey.RData")

set.seed(219)
seeds <- sample(0:10000, nrow(NPGCkey), replace = TRUE)
for(iter in 1:nrow(NPGCkey)) {
  filename <- paste0('NPGC_', iter, '.R')
  fileConn <- file(filename)
  writeLines(c(paste0('source("GC.R")'),
               paste0('load("data.RData")'),
               paste0('numCores <- as.integer(detectCores()) - 1'),
               paste0('options(mc.cores = numCores)'),
               paste0('RNGkind("L', "'", 'Ecuyer-CMRG")'),
               paste0('set.seed(', seeds[iter], ')'),
               paste0('NPGC_', iter, ' <-  mclapply(data, function(zz) npgc(dat = zz$data, type = 0, hdim = ', NPGCkey$N[iter], ', omega = 1, activation = ', NPGCkey$g[iter], ', y_select = 1, z_select = 2, x_select = 3, max_lag = ', NPGCkey$lag[iter], ', m = 199, k = 10, r = 25), mc.preschedule = TRUE, mc.cores = getOption("mc.cores", numCores))'),
               paste0('save(NPGC_', iter, ', file = "NPGC_', iter, '.RData")'),
               paste0('print("COMPLETE")')),
             fileConn)
  close(fileConn)
}
