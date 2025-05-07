library(tidyverse)
library(data.table)
library(stringr)

# ELM Simulations for GAUSS

set.seed(825)
seed <- sample(0:10000, 1, replace = TRUE)
filename <- paste0('GAUSS100.R')
fileConn <- file(filename)
writeLines(c(paste0('source("GC.R")'),
             paste0('load("data.RData")'),
             paste0('numCores <- as.integer(detectCores()) - 1'),
             paste0('options(mc.cores = numCores)'),
             paste0('RNGkind("L', "'", 'Ecuyer-CMRG")'),
             paste0('set.seed(', seed, ')'),
             paste0('GAUSS100 <-  mclapply(data, function(zz) gauss(dat = zz$data, type = 0, hdim = 100, omega = 1, activation = 1, y_select = 1, z_select = 2, x_select = 3, max_lag = 2, m = 199, k = 10, r = 25), mc.preschedule = TRUE, mc.cores = getOption("mc.cores", numCores))'),
             paste0('save(GAUSS100, file = "GAUSS100.RData")'),
             paste0('print("COMPLETE")')),
           fileConn)
close(fileConn)

# ELM Simulations for ZERO

set.seed(829)
seed <- sample(0:10000, 1, replace = TRUE)
filename <- paste0('ZERO100.R')
fileConn <- file(filename)
writeLines(c(paste0('source("GC.R")'),
             paste0('load("data.RData")'),
             paste0('numCores <- as.integer(detectCores()) - 1'),
             paste0('options(mc.cores = numCores)'),
             paste0('RNGkind("L', "'", 'Ecuyer-CMRG")'),
             paste0('set.seed(', seed, ')'),
             paste0('ZERO100 <-  mclapply(data, function(zz) zero(dat = zz$data, type = 0, hdim = 100, omega = 1, activation = 1, y_select = 1, z_select = 2, x_select = 3, max_lag = 2, k = 10, r = 25), mc.preschedule = TRUE, mc.cores = getOption("mc.cores", numCores))'),
             paste0('save(ZERO100, file = "ZERO100.RData")'),
             paste0('print("COMPLETE")')),
           fileConn)
close(fileConn)