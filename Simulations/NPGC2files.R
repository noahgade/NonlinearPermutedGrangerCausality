library(tidyverse)
library(data.table)
library(stringr)

# MLP Simulations for NPGC

set.seed(1002)
seeds <- sample(0:10000, 1800, replace = TRUE)
for(iter in 1:1800) {
  filename <- paste0('NPGC2_', iter, '.R')
  fileConn <- file(filename)
  writeLines(c(paste0('source("GC.R")'),
               paste0('load("data.RData")'),
               paste0('RNGkind("L', "'", 'Ecuyer-CMRG")'),
               paste0('set.seed(', seeds[iter], ')'),
               paste0('NPGC2_', iter, ' <-  npgc(dat = data[[', iter, ']]$data, type = 2, hdim = 100, omega = 1, activation = 0, y_select = 1, z_select = 2, x_select = 3, max_lag = 2, m = 199, k = 10, r = 1)'),
               paste0('save(NPGC2_', iter, ', file = "NPGC2_', iter, '.RData")'),
               paste0('print("COMPLETE")')),
             fileConn)
  close(fileConn)
}
