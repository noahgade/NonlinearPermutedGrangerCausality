library(tidyverse)
library(data.table)
library(stringr)

# MLP Simulations for NPGC
# For each iteration, change 'seed.row' item from 1 to 25
# Builds training instances by using the same files with different set seeds

seed.row <- 1

set.seed(651)
seeds <- matrix(sample(0:10000, 1800 * 25, replace = TRUE), nrow = 25)
for(iter in 1:1800) {
  filename <- paste0('NPGC2_', iter, '.R')
  fileConn <- file(filename)
  writeLines(c(paste0('source("GC.R")'),
               paste0('load("data.RData")'),
               paste0('RNGkind("L', "'", 'Ecuyer-CMRG")'),
               paste0('set.seed(', seeds[seed.row, iter], ')'),
               paste0('NPGC2_', seed.row, '_', iter, ' <- npgc(dat = data[[', iter, ']]$data, type = 2, hdim = 100, omega = 1, activation = 1, y_select = 1, z_select = 2, x_select = 3, max_lag = 2, m = 199, k = 10, r = 1)'),
               paste0('save(NPGC2_', seed.row, '_', iter, ', file = "NPGC2_', seed.row, '_', iter, '.RData")'),
               paste0('print("COMPLETE")')),
             fileConn)
  close(fileConn)
}