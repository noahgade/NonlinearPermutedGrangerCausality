library(tidyverse)
library(data.table)
library(stringr)

# MLP Simulations for LASSO
# For each iteration, change 'seed.row' item from 1 to 25
# Builds training instances by using the same files with different set seeds

seed.row <- 1

penalties <- 1:5
set.seed(846)
seeds <- array(sample(0:10000, 1800 * 5 * 25, replace = TRUE), dim = c(5, 25, 1800))
for(pen in 1:5) {
  for(iter in 1:1800) {
    filename <- paste0('LASSO', pen, '_', iter, '.R')
    fileConn <- file(filename)
    writeLines(c(paste0('source("GC.R")'),
                 paste0('load("data.RData")'),
                 paste0('RNGkind("L', "'", 'Ecuyer-CMRG")'),
                 paste0('set.seed(', seeds[pen, seed.row, iter], ')'),
                 paste0('LASSO', pen, '_', iter, ' <-  lasso(dat = data[[', iter, ']]$data, type = 2, hdim = 100, omega = 1, activation = 1, y_select = 1, z_select = 2, x_select = 3, max_lag = 2, k = 10, r = 1, penalty = ', pen, ', lambdas = 2^seq(-1, 7, 0.1), chooselambda = TRUE)'),
                 paste0('save(LASSO', pen, '_', seed.row, '_', iter, ', file = "LASSO', pen, '_', seed.row, '_', iter, '.RData")'),
                 paste0('print("COMPLETE")')),
               fileConn)
    close(fileConn)
  }
}