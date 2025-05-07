library(tidyverse)
library(data.table)
library(stringr)

# MLP Simulations for LASSO

set.seed(846)
seeds <- array(sample(0:10000, 5 * 1800, replace = TRUE), dim = c(5, 1800))
for(pen in 1:5) {
  for(iter in 1:1800) {
    filename <- paste0('LASSO', pen, '_', iter, '.R')
    fileConn <- file(filename)
    writeLines(c(paste0('source("GC.R")'),
                 paste0('load("data.RData")'),
                 paste0('RNGkind("L', "'", 'Ecuyer-CMRG")'),
                 paste0('set.seed(', seeds[pen, iter], ')'),
                 paste0('LASSO', pen, '_', iter, ' <-  lasso(dat = data[[', iter, ']]$data, type = 2, hdim = 100, omega = 1, activation = 1, y_select = 1, z_select = 2, x_select = 3, max_lag = 2, k = 10, r = 1, penalty = ', pen, ', lambdas = 2^seq(-1, 7, 0.1), chooselambda = TRUE)'),
                 paste0('save(LASSO', pen, '_', iter, ', file = "LASSO', pen, '_', iter, '.RData")'),
                 paste0('print("COMPLETE")')),
               fileConn)
    close(fileConn)
  }
}