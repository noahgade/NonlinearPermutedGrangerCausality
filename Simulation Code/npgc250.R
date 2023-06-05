numCores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) - 1
options(mc.cores = numCores)

library(parallel)
library(tidyverse)
library(devtools)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("NPGC.cpp")
load("simdata250.RData")

data_prep <- function(data, lag = 3) {
  Y <- data[(lag + 1):nrow(data), 1, drop = FALSE]
  Ylag <- lagmatCpp(data[, 1, drop = FALSE], lag)
  Z <- lagmatCpp(data[, 2:4, drop = FALSE], lag)
  X <- lagmatCpp(data[, 5:7, drop = FALSE], lag)
  out <- list(Y = Y, Ylag = Ylag, Z = Z, X = X)
  return(out)
}

method_npgc <- function(data) {
  data1 <- data_prep(data, lag = 3)
  out <- npgc(data1$Y, data1$X, data1$Ylag, data1$Z, n_permutations = 400, n_folds = 5, n_initializations = 50, feature_dimension = 100)
  return(out)
}

# print(sample(0:10000, 1))
# [1] 3971

set.seed(3971)
npgc250 <- mclapply(simdata250, method_npgc, mc.cores = numCores)
save(npgc250, file = "/home/ndg5e/NPGC2/npgc250.RData")



