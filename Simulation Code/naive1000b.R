numCores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) - 1
options(mc.cores = numCores)

library(parallel)
library(tidyverse)
library(devtools)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("NPGC.cpp")
load("simdata1000b.RData")

data_prep <- function(data, lag = 3) {
  Y <- data[(lag + 1):nrow(data), 1, drop = FALSE]
  Ylag <- lagmatCpp(data[, 1, drop = FALSE], lag)
  Z <- lagmatCpp(data[, 2:4, drop = FALSE], lag)
  X <- lagmatCpp(data[, 5:7, drop = FALSE], lag)
  out <- list(Y = Y, Ylag = Ylag, Z = Z, X = X)
  return(out)
}

method_naive <- function(data) {
  data1 <- data_prep(data, lag = 3)
  out <- naive(data1$Y, data1$X, data1$Ylag, data1$Z, n_initializations = 1000, feature_dimension = 100)
  return(out)
}

# print(sample(0:10000, 1))
# [1] 1032

set.seed(1032)
naive1000b <- mclapply(simdata1000b, method_naive, mc.cores = numCores)
save(naive1000b, file = "/home/ndg5e/NPGC3/naive1000b.RData")


