numCores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) - 1
options(mc.cores = numCores)

library(parallel)
library(tidyverse)
library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(R.matlab)
sourceCpp("NPGC.cpp")
datraw <- readMat("dat.mat")[[1]]
dat <- datraw[(datraw[, 1] >= 11 & datraw[, 1] < 15),]

datprep <- function(data, x_stimulus_col, z_stimulus_col) {
  out <- vector("list", 2)
  nr0 <- nrow(data)
  keep <- c(sapply(seq(25, 50, 5) * 4000 / 1000, function(x) -2:2 + x))
  keep_lag <- c(sapply(seq(2, 20, 2) * 4000 / 1000, function(x) -2:2 + x))
  for(col in 1:length(out)) {
    out[[col]]$Y <- data[221:nr0, 4 + col, drop = FALSE]
    out[[col]]$Ylag <- lagmatCpp(data[, 4 + col, drop = FALSE], 85)[(221:nr0 - 85), keep_lag, drop = FALSE]
    out[[col]]$X <- lagmatCpp(data[, x_stimulus_col, drop = FALSE], 220)[, keep, drop = FALSE]
    out[[col]]$Z <- lagmatCpp(data[, z_stimulus_col, drop = FALSE], 220)[, keep, drop = FALSE]
  }
  return(out)
}

# print(sample(0:1000, 1))
# [1] 364
set.seed(364)
test_data <- datprep(dat, x_stimulus_col = 3, z_stimulus_col = 2)

n_permutations <- 400
n_initializations <- 50
feature_dimension <- 200

W1 <- generateW(ncol(test_data[[1]]$Z) + ncol(test_data[[1]]$X) + ncol(test_data[[1]]$Ylag) + 1, feature_dimension, n_initializations)
W2 <- generateW(ncol(test_data[[1]]$Z) + ncol(test_data[[1]]$X) + ncol(test_data[[1]]$Ylag) + 1, feature_dimension, n_initializations)
W3 <- generateW(ncol(test_data[[1]]$Z) + ncol(test_data[[1]]$X) + ncol(test_data[[1]]$Ylag) + 1, feature_dimension, n_initializations)
W4 <- generateW(ncol(test_data[[1]]$Z) + ncol(test_data[[1]]$X) + ncol(test_data[[1]]$Ylag) + 1, feature_dimension, n_initializations)
W5 <- generateW(ncol(test_data[[1]]$Z) + ncol(test_data[[1]]$X) + ncol(test_data[[1]]$Ylag) + 1, feature_dimension, n_initializations)

app_npgc <- function(data) {
  out <- application_npgc(data$Y, data$X, data$Ylag, data$Z, n_permutations, n_initializations, feature_dimension, W1, W2, W3, W4, W5)
  return(out)
}

appresult <- vector("list", 2)
for(item in 1:2) {
  appresult[[item]] <- app_npgc(test_data[[item]])
}

save(appresult, file = "/home/ndg5e/NPGCapp/result32.RData")
