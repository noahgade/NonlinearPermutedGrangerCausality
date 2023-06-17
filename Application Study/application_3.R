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
dat <- datraw[(datraw[, 1] >= 11 & datraw[, 1] < 14.7),]

datprep <- function(data, x_stimulus_col, z_stimulus_col) {
  nr0 <- nrow(data)
  out <- vector("list")
  out$Y <- data[201:nr0, 5, drop = FALSE]
  out$Ylag <- lagmatCpp(data[, 5, drop = FALSE], 40)[(201:nr0 - 40),, drop = FALSE]
  out$X <- lagmatCpp(data[, x_stimulus_col, drop = FALSE], 200)[, 121:160, drop = FALSE]
  out$Z <-  lagmatCpp(data[, z_stimulus_col, drop = FALSE], 200)[, c(121:160, 321:360), drop = FALSE]
  return(out)
}

# print(sample(0:1000, 1))
# [1] 845
set.seed(845)
test_data <- datprep(dat, x_stimulus_col = 4, z_stimulus_col = c(2, 3))

n_permutations <- 400
n_initializations <- 50
feature_dimension <- 250

W1 <- generateW(ncol(test_data$Z) + ncol(test_data$X) + ncol(test_data$Ylag) + 1, feature_dimension, n_initializations)
W2 <- generateW(ncol(test_data$Z) + ncol(test_data$X) + ncol(test_data$Ylag) + 1, feature_dimension, n_initializations)
W3 <- generateW(ncol(test_data$Z) + ncol(test_data$X) + ncol(test_data$Ylag) + 1, feature_dimension, n_initializations)
W4 <- generateW(ncol(test_data$Z) + ncol(test_data$X) + ncol(test_data$Ylag) + 1, feature_dimension, n_initializations)
W5 <- generateW(ncol(test_data$Z) + ncol(test_data$X) + ncol(test_data$Ylag) + 1, feature_dimension, n_initializations)

app_npgc <- function(data) {
  out <- application_npgc(data$Y, data$X, data$Ylag, data$Z, n_permutations, n_initializations, feature_dimension, W1, W2, W3, W4, W5)
  return(out)
}

result3 <- app_npgc(test_data)
save(result3, file = "/home/ndg5e/NPGCapp/result3.RData")