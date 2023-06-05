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
dat <- datraw[(datraw[, 1] >= 11.60 & datraw[, 1] < 12.30),]

datprep <- function(data, stim_lag_delay) {
  out <- vector("list", ncol(data) - 2)
  nr0 <- nrow(data)
  nr1 <- nr0 - (stim_lag_delay + 20) * 4
  for(col in 1:length(out)) {
    out[[col]]$Y <- data[(nr0 - nr1 + 1):nr0, col + 2, drop = FALSE]
    out[[col]]$Ylag <- lagmatCpp(data[, col + 2, drop = FALSE], 80)[(nr0 - 80 - nr1 + 1):(nr0 - 80),, drop = FALSE]
    out[[col]]$X <- lagmatCpp(data[, 2, drop = FALSE], 80)[1:nr1,, drop = FALSE]
  }
  return(out)
}

# print(sample(0:10000, 1))
# [1] 2307
set.seed(2307)
dat1 <- datprep(dat, stim_lag_delay = 30)

n_permutations <- 400
n_initializations <- 50
feature_dimension <- 100

W1 <- generateW(ncol(dat1[[1]]$X) + ncol(dat1[[1]]$Ylag) + 1, feature_dimension, n_initializations)
W2 <- generateW(ncol(dat1[[1]]$X) + ncol(dat1[[1]]$Ylag) + 1, feature_dimension, n_initializations)
W3 <- generateW(ncol(dat1[[1]]$X) + ncol(dat1[[1]]$Ylag) + 1, feature_dimension, n_initializations)
W4 <- generateW(ncol(dat1[[1]]$X) + ncol(dat1[[1]]$Ylag) + 1, feature_dimension, n_initializations)
W5 <- generateW(ncol(dat1[[1]]$X) + ncol(dat1[[1]]$Ylag) + 1, feature_dimension, n_initializations)

app_npgc <- function(data) {
  out <- application_npgc(data$Y, data$X, data$Ylag, n_permutations, n_initializations, feature_dimension, W1, W2, W3, W4, W5)
  return(out)
}

app4 <- mclapply(dat1, app_npgc, mc.cores = numCores)
save(app4, file = "/home/ndg5e/NPGCapp/app4.RData")
