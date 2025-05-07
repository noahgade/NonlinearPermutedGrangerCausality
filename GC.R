library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(Matrix)
library(MASS)
library(parallel)
sourceCpp("GC.cpp")

# Simulation Code for NPGC
npgc <- function(dat, type, hdim, omega, activation, y_select, z_select, x_select, max_lag, m = 199, k = 10, r = 25) {
  lag_set <- max_lag:1
  L <- nrow(dat)
  Y <- dat[(max_lag + 1):L, y_select, drop = FALSE]
  Ylag <- do.call(cbind, lapply(lag_set, function(arg) dat[arg:(L - max_lag + arg - 1), y_select, drop = FALSE]))
  Z <- do.call(cbind, lapply(lag_set, function(arg) dat[arg:(L - max_lag + arg - 1), z_select, drop = FALSE]))
  X <- do.call(cbind, lapply(lag_set, function(arg) dat[arg:(L - max_lag + arg - 1), x_select, drop = FALSE]))
  output <- NPGC(Y = Y, Ylag = Ylag, Z = Z, X = X, Type = type, Omega = omega, M = m, K = k, R = r, HDim = hdim, Activation = activation)
  return(output)
}

# Simulation Code for GAUSS
gauss <- function(dat, type, hdim, omega, activation, y_select, z_select, x_select, max_lag, m = 199, k = 10, r = 25) {
  lag_set <- max_lag:1
  L <- nrow(dat)
  Y <- dat[(max_lag + 1):L, y_select, drop = FALSE]
  Ylag <- do.call(cbind, lapply(lag_set, function(arg) dat[arg:(L - max_lag + arg - 1), y_select, drop = FALSE]))
  Z <- do.call(cbind, lapply(lag_set, function(arg) dat[arg:(L - max_lag + arg - 1), z_select, drop = FALSE]))
  X <- do.call(cbind, lapply(lag_set, function(arg) dat[arg:(L - max_lag + arg - 1), x_select, drop = FALSE]))
  output <- GAUSS(Y = Y, Ylag = Ylag, Z = Z, X = X, Type = type, Omega = omega, M = m, K = k, R = r, HDim = hdim, Activation = activation)
  return(output)
}

# Simulation Code for ZERO
zero <- function(dat, type, hdim, omega, activation, y_select, z_select, x_select, max_lag, k = 10, r = 25) {
  lag_set <- max_lag:1
  L <- nrow(dat)
  Y <- dat[(max_lag + 1):L, y_select, drop = FALSE]
  Ylag <- do.call(cbind, lapply(lag_set, function(arg) dat[arg:(L - max_lag + arg - 1), y_select, drop = FALSE]))
  Z <- do.call(cbind, lapply(lag_set, function(arg) dat[arg:(L - max_lag + arg - 1), z_select, drop = FALSE]))
  X <- do.call(cbind, lapply(lag_set, function(arg) dat[arg:(L - max_lag + arg - 1), x_select, drop = FALSE]))
  output <- ZERO(Y = Y, Ylag = Ylag, Z = Z, X = X, Type = type, Omega = omega, K = k, R = r, HDim = hdim, Activation = activation)
  return(output)
}

# Simulation Code for LASSO
lasso <- function(dat, type, hdim, omega, activation, y_select, z_select, x_select, max_lag, k = 10, r = 25, penalty, lambdas, chooselambda) {
  lag_set <- max_lag:1
  L <- nrow(dat)
  Y <- dat[(max_lag + 1):L, y_select, drop = FALSE]
  Ylag <- do.call(cbind, lapply(lag_set, function(arg) dat[arg:(L - max_lag + arg - 1), y_select, drop = FALSE]))
  Z <- do.call(cbind, lapply(lag_set, function(arg) dat[arg:(L - max_lag + arg - 1), z_select, drop = FALSE]))
  X <- do.call(cbind, lapply(lag_set, function(arg) dat[arg:(L - max_lag + arg - 1), x_select, drop = FALSE]))
  output <- LASSO(Y = Y, Ylag = Ylag, Z = Z, X = X, Type = type, Penalty = penalty, Lambdas = lambdas, Omega = omega, K = k, R = r, HDim = hdim, ChooseLambda = chooselambda, Activation = activation)
  return(output)
}
