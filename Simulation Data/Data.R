library(abind)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("~/Wake Forest University Dropbox/Noah Gade/Research/NPGC/Simulations/Data/Data.cpp")

tarData <- function(L = 1000, Dim = 2, Lag = 2, RegimeL = 10, Causal = NULL) {
  # Generates threshold autoregressive data, N(0, 1) innovations
  # L: length of time series
  # Dim: dimension of time series
  # Lag: lag dependence of time series
  # RegimeL: regime length for threshold, sum of previous RegimeL observations of the first variable
  # Causal: causal structure for VAR (dimension Dim x (Dim * Lag))
  Changes <- NULL
  States <- NULL
  if(is.null(Causal) & is.null(Changes)) {
    Causal <- array(1, dim = c(Dim, Dim * Lag, 2))
  }
  success <- FALSE
  while(!success) {
    Result <- TAR(L, Dim, Lag, RegimeL, Causal)
    maximum <- max(Result$data)
    if(maximum <= 25) {
      success <- TRUE
    }
  }
  return(Result)
}

varData <- function(L = 1000, Dim = 2, Lag = 2, Causal = NULL) {
  # Generates vector autoregressive data, N(0, 1) innovations
  # L: length of time series
  # Dim: dimension of time series
  # Lag: lag dependence of time series
  # RegimeL: regime length for threshold, sum of previous RegimeL observations of the first variable
  # Causal: causal structure for VAR (dimension Dim x (Dim * Lag))
  Changes <- NULL
  States <- NULL
  if(is.null(Causal) & is.null(Changes)) {
    Causal <- matrix(1, nrow = Dim, ncol = Dim * Lag)
  }
  success <- FALSE
  while(!success) {
    Result <- VAR(L, Dim, Lag, Causal)
    maximum <- max(Result$data)
    if(maximum <= 25) {
      success <- TRUE
    }
  }
  return(Result)
}

