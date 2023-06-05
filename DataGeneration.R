rm(list = ls())
setwd("~/Desktop/Research/Rodu/NPGC")
library(tidyverse)
library(purrr)
library(deSolve)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("NPGC.cpp")

threshold <- function(x, cut) {
  return((abs(x) >= cut) * x)
}

l96 <- function(Time, State, Param) {
  p <- length(State)
  state_copy <- State
  states <- rep(State, 3)
  for(i in 1:p) {
    state_copy[i] <- (states[p + i + 1] - states[p + i - 2]) * states[p + i - 1] - states[p + i] + Param
  }
  return(list(state_copy))
}

simLorenz96 <- function(Time, dim) {
  success <- 0
  while(success != 1) {
    FC <- sample(5:20, 1)
    states <- suppressWarnings(unname(deSolve::ode(func = l96, y = rnorm(dim), times = seq(0, (500 + Time) * 0.05, 0.05), parms = FC)))
    if(nrow(states) == (Time + 501)) {
      success <- 1
    }
  }
  
  out <- states[502:(Time + 501), 2:(dim + 1)]
  return(out)
}

simTAR2 <- function(Time, dim) {
  max_eigen1 <- Inf
  rowA1 <- 0
  while(max_eigen1 > 0.8 || rowA1 == 0) {
    A1.1 <- matrix(threshold(runif(dim * dim, -0.5, 0.5), 0.1), nrow = dim)
    A2.1 <- matrix(threshold(runif(dim * dim, -0.5, 0.5), 0.1), nrow = dim)
    A1 <- rbind(cbind(A1.1, A2.1), cbind(diag(dim), matrix(0, dim, dim)))
    rowA1 <- rowSums(A1)
    max_eigen1 <- max(abs(eigen(A1)$values))
  }
  max_eigen2 <- Inf
  rowA2 <- 0
  while(max_eigen2 > 0.8 || rowA2 == 0) {
    A1.2 <- matrix(threshold(runif(dim * dim, -0.5, 0.5), 0.1), nrow = dim)
    A2.2 <- matrix(threshold(runif(dim * dim, -0.5, 0.5), 0.1), nrow = dim)
    A2 <- rbind(cbind(A1.2, A2.2), cbind(diag(dim), matrix(0, dim, dim)))
    rowA2 <- rowSums(A2)
    max_eigen2 <- max(abs(eigen(A2)$values))
  }
  regime <- rep(0, Time)
  states <- matrix(0, nrow = 2 + 500 + Time, ncol = dim)
  for(t in 3:nrow(states)) {
    if(sum(states[t - 2,]) <= 0) {
      newstates <- A1 %*% matrix(c(states[t - 1,], states[t - 2,])) + matrix(c(stats::rnorm(dim, 0, 0.5), rep(0, dim)))
      states[t,] <- newstates[1:dim,]
      regime[t - 2] <- 1
    } else {
      newstates <- A2 %*% matrix(c(states[t - 1,], states[t - 2,])) + matrix(c(stats::rnorm(dim, 0, 0.2), rep(0, dim)))
      states[t,] <- newstates[1:dim,]
      regime[t - 2] <- 2
    }
  }
  out <- states[(2 + 500 + 1):(2 + 500 + Time),]
  return(out)
}

GCsample <- function(GCgroup = 1, dims = c(6, 6), GC = TRUE, NZ = 0) {
  gp1 <- 1:dims[1]
  gp2 <- (dims[1] + 1):sum(dims)
  
  if(GCgroup == 1) {
    Yindx <- sample(gp1, 1)
    
    if(GC) {
      Xindx <- c(sample(gp1[gp1 != Yindx], 2), sample(gp2, 1))
      Zindx <- c(sample(gp1[!(gp1 %in% c(Yindx, Xindx))], NZ), sample(gp2[!(gp2 %in% Xindx)], 3 - NZ))
    } else {
      Xindx <- c(sample(gp1[gp1 != Yindx], 0), sample(gp2, 3))
      Zindx <- c(sample(gp1[!(gp1 %in% c(Yindx, Xindx))], NZ), sample(gp2[!(gp2 %in% Xindx)], 3 - NZ))
    }
    
  } else if(GCgroup == 2) {
    Yindx <- sample(gp2, 1)
    
    if(GC) {
      Xindx <- c(sample(gp2[gp2 != Yindx], 2), sample(gp1, 1))
      Zindx <- c(sample(gp2[!(gp2 %in% c(Yindx, Xindx))], NZ), sample(gp1[!(gp1 %in% Xindx)], 3 - NZ))
    } else {
      Xindx <- c(sample(gp2[gp1 != Yindx], 0), sample(gp1, 3))
      Zindx <- c(sample(gp2[!(gp2 %in% c(Yindx, Xindx))], NZ), sample(gp1[!(gp1 %in% Xindx)], 3 - NZ))
    }
  }
  return(c(Yindx, Zindx, Xindx))
}

gen_sample <- function(nsims) {
  out <- c(replicate(nsims, GCsample(GCgroup = 1, dim = c(6, 6), GC = FALSE, NZ = 0), simplify = FALSE), 
           replicate(nsims, GCsample(GCgroup = 2, dim = c(6, 6), GC = FALSE, NZ = 0), simplify = FALSE),
           replicate(nsims, GCsample(GCgroup = 1, dim = c(6, 6), GC = FALSE, NZ = 2), simplify = FALSE), 
           replicate(nsims, GCsample(GCgroup = 2, dim = c(6, 6), GC = FALSE, NZ = 2), simplify = FALSE),
           replicate(nsims, GCsample(GCgroup = 1, dim = c(6, 6), GC = TRUE, NZ = 0), simplify = FALSE), 
           replicate(nsims, GCsample(GCgroup = 2, dim = c(6, 6), GC = TRUE, NZ = 0), simplify = FALSE),
           replicate(nsims, GCsample(GCgroup = 1, dim = c(6, 6), GC = TRUE, NZ = 2), simplify = FALSE), 
           replicate(nsims, GCsample(GCgroup = 2, dim = c(6, 6), GC = TRUE, NZ = 2), simplify = FALSE))
  return(out)
}

data_sim <- function(Time, dims = c(6, 6), columns) {
  nsims <- length(columns)
  dat <- replicate(nsims, cbind(simLorenz96(Time, dims[1]), simTAR2(Time, dims[2])), simplify = FALSE)
  out_dat <- mapply(function(mat, cols) mat[,cols], dat, columns, SIMPLIFY = FALSE)
  return(out_dat)
}

set.seed(747)
dims = c(6, 6)
max_lag <- 3
nsims <- 200

# Column Selection for Granger Causality
col_select250 <- gen_sample(nsims)
save(col_select250, file = "cs250.RData")
col_select500 <- gen_sample(nsims)
save(col_select500, file = "cs500.RData")
col_select1000 <- gen_sample(nsims)
save(col_select1000, file = "cs1000.RData")

# Data Generation
simdata250 <- data_sim(250 + max_lag, dims, col_select250)
save(simdata250, file = "simdata250.RData")

simdata500 <- data_sim(500 + max_lag, dims, col_select500)
save(simdata500, file = "simdata500.RData")

simdata1000 <- data_sim(1000 + max_lag, dims, col_select1000)
save(simdata1000, file = "simdata1000.RData")

