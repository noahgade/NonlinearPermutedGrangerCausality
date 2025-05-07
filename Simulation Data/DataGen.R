rm(list = ls())
source("~/Wake Forest University Dropbox/Noah Gade/Research/NPGC/Simulations/Data/Data.R")

# Setting parameters for data generation
dim <- 3
lag <- 2
regimelag <- 2

# Key for data generation: dataMMM_LLL_GCN
# MMM: Model, VAR or TAR
# LLL: Length, 250, 500, 1000
# GCN: Granger Causal designation, GC1 = causal at lags 1 & 2, GC2 = causal at lag 1, NC1 = fully disconnected, NC2 = fully disconnected from dim 1, NC3 = not causal for Dims 1 or 2, NC4 = not causal for Dim 1

n <- 50
models <- c("VAR", "TAR")
lens <- c(250, 500, 1000)
gcns <- c("GC1", "GC2", "NC1", "NC2", "NC3", "NC4")
r <- length(models) * length(lens) * length(gcns)
causal <- array(c(matrix(1, nrow = dim, ncol = dim * lag), 
                  matrix(c(rep(1, 6), 1, rep(1, 8), 0, rep(1, 2)), nrow = dim, ncol = dim * lag),
                  matrix(c(rep(1, 2), 0, rep(1, 2), rep(0, 3), rep(1, 3), 0, rep(1, 2), rep(0, 3), 1), nrow = dim, ncol = dim * lag),
                  matrix(c(rep(1, 2), 0, rep(1, 3), 0, rep(1, 4), 0, rep(1, 3), 0, rep(1, 2)), nrow = dim, ncol = dim * lag),
                  matrix(c(rep(1, 6), rep(0, 2), 1, rep(1, 6), rep(0, 2), 1), nrow = dim, ncol = dim * lag),
                  matrix(c(rep(1, 6), 0, rep(1, 8), 0, rep(1, 2)), nrow = dim, ncol = dim * lag)),
                dim = c(dim, dim * lag, 6))
datakey <- data.frame(Model = rep(NA, r * n), Length = rep(NA, r * n), GC = rep(NA, r * n))
data <- vector("list", length = r * n)
set.seed(1158)
seeds <- sample(0:10000, r * n, replace = TRUE)
for(model in 1:length(models)) {
  for(len in 1:length(lens)) {
    for(gc in 1:length(gcns)) {
      for(sim in 1:n) {
        count <- sim + n * (gc - 1) + n * length(gcns) * (len - 1) + n * length(gcns) * length(lens) * (model - 1)
        print(count)
        datakey$Model[count] <- models[model]
        datakey$Length[count] <- lens[len]
        datakey$GC[count] <- gcns[gc]
        set.seed(seeds[count])
        if(model == 1) {
          data[[count]] <- varData(L = lens[len], Dim = dim, Lag = lag, Causal = as.matrix(causal[,, gc]))
        } else {
          data[[count]] <- tarData(L = lens[len], Dim = dim, Lag = lag, RegimeL = regimelag, Causal = array(causal[,, gc], dim = c(dim, dim * lag, 2)))
        }
      }
    }
  }
}
save(datakey, file = "~/Wake Forest University Dropbox/Noah Gade/Research/NPGC/Simulations/Data/datakey.RData")
save(data, file = "~/Wake Forest University Dropbox/Noah Gade/Research/NPGC/Simulations/Data/data.RData")
