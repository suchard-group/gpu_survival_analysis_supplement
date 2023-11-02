library(Matrix)
library(Andromeda)
library(Cyclops)
library(testthat)
library(survival)
library(cmprsk)

GpuDevice <- listGPUDevices()[1]

source("functions/synthetic_experiments.R")

N <- c(100000, 177828, 316228, 562341, 1000000)

# L1 Cox
time_cox <- rep(0, 2*length(N))
for (i in 1:length(N)) {
        time_cox[c(i, i+length(N))] <- runFixedL1(nrows = N[i], modelType = "cox", GpuDevice = GpuDevice)
}
print(time_cox)

# L1 FG
time_fg <- rep(0, 2*length(N))
for (i in 1:length(N)) {
        time_fg[c(i, i+length(N))] <- runFixedL1(nrows = N[i], modelType = "fgr", GpuDevice = GpuDevice)
}
print(time_fg)


