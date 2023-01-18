library(Matrix)
library(Andromeda)
library(Cyclops)
library(testthat)
library(survival)
library(cmprsk)

GpuDevice <- listOpenCLDevices()[1] # "Quadro GV100"
#GpuDevice <- listOpenCLDevices()[2] # "TITAN V"
#GpuDevice <- listOpenCLDevices()[3] # "GeForce RTX 2080"

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


