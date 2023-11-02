library(Matrix)
library(Andromeda)
library(Cyclops)
library(testthat)
library(survival)
library(cmprsk)

GpuDevice <- listGPUDevices()[1]

source("functions/synthetic_experiments.R")

N <- c(100000, 500000, 800000, 1000000)
threadsPool <- c(1L, 2L, 4L, 8L)

# L1 Cox
time_cox_gpu <- rep(0, length(N) * length(threadsPool))
time_cox_cpu <- rep(0, length(N) * length(threadsPool))
for (i in 1:length(N)) {
	times <- runAutoL1(nrows = N[i], modelType = "cox", GpuDevice = GpuDevice, threadsPool = threadsPool)
	s <- (i-1) * length(N) + 1
	time_cox_gpu[s : (s+length(threadsPool)-1)] <- times$gpu
	time_cox_cpu[s : (s+length(threadsPool)-1)] <- times$cpu
}
print(time_cox_gpu)
print(time_cox_cpu)

# L1 FG
time_fg_gpu <- rep(0, length(N) * length(threadsPool))
time_fg_cpu <- rep(0, length(N) * length(threadsPool))
for (i in 1:length(N)) {
    times <- runAutoL1(nrows = N[i], modelType = "fgr", GpuDevice = GpuDevice, threadsPool = threadsPool)
    s <- (i-1) * length(N) + 1
    time_fg_gpu[s : (s+length(threadsPool)-1)] <- times$gpu
    time_fg_cpu[s : (s+length(threadsPool)-1)] <- times$cpu
}
print(time_fg_gpu)
print(time_fg_cpu)


