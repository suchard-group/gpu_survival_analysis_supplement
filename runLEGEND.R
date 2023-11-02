library(Andromeda)
library(Cyclops)

GpuDevice <- listGPUDevices()[1]

source("functions/LEGEND.R")

# Load data
pathToData <- "/home/jianxiao/LEGEND/"
covariateData <- loadAndromeda(paste0(pathToData,"cmd.zip"))
treatmentVarId <- readRDS(paste0(pathToData,"treatmentVarId.rds"))

# Run L1 Cox
resultCox <- runLEGEND(covariateData, treatmentVarId, modelType = "cox", runCPU = FALSE, GpuDevice = GpuDevice)

# Save results
saveRDS(resultCox, "resultCox.rds")
