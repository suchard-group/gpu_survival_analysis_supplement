library(Andromeda)
library(Cyclops)
library(testthat)
library(survival)
library(cmprsk)

runLEGEND <- function(covariateData, treatmentVarId, dp = TRUE, seed = 123, modelType = "cox", runGPU = TRUE, runCPU = TRUE, GpuDevice, threads = 1) {
    
    set.seed(seed)
    
    ######################################
    # Run Experiments
    ######################################

    if (dp) {
        sd = 0.0001
        fp = 64
        tolerance <- 1E-4
    } else {
        sd = 0.1
        fp = 32
        tolerance <- 1E-3
    }

    # Create prior and control for cross-validation
    prior <- createPrior(priorType = "laplace",
                         variance = 1,
                         exclude = treatmentVarId,
                         useCrossValidation = TRUE)
    control <- createControl(cvType = "auto",
                             seed = seed,
                             threads = threads,
                             startingVariance = 0.01,
                             tolerance = 2e-07,
                             cvRepetitions = 10,
                             noiseLevel = "silent")

    # Run L1
    if (runGPU) {
        cyclopsDataG <- convertToCyclopsData(outcomes = covariateData$outcomes,
                                            covariates = covariateData$covariates,
                                            modelType = modelType,
                                            checkRowIds = FALSE,
                                            quiet = TRUE)
        start <- Sys.time()
        fitGPU <- fitCyclopsModel(cyclopsDataG,
                               prior = prior,
                               control = control,
                               computeDevice = GpuDevice)
        delta_g <- Sys.time() - start
        writeLines(paste("GPU sparse analysis took", signif(delta_g,3), attr(delta_g,"units"),
             "(",
             signif(fitGPU$timeFit,3), attr(fitGPU$timeFit,"units"),
             ")"))
             
        ci <- confint(fitGPU, parm = treatmentVarId, includePenalty = TRUE)
    }

    if (runCPU) {
        cyclopsDataC <- convertToCyclopsData(outcomes = covariateData$outcomes,
                                            covariates = covariateData$covariates,
                                            modelType = modelType,
                                            checkRowIds = FALSE,
                                            quiet = TRUE)
        start <- Sys.time()
        fitCPU <- fitCyclopsModel(cyclopsDataC,
                               prior = prior,
                               control = control)
        delta_c <- Sys.time() - start
        writeLines(paste("CPU sparse analysis took", signif(delta_c,3), attr(delta_c,"units"),
             "(",
             signif(fitCPU$timeFit,3), attr(fitCPU$timeFit,"units"),
             ")"))
        
        ci <- confint(fitCPU, parm = treatmentVarId, includePenalty = TRUE)
    }


    ######################################
    # Check results
    ######################################

    # compare GPU v.s. CPU
    if (runGPU && runCPU) {
        expect_equal(coef(fitGPU), coef(fitCPU), tolerance = tolerance)
        expect_equal(logLik(fitGPU), logLik(fitCPU), tolerance = tolerance)
    }

    # return results
    if (runGPU) {
        return(list(fit = fitGPU, ci = ci))
    } else {
        return(list(fit = fitCPU, ci = ci))
    }
}

