library(Matrix)
library(Andromeda)
library(Cyclops)
library(testthat)
library(survival)
library(cmprsk)

simCox <- function(nrows, ncovars = 1000, beta, nstrata = 1, sparseness = 0.95, zeroEffectSizeProp = 0.8, ind = TRUE, seed = 123) {

    set.seed(seed)
 
    effectSizes <- data.frame(covariateId=1:ncovars,rr=exp(beta))

    # generate covariate matrix
    covariates <- rsparsematrix(nrows, ncovars, nnz = nrows*ncovars*(1-sparseness), rand.x = rnorm)
    covariates <- summary(covariates)

    if (ind) {
        covariates$x <- 1 # indicator
    }
    colnames(covariates) <- c("rowId", "covariateId", "covariateValue")

    # calculate sum of exb for each observation
    outcomes <- data.frame(rowId = 1:nrows, stratumId = round(runif(nrows,min=1,max=nstrata)), y=0)
    covariates <- merge(covariates,outcomes[,c("rowId","stratumId")])
    rowId_to_rr <- aggregate(rr ~ rowId, data=merge(covariates,effectSizes), prod)
    outcomes <- merge(outcomes,rowId_to_rr,all.x=TRUE)
    outcomes$rr[is.na(outcomes$rr)] <- 1

    # survival
    strataBackgroundProb <- runif(nstrata,min=0.01,max=0.03)
    outcomes$rate <-  strataBackgroundProb[outcomes$stratumId] * outcomes$rr
    outcomes$timeToOutcome <- 1+round(rexp(n=nrow(outcomes),outcomes$rate))
    outcomes$timeToCensor <- 1+round(runif(n=nrow(outcomes),min=0,max=499))
    outcomes$time <- outcomes$timeToOutcome
    outcomes$time[outcomes$timeToCensor < outcomes$timeToOutcome] <- outcomes$timeToCensor[outcomes$timeToCensor < outcomes$timeToOutcome]
    outcomes$y <- as.integer(outcomes$timeToCensor > outcomes$timeToOutcome)

    # wrap up
    outcomes <- outcomes[order(outcomes$stratumId,outcomes$rowId),]
    covariates <- covariates[order(covariates$stratumId,covariates$rowId,covariates$covariateId),]

    # break ties
    n <- dim(outcomes)[1]
    outcomes$time <- outcomes$time + rnorm(n, mean = 0, sd = 1E-3)

    sparsenessReal <- 1-(nrow(covariates)/(nrows*ncovars))
    writeLines(paste("Simulated sparseness = ",sparsenessReal*100,"%"))

    # andr
    andr <- andromeda(out = outcomes, cov = covariates)
    
    return(andr)
}

simFG <- function(nrows, ncovars = 1000, beta1, nstrata = 1, sparseness = 0.95, zeroEffectSizeProp = 0.8, ind = TRUE, seed = 123) {

    set.seed(seed)

    p = 0.5
    u.min = 0
    u.max = 1
    
    # generate beta
    beta2 <- -beta1

    # generate X
    X <- rsparsematrix(nrows, ncovars, nnz = nrows*ncovars*(1-sparseness), rand.x = rnorm)
    if (ind) {
        X@x <- rep(1, length(X@i)) # indicator
    }

    #- Generate indicator for cause
    pr = (1 - p)^exp(X %*% beta1)
    c.ind <- 1 + rbinom(nrows, 1, prob = pr@x)

    # conditional on cause indicators, we simulate the model.
    ftime <- numeric(nrows)
    eta1 <- X[c.ind == 1, ] %*% beta1 #linear predictor for cause on interest
    eta2 <- X[c.ind == 2, ] %*% beta2 #linear predictor for competing risk

    u1 <- runif(length(eta1@x))
    t1 <- -log(1 - (1 - (1 - u1 * (1 - (1 - p)^exp(eta1)))^(1 / exp(eta1@x))) / p)
    t2 <- rexp(length(eta2@x), rate = exp(eta2@x))
    ci <- runif(nrows, min = u.min, max = u.max) #simulate censoring times

    ftime[c.ind == 1] <- t1
    ftime[c.ind == 2] <- t2
    ftime <- pmin(ftime, ci)
    fstatus <- ifelse(ftime == ci, 0, 1)
    fstatus <- fstatus * c.ind

    # break ties
    tran = log(10^(100) * ftime + 2) +  rnorm(nrows, mean = 0, sd = 0.0001)

    # fg data
    fgDat <- Cyclops:::getFineGrayWeights(tran, fstatus)
    outcomes <- data.frame(rowId = 1:nrows, time = tran, y = fstatus, censorWeights = fgDat$weights)

    # covariates
    covariates <- data.frame(summary(X))
    colnames(covariates) = c("rowId", "covariateId", "covariateValue")

    sparsenessReal <- 1-(nrow(covariates)/(nrows*ncovars))
    writeLines(paste("Simulated sparseness = ",sparsenessReal*100,"%"))

    # andr
    andr <- andromeda(out = outcomes, cov = covariates)

    return(andr)
}

runFixedL1 <- function(nrows, ncovars = 1000, nstrata = 1, sparseness = 0.95, zeroEffectSizeProp = 0.8, effectSizeSd = 1, ind = TRUE, dp = TRUE, seed = 123, modelType = "cox", runGPU = TRUE, runCPU = TRUE, GpuDevice) {
    
    set.seed(seed)

    ######################################
    # Simulate Data
    ######################################

    # initialize beta
    sd <- rep(effectSizeSd, ncovars) * rbinom(ncovars, 1, 1 - zeroEffectSizeProp)
    beta <- rnorm(ncovars,mean=0,sd=sd)

    if (modelType == "cox") {
        andr <- simCox(nrows = nrows, ncovars = ncovars, beta = beta, nstrata = nstrata, sparseness = sparseness, zeroEffectSizeProp = zeroEffectSizeProp, ind = ind, seed = seed)
    } else if (modelType == "fgr") {
        andr <- simFG(nrows = nrows, ncovars = ncovars, beta1 = beta, nstrata = nstrata, sparseness = sparseness, zeroEffectSizeProp = zeroEffectSizeProp, ind = ind, seed = seed)
    }
    
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

    if (runGPU) {
        # sparse gpu
        sparse_gpu_ptr <- convertToCyclopsData(andr$out, andr$cov, modelType = modelType, floatingPoint = fp)
        start <- Sys.time()
        sparse_gpu_fit <- fitCyclopsModel(sparse_gpu_ptr, prior = createPrior("laplace", variance = 1.414), computeDevice = GpuDevice)
        delta_g <- Sys.time() - start
        writeLines(paste("GPU sparse analysis took", signif(delta_g,3), attr(delta_g,"units"),
             "(",
             signif(sparse_gpu_fit$timeFit,3), attr(sparse_gpu_fit$timeFit,"units"),
             ")"))
    }

    if (runCPU) {
        # sparse cpu
        sparse_cpu_ptr <- convertToCyclopsData(andr$out, andr$cov, modelType = modelType, floatingPoint = fp)
        start <- Sys.time()
        sparse_cpu_fit <- fitCyclopsModel(sparse_cpu_ptr, prior = createPrior("laplace", variance = 1.414))
        delta_c <- Sys.time() - start
        writeLines(paste("CPU sparse analysis took", signif(delta_c,3), attr(delta_c,"units"),
             "(",
             signif(sparse_cpu_fit$timeFit,3), attr(sparse_cpu_fit$timeFit,"units"),
             ")"))
    }

    # close Andromeda
    close(andr)

    ######################################
    # Check results
    ######################################

    # compare GPU v.s. CPU
    if (runGPU && runCPU) {
        expect_equal(coef(sparse_gpu_fit), coef(sparse_cpu_fit), tolerance = tolerance)
        expect_equal(logLik(sparse_gpu_fit), logLik(sparse_cpu_fit), tolerance = tolerance)
    }

    # compare with true beta
    if (runGPU && runCPU) {
        mse <- mean((beta - coef(sparse_gpu_fit))^2)
    } else if (runGPU) {
        mse <- mean((beta - coef(sparse_gpu_fit))^2)
    } else if (runCPU) {
        mse <- mean((beta - coef(sparse_cpu_fit))^2)
    }
    writeLines(paste("MSE = ", mse))

    # return time
    if (runGPU && runCPU) {
        return(c(delta_g, delta_c))
    } else if (runGPU) {
        return(delta_g)
    } else if (runCPU) {
        return(delta_c)
    }
}


runAutoL1 <- function(nrows, ncovars = 1000, nstrata = 1, sparseness = 0.95, zeroEffectSizeProp = 0.8, effectSizeSd = 1, ind = TRUE, dp = TRUE, seed = 123, modelType = "cox", runGPU = TRUE, runCPU = TRUE, GpuDevice, threadsPool = c(1L, 2L, 4L, 8L)) {

    set.seed(seed)

    ######################################
    # Simulate Data
    ######################################

    # initialize beta
    sd <- rep(effectSizeSd, ncovars) * rbinom(ncovars, 1, 1 - zeroEffectSizeProp)
    beta <- rnorm(ncovars,mean=0,sd=sd)

    if (modelType == "cox") {
        andr <- simCox(nrows = nrows, ncovars = ncovars, beta = beta, nstrata = nstrata, sparseness = sparseness, zeroEffectSizeProp = zeroEffectSizeProp, ind = ind, seed = seed)
    } else if (modelType == "fgr") {
        andr <- simFG(nrows = nrows, ncovars = ncovars, beta1 = beta, nstrata = nstrata, sparseness = sparseness, zeroEffectSizeProp = zeroEffectSizeProp, ind = ind, seed = seed)
    }

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

    if (runGPU) {
        time_gpu <- rep(0, length(threadsPool))
    }
    if (runCPU) {
        time_cpu <- rep(0, length(threadsPool))
    }

    for (i in 1:length(threadsPool)) {

        # create prior
        prior <- createPrior("laplace", useCrossValidation = TRUE)
        control <- createControl(noiseLevel = "silent",
                                  cvType = "auto",
                                  fold = 10,
                                  cvRepetitions = 10,
                                  startingVariance = 0.01,
                                  tolerance = 2e-07,
                                  threads = threadsPool[i],
                                  seed = seed)

        if (runGPU) {
            # sparse gpu
            sparse_gpu_ptr <- convertToCyclopsData(andr$out, andr$cov, modelType = modelType, floatingPoint = fp)
            start <- Sys.time()
            sparse_gpu_fit <- fitCyclopsModel(sparse_gpu_ptr, prior = prior, control=control, computeDevice = GpuDevice)
            time_gpu[i] <- Sys.time() - start
            writeLines(paste("GPU sparse analysis took", signif(time_gpu[i],3), attr(time_gpu[i],"units"),
                 "(",
                 signif(sparse_gpu_fit$timeFit,3), attr(sparse_gpu_fit$timeFit,"units"),
                 ") using ", threadsPool[i], " threads"))
        }

        if (runCPU) {
            # sparse cpu
            sparse_cpu_ptr <- convertToCyclopsData(andr$out, andr$cov, modelType = modelType, floatingPoint = fp)
            start <- Sys.time()
            sparse_cpu_fit <- fitCyclopsModel(sparse_cpu_ptr, prior = prior, control=control)
            time_cpu[i] <- Sys.time() - start
            writeLines(paste("CPU sparse analysis took", signif(time_cpu[i],3), attr(time_cpu[i],"units"),
                 "(",
                 signif(sparse_cpu_fit$timeFit,3), attr(sparse_cpu_fit$timeFit,"units"),
                 ") using ", threadsPool[i], " threads"))
        }

        ######################################
        # Check results
        ######################################

        # compare GPU v.s. CPU
        if (runGPU && runCPU) {
            expect_equal(coef(sparse_gpu_fit), coef(sparse_cpu_fit), tolerance = tolerance)
            expect_equal(logLik(sparse_gpu_fit), logLik(sparse_cpu_fit), tolerance = tolerance)
        }

        # compare with true beta
        writeLines(paste0(sum(beta == 0), " out of ", ncovars, " true beta are zero"))
        if (runGPU && runCPU) {
            writeLines(paste0(sum(coef(sparse_gpu_fit) == 0), " out of ", ncovars, " fitted coefficients are zero"))
        } else if (runGPU) {
            writeLines(paste0(sum(coef(sparse_gpu_fit) == 0), " out of ", ncovars, " fitted coefficients are zero"))
        } else if (runCPU) {
            writeLines(paste0(sum(coef(sparse_cpu_fit) == 0), " out of ", ncovars, " fitted coefficients are zero"))
        }
    }

    # close Andromeda
    close(andr)

    # return time
    if (runGPU && runCPU) {
        return(list("gpu" = time_gpu, "cpu" = time_cpu))
    } else if (runGPU) {
        return(time_gpu)
    } else if (runCPU) {
        return(time_cpu)
    }
}

