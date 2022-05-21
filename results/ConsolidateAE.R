library(dplyr)
library(survival)

rd <- "/project/flatiron_ucc/programs/kylie/RunMe3"
rootdir <- file.path (rd, "final_results")
iter = 500
J = 10 #number of imputations
mechList <- c("MCAR", "MAR", "MNAR1", "MNAR2")
proportionList <- c(10, 30, 50)

consol = function (i, mdm, pc){
    truth <- read.csv(file.path(rd, "datasets/trueEff", mdm, as.character(pc), "propMiss_trueEffs1.csv"))
    truth = truth[-1]

    fileList = list.files(file.path(rootdir, paste0 ("AE", mdm), as.character(pc)), full.names = TRUE)

    fits = lapply (1:J, function (x){
        data <- read.csv(file.path(rootdir, paste0("AE", mdm), as.character(pc), paste0("result", i, "_num_", x, ".csv")), row.names = 1) # read in instead

        daeFit <- coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + var1 + var2, data = data)
        return (daeFit)
    })
    getCoefs = lapply (fits, function (x){
        vals = x$coefficients
    })
    getSEs = lapply (fits, function (x){
        vals = summary(x)$coefficients[, 3]
    })
    getVars = lapply (getSEs, function (x){
        vals = x^2
    })

    avgCoefs = Reduce('+', getCoefs)/J
    avgSEs = Reduce('+', getSEs)/J
    Ubar = Reduce('+', getVars)/J

    bias = avgCoefs - truth

    coefMat = do.call ('rbind', getCoefs)
    BVar = apply (coefMat, 2, var)
    # Variance formula for MI
    Rubins = Ubar + (1 + (1/J))*BVar
    se = sqrt (Rubins)

    CILower <- avgCoefs - (qnorm(.975) * se)
    CIUpper <- avgCoefs + (qnorm(.975) * se)

    coverage = unlist (lapply (1:length (avgCoefs), function (x){
        tr = as.numeric (as.matrix (truth))
        1*between (tr[x], CILower[x], CIUpper[x])
    }))

    mse = Rubins + bias^2
    return (list ('bias'=bias, 'se'=se, 'coverage'=coverage, 'mse'=mse))
}

final = function (iter, mdm, pc){
    op = lapply (1:iter, function (x) consol (x, mdm, pc))
    allBias = do.call ('rbind', lapply (op, function (x) x$bias))
    allSEs = do.call ('rbind', lapply (op, function (x) x$se))
    allCoverage = do.call ('rbind', lapply (op, function (x) x$coverage))
    allMSE = do.call ('rbind', lapply (op, function (x) x$mse))
    avgBias = apply (allBias, 2, mean)
    avgSE = apply (allSEs, 2, mean)
    avgCoverage = apply (allCoverage, 2, mean)
    avgMSE = apply (allMSE, 2, mean)
    df = cbind (avgBias, avgSE, avgCoverage, avgMSE)
    write.csv (df, file=file.path (rootdir, paste0 ('AEAvgs_', mdm, '_', as.character(pc), '.csv')))
    
    lq = function (x) quantile (x, p=.025)
    uq = function (x) quantile (x, p=.975)
    percLower = apply (allBias, 2, lq)
    percUpper = apply (allBias, 2, uq)
    PercentileAE = cbind (percLower, percUpper)
    write.csv (PercentileAE, file=file.path (rootdir, paste0 ('AEPercentiles_', mdm, '_', as.character(pc), '.csv')))
}

final (iter, 'MCAR', 10)
final (iter, 'MCAR', 30)
final (iter, 'MCAR', 50)
final (iter, 'MAR', 10)
final (iter, 'MAR', 30)
final (iter, 'MAR', 50)
final (iter, 'MNAR1', 10)
final (iter, 'MNAR1', 30)
final (iter, 'MNAR1', 50)
final (iter, 'MNAR2', 10)
final (iter, 'MNAR2', 30)
final (iter, 'MNAR2', 50)

