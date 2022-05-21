library(dplyr)
library(survival)

rd <- "/project/flatiron_ucc/programs/kylie/RunMe3"
rootdir <- file.path (rd, "final_results")
iter = 500

makeTable <- function(iter, mdm='MCAR', pc=10) {
  keep = 1:iter
  fnames = list.files(file.path(rootdir, mdm, as.character (pc)), full.names = TRUE, pattern='result')[keep]
  df <- lapply(fnames, read.csv)

  # For percentile plots
  percOracle = lapply (df, function (x) x$biasComplete)
  percOracle = do.call ('cbind', percOracle)
  percCC = lapply (df, function (x) x$biasExclude)
  percCC = do.call ('cbind', percCC)
  percMICE = lapply (df, function (x) x$biasMICE)
  percMICE = do.call ('cbind', percMICE)
  percForest = lapply (df, function (x) x$biasForest)
  percForest = do.call ('cbind', percForest)

  # Averages
  vars = df[[1]]$X
  dfnum = lapply (df, function (x){
    newx = x[,-1]
    return (newx)
  })
  dfAvg <- Reduce("+", dfnum) / length(dfnum)
  rownames(dfAvg) <- vars
  colnames(dfAvg) = c ('biasOracle', 'biasCC', 'biasMICE', 'biasRF', 'seOracle', 'seCC', 'seMICE', 'seRF', 'coverageOracle', 'coverageCC', 'coverageMICE', 'coverageRF', 'mseOracle', 'mseCC', 'mseMICE', 'mseRF')

  lq = function (x) quantile (x, p=.025)
  uq = function (x) quantile (x, p=.975)
  percOracleLower = apply (percOracle, 1, lq)
  percOracleUpper = apply (percOracle, 1, uq)
  percCCLower = apply (percCC, 1, lq)
  percCCUpper = apply (percCC, 1, uq)
  percMICELower = apply (percMICE, 1, lq)
  percMICEUpper = apply (percMICE, 1, uq)
  percRFLower = apply (percForest, 1, lq)
  percRFUpper = apply (percForest, 1, uq)

  PercentileOracle = cbind (percOracleLower, percOracleUpper)
  PercentileCC = cbind (percCCLower, percCCUpper)
  PercentileMICE = cbind (percMICELower, percMICEUpper)
  PercentileRF = cbind (percRFLower, percRFUpper)
  PercentileDf <- cbind(PercentileOracle, PercentileCC, PercentileMICE, PercentileRF)

  write.csv (PercentileDf, file=file.path (rootdir, paste0 ("Percentiles_", mdm, "_", as.character(pc), '.csv')))
  write.csv (dfAvg, file=file.path (rootdir, paste0 ("Avgs_", mdm, "_", as.character(pc), '.csv')))

  return(list("perc" = PercentileDf, "avg" = dfAvg))
}

MCAR10 <- makeTable(iter, 'MCAR', 10)
MCAR30 <- makeTable(iter, 'MCAR', 30)
MCAR50 <- makeTable(iter, 'MCAR', 50)

MAR10 <- makeTable(iter, 'MAR', 10)
MAR30 <- makeTable(iter, 'MAR', 30)
MAR50 <- makeTable(iter, 'MAR', 50)

MNAR10a <- makeTable(iter, 'MNAR1', 10)
MNAR30a <- makeTable(iter, 'MNAR1', 30)
MNAR50a <- makeTable(iter, 'MNAR1', 50)

MNAR10b <- makeTable(iter, 'MNAR2', 10)
MNAR30b <- makeTable(iter, 'MNAR2', 30)
MNAR50b <- makeTable(iter, 'MNAR2', 50)
