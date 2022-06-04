library(dplyr)
library(survival)

rd <- "/project/flatiron_ucc/programs/kylie/RunMe3"
args = commandArgs(trailingOnly = TRUE)
simple = as.character(args[1])=='simple'
opdir = 'final_results'
if (simple){
    opdir = 'final_simple'
}

iter = 500

makeTable <- function(iter, mdm='MCAR', pc=10) {
  keep = 1:iter
  fnames = list.files(file.path(rd, opdir, mdm, as.character (pc)), full.names = TRUE, pattern='result')[keep]
  df <- lapply(fnames, read.csv)

  # For percentile plots
  percOracle = lapply (df, function (x) x$biasComplete)
  percOracle = do.call ('rbind', percOracle)
  percCC = lapply (df, function (x) x$biasExclude)
  percCC = do.call ('rbind', percCC)
  percMICE = lapply (df, function (x) x$biasMICE)
  percMICE = do.call ('rbind', percMICE)
  percForest = lapply (df, function (x) x$biasForest)
  percForest = do.call ('rbind', percForest)

  # Averages
  vars = df[[1]]$X
  dfnum = lapply (df, function (x){
    newx = x[,-1]
    return (newx)
  })
  dfAvg <- Reduce("+", dfnum) / length(dfnum)
  rownames(dfAvg) <- vars
  colnames(dfAvg) = c ('biasOracle', 'biasCC', 'biasMICE', 'biasRF', 'seOracle', 'seCC', 'seMICE', 'seRF', 'coverageOracle', 'coverageCC', 'coverageMICE', 'coverageRF', 'mseOracle', 'mseCC', 'mseMICE', 'mseRF')

  percOracle = as.data.frame (percOracle)
  colnames (percOracle) = vars
  percOracle$method = 'ORACLE'
  percCC = as.data.frame (percCC)
  colnames (percCC) = vars
  percCC$method = 'CC'
  percMICE = as.data.frame (percMICE)
  colnames (percMICE) = vars
  percMICE$method = 'MICE'
  percRF = as.data.frame (percForest)
  colnames (percRF) = vars
  percRF$method = 'RF'
  PercentileDf <- rbind(percOracle, percCC, percMICE, percRF)

  write.csv (PercentileDf, file=file.path (rd, opdir, paste0 ("Percentiles_", mdm, "_", as.character(pc), '.csv')))
  write.csv (dfAvg, file=file.path (rd, opdir, paste0 ("Avgs_", mdm, "_", as.character(pc), '.csv')))

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
