library(mice)
library(CALIBERrfimpute)
library(survival)
library(dplyr)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
mech = as.character(args[2])
simple = as.character(args[3])=='TRUE'

rootdir <- "/project/flatiron_ucc/programs/kylie/RunMe4"
if (simple){
  rd = file.path (rootdir, 'final_simple')
} else{
  rd = file.path (rootdir, 'final_results')
}
system (paste0 ('mkdir ', rd))
system (paste0 ('mkdir ', file.path (rd, mech)))

proportionList <- c(10, 30, 50)
set.seed(100)
for (k in 1:length(proportionList)) {
  propName <- proportionList[k]
  system (paste0 ('mkdir ', file.path (rd, mech, propName)))
  print(propName)
  i <- as.numeric(args[1])
  mData <- read.csv(file.path(rootdir, "datasets/mDats", mech, propName, paste0("mData", i, ".csv")))
  cData <- read.csv(file.path(rootdir, "datasets/cDats", mech, propName, paste0("cData", i, ".csv")))
  cData <- cData[, -1]
  if (simple){
    truth <- read.csv(file.path(rootdir, "datasets/trueEff_simple", mech, propName, "propMiss_trueEffs1.csv"))
  } else{
    truth <- read.csv(file.path(rootdir, "datasets/trueEff", mech, propName, "propMiss_trueEffs1.csv"))
  }
  propMiss = truth[1]
  truth = as.numeric (as.matrix (truth[-1]))
#  trueCoefs <- truth[names (truth)%in%c("TREAT", "b.ecogvalue", "var1", "var2")]

  mCat <- mData
  mCat$race <- rep(NA, nrow(mData))
  mCat$race[which(mData$reth_black == 1)] <- "black"
  mCat$race[which(mData$reth_hisp == 1)] <- "hisp"
  mCat$race[which(mData$reth_oth == 1)] <- "oth"
  mCat$race[which(mData$reth_black == 0 & mData$reth_hisp == 0 & mData$reth_oth == 0)] <- "white"
  mCat$race <- factor(mCat$race)

  mCat$site <- rep(NA, nrow(mData))
  mCat$site[which(mData$site_ureter == 1)] <- "ureter"
  mCat$site[which(mData$site_renal == 1)] <- "renal"
  mCat$site[which(mData$site_urethra == 1)] <- "urethra"
  mCat$site[which(mData$site_ureter == 0 & mData$site_renal == 0 & mData$site_urethra == 0)] <- "bladder"
  mCat$site <- factor(mCat$site)

  mData$reth_black <- NULL
  mData$reth_hisp <- NULL
  mData$reth_oth <- NULL
  mData$race <- mCat$race

  mData$site_ureter <- NULL
  mData$site_renal <- NULL
  mData$site_urethra <- NULL
  mData$site <- mCat$site

  imp <- mice(mData, maxit = 0)
  predM <- imp$predictorMatrix

  # Setting values of variables I'd like to leave out to 0 in the predictor matrix
  predM[, c("time")] <- 0

  imp2 <- mice(mData,
    m = 10, maxit = 5,
    predictorMatrix = predM,
    print = FALSE
  )

  anesimp_long <- mice::complete(imp2, action = "long", include = TRUE)

  reth_black <- 1 * (anesimp_long$race == "black")
  reth_hisp <- 1 * (anesimp_long$race == "hisp")
  reth_oth <- 1 * (anesimp_long$race == "oth")
  site_ureter <- 1 * (anesimp_long$site == "ureter")
  site_renal <- 1 * (anesimp_long$site == "renal")
  site_urethra <- 1 * (anesimp_long$site == "urethra")
  anesimp_long$race <- NULL
  anesimp_long$site <- NULL

  anesimp_long <- add_column(anesimp_long, reth_black = reth_black, .after = "genderf")
  anesimp_long <- add_column(anesimp_long, reth_hisp = reth_hisp, .after = "reth_black")
  anesimp_long <- add_column(anesimp_long, reth_oth = reth_oth, .after = "reth_hisp")

  anesimp_long <- add_column(anesimp_long, site_ureter = site_ureter, .after = "practypec")
  anesimp_long <- add_column(anesimp_long, site_renal = site_renal, .after = "site_ureter")
  anesimp_long <- add_column(anesimp_long, site_urethra = site_urethra, .after = "site_renal")

  write.csv(anesimp_long, file.path(rd, mech, propName, paste0("mice_", i, ".csv")))
  anesimp_long_mids <- as.mids(anesimp_long)

  if (simple){
    fitimp1 <- with(
      anesimp_long_mids,
      coxph(Surv(time, event) ~ treat + b.ecogvalue)
    )
  } else{
    fitimp1 <- with(
      anesimp_long_mids,
      coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + var1 + var2)
    )
  }

  pool <- summary(pool(fitimp1))

  terms = pool$term
  biasMICE <- as.numeric (pool$estimate - truth)
  SEMICE <- pool$std.error
  VarMICE <- SEMICE^2
  CIMICELower <- pool$estimate - (qnorm(.975) * (pool$std.error))
  CIMICEUpper <- pool$estimate + (qnorm(.975) * (pool$std.error))
  coverageMICE <- unlist (lapply (1:length (truth), function (x) 1 * between(truth[x], CIMICELower[x], CIMICEUpper[x])))
  mseMICE <- VarMICE + biasMICE^2

  print("RF imputation")

  mData <- read.csv(file.path(rootdir, "datasets/mDats", mech, propName, paste0("mData", i, ".csv")))

  mCat <- mData
  mCat$race <- rep(NA, nrow(mData))
  mCat$race[which(mData$reth_black == 1)] <- "black"
  mCat$race[which(mData$reth_hisp == 1)] <- "hisp"
  mCat$race[which(mData$reth_oth == 1)] <- "oth"
  mCat$race[which(mData$reth_black == 0 & mData$reth_hisp == 0 & mData$reth_oth == 0)] <- "white"
  mCat$race <- factor(mCat$race)

  mCat$site <- rep(NA, nrow(mData))
  mCat$site[which(mData$site_ureter == 1)] <- "ureter"
  mCat$site[which(mData$site_renal == 1)] <- "renal"
  mCat$site[which(mData$site_urethra == 1)] <- "urethra"
  mCat$site[which(mData$site_ureter == 0 & mData$site_renal == 0 & mData$site_urethra == 0)] <- "bladder"
  mCat$site <- factor(mCat$site)

  mData$reth_black <- NULL
  mData$reth_hisp <- NULL
  mData$reth_oth <- NULL
  mData$race <- mCat$race

  mData$site_ureter <- NULL
  mData$site_renal <- NULL
  mData$site_urethra <- NULL
  mData$site <- mCat$site

  mData$event <- as.numeric(mData$event)
  mData$genderf <- as.factor(mData$genderf)
  mData$practypec <- as.factor(mData$practypec)
  mData$smokey <- as.factor(mData$smokey)
  mData$dgradeh <- as.factor(mData$dgradeh)
  mData$surgery <- as.factor(mData$surgery)
  mData$treat <- as.factor(mData$treat)
  mData$b.ecogvalue <- as.numeric(mData$b.ecogvalue)
  mData$var1 <- as.numeric(mData$var1)
  mData$var2 <- as.numeric(mData$var2)
  mData$race <- as.factor(mData$race)
  mData$site <- as.factor(mData$site)
  keepTime = rep (mData$time, 11)
  mData$time = NULL

  forest <- mice(mData, method = c("rfcat", "rfcont", "rfcat", "rfcont", "rfcat", "rfcat", "rfcat", "rfcat", "rfcat", "rfcont", "rfcont", "rfcont", "rfcat", "rfcat"), m = 10, maxit = 5)

  anesimp_long <- mice::complete(forest, action = "long", include = TRUE)
  anesimp_long$b.ecogvalue <- round(anesimp_long$b.ecogvalue)
  anesimp_long <- add_column(anesimp_long, time = keepTime, .after = "hazard")

  reth_black <- 1 * (anesimp_long$race == "black")
  reth_hisp <- 1 * (anesimp_long$race == "hisp")
  reth_oth <- 1 * (anesimp_long$race == "oth")
  site_ureter <- 1 * (anesimp_long$site == "ureter")
  site_renal <- 1 * (anesimp_long$site == "renal")
  site_urethra <- 1 * (anesimp_long$site == "urethra")
  anesimp_long$race <- NULL
  anesimp_long$site <- NULL

  anesimp_long <- add_column(anesimp_long, reth_black = reth_black, .after = "genderf")
  anesimp_long <- add_column(anesimp_long, reth_hisp = reth_hisp, .after = "reth_black")
  anesimp_long <- add_column(anesimp_long, reth_oth = reth_oth, .after = "reth_hisp")

  anesimp_long <- add_column(anesimp_long, site_ureter = site_ureter, .after = "practypec")
  anesimp_long <- add_column(anesimp_long, site_renal = site_renal, .after = "site_ureter")
  anesimp_long <- add_column(anesimp_long, site_urethra = site_urethra, .after = "site_renal")

  write.csv(anesimp_long, file.path(rd, mech, propName, paste0("rf_", i, ".csv")))

  anesimp_long_mids <- as.mids(anesimp_long)
  if (simple){
    fitimp1 <- with(
      anesimp_long_mids,
      coxph(Surv(time, event) ~ treat + b.ecogvalue)
    )
  } else{
    fitimp1 <- with(
      anesimp_long_mids,
      coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + var1 + var2)
    )
  }

  pool <- summary(pool(fitimp1))
  terms = pool$term
  biasForest <- as.numeric (pool$estimate - truth)
  SEForest <- pool$std.error
  VarForest <- SEForest^2
  CIForestLower <- pool$estimate - (qnorm(.975) * (pool$std.error))
  CIForestUpper <- pool$estimate + (qnorm(.975) * (pool$std.error))
  coverageForest <- unlist (lapply (1:length (truth), function (x) 1 * between(truth[x], CIForestLower[x], CIForestUpper[x])))
  mseForest <- VarForest + biasForest^2

  print("Oracle")

  # compare complete dataset to OS1 estimates
  if (simple){
    CompleteFit <- coxph(Surv(time, event) ~ treat + b.ecogvalue, data = cData)
  } else{
    CompleteFit <- coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + var1 + var2, data = cData)
  }

  biasComplete <- as.numeric (CompleteFit$coefficients - truth)
  SEComplete <- summary(CompleteFit)$coefficients[,3]
  VarComplete <- SEComplete^2
  CICompleteLower <- CompleteFit$coefficients - (qnorm(.975) * SEComplete)
  CICompleteUpper <- CompleteFit$coefficients + (qnorm(.975) * SEComplete)
  coverageComplete <- unlist (lapply (1:length (truth), function (x) 1 * between(truth[x], CICompleteLower[x], CICompleteUpper[x])))
  mseComplete <- VarComplete + biasComplete^2

  print("Complete Case")

  # compare exclude missingness dataset to OS1 estimates
  mData <- read.csv(file.path(rootdir, "datasets/mDats", mech, propName, paste0("mData", i, ".csv")))

  ExcludeData <- mData[complete.cases(mData), ]

  if (simple){
    ExcludeFit <- coxph(Surv(time, event) ~ treat + b.ecogvalue, data = ExcludeData)
  } else{
    ExcludeFit <- coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + var1 + var2, data = ExcludeData)
  }

  biasExclude <- as.numeric (ExcludeFit$coefficients - truth)
  SEExclude <- summary(ExcludeFit)$coefficients[,3]
  VarExclude <- SEExclude^2
  CIExcludeLower <- ExcludeFit$coefficients - (qnorm(.975) * SEExclude)
  CIExcludeUpper <- ExcludeFit$coefficients + (qnorm(.975) * SEExclude)
  coverageExclude <- unlist (lapply (1:length (truth), function (x) 1 * between(truth[x], CIExcludeLower[x], CIExcludeUpper[x])))
  mseExclude <- VarExclude + biasExclude^2

  bias <- cbind (biasComplete, biasExclude, biasMICE, biasForest)
  se <- cbind (SEComplete, SEExclude, SEMICE, SEForest)
  mse <- cbind (mseComplete, mseExclude, mseMICE, mseForest)
  coverage <- cbind(coverageComplete, coverageExclude, coverageMICE, coverageForest)

  results <- cbind(bias, se, coverage, mse)
  write.csv(results, file.path(rd, mech, propName, paste0("result", i, ".csv")))
}
