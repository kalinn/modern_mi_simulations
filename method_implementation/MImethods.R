library(mice)
library(CALIBERrfimpute)
library(survival)
library(dplyr)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
mech = as.character(args[2])

rootdir <- "/project/flatiron_ucc/programs/kylie/RunMe2"
system (paste0 ('mkdir ', file.path (rootdir, 'final_results')))

proportionList <- c(10, 30, 50)
set.seed(100)
for (k in 1:length(proportionList)) {
  propName <- proportionList[k]
  print(propName)
  i <- as.numeric(args[1])
  mData <- read.csv(file.path(rootdir, "datasets/mDats", mech, propName, paste0("mData", i, ".csv")))
  cData <- read.csv(file.path(rootdir, "datasets/cDats", mech, propName, paste0("cData", i, ".csv")))
  cData <- cData[, -1]
  truth <- read.csv(file.path(rootdir, "datasets/trueEff", mech, propName, "propMiss_trueEffs1.csv"))
  trueCoefs <- truth[names (truth)%in%c("TREAT", "b.ecogvalue", "var1", "var2")]

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

  anesimp_long_mids <- as.mids(anesimp_long)

  fitimp1 <- with(
    anesimp_long_mids,
    coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + var1 + var2)
  )

  pool <- summary(pool(fitimp1))

  TreatRow <- which(pool$term == "treat")
  EcogRow <- which(pool$term == "b.ecogvalue")
  Var1Row <- which(pool$term == "var1")
  Var2Row <- which(pool$term == "var2")

  MICEtreat <- pool$estimate[TreatRow]
  MICEEcog <- pool$estimate[EcogRow]
  MICEVar1 <- pool$estimate[Var1Row]
  MICEVar2 <- pool$estimate[Var2Row]

  biasTreatMICE <- as.numeric(MICEtreat - trueCoefs[1])
  biasEcogMICE <- as.numeric(MICEEcog - trueCoefs[2])
  biasVar1MICE <- as.numeric(MICEVar1 - trueCoefs[3])
  biasVar2MICE <- as.numeric(MICEVar2 - trueCoefs[4])

  VarTreatMICE <- (pool$std.error[TreatRow])^2
  VarEcogMICE <- (pool$std.error[EcogRow])^2
  VarVar1MICE <- (pool$std.error[Var1Row])^2
  VarVar2MICE <- (pool$std.error[Var2Row])^2

  CIMICELowerTreat <- MICEtreat - (qnorm(.975) * (pool$std.error[TreatRow]))
  CIMICEUpperTreat <- MICEtreat + (qnorm(.975) * (pool$std.error[TreatRow]))

  CIMICELowerEcog <- MICEEcog - (qnorm(.975) * (pool$std.error[EcogRow]))
  CIMICEUpperEcog <- MICEEcog + (qnorm(.975) * (pool$std.error[EcogRow]))

  CIMICELowerVar1 <- MICEVar1 - (qnorm(.975) * (pool$std.error[Var1Row]))
  CIMICEUpperVar1 <- MICEVar1 + (qnorm(.975) * (pool$std.error[Var1Row]))

  CIMICELowerVar2 <- MICEVar2 - (qnorm(.975) * (pool$std.error[Var2Row]))
  CIMICEUpperVar2 <- MICEVar2 + (qnorm(.975) * (pool$std.error[Var2Row]))

  covereageMICETreat <- 1 * between(c(trueCoefs[1]), CIMICELowerTreat, CIMICEUpperTreat)
  covereageMICEEcog <- 1 * between(c(trueCoefs[2]), CIMICELowerEcog, CIMICEUpperEcog)
  covereageMICEVar1 <- 1 * between(c(trueCoefs[3]), CIMICELowerVar1, CIMICEUpperVar1)
  covereageMICEVar2 <- 1 * between(c(trueCoefs[4]), CIMICELowerVar2, CIMICEUpperVar2)

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

  forest <- mice(mData, method = c("rfcat", "rfcont", "rfcont", "rfcat", "rfcont", "rfcat", "rfcat", "rfcat", "rfcat", "rfcat", "rfcont", "rfcont", "rfcont", "rfcat", "rfcat"), m = 10, maxit = 5)

  anesimp_long <- mice::complete(forest, action = "long", include = TRUE)
  anesimp_long$b.ecogvalue <- round(anesimp_long$b.ecogvalue)

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

  anesimp_long_mids <- as.mids(anesimp_long)
  fitimp1 <- with(
    anesimp_long_mids,
    coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + var1 + var2)
  )

  pool <- summary(pool(fitimp1))

  TreatRow <- which(pool$term == "treat1")
  EcogRow <- which(pool$term == "b.ecogvalue")
  Var1Row <- which(pool$term == "var1")
  Var2Row <- which(pool$term == "var2")

  Foresttreat <- pool$estimate[TreatRow]
  ForestEcog <- pool$estimate[EcogRow]
  ForestVar1 <- pool$estimate[Var1Row]
  ForestVar2 <- pool$estimate[Var2Row]

  biasTreatForest <- as.numeric(Foresttreat - trueCoefs[1])
  biasEcogForest <- as.numeric(ForestEcog - trueCoefs[2])
  biasVar1Forest <- as.numeric(ForestVar1 - trueCoefs[3])
  biasVar2Forest <- as.numeric(ForestVar2 - trueCoefs[4])

  VarTreatForest <- (pool$std.error[TreatRow])^2
  VarEcogForest <- (pool$std.error[EcogRow])^2
  VarVar1Forest <- (pool$std.error[Var1Row])^2
  VarVar2Forest <- (pool$std.error[Var2Row])^2

  CIForestLowerTreat <- Foresttreat - (qnorm(.975) * (pool$std.error[TreatRow]))
  CIForestUpperTreat <- Foresttreat + (qnorm(.975) * (pool$std.error[TreatRow]))

  CIForestLowerEcog <- ForestEcog - (qnorm(.975) * (pool$std.error[EcogRow]))
  CIForestUpperEcog <- ForestEcog + (qnorm(.975) * (pool$std.error[EcogRow]))

  CIForestLowerVar1 <- ForestVar1 - (qnorm(.975) * (pool$std.error[Var1Row]))
  CIForestUpperVar1 <- ForestVar1 + (qnorm(.975) * (pool$std.error[Var1Row]))

  CIForestLowerVar2 <- ForestVar2 - (qnorm(.975) * (pool$std.error[Var2Row]))
  CIForestUpperVar2 <- ForestVar2 + (qnorm(.975) * (pool$std.error[Var2Row]))

  covereageForestTreat <- 1 * between(c(trueCoefs[1]), CIForestLowerTreat, CIForestUpperTreat)
  covereageForestEcog <- 1 * between(c(trueCoefs[2]), CIForestLowerEcog, CIForestUpperEcog)
  covereageForestVar1 <- 1 * between(c(trueCoefs[3]), CIForestLowerVar1, CIForestUpperVar1)
  covereageForestVar2 <- 1 * between(c(trueCoefs[4]), CIForestLowerVar2, CIForestUpperVar2)

  print("Oracle")

  # compare complete dataset to OS1 estimates
  CompleteFit <- coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + var1 + var2, data = cData)

  TreatRow <- which(names(CompleteFit$coefficients) == "treat")
  EcogRow <- which(names(CompleteFit$coefficients) == "b.ecogvalue")
  Var1Row <- which(names(CompleteFit$coefficients) == "var1")
  Var2Row <- which(names(CompleteFit$coefficients) == "var2")

  Completetreat <- CompleteFit$coefficients[TreatRow]
  CompleteEcog <- CompleteFit$coefficients[EcogRow]
  CompleteVar1 <- CompleteFit$coefficients[Var1Row]
  CompleteVar2 <- CompleteFit$coefficients[Var2Row]

  CICompleteLowerTreat <- Completetreat - (qnorm(.975) * (summary(CompleteFit)$coefficients[TreatRow, 3]))
  CICompleteUpperTreat <- Completetreat + (qnorm(.975) * (summary(CompleteFit)$coefficients[TreatRow, 3]))

  CICompleteLowerEcog <- CompleteEcog - (qnorm(.975) * (summary(CompleteFit)$coefficients[EcogRow, 3]))
  CICompleteUpperEcog <- CompleteEcog + (qnorm(.975) * (summary(CompleteFit)$coefficients[EcogRow, 3]))

  CICompleteLowerVar1 <- CompleteVar1 - (qnorm(.975) * (summary(CompleteFit)$coefficients[Var1Row, 3]))
  CICompleteUpperVar1 <- CompleteVar1 + (qnorm(.975) * (summary(CompleteFit)$coefficients[Var1Row, 3]))

  CICompleteLowerVar2 <- CompleteVar2 - (qnorm(.975) * (summary(CompleteFit)$coefficients[Var2Row, 3]))
  CICompleteUpperVar2 <- CompleteVar2 + (qnorm(.975) * (summary(CompleteFit)$coefficients[Var2Row, 3]))

  covereageCompleteTreat <- 1 * between(c(trueCoefs[1]), CICompleteLowerTreat, CICompleteUpperTreat)
  covereageCompleteEcog <- 1 * between(c(trueCoefs[2]), CICompleteLowerEcog, CICompleteUpperEcog)
  covereageCompleteVar1 <- 1 * between(c(trueCoefs[3]), CICompleteLowerVar1, CICompleteUpperVar1)
  covereageCompleteVar2 <- 1 * between(c(trueCoefs[4]), CICompleteLowerVar2, CICompleteUpperVar2)

  biasTreatComplete <- as.numeric(Completetreat - trueCoefs[1])
  biasEcogComplete <- as.numeric(CompleteEcog - trueCoefs[2])
  biasVar1Complete <- as.numeric(CompleteVar1 - trueCoefs[3])
  biasVar2Complete <- as.numeric(CompleteVar2 - trueCoefs[4])

  VarTreatComplete <- (summary(CompleteFit)$coefficients[TreatRow, 3])^2
  VarEcogComplete <- (summary(CompleteFit)$coefficients[EcogRow, 3])^2
  VarVar1Complete <- (summary(CompleteFit)$coefficients[Var1Row, 3])^2
  VarVar2Complete <- (summary(CompleteFit)$coefficients[Var2Row, 3])^2

  print("Complete Case")

  # compare exclude missingness dataset to OS1 estimates
  mData <- read.csv(file.path(rootdir, "datasets/mDats", mech, propName, paste0("mData", i, ".csv")))

  ExcludeData <- mData[complete.cases(mData), ]

  ExcludeFit <- coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + var1 + var2, data = ExcludeData)

  TreatRow <- which(names(ExcludeFit$coefficients) == "treat")
  EcogRow <- which(names(ExcludeFit$coefficients) == "b.ecogvalue")
  Var1Row <- which(names(ExcludeFit$coefficients) == "var1")
  Var2Row <- which(names(ExcludeFit$coefficients) == "var2")

  Excludetreat <- ExcludeFit$coefficients[TreatRow]
  ExcludeEcog <- ExcludeFit$coefficients[EcogRow]
  ExcludeVar1 <- ExcludeFit$coefficients[Var1Row]
  ExcludeVar2 <- ExcludeFit$coefficients[Var2Row]

  biasTreatExclude <- as.numeric(Excludetreat - trueCoefs[1])
  biasEcogExclude <- as.numeric(ExcludeEcog - trueCoefs[2])
  biasVar1Exclude <- as.numeric(ExcludeVar1 - trueCoefs[3])
  biasVar2Exclude <- as.numeric(ExcludeVar2 - trueCoefs[4])

  VarTreatExclude <- (summary(ExcludeFit)$coefficients[TreatRow, 3])^2
  VarEcogExclude <- (summary(ExcludeFit)$coefficients[EcogRow, 3])^2
  VarVar1Exclude <- (summary(ExcludeFit)$coefficients[Var1Row, 3])^2
  VarVar2Exclude <- (summary(ExcludeFit)$coefficients[Var2Row, 3])^2

  CIExcludeLowerTreat <- Excludetreat - (qnorm(.975) * (summary(ExcludeFit)$coefficients[TreatRow, 3]))
  CIExcludeUpperTreat <- Excludetreat + (qnorm(.975) * (summary(ExcludeFit)$coefficients[TreatRow, 3]))

  CIExcludeLowerEcog <- ExcludeEcog - (qnorm(.975) * (summary(ExcludeFit)$coefficients[EcogRow, 3]))
  CIExcludeUpperEcog <- ExcludeEcog + (qnorm(.975) * (summary(ExcludeFit)$coefficients[EcogRow, 3]))

  CIExcludeLowerVar1 <- ExcludeVar1 - (qnorm(.975) * (summary(ExcludeFit)$coefficients[Var1Row, 3]))
  CIExcludeUpperVar1 <- ExcludeVar1 + (qnorm(.975) * (summary(ExcludeFit)$coefficients[Var1Row, 3]))

  CIExcludeLowerVar2 <- ExcludeVar2 - (qnorm(.975) * (summary(ExcludeFit)$coefficients[Var2Row, 3]))
  CIExcludeUpperVar2 <- ExcludeVar2 + (qnorm(.975) * (summary(ExcludeFit)$coefficients[Var2Row, 3]))

  covereageExcludeTreat <- 1 * between(c(trueCoefs[1]), CIExcludeLowerTreat, CIExcludeUpperTreat)
  covereageExcludeEcog <- 1 * between(c(trueCoefs[2]), CIExcludeLowerEcog, CIExcludeUpperEcog)
  covereageExcludeVar1 <- 1 * between(c(trueCoefs[3]), CIExcludeLowerVar1, CIExcludeUpperVar1)
  covereageExcludeVar2 <- 1 * between(c(trueCoefs[4]), CIExcludeLowerVar2, CIExcludeUpperVar2)
  
  mseOracleTreat <- VarTreatComplete + biasTreatComplete^2
  mseOracleEcog <- VarEcogComplete + biasEcogComplete^2
  mseOracleVar1 <- VarVar1Complete + biasVar1Complete^2
  mseOracleVar2 <- VarVar2Complete + biasVar2Complete^2
  mseCCTreat <- VarTreatExclude + biasTreatExclude^2
  mseCCEcog <- VarEcogExclude + biasEcogExclude^2
  mseCCVar1 <- VarVar1Exclude + biasVar1Exclude^2
  mseCCVar2 <- VarVar2Exclude + biasVar2Exclude^2

  mseMiceTreat <- VarTreatMICE + biasTreatMICE^2
  mseMiceEcog <- VarEcogMICE + biasEcogMICE^2
  mseMiceVar1 <- VarVar1MICE + biasVar1MICE^2
  mseMiceVar2 <- VarVar2MICE + biasVar2MICE^2
  mseForestTreat <- VarTreatForest + biasTreatForest^2
  mseForestEcog <- VarEcogForest + biasEcogForest^2
  mseForestVar1 <- VarVar1Forest + biasVar1Forest^2
  mseForestVar2 <- VarVar2Forest + biasVar2Forest^2

  names <- c("oracleTreat", "oracleEcog", "oracleVar1", "oracleVar2", "complete-caseTreat", "complete-caseEcog", "complete-caseVar1", "complete-caseVar2", "MiceTreat", "MiceEcog", "MiceVar1", "MiceVar2", "forestTreat", "forestEcog", "forestVar1", "forestVar2")

  bias <- c(biasTreatComplete, biasEcogComplete, biasVar1Complete, biasVar2Complete, biasTreatExclude, biasEcogExclude, biasVar1Exclude, biasVar2Exclude, biasTreatMICE, biasEcogMICE, biasVar1MICE, biasVar2MICE, biasTreatForest, biasEcogForest, biasVar1Forest, biasVar2Forest)

  se <- c(sqrt(VarTreatComplete), sqrt(VarEcogComplete), sqrt(VarVar1Complete), sqrt(VarVar2Complete), sqrt(VarTreatExclude), sqrt(VarEcogExclude), sqrt(VarVar1Exclude), sqrt(VarVar2Exclude), sqrt(VarTreatMICE), sqrt(VarEcogMICE), sqrt(VarVar1MICE), sqrt(VarVar2MICE), sqrt(VarTreatForest), sqrt(VarEcogForest), sqrt(VarVar1Forest), sqrt(VarVar2Forest))

  coverage <- c(covereageCompleteTreat, covereageCompleteEcog, covereageCompleteVar1, covereageCompleteVar2, covereageExcludeTreat, covereageExcludeEcog, covereageExcludeVar1, covereageExcludeVar2, covereageMICETreat, covereageMICEEcog, covereageMICEVar1, covereageMICEVar2, covereageForestTreat, covereageForestEcog, covereageForestVar1, covereageForestVar2)

  MSE <- c(mseOracleTreat, mseOracleEcog, mseOracleVar1, mseOracleVar2, mseCCTreat, mseCCEcog, mseCCVar1, mseCCVar2, mseMiceTreat, mseMiceEcog, mseMiceVar1, mseMiceVar2, mseForestTreat, mseForestEcog, mseForestVar1, mseForestVar2)

  results <- t(cbind(bias, se, coverage, MSE))
  colnames(results) <- names
  wdir <- file.path(rootdir, "final_results")
  system(paste0("mkdir ", file.path(wdir, mech)))
  system(paste0("mkdir ", file.path(wdir, mech, propName)))
  write.csv(results, file.path(rootdir, "final_results", mech, propName, paste0("result", i, ".csv")))
}
