rootdir <- "/project/flatiron_ucc/programs/kylie/RunMe"
library(mice)
library(CALIBERrfimpute)
library(survival)
library(dplyr)
library(tidyverse)

mechList <- c("MCAR", "MAR", "MNAR")
proportionList <- c(10, 30, 50)
set.seed(100)
for (w in 1:3) {
  mech <- mechList[w]
  for (k in 1:3) {
    propName <- proportionList[k]
    args <- commandArgs(trailingOnly = F)
    i <- args[6]
    mData <- read.csv(file.path(rootdir, "datasets/mDats", mech, propName, paste0("mData", i, ".csv")))
    cData <- read.csv(file.path(rootdir, "datasets/cDats", mech, propName, paste0("cData", i, ".csv")))
    cData <- cData[, -1]
    truth <- read.csv(file.path(rootdir, "datasets/trueEff", mech, propName, "trueEff_propMiss1.csv"))
    trueTreat <- truth[1]
    trueEcog <- truth[2]

    mCat = mData
    mCat$race = rep (NA, nrow (mData))
    mCat$race[which (mData$reth_black == 1)] = 'black'
    mCat$race[which (mData$reth_hisp == 1)] = 'hisp'
    mCat$race[which (mData$reth_oth == 1)] = 'oth'
    mCat$race[which (mData$reth_black == 0 & mData$reth_hisp == 0 & mData$reth_oth == 0)] = 'white'
    mCat$race = factor (mCat$race)

    mCat$site = rep (NA, nrow (mData))
    mCat$site[which (mData$site_ureter == 1)] = 'ureter'
    mCat$site[which (mData$site_renal == 1)] = 'renal'
    mCat$site[which (mData$site_urethra == 1)] = 'urethra'
    mCat$site[which (mData$site_ureter == 0 & mData$site_renal == 0 & mData$site_urethra == 0)] = 'bladder'
    mCat$site = factor (mCat$site)

    mData$reth_black = NULL
    mData$reth_hisp = NULL
    mData$reth_oth = NULL
    mData$race = mCat$race

    mData$site_ureter = NULL
    mData$site_renal = NULL
    mData$site_urethra = NULL
    mData$site = mCat$site

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

    reth_black = 1*(anesimp_long$race=='black')
    reth_hisp = 1*(anesimp_long$race=='hisp')
    reth_oth = 1*(anesimp_long$race=='oth')
    site_ureter = 1*(anesimp_long$site=='ureter')
    site_renal = 1*(anesimp_long$site=='renal')
    site_urethra = 1*(anesimp_long$site=='urethra')
    anesimp_long$race = NULL
    anesimp_long$site = NULL
    
    anesimp_long <- add_column(anesimp_long, reth_black = reth_black, .after = "genderf")
    anesimp_long <- add_column(anesimp_long, reth_hisp = reth_hisp, .after = "reth_black")
    anesimp_long <- add_column(anesimp_long, reth_oth = reth_oth, .after = "reth_hisp")

    anesimp_long <- add_column(anesimp_long, site_ureter = site_ureter, .after = "practypec")
    anesimp_long <- add_column(anesimp_long, site_renal = site_renal, .after = "site_ureter")
    anesimp_long <- add_column(anesimp_long, site_urethra = site_urethra, .after = "site_renal")

    anesimp_long_mids <- as.mids(anesimp_long)

    fitimp1 <- with(
      anesimp_long_mids,
      coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age)
    )

    pool <- summary(pool(fitimp1))

    TreatRow <- which(pool$term == "treat")
    EcogRow <- which(pool$term == "b.ecogvalue")

    MICEtreat <- pool$estimate[TreatRow]
    MICEEcog <- pool$estimate[EcogRow]

    biasTreatMICE <- as.numeric(MICEtreat - trueTreat)
    biasEcogMICE <- as.numeric(MICEEcog - trueEcog)

    VarTreatMICE <- (pool$std.error[TreatRow])^2
    VarEcogMICE <- (pool$std.error[EcogRow])^2

    CIMICELowerTreat <- MICEtreat - (qnorm(.975) * (pool$std.error[TreatRow]))
    CIMICEUpperTreat <- MICEtreat + (qnorm(.975) * (pool$std.error[TreatRow]))

    CIMICELowerEcog <- MICEEcog - (qnorm(.975) * (pool$std.error[EcogRow]))
    CIMICEUpperEcog <- MICEEcog + (qnorm(.975) * (pool$std.error[EcogRow]))

    covereageMICETreat <- 1*between(c (trueTreat), CIMICELowerTreat, CIMICEUpperTreat)
    covereageMICEEcog <- 1*between(c (trueEcog), CIMICELowerEcog, CIMICEUpperEcog)

    mData <- read.csv(file.path(rootdir, "datasets/mDats", mech, propName, paste0("mData", i, ".csv")))

    mCat = mData
    mCat$race = rep (NA, nrow (mData))
    mCat$race[which (mData$reth_black == 1)] = 'black'
    mCat$race[which (mData$reth_hisp == 1)] = 'hisp'
    mCat$race[which (mData$reth_oth == 1)] = 'oth'
    mCat$race[which (mData$reth_black == 0 & mData$reth_hisp == 0 & mData$reth_oth == 0)] = 'white'
    mCat$race = factor (mCat$race)

    mCat$site = rep (NA, nrow (mData))
    mCat$site[which (mData$site_ureter == 1)] = 'ureter'
    mCat$site[which (mData$site_renal == 1)] = 'renal'
    mCat$site[which (mData$site_urethra == 1)] = 'urethra'
    mCat$site[which (mData$site_ureter == 0 & mData$site_renal == 0 & mData$site_urethra == 0)] = 'bladder'
    mCat$site = factor (mCat$site)

    mData$reth_black = NULL
    mData$reth_hisp = NULL
    mData$reth_oth = NULL
    mData$race = mCat$race

    mData$site_ureter = NULL
    mData$site_renal = NULL
    mData$site_urethra = NULL
    mData$site = mCat$site
    
    mData$event <- as.numeric(mData$event)
    mData$genderf <- as.factor(mData$genderf)
    mData$practypec <- as.factor(mData$practypec)
    mData$smokey <- as.factor(mData$smokey)
    mData$dgradeh <- as.factor(mData$dgradeh)
    mData$surgery <- as.factor(mData$surgery)
    mData$treat <- as.factor(mData$treat)
    mData$b.ecogvalue <- as.numeric(mData$b.ecogvalue)
    mData$race <- as.factor(mData$race)
    mData$site <- as.factor(mData$site)

    forest <- mice(mData, method = c("rfcat", "rfcont", "rfcont", "rfcat", "rfcont", "rfcat", "rfcat", "rfcat", "rfcat", "rfcat", "rfcont", "rfcat", "rfcat"), m = 10, maxit = 5)

    anesimp_long <- mice::complete(forest, action = "long", include = TRUE)

    reth_black = 1*(anesimp_long$race=='black')
    reth_hisp = 1*(anesimp_long$race=='hisp')
    reth_oth = 1*(anesimp_long$race=='oth')
    site_ureter = 1*(anesimp_long$site=='ureter')
    site_renal = 1*(anesimp_long$site=='renal')
    site_urethra = 1*(anesimp_long$site=='urethra')
    anesimp_long$race = NULL
    anesimp_long$site = NULL
    
    anesimp_long <- add_column(anesimp_long, reth_black = reth_black, .after = "genderf")
    anesimp_long <- add_column(anesimp_long, reth_hisp = reth_hisp, .after = "reth_black")
    anesimp_long <- add_column(anesimp_long, reth_oth = reth_oth, .after = "reth_hisp")

    anesimp_long <- add_column(anesimp_long, site_ureter = site_ureter, .after = "practypec")
    anesimp_long <- add_column(anesimp_long, site_renal = site_renal, .after = "site_ureter")
    anesimp_long <- add_column(anesimp_long, site_urethra = site_urethra, .after = "site_renal")

    anesimp_long_mids <- as.mids(anesimp_long)
    fitimp1 <- with(
      anesimp_long_mids,
      coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age)
    )

    pool <- summary(pool(fitimp1))

    TreatRow <- which(pool$term == "treat1")
    EcogRow <- which(pool$term == "b.ecogvalue")

    forestTreat <- pool$estimate[TreatRow]
    forestEcog <- pool$estimate[EcogRow]

    biasTreatForest <- as.numeric(forestTreat - trueTreat)
    biasEcogForest <- as.numeric(forestEcog - trueEcog)

    VarTreatforest <- (pool$std.error[TreatRow])^2
    VarEcogforest <- (pool$std.error[EcogRow])^2

    CIForestLowerTreat <- forestTreat - (qnorm(.975) * (pool$std.error[TreatRow]))
    CIForestUpperTreat <- forestTreat + (qnorm(.975) * (pool$std.error[TreatRow]))

    CIForestLowerEcog <- forestEcog - (qnorm(.975) * (pool$std.error[EcogRow]))
    CIForestUpperEcog <- forestEcog + (qnorm(.975) * (pool$std.error[EcogRow]))

    covereageForestTreat <- 1*between(c (trueTreat), CIForestLowerTreat, CIForestUpperTreat)
    covereageForestEcog <- 1*between(c (trueEcog), CIForestLowerEcog, CIForestUpperEcog)

    # compare complete dataset to OS1 estimates
    CompleteFit <- coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age, data = cData)

    TreatRow <- which(names (CompleteFit$coefficients) == "treat")
    EcogRow <- which(names (CompleteFit$coefficients) == "b.ecogvalue")

    Completetreat <- CompleteFit$coefficients[TreatRow]
    CompleteEcog <- CompleteFit$coefficients[EcogRow]

    CICompleteLowerTreat <- Completetreat - (qnorm(.975) * (summary(CompleteFit)$coefficients[TreatRow, 3]))
    CICompleteUpperTreat <- Completetreat + (qnorm(.975) * (summary(CompleteFit)$coefficients[TreatRow, 3]))

    CICompleteLowerEcog <- CompleteEcog - (qnorm(.975) * (summary(CompleteFit)$coefficients[EcogRow, 3]))
    CICompleteUpperEcog <- CompleteEcog + (qnorm(.975) * (summary(CompleteFit)$coefficients[EcogRow, 3]))

    covereageCompleteTreat <- 1*between(c (trueTreat), CICompleteLowerTreat, CICompleteUpperTreat)
    covereageCompleteEcog <- 1*between(c (trueEcog), CICompleteLowerEcog, CICompleteUpperEcog)

    biasTreatComplete <- as.numeric(Completetreat - trueTreat)
    biasEcogComplete <- as.numeric(CompleteEcog - trueEcog)

    VarTreatComplete <- (summary(CompleteFit)$coefficients[TreatRow, 3])^2
    VarEcogComplete <- (summary(CompleteFit)$coefficients[EcogRow, 3])^2

    # compare exclude missingness dataset to OS1 estimates
    mData <- read.csv(file.path(rootdir, "datasets/mDats", mech, propName, paste0("mData", i, ".csv")))

    ExcludeData <- mData[complete.cases(mData), ]

    ExcludeFit <- coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age, data = ExcludeData)

    TreatRow <- which(names (ExcludeFit$coefficients) == "treat")
    EcogRow <- which(names (ExcludeFit$coefficients) == "b.ecogvalue")

    Excludetreat <- ExcludeFit$coefficients[TreatRow]
    ExcludeEcog <- ExcludeFit$coefficients[EcogRow]

    biasTreatExclude <- as.numeric(Excludetreat - trueTreat)
    biasEcogExclude <- as.numeric(ExcludeEcog - trueEcog)

    VarTreatExclude <- (summary(ExcludeFit)$coefficients[TreatRow, 3])^2
    VarEcogExclude <- (summary(ExcludeFit)$coefficients[EcogRow, 3])^2

    CIExcludeLowerTreat <- ExcludeFit$coefficients[TreatRow] - (qnorm(.975) * (summary(ExcludeFit)$coefficients[TreatRow, 3]))
    CIExcludeUpperTreat <- ExcludeFit$coefficients[TreatRow] + (qnorm(.975) * (summary(ExcludeFit)$coefficients[TreatRow, 3]))

    CIExcludeLowerEcog <- ExcludeFit$coefficients[EcogRow] - (qnorm(.975) * (summary(ExcludeFit)$coefficients[EcogRow, 3]))
    CIExcludeUpperEcog <- ExcludeFit$coefficients[EcogRow] + (qnorm(.975) * (summary(ExcludeFit)$coefficients[EcogRow, 3]))

    covereageExcludeTreat <- 1*between(c (trueTreat), CIExcludeLowerTreat, CIExcludeUpperTreat)
    covereageExcludeEcog <- 1*between(c (trueEcog), CIExcludeLowerEcog, CIExcludeUpperEcog)

    mseOracleTreat <- VarTreatComplete + biasTreatComplete^2
    mseOracleEcog <- VarEcogComplete + biasEcogComplete^2
    mseCCTreat <- VarTreatExclude + biasTreatExclude^2
    mseCCEcog <- VarEcogExclude + biasEcogExclude^2
    mseMiceTreat <- VarTreatMICE + biasTreatMICE^2
    mseMiceEcog <- VarEcogMICE + biasEcogMICE^2
    mseForestTreat <- VarTreatforest + biasTreatForest^2
    mseForestEcog <- VarEcogforest + biasEcogForest^2

    names <- c("oracleTreat", "oracleEcog", "complete-caseTreat", "complete-caseEcog", "MiceTreat", "MiceEcog", "forestTreat", "forestEcog")
    bias <- c(biasTreatComplete, biasEcogComplete, biasTreatExclude, biasEcogExclude, biasTreatMICE, biasEcogMICE, biasTreatForest, biasEcogForest)
    se <- c(sqrt(VarTreatComplete), sqrt(VarEcogComplete), sqrt(VarTreatExclude), sqrt(VarEcogExclude), sqrt(VarTreatMICE), sqrt(VarEcogMICE), sqrt(VarTreatforest), sqrt(VarEcogforest))
    coverage <- c(covereageCompleteTreat, covereageCompleteEcog, covereageExcludeTreat, covereageExcludeEcog, covereageMICETreat, covereageMICEEcog, covereageForestTreat, covereageForestEcog)
    MSE <- c(mseOracleTreat, mseOracleEcog, mseCCTreat, mseCCEcog, mseMiceTreat, mseMiceEcog, mseForestTreat, mseForestEcog)
    results <- t(cbind(bias, se, coverage, MSE))
    colnames(results) <- names
    write.csv(results, file.path(rootdir, "final_results", mech, propName, paste0("result", i, ".csv")))
  }
}