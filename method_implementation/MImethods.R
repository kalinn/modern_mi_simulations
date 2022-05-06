library(mice)
library(CALIBERrfimpute)
library(survival)
library(dplyr)
library(tidyverse)

rootdir <- "/project/flatiron_ucc/programs/kylie/RunMe2"
system (paste0 ('mkdir ', file.path (rootdir, 'final_results')))

mechList <- c("MCAR", "MAR", "MNAR1", "MNAR2")
proportionList <- c(10, 30, 50)
set.seed(100)
for (w in 1:length (mechList)) {
  mech <- mechList[w]
  print (mech)
  for (k in 1:length (proportionList)) {
    propName <- proportionList[k]
    print (propName)
    args <- commandArgs(trailingOnly = TRUE)
    i <- as.numeric (args[1])
    mData <- read.csv(file.path(rootdir, "datasets/mDats", mech, propName, paste0("mData", i, ".csv")))
    cData <- read.csv(file.path(rootdir, "datasets/cDats", mech, propName, paste0("cData", i, ".csv")))
    cData <- cData[, -1]
    truth <- read.csv(file.path(rootdir, "datasets/trueEff", mech, propName, "propMiss_trueEffs1.csv"))$x
    tnames = read.csv (file.path(rootdir, "datasets/names.csv"))
    trueTreat <- truth[which (tnames$x=='TREAT')]
    trueEcog <- truth[which (tnames$x=='b.ecogvalue')]
    trueNewVar <- truth[which (tnames$x=='newVar')]

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
      coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + newVar)
    )

    pool <- summary(pool(fitimp1))

    TreatRow <- which(pool$term == "treat")
    EcogRow <- which(pool$term == "b.ecogvalue")
    NewVarRow <- which(pool$term == "newVar")

    MICEtreat <- pool$estimate[TreatRow]
    MICEEcog <- pool$estimate[EcogRow]
    MICENewVar <- pool$estimate[NewVarRow]

    biasTreatMICE <- as.numeric(MICEtreat - trueTreat)
    biasEcogMICE <- as.numeric(MICEEcog - trueEcog)
    biasNewVarMICE <- as.numeric(MICENewVar - trueNewVar)

    VarTreatMICE <- (pool$std.error[TreatRow])^2
    VarEcogMICE <- (pool$std.error[EcogRow])^2
    VarNewVarMICE <- (pool$std.error[NewVarRow])^2

    CIMICELowerTreat <- MICEtreat - (qnorm(.975) * (pool$std.error[TreatRow]))
    CIMICEUpperTreat <- MICEtreat + (qnorm(.975) * (pool$std.error[TreatRow]))

    CIMICELowerEcog <- MICEEcog - (qnorm(.975) * (pool$std.error[EcogRow]))
    CIMICEUpperEcog <- MICEEcog + (qnorm(.975) * (pool$std.error[EcogRow]))

    CIMICELowerNewVar <- MICENewVar - (qnorm(.975) * (pool$std.error[NewVarRow]))
    CIMICEUpperNewVar <- MICENewVar + (qnorm(.975) * (pool$std.error[NewVarRow]))

    covereageMICETreat <- 1*between(c (trueTreat), CIMICELowerTreat, CIMICEUpperTreat)
    covereageMICEEcog <- 1*between(c (trueEcog), CIMICELowerEcog, CIMICEUpperEcog)
    covereageMICENewVar <- 1*between(c (trueNewVar), CIMICELowerNewVar, CIMICEUpperNewVar)

    print ("RF imputation")

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
    mData$newVar <- as.numeric(mData$newVar)
    mData$race <- as.factor(mData$race)
    mData$site <- as.factor(mData$site)

    forest <- mice(mData, method = c("rfcat", "rfcont", "rfcont", "rfcat", "rfcont", "rfcat", "rfcat", "rfcat", "rfcat", "rfcat", "rfcont", "rfcont", "rfcat", "rfcat"), m = 10, maxit = 5)

    anesimp_long <- mice::complete(forest, action = "long", include = TRUE)
    anesimp_long$b.ecogvalue = round (anesimp_long$b.ecogvalue)

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
      coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + newVar)
    )

    pool <- summary(pool(fitimp1))

    TreatRow <- which(pool$term == "treat1")
    EcogRow <- which(pool$term == "b.ecogvalue")
    NewVarRow <- which(pool$term == "newVar")

    forestTreat <- pool$estimate[TreatRow]
    forestEcog <- pool$estimate[EcogRow]
    forestNewVar <- pool$estimate[NewVarRow]

    biasTreatForest <- as.numeric(forestTreat - trueTreat)
    biasEcogForest <- as.numeric(forestEcog - trueEcog)
    biasNewVarForest <- as.numeric(forestNewVar - trueNewVar)

    VarTreatForest <- (pool$std.error[TreatRow])^2
    VarEcogForest <- (pool$std.error[EcogRow])^2
    VarNewVarForest <- (pool$std.error[NewVarRow])^2

    CIForestLowerTreat <- forestTreat - (qnorm(.975) * (pool$std.error[TreatRow]))
    CIForestUpperTreat <- forestTreat + (qnorm(.975) * (pool$std.error[TreatRow]))

    CIForestLowerEcog <- forestEcog - (qnorm(.975) * (pool$std.error[EcogRow]))
    CIForestUpperEcog <- forestEcog + (qnorm(.975) * (pool$std.error[EcogRow]))

    CIForestLowerNewVar <- forestNewVar - (qnorm(.975) * (pool$std.error[NewVarRow]))
    CIForestUpperNewVar <- forestNewVar + (qnorm(.975) * (pool$std.error[NewVarRow]))

    covereageForestTreat <- 1*between(c (trueTreat), CIForestLowerTreat, CIForestUpperTreat)
    covereageForestEcog <- 1*between(c (trueEcog), CIForestLowerEcog, CIForestUpperEcog)
    covereageForestNewVar <- 1*between(c (trueNewVar), CIForestLowerNewVar, CIForestUpperNewVar)

    print ("Oracle")

    # compare complete dataset to OS1 estimates
    CompleteFit <- coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + newVar, data = cData)

    TreatRow <- which(names (CompleteFit$coefficients) == "treat")
    EcogRow <- which(names (CompleteFit$coefficients) == "b.ecogvalue")
    NewVarRow <- which(names (CompleteFit$coefficients) == "newVar")

    Completetreat <- CompleteFit$coefficients[TreatRow]
    CompleteEcog <- CompleteFit$coefficients[EcogRow]
    CompleteNewVar <- CompleteFit$coefficients[NewVarRow]

    CICompleteLowerTreat <- Completetreat - (qnorm(.975) * (summary(CompleteFit)$coefficients[TreatRow, 3]))
    CICompleteUpperTreat <- Completetreat + (qnorm(.975) * (summary(CompleteFit)$coefficients[TreatRow, 3]))

    CICompleteLowerEcog <- CompleteEcog - (qnorm(.975) * (summary(CompleteFit)$coefficients[EcogRow, 3]))
    CICompleteUpperEcog <- CompleteEcog + (qnorm(.975) * (summary(CompleteFit)$coefficients[EcogRow, 3]))

    CICompleteLowerNewVar <- CompleteNewVar - (qnorm(.975) * (summary(CompleteFit)$coefficients[NewVarRow, 3]))
    CICompleteUpperNewVar <- CompleteNewVar + (qnorm(.975) * (summary(CompleteFit)$coefficients[NewVarRow, 3]))

    covereageCompleteTreat <- 1*between(c (trueTreat), CICompleteLowerTreat, CICompleteUpperTreat)
    covereageCompleteEcog <- 1*between(c (trueEcog), CICompleteLowerEcog, CICompleteUpperEcog)
    covereageCompleteNewVar <- 1*between(c (trueNewVar), CICompleteLowerNewVar, CICompleteUpperNewVar)

    biasTreatComplete <- as.numeric(Completetreat - trueTreat)
    biasEcogComplete <- as.numeric(CompleteEcog - trueEcog)
    biasNewVarComplete <- as.numeric(CompleteNewVar - trueNewVar)

    VarTreatComplete <- (summary(CompleteFit)$coefficients[TreatRow, 3])^2
    VarEcogComplete <- (summary(CompleteFit)$coefficients[EcogRow, 3])^2
    VarNewVarComplete <- (summary(CompleteFit)$coefficients[NewVarRow, 3])^2

    print ("Complete Case")

    # compare exclude missingness dataset to OS1 estimates
    mData <- read.csv(file.path(rootdir, "datasets/mDats", mech, propName, paste0("mData", i, ".csv")))

    ExcludeData <- mData[complete.cases(mData), ]

    ExcludeFit <- coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + newVar, data = ExcludeData)

    TreatRow <- which(names (ExcludeFit$coefficients) == "treat")
    EcogRow <- which(names (ExcludeFit$coefficients) == "b.ecogvalue")
    NewVarRow <- which(names (ExcludeFit$coefficients) == "newVar")

    Excludetreat <- ExcludeFit$coefficients[TreatRow]
    ExcludeEcog <- ExcludeFit$coefficients[EcogRow]
    ExcludeNewVar <- ExcludeFit$coefficients[NewVarRow]

    biasTreatExclude <- as.numeric(Excludetreat - trueTreat)
    biasEcogExclude <- as.numeric(ExcludeEcog - trueEcog)
    biasNewVarExclude <- as.numeric(ExcludeNewVar - trueNewVar)

    VarTreatExclude <- (summary(ExcludeFit)$coefficients[TreatRow, 3])^2
    VarEcogExclude <- (summary(ExcludeFit)$coefficients[EcogRow, 3])^2
    VarNewVarExclude <- (summary(ExcludeFit)$coefficients[NewVarRow, 3])^2

    CIExcludeLowerTreat <- ExcludeFit$coefficients[TreatRow] - (qnorm(.975) * (summary(ExcludeFit)$coefficients[TreatRow, 3]))
    CIExcludeUpperTreat <- ExcludeFit$coefficients[TreatRow] + (qnorm(.975) * (summary(ExcludeFit)$coefficients[TreatRow, 3]))

    CIExcludeLowerEcog <- ExcludeFit$coefficients[EcogRow] - (qnorm(.975) * (summary(ExcludeFit)$coefficients[EcogRow, 3]))
    CIExcludeUpperEcog <- ExcludeFit$coefficients[EcogRow] + (qnorm(.975) * (summary(ExcludeFit)$coefficients[EcogRow, 3]))
    
    CIExcludeLowerNewVar <- ExcludeFit$coefficients[NewVarRow] - (qnorm(.975) * (summary(ExcludeFit)$coefficients[NewVarRow, 3]))
    CIExcludeUpperNewVar <- ExcludeFit$coefficients[NewVarRow] + (qnorm(.975) * (summary(ExcludeFit)$coefficients[NewVarRow, 3]))

    covereageExcludeTreat <- 1*between(c (trueTreat), CIExcludeLowerTreat, CIExcludeUpperTreat)
    covereageExcludeEcog <- 1*between(c (trueEcog), CIExcludeLowerEcog, CIExcludeUpperEcog)
    covereageExcludeNewVar <- 1*between(c (trueNewVar), CIExcludeLowerNewVar, CIExcludeUpperNewVar)

    mseOracleTreat <- VarTreatComplete + biasTreatComplete^2
    mseOracleEcog <- VarEcogComplete + biasEcogComplete^2
    mseOracleNewVar <- VarNewVarComplete + biasNewVarComplete^2
    mseCCTreat <- VarTreatExclude + biasTreatExclude^2
    mseCCEcog <- VarEcogExclude + biasEcogExclude^2
    mseCCNewVar <- VarNewVarExclude + biasNewVarExclude^2

    mseMiceTreat <- VarTreatMICE + biasTreatMICE^2
    mseMiceEcog <- VarEcogMICE + biasEcogMICE^2
    mseMiceNewVar <- VarNewVarMICE + biasNewVarMICE^2
    mseForestTreat <- VarTreatForest + biasTreatForest^2
    mseForestEcog <- VarEcogForest + biasEcogForest^2
    mseForestNewVar <- VarNewVarForest + biasNewVarForest^2

    names <- c("oracleTreat", "oracleEcog", "oracleNewVar", "complete-caseTreat", "complete-caseEcog", "complete-caseNewVar", "MiceTreat", "MiceEcog", "MiceNewVar", "forestTreat", "forestEcog", "forestNewVar")
    bias <- c(biasTreatComplete, biasEcogComplete, biasNewVarComplete, biasTreatExclude, biasEcogExclude, biasNewVarExclude, biasTreatMICE, biasEcogMICE, biasNewVarMICE, biasTreatForest, biasEcogForest, biasNewVarForest)
    se <- c(sqrt(VarTreatComplete), sqrt(VarEcogComplete), sqrt(VarNewVarComplete), sqrt(VarTreatExclude), sqrt(VarEcogExclude), sqrt(VarNewVarExclude), sqrt(VarTreatMICE), sqrt(VarEcogMICE), sqrt(VarNewVarMICE), sqrt(VarTreatForest), sqrt(VarEcogForest), sqrt(VarNewVarForest))
    coverage <- c(covereageCompleteTreat, covereageCompleteEcog, covereageCompleteNewVar, covereageExcludeTreat, covereageExcludeEcog, covereageExcludeNewVar, covereageMICETreat, covereageMICEEcog, covereageMICENewVar, covereageForestTreat, covereageForestEcog, covereageForestNewVar)
    MSE <- c(mseOracleTreat, mseOracleEcog, mseOracleNewVar, mseCCTreat, mseCCEcog, mseCCNewVar, mseMiceTreat, mseMiceEcog, mseMiceNewVar, mseForestTreat, mseForestEcog, mseForestNewVar)
    results <- t(cbind(bias, se, coverage, MSE))
    colnames(results) <- names
    wdir = file.path (rootdir, 'final_results')
    system (paste0 ('mkdir ', file.path (wdir, mech)))
    system (paste0 ('mkdir ', file.path (wdir, mech, propName)))
    write.csv(results, file.path(rootdir, "final_results", mech, propName, paste0("result", i, ".csv")))
  }
}