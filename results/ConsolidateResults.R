library(dplyr)
library(survival)

rootdir <- "/project/flatiron_ucc/programs/kylie/RunMe/final_results"
filesMCAR10 <- list.files(file.path(rootdir, "MCAR/10"), full.names = TRUE)
filesMCAR30 <- list.files(file.path(rootdir, "MCAR/30"), full.names = TRUE)
filesMCAR50 <- list.files(file.path(rootdir, "MCAR/50"), full.names = TRUE)
filesMAR10 <- list.files(file.path(rootdir, "MAR/10"), full.names = TRUE)
filesMAR30 <- list.files(file.path(rootdir, "MAR/30"), full.names = TRUE)
filesMAR50 <- list.files(file.path(rootdir, "MAR/50"), full.names = TRUE)
filesMNAR10 <- list.files(file.path(rootdir, "MNAR/10"), full.names = TRUE)
filesMNAR30 <- list.files(file.path(rootdir, "MNAR/30"), full.names = TRUE)
filesMNAR50 <- list.files(file.path(rootdir, "MNAR/50"), full.names = TRUE)

makeTable <- function(data) {
  df <- lapply(data, read.csv)
  df <- lapply(df, function(x) x <- x[, -1])
  dfAvg <- Reduce("+", df) / length(df)
  rownames(dfAvg) <- c("bias", "se", "coverage", "MSE")
  biasList <- (lapply(df, function(x) x <- x[1, ]))
  biasTable <- bind_rows(biasList, .id = "column_label")
  biasTable <- biasTable[, -1]
  PercentileOracleTreat <- data.frame(quantile(biasTable[, 1], c(0.025, 0.975)))
  PercentileOracleEcog <- data.frame(quantile(biasTable[, 2], c(0.025, 0.975)))
  PercentileCCTreat <- data.frame(quantile(biasTable[, 3], c(0.025, 0.975)))
  PercentileCCEcog <- data.frame(quantile(biasTable[, 4], c(0.025, 0.975)))
  PercentileMICETreat <- data.frame(quantile(biasTable[, 5], c(0.025, 0.975)))
  PercentileMICEEcog <- data.frame(quantile(biasTable[, 6], c(0.025, 0.975)))
  PercentileForestTreat <- data.frame(quantile(biasTable[, 7], c(0.025, 0.975)))
  PercentileForestEcog <- data.frame(quantile(biasTable[, 8], c(0.025, 0.975)))
  colnames(PercentileOracleTreat) <- "Percentile"
  colnames(PercentileOracleEcog) <- "Percentile"
  colnames(PercentileCCTreat) <- "Percentile"
  colnames(PercentileCCEcog) <- "Percentile"
  colnames(PercentileMICETreat) <- "Percentile"
  colnames(PercentileMICEEcog) <- "Percentile"
  colnames(PercentileForestTreat) <- "Percentile"
  colnames(PercentileForestEcog) <- "Percentile"
  PercentileTreatDf <- rbind(PercentileOracleTreat, PercentileCCTreat, PercentileMICETreat, PercentileForestTreat)
  rownames(PercentileTreatDf) <- c("OracleLower", "OracleUpper", "CCLower", "CCUpper", "MICELower", "MICEUpper", "ForestLower", "ForestUpper")
  PercentileEcogDf <- rbind(PercentileOracleEcog, PercentileCCEcog, PercentileMICEEcog, PercentileForestEcog)
  rownames(PercentileEcogDf) <- c("OracleLower", "OracleUpper", "CCLower", "CCUpper", "MICELower", "MICEUpper", "ForestLower", "ForestUpper")
  return(list("treat" = PercentileTreatDf, "ecog" = PercentileEcogDf, "avg" = dfAvg))
}

MCAR10 <- makeTable(filesMCAR10)
PercentileTreatMCAR10 <- MCAR10$treat
PercentileEcogMCAR10 <- MCAR10$ecog

MCAR30 <- makeTable(filesMCAR30)
PercentileTreatMCAR30 <- MCAR30$treat
PercentileEcogMCAR30 <- MCAR30$ecog

MCAR50 <- makeTable(filesMCAR50)
PercentileTreatMCAR50 <- MCAR50$treat
PercentileEcogMCAR50 <- MCAR50$ecog

MAR10 <- makeTable(filesMAR10)
PercentileTreatMAR10 <- MAR10$treat
PercentileEcogMAR10 <- MAR10$ecog

MAR30 <- makeTable(filesMAR30)
PercentileTreatMAR30 <- MAR30$treat
PercentileEcogMAR30 <- MAR30$ecog

MAR50 <- makeTable(filesMAR50)
PercentileTreatMAR50 <- MAR50$treat
PercentileEcogMAR50 <- MAR50$ecog

MNAR10 <- makeTable(filesMNAR10)
PercentileTreatMNAR10 <- MNAR10$treat
PercentileEcogMNAR10 <- MNAR10$ecog

MNAR30 <- makeTable(filesMNAR30)
PercentileTreatMNAR30 <- MNAR30$treat
PercentileEcogMNAR30 <- MNAR30$ecog

MNAR50 <- makeTable(filesMNAR50)
PercentileTreatMNAR50 <- MNAR50$treat
PercentileEcogMNAR50 <- MNAR50$ecog


PercentileTreatMCAR <- cbind(PercentileTreatMCAR10, PercentileTreatMCAR30, PercentileTreatMCAR50)
colnames(PercentileTreatMCAR) <- c("10", "30", "50")
PercentileEcogMCAR <- cbind(PercentileEcogMCAR10, PercentileEcogMCAR30, PercentileEcogMCAR50)
colnames(PercentileEcogMCAR) <- c("10", "30", "50")
PercentileTreatMAR <- cbind(PercentileTreatMAR10, PercentileTreatMAR30, PercentileTreatMAR50)
colnames(PercentileTreatMAR) <- c("10", "30", "50")
PercentileEcogMAR <- cbind(PercentileEcogMAR10, PercentileEcogMAR30, PercentileEcogMAR50)
colnames(PercentileEcogMAR) <- c("10", "30", "50")
PercentileTreatMNAR <- cbind(PercentileTreatMNAR10, PercentileTreatMNAR30, PercentileTreatMNAR50)
colnames(PercentileTreatMNAR) <- c("10", "30", "50")
PercentileEcogMNAR <- cbind(PercentileEcogMNAR10, PercentileEcogMNAR30, PercentileEcogMNAR50)
colnames(PercentileEcogMNAR) <- c("10", "30", "50")

# fixed percentile treat should be for MCAR, MAR, and MNAR
write.csv(PercentileTreatMCAR, file.path(rootdir, "PercentileTreatMCAR10TreeV2.csv"))
write.csv(PercentileEcogMCAR, file.path(rootdir, "PercentileEcogMCAR10TreeV2.csv"))
write.csv(PercentileTreatMAR, file.path(rootdir, "PercentileTreatMAR10TreeV2.csv"))
write.csv(PercentileEcogMAR, file.path(rootdir, "PercentileEcogMAR10TreeV2.csv"))
write.csv(PercentileTreatMNAR, file.path(rootdir, "PercentileTreatMNAR10TreeV2.csv"))
write.csv(PercentileEcogMNAR, file.path(rootdir, "PercentileEcogMNAR10TreeV2.csv"))

MCAR10Avg <- MCAR10$avg
MCAR30Avg <- MCAR30$avg
MCAR50Avg <- MCAR50$avg
MAR10Avg <- MAR10$avg
MAR30Avg <- MAR30$avg
MAR50Avg <- MAR50$avg
MNAR10Avg <- MNAR10$avg
MNAR30Avg <- MNAR30$avg
MNAR50Avg <- MNAR50$avg

write.csv(MCAR10Avg, file.path(rootdir, "MCAR10Avg10TreeV2.csv"))
write.csv(MCAR30Avg, file.path(rootdir, "MCAR30Avg10TreeV2.csv"))
write.csv(MCAR50Avg, file.path(rootdir, "MCAR50Avg10TreeV2.csv"))
write.csv(MAR10Avg, file.path(rootdir, "MAR10Avg10TreeV2.csv"))
write.csv(MAR30Avg, file.path(rootdir, "MAR30Avg10TreeV2.csv"))
write.csv(MAR50Avg, file.path(rootdir, "MAR50Avg10TreeV2.csv"))
write.csv(MNAR10Avg, file.path(rootdir, "MNAR10Avg10TreeV2.csv"))
write.csv(MNAR30Avg, file.path(rootdir, "MNAR30Avg10TreeV2.csv"))
write.csv(MNAR50Avg, file.path(rootdir, "MNAR50Avg10TreeV2.csv"))

# Autoencoder
filesMCAR10 <- list.files(file.path(rootdir, "AEMCAR/10"), full.names = TRUE)
filesMCAR30 <- list.files(file.path(rootdir, "AEMCAR/30"), full.names = TRUE)
filesMCAR50 <- list.files(file.path(rootdir, "AEMCAR/50"), full.names = TRUE)
filesMAR10 <- list.files(file.path(rootdir, "AEMAR/10"), full.names = TRUE)
filesMAR30 <- list.files(file.path(rootdir, "AEMAR/30"), full.names = TRUE)
filesMAR50 <- list.files(file.path(rootdir, "AEMAR/50"), full.names = TRUE)
filesMNAR10 <- list.files(file.path(rootdir, "AEMNAR/10"), full.names = TRUE)
filesMNAR30 <- list.files(file.path(rootdir, "AEMNAR/30"), full.names = TRUE)
filesMNAR50 <- list.files(file.path(rootdir, "AEMNAR/50"), full.names = TRUE)

iter <- 1000
mechList <- c("MCAR", "MAR", "MNAR")
proportionList <- c(10, 30, 50)
# J is number of imputations
J <- 10

for (w in 1:3) {
  for (y in 1:3) {
    mech <- mechList[w]
    propName <- proportionList[y]
    truth <- read.csv(file.path("/project/flatiron_ucc/programs/kylie/RunMe/datasets/trueEff", mech, propName, "trueEff_propMiss1.csv"))
    trueTreat <- truth[1]
    trueEcog <- truth[2]
    fileListName <- paste0("files", mech, propName)
    fileCur <- get(fileListName)

    biasTreatAE <- c()
    biasEcogAE <- c()
    VarTreatAE <- c()
    VarEcogAE <- c()
    covereageAETreat <- c()
    covereageAEEcog <- c()
    mseAETreat <- c()
    mseAEEcog <- c()

    for (i in 1:iter) {
      print (mech)
      print (propName)
      print(i)
      daeTreatList <- c()
      daeEcogList <- c()
      VarTreatdaeList <- c()
      VarEcogdaeList <- c()
      for (j in 1:J) {
        print (j)
        data <- read.csv(file.path("/project/flatiron_ucc/programs/kylie/RunMe/final_results", paste0("AE", mech), propName, paste0("result", i, "_num_", j, ".csv")), row.names = 1) # read in instead

        daeFit <- coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age, data = data)
        treatRow <- which(names(daeFit$coefficients) == "treat")
        ecogRow <- which(names(daeFit$coefficients) == "b.ecogvalue")
        daeTreatList[j] <- daeFit$coefficients[treatRow]
        daeEcogList[j] <- daeFit$coefficients[ecogRow]
        VarTreatdaeList[j] <- (summary(daeFit)$coefficients[treatRow, 3])^2
        VarEcogdaeList[j] <- (summary(daeFit)$coefficients[ecogRow, 3])^2
      }

      # pool them
      daetreat <- mean(daeTreatList)
      daeEcog <- mean(daeEcogList)

      biasTreatAE[i] <- as.numeric(daetreat - trueTreat)
      biasEcogAE[i] <- as.numeric(daeEcog - trueEcog)

      UbarTreat <- mean(VarTreatdaeList)
      UbarEcog <- mean(VarEcogdaeList)

      BTreat <- var(daeTreatList)
      BEcog <- var(daeEcogList)

      # Variance formula for MI
      VarTreatAE[i] <- UbarTreat + (1 + (1 / J)) * BTreat
      VarEcogAE[i] <- UbarEcog + (1 + (1 / J)) * BEcog

      CIAELowerTreat <- daetreat - (qnorm(.975) * sqrt(VarTreatAE[i]))
      CIAEUpperTreat <- daetreat + (qnorm(.975) * sqrt(VarTreatAE[i]))
      CIAELowerEcog <- daeEcog - (qnorm(.975) * sqrt(VarEcogAE[i]))
      CIAEUpperEcog <- daeEcog + (qnorm(.975) * sqrt(VarEcogAE[i]))

      covereageAETreat[i] <- 1 * between(c(trueTreat), CIAELowerTreat, CIAEUpperTreat)
      covereageAEEcog[i] <- 1 * between(c(trueEcog), CIAELowerEcog, CIAEUpperEcog)

      mseAETreat[i] <- VarTreatAE[i] + biasTreatAE[i]^2
      mseAEEcog[i] <- VarEcogAE[i] + biasEcogAE[i]^2
    }

    PercentileTreatAE <- quantile(biasTreatAE, c(0.025, 0.975))
    PercentileEcogAE <- quantile(biasEcogAE, c(0.025, 0.975))
    biasTreatAvg <- mean(biasTreatAE)
    biasEcogAvg <- mean(biasEcogAE)
    VarTreatAvg <- mean(VarTreatAE)
    VarEcogAvg <- mean(VarEcogAE)
    coverForestTreat <- mean(covereageAETreat)
    coverForestEcog <- mean(covereageAEEcog)
    mseAETreatAvg <- mean(mseAETreat)
    mseAEEcogAvg <- mean(mseAEEcog)
    results <- rbind(cbind(biasTreatAvg, biasEcogAvg), cbind(sqrt(VarTreatAvg), sqrt(VarEcogAvg)), cbind(coverForestTreat, coverForestEcog), cbind(mseAETreatAvg, mseAEEcogAvg))
    rownames(results) <- c("Bias", "se", "coverage", "MSE")
    colnames(results) <- c("AETreat", "AEEcog")
    percentileTreat <- PercentileTreatAE
    percentileEcog <- PercentileEcogAE
    percentile <- cbind(percentileTreat, percentileEcog)
    write.csv(results, file.path(rootdir, paste0("AEresult10Tree", mech, propName, ".csv")))
    write.csv(percentile, file.path(rootdir, paste0("AEpercentile10Tree", mech, propName, ".csv")))
  }
}