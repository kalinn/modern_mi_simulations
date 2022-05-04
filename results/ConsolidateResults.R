library(dplyr)
library(survival)

rootdir <- "/project/flatiron_ucc/programs/kylie/RunMe2/final_results"
iter = 1000

filesMCAR10 <- list.files(file.path(rootdir, "MCAR/10"), full.names = TRUE)
filesMCAR30 <- list.files(file.path(rootdir, "MCAR/30"), full.names = TRUE)
filesMCAR50 <- list.files(file.path(rootdir, "MCAR/50"), full.names = TRUE)
filesMAR10 <- list.files(file.path(rootdir, "MAR/10"), full.names = TRUE)
filesMAR30 <- list.files(file.path(rootdir, "MAR/30"), full.names = TRUE)
filesMAR50 <- list.files(file.path(rootdir, "MAR/50"), full.names = TRUE)
filesMNAR10a <- list.files(file.path(rootdir, "MNAR1/10"), full.names = TRUE)
filesMNAR30a <- list.files(file.path(rootdir, "MNAR1/30"), full.names = TRUE)
filesMNAR50a <- list.files(file.path(rootdir, "MNAR1/50"), full.names = TRUE)
filesMNAR10b <- list.files(file.path(rootdir, "MNAR2/10"), full.names = TRUE)
filesMNAR30b <- list.files(file.path(rootdir, "MNAR2/30"), full.names = TRUE)
filesMNAR50b <- list.files(file.path(rootdir, "MNAR2/50"), full.names = TRUE)

makeTable <- function(data) {
  df <- lapply(data, read.csv)
  df <- lapply(df, function(x) x <- x[, -1])
  dfAvg <- Reduce("+", df) / length(df)
  rownames(dfAvg) <- c("bias", "se", "coverage", "MSE")
  biasList <- (lapply(df, function(x) x <- x[1, ]))
  biasTable <- bind_rows(biasList, .id = "column_label")
  biasTable <- biasTable[, -1]
#  PercentileOracleTreat <- data.frame(quantile(biasTable[, 1], c(0.025, 0.975)))
#  PercentileOracleEcog <- data.frame(quantile(biasTable[, 2], c(0.025, 0.975)))
#  PercentileOracleNewVar <- data.frame(quantile(biasTable[, 3], c(0.025, 0.975)))
#  PercentileCCTreat <- data.frame(quantile(biasTable[, 4], c(0.025, 0.975)))
#  PercentileCCEcog <- data.frame(quantile(biasTable[, 5], c(0.025, 0.975)))
#  PercentileCCNewVar <- data.frame(quantile(biasTable[, 6], c(0.025, 0.975)))
#  PercentileMICETreat <- data.frame(quantile(biasTable[, 7], c(0.025, 0.975)))
#  PercentileMICEEcog <- data.frame(quantile(biasTable[, 8], c(0.025, 0.975)))
#  PercentileMICENewVar <- data.frame(quantile(biasTable[, 9], c(0.025, 0.975)))
#  PercentileForestTreat <- data.frame(quantile(biasTable[, 10], c(0.025, 0.975)))
#  PercentileForestEcog <- data.frame(quantile(biasTable[, 11], c(0.025, 0.975)))
#  PercentileForestNewVar <- data.frame(quantile(biasTable[, 12], c(0.025, 0.975)))

  # Asymptotic CIs:
  PercentileOracleTreat <- data.frame(c (mean(biasTable[, 1]) - 1.96*sd (biasTable[, 1])/sqrt (iter), mean(biasTable[, 1]) + 1.96*sd (biasTable[, 1])/sqrt (iter)))
  PercentileOracleEcog <- data.frame(c (mean(biasTable[, 2]) - 1.96*sd (biasTable[, 2])/sqrt (iter), mean(biasTable[, 2]) + 1.96*sd (biasTable[, 2])/sqrt (iter)))
  PercentileOracleNewVar <- data.frame(c (mean(biasTable[, 3]) - 1.96*sd (biasTable[, 3])/sqrt (iter), mean(biasTable[, 3]) + 1.96*sd (biasTable[, 3])/sqrt (iter)))
  PercentileCCTreat <- data.frame(c (mean(biasTable[, 4]) - 1.96*sd (biasTable[, 4])/sqrt (iter), mean(biasTable[, 4]) + 1.96*sd (biasTable[, 4])/sqrt (iter)))
  PercentileCCEcog <- data.frame(c(mean(biasTable[, 5]) - 1.96*sd (biasTable[, 5])/sqrt (iter), mean(biasTable[, 5]) + 1.96*sd (biasTable[, 5])/sqrt (iter)))
  PercentileCCNewVar <- data.frame(c(mean(biasTable[, 6]) - 1.96*sd (biasTable[, 6])/sqrt (iter), mean(biasTable[, 6]) + 1.96*sd (biasTable[, 6])/sqrt (iter)))
  PercentileMICETreat <- data.frame(c(mean(biasTable[, 7]) - 1.96*sd (biasTable[, 7])/sqrt (iter), mean(biasTable[, 7]) + 1.96*sd (biasTable[, 7])/sqrt (iter)))
  PercentileMICEEcog <- data.frame(c(mean(biasTable[, 8]) - 1.96*sd (biasTable[, 8])/sqrt (iter), mean(biasTable[, 8]) + 1.96*sd (biasTable[, 8])/sqrt (iter)))
  PercentileMICENewVar <- data.frame(c(mean(biasTable[, 9]) - 1.96*sd (biasTable[, 9])/sqrt (iter), mean(biasTable[, 9]) + 1.96*sd (biasTable[, 9])/sqrt (iter)))
  PercentileForestTreat <- data.frame(c(mean(biasTable[, 10]) - 1.96*sd (biasTable[, 10])/sqrt (iter), mean(biasTable[, 10]) + 1.96*sd (biasTable[, 10])/sqrt (iter)))
  PercentileForestEcog <- data.frame(c(mean(biasTable[, 11]) - 1.96*sd (biasTable[, 11])/sqrt (iter), mean(biasTable[, 11]) + 1.96*sd (biasTable[, 11])/sqrt (iter)))
  PercentileForestNewVar <- data.frame(c(mean(biasTable[, 12]) - 1.96*sd (biasTable[, 12])/sqrt (iter), mean(biasTable[, 12]) + 1.96*sd (biasTable[, 12])/sqrt (iter)))
  colnames(PercentileOracleTreat) <- ""
  colnames(PercentileOracleEcog) <- ""
  colnames(PercentileOracleNewVar) <- ""
  colnames(PercentileCCTreat) <- ""
  colnames(PercentileCCEcog) <- ""
  colnames(PercentileCCNewVar) <- ""
  colnames(PercentileMICETreat) <- ""
  colnames(PercentileMICEEcog) <- ""
  colnames(PercentileMICENewVar) <- ""
  colnames(PercentileForestTreat) <- ""
  colnames(PercentileForestEcog) <- ""
  colnames(PercentileForestNewVar) <- ""

  PercentileTreatDf <- rbind(PercentileOracleTreat, PercentileCCTreat, PercentileMICETreat, PercentileForestTreat)
  rownames(PercentileTreatDf) <- c("OracleLower", "OracleUpper", "CCLower", "CCUpper", "MICELower", "MICEUpper", "ForestLower", "ForestUpper")
  PercentileEcogDf <- rbind(PercentileOracleEcog, PercentileCCEcog, PercentileMICEEcog, PercentileForestEcog)
  rownames(PercentileEcogDf) <- c("OracleLower", "OracleUpper", "CCLower", "CCUpper", "MICELower", "MICEUpper", "ForestLower", "ForestUpper")
  PercentileNewVarDf <- rbind(PercentileOracleNewVar, PercentileCCNewVar, PercentileMICENewVar, PercentileForestNewVar)
  rownames(PercentileNewVarDf) <- c("OracleLower", "OracleUpper", "CCLower", "CCUpper", "MICELower", "MICEUpper", "ForestLower", "ForestUpper")
  return(list("treat" = PercentileTreatDf, "ecog" = PercentileEcogDf, "newVar" = PercentileNewVarDf, "avg" = dfAvg))
}

MCAR10 <- makeTable(filesMCAR10)
PercentileTreatMCAR10 <- MCAR10$treat
PercentileEcogMCAR10 <- MCAR10$ecog
PercentileNewVarMCAR10 <- MCAR10$newVar

MCAR30 <- makeTable(filesMCAR30)
PercentileTreatMCAR30 <- MCAR30$treat
PercentileEcogMCAR30 <- MCAR30$ecog
PercentileNewVarMCAR30 <- MCAR30$newVar

MCAR50 <- makeTable(filesMCAR50)
PercentileTreatMCAR50 <- MCAR50$treat
PercentileEcogMCAR50 <- MCAR50$ecog
PercentileNewVarMCAR50 <- MCAR50$newVar

MAR10 <- makeTable(filesMAR10)
PercentileTreatMAR10 <- MAR10$treat
PercentileEcogMAR10 <- MAR10$ecog
PercentileNewVarMAR10 <- MAR10$newVar

MAR30 <- makeTable(filesMAR30)
PercentileTreatMAR30 <- MAR30$treat
PercentileEcogMAR30 <- MAR30$ecog
PercentileNewVarMAR30 <- MAR30$newVar

MAR50 <- makeTable(filesMAR50)
PercentileTreatMAR50 <- MAR50$treat
PercentileEcogMAR50 <- MAR50$ecog
PercentileNewVarMAR50 <- MAR50$newVar

MNAR10a <- makeTable(filesMNAR10a)
PercentileTreatMNAR10a <- MNAR10a$treat
PercentileEcogMNAR10a <- MNAR10a$ecog
PercentileNewVarMNAR10a <- MNAR10a$newVar

MNAR30a <- makeTable(filesMNAR30a)
PercentileTreatMNAR30a <- MNAR30a$treat
PercentileEcogMNAR30a <- MNAR30a$ecog
PercentileNewVarMNAR30a <- MNAR30a$newVar

MNAR50a <- makeTable(filesMNAR50a)
PercentileTreatMNAR50a <- MNAR50a$treat
PercentileEcogMNAR50a <- MNAR50a$ecog
PercentileNewVarMNAR50a <- MNAR50a$newVar

MNAR10b <- makeTable(filesMNAR10b)
PercentileTreatMNAR10b <- MNAR10b$treat
PercentileEcogMNAR10b <- MNAR10b$ecog
PercentileNewVarMNAR10b <- MNAR10b$newVar

MNAR30b <- makeTable(filesMNAR30b)
PercentileTreatMNAR30b <- MNAR30b$treat
PercentileEcogMNAR30b <- MNAR30b$ecog
PercentileNewVarMNAR30b <- MNAR30b$newVar

MNAR50b <- makeTable(filesMNAR50b)
PercentileTreatMNAR50b <- MNAR50b$treat
PercentileEcogMNAR50b <- MNAR50b$ecog
PercentileNewVarMNAR50b <- MNAR50b$newVar


PercentileTreatMCAR <- cbind(PercentileTreatMCAR10, PercentileTreatMCAR30, PercentileTreatMCAR50)
colnames(PercentileTreatMCAR) <- c("10", "30", "50")
PercentileEcogMCAR <- cbind(PercentileEcogMCAR10, PercentileEcogMCAR30, PercentileEcogMCAR50)
colnames(PercentileEcogMCAR) <- c("10", "30", "50")
PercentileNewVarMCAR <- cbind(PercentileNewVarMCAR10, PercentileNewVarMCAR30, PercentileNewVarMCAR50)
colnames(PercentileNewVarMCAR) <- c("10", "30", "50")

PercentileTreatMAR <- cbind(PercentileTreatMAR10, PercentileTreatMAR30, PercentileTreatMAR50)
colnames(PercentileTreatMAR) <- c("10", "30", "50")
PercentileEcogMAR <- cbind(PercentileEcogMAR10, PercentileEcogMAR30, PercentileEcogMAR50)
colnames(PercentileEcogMAR) <- c("10", "30", "50")
PercentileNewVarMAR <- cbind(PercentileNewVarMAR10, PercentileNewVarMAR30, PercentileNewVarMAR50)
colnames(PercentileNewVarMAR) <- c("10", "30", "50")

PercentileTreatMNARa <- cbind(PercentileTreatMNAR10a, PercentileTreatMNAR30a, PercentileTreatMNAR50a)
colnames(PercentileTreatMNARa) <- c("10", "30", "50")
PercentileEcogMNARa <- cbind(PercentileEcogMNAR10a, PercentileEcogMNAR30a, PercentileEcogMNAR50a)
colnames(PercentileEcogMNARa) <- c("10", "30", "50")
PercentileNewVarMNARa <- cbind(PercentileNewVarMNAR10a, PercentileNewVarMNAR30a, PercentileNewVarMNAR50a)
colnames(PercentileNewVarMNARa) <- c("10", "30", "50")

PercentileTreatMNARb <- cbind(PercentileTreatMNAR10b, PercentileTreatMNAR30b, PercentileTreatMNAR50b)
colnames(PercentileTreatMNARb) <- c("10", "30", "50")
PercentileEcogMNARb <- cbind(PercentileEcogMNAR10b, PercentileEcogMNAR30b, PercentileEcogMNAR50b)
colnames(PercentileEcogMNARb) <- c("10", "30", "50")
PercentileNewVarMNARb <- cbind(PercentileNewVarMNAR10b, PercentileNewVarMNAR30b, PercentileNewVarMNAR50b)
colnames(PercentileNewVarMNARb) <- c("10", "30", "50")

# fixed percentile treat should be for MCAR, MAR, and MNAR
write.csv(PercentileTreatMCAR, file.path(rootdir, "PercentileTreatMCAR10TreeV2.csv"))
write.csv(PercentileEcogMCAR, file.path(rootdir, "PercentileEcogMCAR10TreeV2.csv"))
write.csv(PercentileNewVarMCAR, file.path(rootdir, "PercentileNewVarMCAR10TreeV2.csv"))
write.csv(PercentileTreatMAR, file.path(rootdir, "PercentileTreatMAR10TreeV2.csv"))
write.csv(PercentileEcogMAR, file.path(rootdir, "PercentileEcogMAR10TreeV2.csv"))
write.csv(PercentileNewVarMAR, file.path(rootdir, "PercentileNewVarMAR10TreeV2.csv"))
write.csv(PercentileTreatMNARa, file.path(rootdir, "PercentileTreatMNAR10aTreeV2.csv"))
write.csv(PercentileEcogMNARa, file.path(rootdir, "PercentileEcogMNAR10aTreeV2.csv"))
write.csv(PercentileNewVarMNARa, file.path(rootdir, "PercentileNewVarMNAR10aTreeV2.csv"))
write.csv(PercentileTreatMNARb, file.path(rootdir, "PercentileTreatMNAR10bTreeV2.csv"))
write.csv(PercentileEcogMNARb, file.path(rootdir, "PercentileEcogMNAR10bTreeV2.csv"))
write.csv(PercentileNewVarMNARb, file.path(rootdir, "PercentileNewVarMNAR10bTreeV2.csv"))

MCAR10Avg <- MCAR10$avg
MCAR30Avg <- MCAR30$avg
MCAR50Avg <- MCAR50$avg
MAR10Avg <- MAR10$avg
MAR30Avg <- MAR30$avg
MAR50Avg <- MAR50$avg
MNAR10aAvg <- MNAR10a$avg
MNAR30aAvg <- MNAR30a$avg
MNAR50aAvg <- MNAR50a$avg
MNAR10bAvg <- MNAR10b$avg
MNAR30bAvg <- MNAR30b$avg
MNAR50bAvg <- MNAR50b$avg

write.csv(MCAR10Avg, file.path(rootdir, "MCAR10Avg10TreeV2.csv"))
write.csv(MCAR30Avg, file.path(rootdir, "MCAR30Avg10TreeV2.csv"))
write.csv(MCAR50Avg, file.path(rootdir, "MCAR50Avg10TreeV2.csv"))
write.csv(MAR10Avg, file.path(rootdir, "MAR10Avg10TreeV2.csv"))
write.csv(MAR30Avg, file.path(rootdir, "MAR30Avg10TreeV2.csv"))
write.csv(MAR50Avg, file.path(rootdir, "MAR50Avg10TreeV2.csv"))
write.csv(MNAR10aAvg, file.path(rootdir, "MNAR10aAvg10TreeV2.csv"))
write.csv(MNAR30aAvg, file.path(rootdir, "MNAR30aAvg10TreeV2.csv"))
write.csv(MNAR50aAvg, file.path(rootdir, "MNAR50aAvg10TreeV2.csv"))
write.csv(MNAR10bAvg, file.path(rootdir, "MNAR10bAvg10TreeV2.csv"))
write.csv(MNAR30bAvg, file.path(rootdir, "MNAR30bAvg10TreeV2.csv"))
write.csv(MNAR50bAvg, file.path(rootdir, "MNAR50bAvg10TreeV2.csv"))


# Autoencoder
filesMCAR10 <- list.files(file.path(rootdir, "AEMCAR/10"), full.names = TRUE)
filesMCAR30 <- list.files(file.path(rootdir, "AEMCAR/30"), full.names = TRUE)
filesMCAR50 <- list.files(file.path(rootdir, "AEMCAR/50"), full.names = TRUE)
filesMAR10 <- list.files(file.path(rootdir, "AEMAR/10"), full.names = TRUE)
filesMAR30 <- list.files(file.path(rootdir, "AEMAR/30"), full.names = TRUE)
filesMAR50 <- list.files(file.path(rootdir, "AEMAR/50"), full.names = TRUE)
filesMNAR10a <- list.files(file.path(rootdir, "AEMNAR1/10"), full.names = TRUE)
filesMNAR30a <- list.files(file.path(rootdir, "AEMNAR1/30"), full.names = TRUE)
filesMNAR50a <- list.files(file.path(rootdir, "AEMNAR1/50"), full.names = TRUE)
filesMNAR10b <- list.files(file.path(rootdir, "AEMNAR2/10"), full.names = TRUE)
filesMNAR30b <- list.files(file.path(rootdir, "AEMNAR2/30"), full.names = TRUE)
filesMNAR50b <- list.files(file.path(rootdir, "AEMNAR2/50"), full.names = TRUE)

iter <- 1000
mechList <- c("MCAR", "MAR", "MNAR1", "MNAR2")
proportionList <- c(10, 30, 50)
# J is number of imputations
J <- 10

for (w in 1:length(mechList)) {
  mech <- mechList[w]
  print (mech)
  for (y in 1:length(proportionList)) {
    propName <- proportionList[y]
    print (propName)
    truth <- read.csv(file.path("/project/flatiron_ucc/programs/kylie/RunMe2/datasets/trueEff", mech, propName, "propMiss_trueEffs1.csv"))[[1]]
    trueTreat <- truth[1]
    trueEcog <- truth[2]
    trueNewVar <- truth[3]
    fileListName <- paste0("files", mech, propName)
    fileCur <- get(fileListName)

    biasTreatAE <- c()
    biasEcogAE <- c()
    biasNewVarAE <- c()
    VarTreatAE <- c()
    VarEcogAE <- c()
    VarNewVarAE <- c()
    covereageAETreat <- c()
    covereageAEEcog <- c()
    covereageAENewVar <- c()
    mseAETreat <- c()
    mseAEEcog <- c()
    mseAENewVar <- c()

    for (i in 1:iter) {
      print(i)
      daeTreatList <- c()
      daeEcogList <- c()
      daeNewVarList <- c()
      VarTreatdaeList <- c()
      VarEcogdaeList <- c()
      VarNewVardaeList <- c()
      for (j in 1:J) {
        data <- read.csv(file.path("/project/flatiron_ucc/programs/kylie/RunMe2/final_results", paste0("AE", mech), propName, paste0("result", i, "_num_", j, ".csv")), row.names = 1) # read in instead

        daeFit <- coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + newVar, data = data)
        treatRow <- which(names(daeFit$coefficients) == "treat")
        ecogRow <- which(names(daeFit$coefficients) == "b.ecogvalue")
        newVarRow <- which(names(daeFit$coefficients) == "newVar")
        daeTreatList[j] <- daeFit$coefficients[treatRow]
        daeEcogList[j] <- daeFit$coefficients[ecogRow]
        daeNewVarList[j] <- daeFit$coefficients[newVarRow]
        VarTreatdaeList[j] <- (summary(daeFit)$coefficients[treatRow, 3])^2
        VarEcogdaeList[j] <- (summary(daeFit)$coefficients[ecogRow, 3])^2
        VarNewVardaeList[j] <- (summary(daeFit)$coefficients[newVarRow, 3])^2
      }

      # pool them
      daetreat <- mean(daeTreatList)
      daeEcog <- mean(daeEcogList)
      daeNewVar <- mean(daeNewVarList)

      biasTreatAE[i] <- as.numeric(daetreat - trueTreat)
      biasEcogAE[i] <- as.numeric(daeEcog - trueEcog)
      biasNewVarAE[i] <- as.numeric(daeNewVar - trueNewVar)

      UbarTreat <- mean(VarTreatdaeList)
      UbarEcog <- mean(VarEcogdaeList)
      UbarNewVar <- mean(VarNewVardaeList)

      BTreat <- var(daeTreatList)
      BEcog <- var(daeEcogList)
      BNewVar <- var(daeNewVarList)

      # Variance formula for MI
      VarTreatAE[i] <- UbarTreat + (1 + (1 / J)) * BTreat
      VarEcogAE[i] <- UbarEcog + (1 + (1 / J)) * BEcog
      VarNewVarAE[i] <- UbarNewVar + (1 + (1 / J)) * BNewVar

      CIAELowerTreat <- daetreat - (qnorm(.975) * sqrt(VarTreatAE[i]))
      CIAEUpperTreat <- daetreat + (qnorm(.975) * sqrt(VarTreatAE[i]))
      CIAELowerEcog <- daeEcog - (qnorm(.975) * sqrt(VarEcogAE[i]))
      CIAEUpperEcog <- daeEcog + (qnorm(.975) * sqrt(VarEcogAE[i]))
      CIAELowerNewVar <- daeNewVar - (qnorm(.975) * sqrt(VarNewVarAE[i]))
      CIAEUpperNewVar <- daeNewVar + (qnorm(.975) * sqrt(VarNewVarAE[i]))

      covereageAETreat[i] <- 1 * between(c(trueTreat), CIAELowerTreat, CIAEUpperTreat)
      covereageAEEcog[i] <- 1 * between(c(trueEcog), CIAELowerEcog, CIAEUpperEcog)
      covereageAENewVar[i] <- 1 * between(c(trueNewVar), CIAELowerNewVar, CIAEUpperNewVar)

      mseAETreat[i] <- VarTreatAE[i] + biasTreatAE[i]^2
      mseAEEcog[i] <- VarEcogAE[i] + biasEcogAE[i]^2
      mseAENewVar[i] <- VarNewVarAE[i] + biasNewVarAE[i]^2
    }

    PercentileTreatAE <- data.frame (c(mean (biasTreatAE)-1.96*sd(biasTreatAE)/sqrt (iter),mean (biasTreatAE)+1.96*sd(biasTreatAE)/sqrt (iter)))
    PercentileEcogAE <- data.frame (c(mean (biasEcogAE)-1.96*sd(biasEcogAE)/sqrt (iter),mean (biasEcogAE)+1.96*sd(biasEcogAE)/sqrt (iter)))
    PercentileNewVarAE <- data.frame (c(mean (biasNewVarAE)-1.96*sd(biasNewVarAE)/sqrt (iter),mean (biasNewVarAE)+1.96*sd(biasNewVarAE)/sqrt (iter)))
    biasTreatAvg <- mean(biasTreatAE)
    biasEcogAvg <- mean(biasEcogAE)
    biasNewVarAvg <- mean(biasNewVarAE)
    VarTreatAvg <- mean(VarTreatAE)
    VarEcogAvg <- mean(VarEcogAE)
    VarNewVarAvg <- mean(VarNewVarAE)
    coverForestTreat <- mean(covereageAETreat)
    coverForestEcog <- mean(covereageAEEcog)
    coverForestNewVar <- mean(covereageAENewVar)
    mseAETreatAvg <- mean(mseAETreat)
    mseAEEcogAvg <- mean(mseAEEcog)
    mseAENewVarAvg <- mean(mseAENewVar)
    results <- rbind(cbind(biasTreatAvg, biasEcogAvg, biasNewVarAvg), cbind(sqrt(VarTreatAvg), sqrt(VarEcogAvg), sqrt(VarNewVarAvg)), cbind(coverForestTreat, coverForestEcog, coverForestNewVar), cbind(mseAETreatAvg, mseAEEcogAvg, mseAENewVarAvg))
    rownames(results) <- c("Bias", "se", "coverage", "MSE")
    colnames(results) <- c("AETreat", "AEEcog", "AENewVar")
    percentileTreat <- PercentileTreatAE
    percentileEcog <- PercentileEcogAE
    percentileNewVar <- PercentileNewVarAE
    percentile <- cbind(percentileTreat, percentileEcog, percentileNewVar)
    write.csv(results, file.path(rootdir, paste0("AEresult10Tree", mech, propName, ".csv")))
    write.csv(percentile, file.path(rootdir, paste0("AEpercentile10Tree", mech, propName, ".csv")))
  }
}