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
  df <- lapply(df, function(x){
    newx = x
    newx$X <- NULL
    return (newx)
  })
  dfAvg <- Reduce("+", df) / length(df)
  rownames(dfAvg) <- c("bias", "se", "coverage", "MSE")
  biasList <- (lapply(df, function(x) x <- x[1, ]))
  biasTable <- bind_rows(biasList, .id = "column_label")
  biasTable <- biasTable[, -1]

  PercentileOracleTreat <- data.frame(quantile(biasTable[, 1], c(0.025, 0.975)))
  PercentileOracleEcog <- data.frame(quantile(biasTable[, 2], c(0.025, 0.975)))
  PercentileOracleVar1 <- data.frame(quantile(biasTable[, 3], c(0.025, 0.975)))
  PercentileOracleVar2 <- data.frame(quantile(biasTable[, 4], c(0.025, 0.975)))

  PercentileCCTreat <- data.frame(quantile(biasTable[, 5], c(0.025, 0.975)))
  PercentileCCEcog <- data.frame(quantile(biasTable[, 6], c(0.025, 0.975)))
  PercentileCCVar1 <- data.frame(quantile(biasTable[, 7], c(0.025, 0.975)))
  PercentileCCVar2 <- data.frame(quantile(biasTable[, 8], c(0.025, 0.975)))

  PercentileMICETreat <- data.frame(quantile(biasTable[, 9], c(0.025, 0.975)))
  PercentileMICEEcog <- data.frame(quantile(biasTable[, 10], c(0.025, 0.975)))
  PercentileMICEVar1 <- data.frame(quantile(biasTable[, 11], c(0.025, 0.975)))
  PercentileMICEVar2 <- data.frame(quantile(biasTable[, 12], c(0.025, 0.975)))

  PercentileForestTreat <- data.frame(quantile(biasTable[, 13], c(0.025, 0.975)))
  PercentileForestEcog <- data.frame(quantile(biasTable[, 14], c(0.025, 0.975)))
  PercentileForestVar1 <- data.frame(quantile(biasTable[, 15], c(0.025, 0.975)))
  PercentileForestVar2 <- data.frame(quantile(biasTable[, 16], c(0.025, 0.975)))
  
  colnames(PercentileOracleTreat) <- ""
  colnames(PercentileOracleEcog) <- ""
  colnames(PercentileOracleVar1) <- ""
  colnames(PercentileOracleVar2) <- ""
  colnames(PercentileCCTreat) <- ""
  colnames(PercentileCCEcog) <- ""
  colnames(PercentileCCVar1) <- ""
  colnames(PercentileCCVar2) <- ""
  colnames(PercentileMICETreat) <- ""
  colnames(PercentileMICEEcog) <- ""
  colnames(PercentileMICEVar1) <- ""
  colnames(PercentileMICEVar2) <- ""
  colnames(PercentileForestTreat) <- ""
  colnames(PercentileForestEcog) <- ""
  colnames(PercentileForestVar1) <- ""
  colnames(PercentileForestVar2) <- ""

  PercentileTreatDf <- rbind(PercentileOracleTreat, PercentileCCTreat, PercentileMICETreat, PercentileForestTreat)
  rownames(PercentileTreatDf) <- c("OracleLower", "OracleUpper", "CCLower", "CCUpper", "MICELower", "MICEUpper", "ForestLower", "ForestUpper")

  PercentileEcogDf <- rbind(PercentileOracleEcog, PercentileCCEcog, PercentileMICEEcog, PercentileForestEcog)
  rownames(PercentileEcogDf) <- c("OracleLower", "OracleUpper", "CCLower", "CCUpper", "MICELower", "MICEUpper", "ForestLower", "ForestUpper")

  PercentileVar1Df <- rbind(PercentileOracleVar1, PercentileCCVar1, PercentileMICEVar1, PercentileForestVar1)
  rownames(PercentileVar1Df) <- c("OracleLower", "OracleUpper", "CCLower", "CCUpper", "MICELower", "MICEUpper", "ForestLower", "ForestUpper")

  PercentileVar2Df <- rbind(PercentileOracleVar2, PercentileCCVar2, PercentileMICEVar2, PercentileForestVar2)
  rownames(PercentileVar2Df) <- c("OracleLower", "OracleUpper", "CCLower", "CCUpper", "MICELower", "MICEUpper", "ForestLower", "ForestUpper")

  return(list("treat" = PercentileTreatDf, "ecog" = PercentileEcogDf, "Var1" = PercentileVar1Df, "Var2" = PercentileVar2Df, "avg" = dfAvg))
}

first = unlist (lapply (filesMCAR10, function (x) strsplit (x, '/')[[1]][10]))
second = unlist (lapply (first, function (x) strsplit (x, '.csv')[[1]]))
third = unlist (lapply (second, function (x) strsplit (x, 'result')[[1]][2]))
keep = which (third%in%as.character (1:iter)==TRUE)

MCAR10 <- makeTable(filesMCAR10[keep])
PercentileTreatMCAR10 <- MCAR10$treat
PercentileEcogMCAR10 <- MCAR10$ecog
PercentileVar1MCAR10 <- MCAR10$Var1
PercentileVar2MCAR10 <- MCAR10$Var2

MCAR30 <- makeTable(filesMCAR30[keep])
PercentileTreatMCAR30 <- MCAR30$treat
PercentileEcogMCAR30 <- MCAR30$ecog
PercentileVar1MCAR30 <- MCAR30$Var1
PercentileVar2MCAR30 <- MCAR30$Var2

MCAR50 <- makeTable(filesMCAR50[keep])
PercentileTreatMCAR50 <- MCAR50$treat
PercentileEcogMCAR50 <- MCAR50$ecog
PercentileVar1MCAR50 <- MCAR50$Var1
PercentileVar2MCAR50 <- MCAR50$Var2

MAR10 <- makeTable(filesMAR10[keep])
PercentileTreatMAR10 <- MAR10$treat
PercentileEcogMAR10 <- MAR10$ecog
PercentileVar1MAR10 <- MAR10$Var1
PercentileVar2MAR10 <- MAR10$Var2

MAR30 <- makeTable(filesMAR30[keep])
PercentileTreatMAR30 <- MAR30$treat
PercentileEcogMAR30 <- MAR30$ecog
PercentileVar1MAR30 <- MAR30$Var1
PercentileVar2MAR30 <- MAR30$Var2

MAR50 <- makeTable(filesMAR50[keep])
PercentileTreatMAR50 <- MAR50$treat
PercentileEcogMAR50 <- MAR50$ecog
PercentileVar1MAR50 <- MAR50$Var1
PercentileVar2MAR50 <- MAR50$Var2

MNAR10a <- makeTable(filesMNAR10a[keep])
PercentileTreatMNAR10a <- MNAR10a$treat
PercentileEcogMNAR10a <- MNAR10a$ecog
PercentileVar1MNAR10a <- MNAR10a$Var1
PercentileVar2MNAR10a <- MNAR10a$Var2

MNAR30a <- makeTable(filesMNAR30a[keep])
PercentileTreatMNAR30a <- MNAR30a$treat
PercentileEcogMNAR30a <- MNAR30a$ecog
PercentileVar1MNAR30a <- MNAR30a$Var1
PercentileVar2MNAR30a <- MNAR30a$Var2

MNAR50a <- makeTable(filesMNAR50a[keep])
PercentileTreatMNAR50a <- MNAR50a$treat
PercentileEcogMNAR50a <- MNAR50a$ecog
PercentileVar1MNAR50a <- MNAR50a$Var1
PercentileVar2MNAR50a <- MNAR50a$Var2

MNAR10b <- makeTable(filesMNAR10b[keep])
PercentileTreatMNAR10b <- MNAR10b$treat
PercentileEcogMNAR10b <- MNAR10b$ecog
PercentileVar1MNAR10b <- MNAR10b$Var1
PercentileVar2MNAR10b <- MNAR10b$Var2

MNAR30b <- makeTable(filesMNAR30b[keep])
PercentileTreatMNAR30b <- MNAR30b$treat
PercentileEcogMNAR30b <- MNAR30b$ecog
PercentileVar1MNAR30b <- MNAR30b$Var1
PercentileVar2MNAR30b <- MNAR30b$Var2

MNAR50b <- makeTable(filesMNAR50b[keep])
PercentileTreatMNAR50b <- MNAR50b$treat
PercentileEcogMNAR50b <- MNAR50b$ecog
PercentileVar1MNAR50b <- MNAR50b$Var1
PercentileVar2MNAR50b <- MNAR50b$Var2

PercentileTreatMCAR <- cbind(PercentileTreatMCAR10, PercentileTreatMCAR30, PercentileTreatMCAR50)
colnames(PercentileTreatMCAR) <- c("10", "30", "50")
PercentileEcogMCAR <- cbind(PercentileEcogMCAR10, PercentileEcogMCAR30, PercentileEcogMCAR50)
colnames(PercentileEcogMCAR) <- c("10", "30", "50")
PercentileVar1MCAR <- cbind(PercentileVar1MCAR10, PercentileVar1MCAR30, PercentileVar1MCAR50)
colnames(PercentileVar1MCAR) <- c("10", "30", "50")
PercentileVar2MCAR <- cbind(PercentileVar2MCAR10, PercentileVar2MCAR30, PercentileVar2MCAR50)
colnames(PercentileVar2MCAR) <- c("10", "30", "50")

PercentileTreatMAR <- cbind(PercentileTreatMAR10, PercentileTreatMAR30, PercentileTreatMAR50)
colnames(PercentileTreatMAR) <- c("10", "30", "50")
PercentileEcogMAR <- cbind(PercentileEcogMAR10, PercentileEcogMAR30, PercentileEcogMAR50)
colnames(PercentileEcogMAR) <- c("10", "30", "50")
PercentileVar1MAR <- cbind(PercentileVar1MAR10, PercentileVar1MAR30, PercentileVar1MAR50)
colnames(PercentileVar1MAR) <- c("10", "30", "50")
PercentileVar2MAR <- cbind(PercentileVar2MAR10, PercentileVar2MAR30, PercentileVar2MAR50)
colnames(PercentileVar2MAR) <- c("10", "30", "50")

PercentileTreatMNARa <- cbind(PercentileTreatMNAR10a, PercentileTreatMNAR30a, PercentileTreatMNAR50a)
colnames(PercentileTreatMNARa) <- c("10", "30", "50")
PercentileEcogMNARa <- cbind(PercentileEcogMNAR10a, PercentileEcogMNAR30a, PercentileEcogMNAR50a)
colnames(PercentileEcogMNARa) <- c("10", "30", "50")
PercentileVar1MNARa <- cbind(PercentileVar1MNAR10a, PercentileVar1MNAR30a, PercentileVar1MNAR50a)
colnames(PercentileVar1MNARa) <- c("10", "30", "50")
PercentileVar2MNARa <- cbind(PercentileVar2MNAR10a, PercentileVar2MNAR30a, PercentileVar2MNAR50a)
colnames(PercentileVar2MNARa) <- c("10", "30", "50")

PercentileTreatMNARb <- cbind(PercentileTreatMNAR10b, PercentileTreatMNAR30b, PercentileTreatMNAR50b)
colnames(PercentileTreatMNARb) <- c("10", "30", "50")
PercentileEcogMNARb <- cbind(PercentileEcogMNAR10b, PercentileEcogMNAR30b, PercentileEcogMNAR50b)
colnames(PercentileEcogMNARb) <- c("10", "30", "50")
PercentileVar1MNARb <- cbind(PercentileVar1MNAR10b, PercentileVar1MNAR30b, PercentileVar1MNAR50b)
colnames(PercentileVar1MNARb) <- c("10", "30", "50")
PercentileVar2MNARb <- cbind(PercentileVar2MNAR10b, PercentileVar2MNAR30b, PercentileVar2MNAR50b)
colnames(PercentileVar2MNARb) <- c("10", "30", "50")

# fixed percentile treat should be for MCAR, MAR, and MNAR
write.csv(PercentileTreatMCAR, file.path(rootdir, "PercentileTreatMCAR10TreeV2.csv"))
write.csv(PercentileEcogMCAR, file.path(rootdir, "PercentileEcogMCAR10TreeV2.csv"))
write.csv(PercentileVar1MCAR, file.path(rootdir, "PercentileVar1MCAR10TreeV2.csv"))
write.csv(PercentileVar2MCAR, file.path(rootdir, "PercentileVar2MCAR10TreeV2.csv"))

write.csv(PercentileTreatMAR, file.path(rootdir, "PercentileTreatMAR10TreeV2.csv"))
write.csv(PercentileEcogMAR, file.path(rootdir, "PercentileEcogMAR10TreeV2.csv"))
write.csv(PercentileVar1MAR, file.path(rootdir, "PercentileVar1MAR10TreeV2.csv"))
write.csv(PercentileVar2MAR, file.path(rootdir, "PercentileVar2MAR10TreeV2.csv"))

write.csv(PercentileTreatMNARa, file.path(rootdir, "PercentileTreatMNAR10aTreeV2.csv"))
write.csv(PercentileEcogMNARa, file.path(rootdir, "PercentileEcogMNAR10aTreeV2.csv"))
write.csv(PercentileVar1MNARa, file.path(rootdir, "PercentileVar1MNAR10aTreeV2.csv"))
write.csv(PercentileVar2MNARa, file.path(rootdir, "PercentileVar2MNAR10aTreeV2.csv"))

write.csv(PercentileTreatMNARb, file.path(rootdir, "PercentileTreatMNAR10bTreeV2.csv"))
write.csv(PercentileEcogMNARb, file.path(rootdir, "PercentileEcogMNAR10bTreeV2.csv"))
write.csv(PercentileVar1MNARb, file.path(rootdir, "PercentileVar1MNAR10bTreeV2.csv"))
write.csv(PercentileVar2MNARb, file.path(rootdir, "PercentileVar2MNAR10bTreeV2.csv"))

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
    truth <- read.csv(file.path("/project/flatiron_ucc/programs/kylie/RunMe2/datasets/trueEff", mech, propName, "propMiss_trueEffs1.csv"))
    trueCoefs <- truth[names (truth)%in%c ("TREAT", "b.ecogvalue", "var1", "var2")]
    if (mech=='MCAR'){
      fileListName <- paste0("filesMCAR", propName)
    }
    if (mech=='MAR'){
      fileListName <- paste0("filesMAR", propName)
    }
    if (mech=='MNAR1'){
      fileListName <- paste0("filesMNAR", propName, "a")
    }
    if (mech=='MNAR2'){
      fileListName <- paste0("filesMNAR", propName, "b")
    }
    fileCur <- get(fileListName)

    biasTreatAE <- c()
    biasEcogAE <- c()
    biasVar1AE <- c()
    biasVar2AE <- c()
    VarTreatAE <- c()
    VarEcogAE <- c()
    VarVar1AE <- c()
    VarVar2AE <- c()
    covereageAETreat <- c()
    covereageAEEcog <- c()
    covereageAEVar1 <- c()
    covereageAEVar2 <- c()
    mseAETreat <- c()
    mseAEEcog <- c()
    mseAEVar1 <- c()
    mseAEVar2 <- c()

    for (i in 1:iter) {
      print(i)
      daeTreatList <- c()
      daeEcogList <- c()
      daeVar1List <- c()
      daeVar2List <- c()
      VarTreatdaeList <- c()
      VarEcogdaeList <- c()
      VarVar1daeList <- c()
      VarVar2daeList <- c()
      for (j in 1:J) {
        data <- read.csv(file.path("/project/flatiron_ucc/programs/kylie/RunMe2/final_results", paste0("AE", mech), propName, paste0("result", i, "_num_", j, ".csv")), row.names = 1) # read in instead

        daeFit <- coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + var1 + var2, data = data)
        treatRow <- which(names(daeFit$coefficients) == "treat")
        ecogRow <- which(names(daeFit$coefficients) == "b.ecogvalue")
        Var1Row <- which(names(daeFit$coefficients) == "var1")
        Var2Row <- which(names(daeFit$coefficients) == "var2")
        daeTreatList[j] <- daeFit$coefficients[treatRow]
        daeEcogList[j] <- daeFit$coefficients[ecogRow]
        daeVar1List[j] <- daeFit$coefficients[Var1Row]
        daeVar2List[j] <- daeFit$coefficients[Var2Row]
        VarTreatdaeList[j] <- (summary(daeFit)$coefficients[treatRow, 3])^2
        VarEcogdaeList[j] <- (summary(daeFit)$coefficients[ecogRow, 3])^2
        VarVar1daeList[j] <- (summary(daeFit)$coefficients[Var1Row, 3])^2
        VarVar2daeList[j] <- (summary(daeFit)$coefficients[Var2Row, 3])^2
      }

      # pool them
      daetreat <- mean(daeTreatList)
      daeEcog <- mean(daeEcogList)
      daeVar1 <- mean(daeVar1List)
      daeVar2 <- mean(daeVar2List)

      biasTreatAE[i] <- as.numeric(daetreat - trueCoefs[1])
      biasEcogAE[i] <- as.numeric(daeEcog - trueCoefs[2])
      biasVar1AE[i] <- as.numeric(daeVar1 - trueCoefs[3])
      biasVar2AE[i] <- as.numeric(daeVar2 - trueCoefs[4])

      UbarTreat <- mean(VarTreatdaeList)
      UbarEcog <- mean(VarEcogdaeList)
      UbarVar1 <- mean(VarVar1daeList)
      UbarVar2 <- mean(VarVar2daeList)

      BTreat <- var(daeTreatList)
      BEcog <- var(daeEcogList)
      BVar1 <- var(daeVar1List)
      BVar2 <- var(daeVar2List)

      # Variance formula for MI
      VarTreatAE[i] <- UbarTreat + (1 + (1 / J)) * BTreat
      VarEcogAE[i] <- UbarEcog + (1 + (1 / J)) * BEcog
      VarVar1AE[i] <- UbarVar1 + (1 + (1 / J)) * BVar1
      VarVar2AE[i] <- UbarVar2 + (1 + (1 / J)) * BVar2

      CIAELowerTreat <- daetreat - (qnorm(.975) * sqrt(VarTreatAE[i]))
      CIAEUpperTreat <- daetreat + (qnorm(.975) * sqrt(VarTreatAE[i]))
      CIAELowerEcog <- daeEcog - (qnorm(.975) * sqrt(VarEcogAE[i]))
      CIAEUpperEcog <- daeEcog + (qnorm(.975) * sqrt(VarEcogAE[i]))
      CIAELowerVar1 <- daeVar1 - (qnorm(.975) * sqrt(VarVar1AE[i]))
      CIAEUpperVar1 <- daeVar1 + (qnorm(.975) * sqrt(VarVar1AE[i]))
      CIAELowerVar2 <- daeVar2 - (qnorm(.975) * sqrt(VarVar2AE[i]))
      CIAEUpperVar2 <- daeVar2 + (qnorm(.975) * sqrt(VarVar2AE[i]))

      covereageAETreat[i] <- 1 * between(c(trueCoefs[1]), CIAELowerTreat, CIAEUpperTreat)
      covereageAEEcog[i] <- 1 * between(c(trueCoefs[2]), CIAELowerEcog, CIAEUpperEcog)
      covereageAEVar1[i] <- 1 * between(c(trueCoefs[3]), CIAELowerVar1, CIAEUpperVar1)
      covereageAEVar2[i] <- 1 * between(c(trueCoefs[4]), CIAELowerVar2, CIAEUpperVar2)

      mseAETreat[i] <- VarTreatAE[i] + biasTreatAE[i]^2
      mseAEEcog[i] <- VarEcogAE[i] + biasEcogAE[i]^2
      mseAEVar1[i] <- VarVar1AE[i] + biasVar1AE[i]^2
      mseAEVar2[i] <- VarVar2AE[i] + biasVar2AE[i]^2
    }

    PercentileTreatAE <- data.frame (quantile (biasTreatAE, c(0.025, 0.975)))
    PercentileEcogAE <- data.frame (quantile (biasEcogAE, c(0.025, 0.975)))
    PercentileVar1AE <- data.frame (quantile (biasVar1AE, c(0.025, 0.975)))
    PercentileVar2AE <- data.frame (quantile (biasVar2AE, c(0.025, 0.975)))
    
    biasTreatAvg <- mean(biasTreatAE)
    biasEcogAvg <- mean(biasEcogAE)
    biasVar1Avg <- mean(biasVar1AE)
    biasVar2Avg <- mean(biasVar2AE)

    VarTreatAvg <- mean(VarTreatAE)
    VarEcogAvg <- mean(VarEcogAE)
    VarVar1Avg <- mean(VarVar1AE)
    VarVar2Avg <- mean(VarVar2AE)

    coverAETreat <- mean(covereageAETreat)
    coverAEEcog <- mean(covereageAEEcog)
    coverAEVar1 <- mean(covereageAEVar1)
    coverAEVar2 <- mean(covereageAEVar2)

    mseAETreatAvg <- mean(mseAETreat)
    mseAEEcogAvg <- mean(mseAEEcog)
    mseAEVar1Avg <- mean(mseAEVar1)
    mseAEVar2Avg <- mean(mseAEVar2)

    results <- rbind(cbind(biasTreatAvg, biasEcogAvg, biasVar1Avg, biasVar2Avg), cbind(sqrt(VarTreatAvg), sqrt(VarEcogAvg), sqrt(VarVar1Avg), sqrt(VarVar2Avg)), cbind(coverAETreat, coverAEEcog, coverAEVar1, coverAEVar2), cbind(mseAETreatAvg, mseAEEcogAvg, mseAEVar1Avg, mseAEVar2Avg))
    rownames(results) <- c("Bias", "se", "coverage", "MSE")
    colnames(results) <- c("AETreat", "AEEcog", "AEVar1", "AEVar2")
    percentileTreat <- PercentileTreatAE
    percentileEcog <- PercentileEcogAE
    percentileVar1 <- PercentileVar1AE
    percentileVar2 <- PercentileVar2AE
    percentile <- cbind(percentileTreat, percentileEcog, percentileVar1, percentileVar2)
    write.csv(results, file.path(rootdir, paste0("AEresult10Tree", mech, propName, ".csv")))
    write.csv(percentile, file.path(rootdir, paste0("AEpercentile10Tree", mech, propName, ".csv")))
  }
}