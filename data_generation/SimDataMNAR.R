library(Plasmode)
library(survival)
library(dplyr)
library(mice)
library(tidyverse)

rootdir <- "/project/flatiron_ucc/programs/kylie/RunMe"
system(paste0("mkdir ", file.path(rootdir, "datasets")))
system(paste0("mkdir ", file.path(rootdir, "datasets", "mDats")))
system(paste0("mkdir ", file.path(rootdir, "datasets", "cDats")))
system (paste0 ('mkdir ', file.path (rootdir, 'datasets', 'trueEff')))
system(paste0("mkdir ", file.path(rootdir, "datasets", "mDats", "MNAR")))
system(paste0("mkdir ", file.path(rootdir, "datasets", "cDats", "MNAR")))
system (paste0 ('mkdir ', file.path (rootdir, 'datasets', 'trueEff', 'MNAR')))

proportionList <- c(10, 30, 50)
system (paste0 ('mkdir ', file.path (rootdir, 'datasets', 'trueEff', 'MNAR', proportionList[1])))
system (paste0 ('mkdir ', file.path (rootdir, 'datasets', 'trueEff', 'MNAR', proportionList[2])))
system (paste0 ('mkdir ', file.path (rootdir, 'datasets', 'trueEff', 'MNAR', proportionList[3])))
iter <- 1000
effOR <- 0.5
# Optional stages for debuging / descriptive stats
assess.missing <- F
# set number of Plasmode simulations
nsim <- 1
# Magic intercepts for missingness proportions
b0 <- c(-21.40, -19.26, -17.83)

for (w in 1:3) {
  set.seed(123)

  ################################################################################
  # Step 1: load data and select variables for analysis
  # Outcome is overall survival using treatment start -> death
  # Primary exposure is immunotherapy vs chemotherapy
  # (should we do carboplatin vs immunotherapy?).
  # Covariables should be selected to build structure for the simulated exposure
  # + outcome and to be included in the propensity score adjustment.
  ################################################################################

  sdat <- read.csv(file.path("/project/flatiron_ucc/programs/kylie/plasmode_dan/deident_comorb_flatiron_ucc_2020-11-03.csv")) %>%
    # what variables are we using? Other project recommends albumin, hemoglobin, number of mets, and BMI
    dplyr::mutate(
      age = baseline_year - birthyear,
      # Inverse of dead for plasmode sim
      ndead = 1 - dead,
      # Set treatment class to immuno vs non-immuno
      treat = ifelse(trtclass1 == "immuno", 1, 0),
      cmonth = lastcontactmonth + 1,
      # Create dummy vars
      reth_black = ifelse(race.ethnicity == "Black", 1, 0),
      reth_hisp = ifelse(race.ethnicity == "Hispanic or Latino", 1, 0),
      reth_oth = ifelse(race.ethnicity == "Other", 1, 0),
      genderf = ifelse(gender == "F", 1, 0),
      practypec = ifelse(practicetype == "COMMUNITY", 1, 0),
      smokey = ifelse(smokingstatus == "History of smoking", 1, 0),
      dgradeh = ifelse(diseasegrade == "High grade (G2/G3/G4)", 1, 0),
      surgery = ifelse(surgery == T, 1, 0),
      site_ureter = ifelse(primarysite == "Ureter", 1, 0),
      site_renal = ifelse(primarysite == "Renal Pelvis", 1, 0),
      site_urethra = ifelse(primarysite == "Urethra", 1, 0),
      # Change unknown to NA
      race.ethnicity = ifelse(race.ethnicity == "Unknown", NA, race.ethnicity),
      diseasegrade = ifelse(diseasegrade == "Unknown/not documented", NA, diseasegrade),
      smokingstatus = ifelse(smokingstatus == "Unknown", NA, smokingstatus),
      groupstage = ifelse(groupstage == "Unknown/not documented", NA, groupstage),
    ) %>%
    dplyr::select(
      patientid, genderf, reth_black, reth_hisp, reth_oth, age, practypec, site_ureter, site_renal, site_urethra, smokey, dead, ndead, cmonth, dgradeh, surgery, treat, b.ecogvalue
      # elixhauser score vars
      # chf, carit, valv, pcd, pvd, hypunc, hypc, para, ond, cpd,
      # diabunc, diabc, hypothy, rf, ld, pud, aids, lymph, metacanc,
      # solidtum, rheumd, coag, obes, wloss, fed, blane, dane, alcohol,
      # drug, psycho, depre
    )

  # Remove incomplete data. N = 3858 (Previous N = 6102)
  sdat.c <- sdat %>% na.omit()
  # Round the few half scores up to next integer
  sdat.c$b.ecogvalue <- ceiling(sdat.c$b.ecogvalue)

  # Check NA proportions
  # ECOG is the main reason for missing data.
  if (assess.missing) {
    natab <- function(x) {
      if (anyNA(x)) {
        table(is.na(x))
      } else {
        (NA)
      }
    }
    na.omit.list <- function(y) {
      return(y[!sapply(y, function(x) all(is.na(x)))])
    }
    apply(sdat, 2, natab) %>% na.omit.list()
    # Missing data visualization (generate and output to file)
    # 1778 missing due to ECOG alone, 384 missing due to disease grade alone,
    # 264 due to race/ethnicity alone, 465 missing due to ECOG with others,
    # # 34 due to race/eth + disease grade combo, 1 missing due to gender
    # p1 <- gg_miss_upset(sdat)
    # tiff("U:/projects/FLATIRON missing data sim/missing-data-pattern.tiff", units="in",
    #      width=7, height=6, res=300)
    # p1
    # dev.off()
    # Compare event rate in base data to incomplete data
    brate.b <- pyears(Surv(sdat$cmonth, sdat$dead) ~ 1, scale = 12)
    brate.cc <- pyears(Surv(sdat.c$cmonth, sdat.c$dead) ~ 1, scale = 12)
    # Check exposure prevalence in base data and incomplete data
    table(sdat$treat) %>% prop.table()
    table(sdat.c$treat) %>% prop.table()
  }

  # Step 2: Estimate associations with outcome and censoring
  # One model for the hazard of the outcome event
  # One model for the hazard of censoring (this is simply the reverse of the
  # outcome, aka 1-Y)

  # Note: the Plasmode package actually selects the second independent variable in
  # the formula for the exposure, rather than the first if you are having the
  # package run the model. If you run it yourself it gets it right.

  # Build outcome hazard
  os1 <- coxph(Surv(sdat.c$cmonth, sdat.c$dead) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age, data = sdat.c, x = T)

  # Censoring hazard
  oc1 <- coxph(Surv(sdat.c$cmonth, sdat.c$ndead) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age, data = sdat.c, x = T)

  # KAL: WHY ARE THESE THE SAME AS ABOVE?
  # Commenting out for now
  # Get true estimates
  #  os1True <- coxph(Surv(sdat.c$cmonth, sdat.c$dead) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age, data = sdat.c, x = T)

  #  oc1True <- coxph(Surv(sdat.c$cmonth, sdat.c$ndead) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age, data = sdat.c, x = T)

  # Simulate dichotomous variable zeta with an exposure probability of 20%
  # THE FOLLOWING WAS KYLIE'S CODE:
  #  time <- ifelse(sdat.c$dead == 0, 111, sdat.c$cmonth)

#  sdat.c$zeta <- rnorm(n = nrow(sdat.c), mean = sdat.c$b.ecogvalue)
  sdat.c$zeta <- rnorm(n = nrow(sdat.c), mean = sdat.c$b.ecogvalue)


  sor <- PlasmodeSur(
    objectOut = os1,
    objectCen = oc1,
    idVar = sdat.c$patientid,
    effectOR = effOR,
    nsim = nsim,
    size = 1000000
  )
  ### Generate missingness in simulated data
  # Remove unneeded variables from primary data
  sdat.m <- sdat.c %>% select(-ndead)
  # initialize empty storage list
  s.list <- list()
  # Pull simulation data
  wdat <- sor[["Sim_Data"]] %>%
    select(paste("ID", nsim, sep = ""), paste("EVENT", nsim, sep = ""), paste("TIME", nsim, sep = ""))

  names(wdat) <- c("patientid", "event", "time")

  # Attach full data to simulation i
  id <- wdat$patientid
  wdat <- dplyr::left_join(wdat, sdat.m, by = "patientid") %>%
    select(-patientid) %>%
    mutate(event = ifelse(event == T, 1, 0))
  dealcol <- which(colnames(wdat) == "dead")
  cmonthcol <- which(colnames(wdat) == "cmonth")
  wdat <- wdat[, -c(dealcol, cmonthcol)]
  haz <- nelsonaalen(wdat, time, event)
  wdat <- add_column(wdat, hazard = haz, .after = "event")
  s.list[[paste("FullSim", nsim, sep = "")]] <- cbind(id, wdat)

  cData <- s.list$FullSim1
  TrueFit <- coxph(Surv(time, event) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age, data = cData)

  TreatRow <- which(names(TrueFit$coefficients) == "treat")
  EcogRow <- which(names(TrueFit$coefficients) == "b.ecogvalue")
  trueTreat <- TrueFit$coefficients[TreatRow]
  trueEcog <- TrueFit$coefficients[EcogRow]

  for (i in 1:iter) {
    sor <- PlasmodeSur(
      objectOut = os1,
      objectCen = oc1,
      idVar = sdat.c$patientid,
      effectOR = 0.5,
      nsim = nsim,
      size = nrow(sdat.c)
    )

    ### Generate missingness in simulated data
    # Remove unneeded variables from primary data
    sdat.m <- sdat.c %>% select(-ndead)
    # initialize empty storage list
    s.list <- list()
    # Pull simulation data
    wdat <- sor[["Sim_Data"]] %>%
      select(paste("ID", nsim, sep = ""), paste("EVENT", nsim, sep = ""), paste("TIME", nsim, sep = ""))

    names(wdat) <- c("patientid", "event", "time")

    # Attach full data to simulation i
    wdat <- dplyr::left_join(wdat, sdat.m, by = "patientid") # %>% select(-patientid)
    wdat$event <- wdat$event * 1
    wdat$dead <- NULL
    wdat$cmonth <- NULL

    haz <- nelsonaalen(wdat, time, event)
    wdat <- add_column(wdat, hazard = haz, .after = "event")
    id <- wdat$patientid
    wdat$patientid <- NULL

    X <- wdat
    X$hazard <- NULL
    X <- as.matrix(X)
    # Note ECOG is second to last column in X
    # Note zeta is last column in X
    miss.coef <- c(rep(log(1.1), 15), log(5), log(5))
    # Computes columnwise prob of missing to reach
    # Pr(missing at least 1) = 10/30/50%
    # 9 variables including ECOG
    nvars <- 9
    multiplier <- 1 - exp(log(1 - proportionList[w] / 100) / nvars)

    # Add MNAR missingness to ECOG
    expit <- function(x) exp(x) / (1 + exp(x))
    logit <- function(x) log(x / (1 - x))
    beta0 <- b0[w]
    p <- expit(beta0 + X %*% miss.coef)
    miss.ind <- rbinom(nrow(X), 1, p)
    if (FALSE) {
      J <- 1000
      mmi <- rep(NA, J)
      for (j in 1:J) {
        miss.ind <- rbinom(nrow(X), 1, p)
        mmi[j] <- mean(miss.ind)
      }
      summary(mmi)
    }

    NAEcog <- wdat$b.ecogvalue
    NAEcog[miss.ind == 1] <- NA

    # 3 race and 3 site dummies are considered as blocks
    # NOTE: made decision to NOT induce missing in treat

    # Add MCAR missingness
    wdatNA <- select(wdat, -event, -hazard, -time, -reth_black, -reth_hisp, -reth_oth, -site_ureter, -site_renal, -site_urethra, -treat, -b.ecogvalue)

    wdatNA <- apply(wdatNA, 2, function(x) {
      x[sample(c(1:nrow(wdatNA)), floor(nrow(wdatNA) * multiplier))] <- NA
      x
    })
    wdatNA <- data.frame(wdatNA)
    wdatNA <- add_column(wdatNA, treat = wdat$treat, .after = "surgery")
    wdatNA <- add_column(wdatNA, event = wdat$event, .before = "genderf")
    wdatNA <- add_column(wdatNA, hazard = wdat$hazard, .after = "event")
    wdatNA <- add_column(wdatNA, time = wdat$time, .after = "hazard")

    # Induce MCAR missingness in race
    # Treat 3 dummies as 1
    raceNA <- select (wdat, reth_black, reth_hisp, reth_oth)
    s = sample(c(1:nrow(raceNA)), floor(nrow(raceNA) * multiplier))
    raceNA <- apply(raceNA, 2, function(x) {
      x[s] <- NA
      x
    })    
    raceNA <- data.frame(raceNA)
    wdatNA <- add_column(wdatNA, reth_black = raceNA$reth_black, .after = "genderf")
    wdatNA <- add_column(wdatNA, reth_hisp = raceNA$reth_hisp, .after = "reth_black")
    wdatNA <- add_column(wdatNA, reth_oth = raceNA$reth_oth, .after = "reth_hisp")

    # Induce MCAR missingness in site
    # Treat 3 dummies as 1
    siteNA <- select (wdat, site_ureter, site_renal, site_urethra)
    s = sample(c(1:nrow(siteNA)), floor(nrow(siteNA) * multiplier))
    siteNA <- apply(siteNA, 2, function(x) {
      x[s] <- NA
      x
    })    
    siteNA <- data.frame(siteNA)
    wdatNA <- add_column(wdatNA, site_ureter = siteNA$site_ureter, .after = "practypec")
    wdatNA <- add_column(wdatNA, site_renal = siteNA$site_renal, .after = "site_ureter")
    wdatNA <- add_column(wdatNA, site_urethra = siteNA$site_urethra, .after = "site_renal")

    wdatNA <- add_column(wdatNA, b.ecogvalue = NAEcog, .after = 'treat')

    nas <- wdatNA
    nas$reth_hisp <- NULL
    nas$reth_oth <- NULL
    nas$site_renal <- NULL
    nas$site_urethra <- NULL
    prop <- mean(apply(nas, 1, function(x) any(is.na(x))))

    s.list[[paste("MissingSim", nsim, sep = "")]] <- cbind(id, wdatNA)
    s.list[[paste("FullSim", nsim, sep = "")]] <- cbind(id, wdat)

    # True effects
    trueTreat <- sor$TrueOutBeta[which(names(sor$TrueOutBeta) == "TREAT")]
    trueEcog <- sor$TrueOutBeta[which(names(sor$TrueOutBeta) == "b.ecogvalue")]

    # MNAR data
    mData <- s.list$MissingSim1
    cData <- s.list$FullSim1

    zetaCol <- which(colnames(mData) == "zeta")
    train <- mData[, -c(1, zetaCol)]
    filename <- paste0("mData", i, ".csv")
    Cfilename <- paste0("cData", i, ".csv")
    trueEffect <- cbind(trueTreat, trueEcog, prop)
    colnames(trueEffect)[3] <- "propMissAtLeast1Cov"
    effname <- paste0("trueEff_propMiss", i, ".csv")

    propName <- proportionList[w]
    system(paste0("mkdir ", file.path(rootdir, "datasets", "mDats", "MNAR", propName)))
    system(paste0("mkdir ", file.path(rootdir, "datasets", "cDats", "MNAR", propName)))
    write.csv(train, file.path(rootdir, "datasets/mDats/MNAR", propName, filename), row.names = F)
    write.csv(cData, file.path(rootdir, "datasets/cDats/MNAR", propName, Cfilename), row.names = F)
    write.csv(trueEffect, file.path(rootdir, "datasets/trueEff/MNAR", propName, effname), row.names = F)
  }
}