library(Plasmode)
library(survival)
library(dplyr)
library(mice)
library(tidyverse)

rootdir <- "/project/flatiron_ucc/programs/kylie/RunMe2"
system (paste0 ('mkdir ', file.path (rootdir, 'datasets')))
system (paste0 ('mkdir ', file.path (rootdir, 'datasets', 'mDats')))
system (paste0 ('mkdir ', file.path (rootdir, 'datasets', 'cDats')))
system (paste0 ('mkdir ', file.path (rootdir, 'datasets', 'trueEff')))
system (paste0 ('mkdir ', file.path (rootdir, 'datasets', 'mDats', 'MCAR')))
system (paste0 ('mkdir ', file.path (rootdir, 'datasets', 'cDats', 'MCAR')))
system (paste0 ('mkdir ', file.path (rootdir, 'datasets', 'trueEff', 'MCAR')))

# These are for row-wise missing covariates
# other than ECOG
proportionList <- c(10, 30, 50)
# These are for column-wise missing in ECOG
ecogList <- c(10, 20, 30)

system (paste0 ('mkdir ', file.path (rootdir, 'datasets', 'trueEff', 'MCAR', proportionList[1])))
system (paste0 ('mkdir ', file.path (rootdir, 'datasets', 'trueEff', 'MCAR', proportionList[2])))
system (paste0 ('mkdir ', file.path (rootdir, 'datasets', 'trueEff', 'MCAR', proportionList[3])))
iter <- 1000
effOR <- 0.5

# set number of Plasmode simulations
nsim <- 1

# w will loop through proportions of missing
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

  # Here we generate a new variables for the outcome model that are nonlinear and involve interactions
  # We will include these terms in the outcome model so it is correctly specified, but imputation models will be misspecified by MICE
  sdat.c$var1 <- (sdat.c$age - 73)^2 + rnorm (nrow (sdat.c), sd=100)
  sdat.c$var2 <- (sdat.c$age - 73)*sdat.c$smokey - (sdat.c$age - 73)*sdat.c$surgery - 5*(sdat.c$age - 73)*sdat.c$smokey*sdat.c$surgery + rnorm (nrow (sdat.c), sd=10)

  if (FALSE){
    cor (cbind (sdat.c$age, sdat.c$smokey, sdat.c$surgery, sdat.c$var1, sdat.c$var2), method='spearman')

    age2 = sdat.c$age^2
    a = lm(sdat.c$var1~sdat.c$age + age2)
    b = lm(var1~age,data=sdat.c)
    summary(a)
    summary(b)
    sum ((a$fitted.values - sdat.c$var1)^2)
    sum ((b$fitted.values - sdat.c$var1)^2)

    a = lm(var2~age*smokey*surgery,data=sdat.c)
    b = lm(var2~age+smokey+surgery,data=sdat.c)
    summary(a)
    summary(b)
    sum ((a$fitted.values - sdat.c$var2)^2)
    sum ((b$fitted.values - sdat.c$var2)^2)
  }

  # Step 2: Estimate associations with outcome and censoring
  # One model for the hazard of the outcome event
  # One model for the hazard of censoring (this is simply the reverse of the
  # outcome, aka 1-Y)

  # Note: the Plasmode package actually selects the second independent variable in
  # the formula for the exposure, rather than the first if you are having the
  # package run the model. If you run it yourself it gets it right.

  # Build outcome hazard with newVar
  os1 <- coxph(Surv(sdat.c$cmonth, sdat.c$dead) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + var1 + var2, data = sdat.c, x = TRUE)

  # Censoring hazard with newVar
  oc1 <- coxph(Surv(sdat.c$cmonth, sdat.c$ndead) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + var1 + var2, data = sdat.c, x = TRUE)

  genDat <- function (k){
    # Increase effect of vars on outcome using MMOut
    outEff <- rep(1, length(coef(os1)))
    lc = length (outEff)
    outEff[lc-1] <- 40
    outEff[lc] <- 20
    sor <- PlasmodeSur(
      objectOut = os1,
      objectCen = oc1,
      idVar = sdat.c$patientid,
      effectOR = effOR,
      MMOut = outEff,
      nsim = nsim,
      size = nrow(sdat.c)
    )

    if (FALSE) {
      # Check that associations reflect real
      # associations and altered associations
      # for treat and new variables
      ps = matrix (NA, nrow=50, ncol=2)
      for (l in 1:50) {
        print (l)
        sor <- PlasmodeSur(
          objectOut = os1,
          objectCen = oc1,
          idVar = sdat.c$patientid,
          effectOR = effOR,
          MMOut = outEff,
          nsim = nsim,
          size = nrow(sdat.c)
        )

        colnames(sor$Sim_Data)[1] <- "patientid"
        testdf <- dplyr::left_join(sor$Sim_Data, sdat.c, by = "patientid")
        testdf$dead <- sor$Sim_Data$EVENT1 * 1
        testdf$cmonth <- sor$Sim_Data$TIME1
        os1test <- coxph(Surv(testdf$cmonth, testdf$dead) ~ treat + genderf + reth_black + reth_hisp + reth_oth + practypec + b.ecogvalue + smokey + dgradeh + surgery + site_ureter + site_renal + site_urethra + age + var1 + var2, data = testdf, x = TRUE)
        ps[l,] = summary(os1test)$coefficients[(lc-1):lc,5]
      }
      apply (ps, 2, function (x) mean (x<=0.05))
    }

    ### Generate missingness in simulated data
    # Remove unneeded variables from primary data
    sdat.m <- sdat.c %>% select(-ndead)

    # Attach full data to simulation i
    colnames(sor$Sim_Data)[1] = "patientid"
    wdat <- dplyr::left_join(sor$Sim_Data, sdat.m, by = "patientid") 
    wdat$event <- wdat$EVENT1 * 1
    wdat$EVENT1 = NULL
    wdat$time = wdat$TIME1
    wdat$TIME1 = NULL
    wdat$dead <- NULL
    wdat$cmonth <- NULL

    haz <- nelsonaalen(wdat, time, event)
    wdat <- add_column(wdat, hazard = haz, .after = "event")
    id <- wdat$patientid
    wdat$patientid <- NULL

    # Computes columnwise prob of missing to reach
    # Pr(missing at least 1) = 10/30/50%
    # No missingness in txt
    # ECOG treated separately, columnwise
    # Subtract 4 for race and site categories
    # Subtract 3 to exclude event, time, hazard
    # Subtract 2 for treat and ECOG
    nvars = ncol (wdat) - 4 - 3 - 2
    multiplier <- 1 - exp(log(1 - proportionList[w] / 100) / nvars)

    # Add MCAR missingness to ECOG
    miss.ind <- rbinom(wdat, 1, ecogList[w]/100)
    NAEcog <- wdat$b.ecogvalue
    NAEcog[miss.ind == 1] <- NA

    # Add MCAR missingness to rest
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
    raceNA <- select(wdat, reth_black, reth_hisp, reth_oth)
    s <- sample(c(1:nrow(raceNA)), floor(nrow(raceNA) * multiplier))
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
    siteNA <- select(wdat, site_ureter, site_renal, site_urethra)
    s <- sample(c(1:nrow(siteNA)), floor(nrow(siteNA) * multiplier))
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
    # Overall row-wise missing INCLUDING
    # ECOG and newVar
    prop <- mean(apply(nas, 1, function(x) any(is.na(x))))

    s.list = list ()
    s.list[[paste("MissingSim", nsim, sep = "")]] <- cbind(id, wdatNA)
    s.list[[paste("FullSim", nsim, sep = "")]] <- cbind(id, wdat)
    s.list[[paste("trueBeta", nsim, sep = "")]] <- sor$TrueOutBeta
    s.list[[paste("propMiss", nsim, sep = "")]] <- prop

    return (s.list)
  }

  for (i in 1:iter) {
    keep = FALSE
    while (!keep){
      count = i
      print (count)
      op = genDat(count)
      # MCAR data
      mData <- op$MissingSim1
      cData <- op$FullSim1
      # True effects
      trueBeta <- op$trueBeta1
      # True rowwise missing rate
      prop <- op$propMiss1
      keep = !any (apply (mData, 1, function (x) all (is.na (x))))
      count = count + 1000
    }
    train <- mData[, -1]
    filename <- paste0("mData", i, ".csv")
    Cfilename <- paste0("cData", i, ".csv")
    trueEffect <- c(prop, trueBeta)
    names(trueEffect)[1] <- "propMissAtLeast1Cov"
    effname <- paste0("propMiss_trueEffs", i, ".csv")

    propName <- proportionList[w]
    system(paste0("mkdir ", file.path(rootdir, "datasets", "mDats", "MCAR", propName)))
    system(paste0("mkdir ", file.path(rootdir, "datasets", "cDats", "MCAR", propName)))
    write.csv(train, file.path(rootdir, "datasets/mDats/MCAR", propName, filename), row.names = F)
    write.csv(cData, file.path(rootdir, "datasets/cDats/MCAR", propName, Cfilename), row.names = F)
    write.csv(t (trueEffect), file.path(rootdir, "datasets/trueEff/MCAR", propName, effname), row.names = F)
  }
}