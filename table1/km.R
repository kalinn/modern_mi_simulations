library(Plasmode)
library(dplyr)
library(Survival)

rootdir <- "/project/flatiron_ucc/programs/kylie/RunMe"
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


completeCaseVersion = TRUE
if (completeCaseVersion){
    sdat.c <- sdat %>% na.omit()
} else{
    sdat.c <- sdat
}

f1 <- survfit(Surv(cmonth, dead) ~ 1, data = sdat.c)

