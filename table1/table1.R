library(Plasmode)
library(dplyr)
library(kableExtra)

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
# Round the few half scores up to next integer
sdat.c$b.ecogvalue <- ceiling(sdat.c$b.ecogvalue)

# Reconstruct categorical variables
sdat.c$genderf = factor (sdat.c$genderf)
levels (sdat.c$genderf) = c ('Male', 'Female')

sdat.c$reth_hisp[sdat.c$reth_hisp==1] = 2
sdat.c$reth_oth[sdat.c$reth_oth==1] = 3

sdat.c$race.ethnicity = sdat.c$reth_black + sdat.c$reth_hisp + sdat.c$reth_oth
sdat.c$race.ethnicity = factor (sdat.c$race.ethnicity)
levels (sdat.c$race.ethnicity) = c ('White', 'Black', 'Hispanic or Latino', 'Other')
sdat.c$reth_black = NULL
sdat.c$reth_hisp = NULL
sdat.c$reth_oth = NULL

sdat.c$practypec = factor (sdat.c$practypec)
levels (sdat.c$practypec) = c ('Not Community', 'Community')

sdat.c$site_renal[sdat.c$site_renal==1] = 2
sdat.c$site_urethra[sdat.c$site_urethra==1] = 3
sdat.c$primarysite = sdat.c$site_ureter + sdat.c$site_renal + sdat.c$site_urethra
sdat.c$primarysite = factor (sdat.c$primarysite)
levels (sdat.c$primarysite) = c ('Bladder', 'Ureter', 'Renal', 'Urethra')
sdat.c$site_ureter = NULL
sdat.c$site_renal = NULL
sdat.c$site_urethra = NULL

sdat.c$smokey = factor (sdat.c$smokey)
levels (sdat.c$smokey) = c ('No history of smoking', 'History of smoking')

sdat.c$dgradeh = factor (sdat.c$dgradeh)
levels (sdat.c$dgradeh) = c ('Low grade', 'High grade (G2/G3/G4)')

sdat.c$surgery = factor (sdat.c$surgery)
levels (sdat.c$surgery) = c ('No surgery', 'Surgery')

sdat.c$treat = factor (sdat.c$treat)
levels (sdat.c$treat) = c ('Chemotherapy', 'Immunotherapy')


n = group_by(sdat.c, treat) %>% summarize (., 'N'=n())
N = as.integer (n$N)
#name1 = paste0 ('Chemotherapy (N=', N[1], ')')
#name2 = paste0 ('Immunotherapy (N=', N[2], ')')
name1 = paste0 ('Chemotherapy')
name2 = paste0 ('Immunotherapy')

age = group_by(sdat.c, treat) %>% summarize (., 'Age'=mean (age))
ageSD = group_by(sdat.c, treat) %>% summarize (., 'Age'=sd (age))

#ecog = group_by(sdat.c, treat) %>% summarize (., 'ECOG PS'=mean (b.ecogvalue))
#ecogSD = group_by(sdat.c, treat) %>% summarize (., 'ECOG PS'=sd (b.ecogvalue))

ageChar = mapply (function (x, y){
    mu = round (x, 0)
    sig = round (y, 1)
    paste0 (mu, ' (', sig, ')')
}, age$Age, ageSD$Age)

#ecogChar = mapply (function (x, y){
#    mu = round (x, 2)
#    sig = round (y, 1)
#    paste0 (mu, ' (', sig, ')')
#}, ecog$`ECOG PS`, ecogSD$`ECOG PS`)

#contDf = rbind (ageChar, ecogChar)
contDf = as.data.frame (t (ageChar))
colnames (contDf) = c (name1, name2)
rownames (contDf) = c ('Age (mean, SD)')

ecogCat = group_by (sdat.c, treat) %>% summarise(., 'Missing'=100*mean (is.na (b.ecogvalue)), 'ECOG0'=100*sum (b.ecogvalue==0, na.rm=TRUE)/n(), 'ECOG1'=100*sum (b.ecogvalue==1, na.rm=TRUE)/n(), 'ECOG2'=100*sum (b.ecogvalue==2, na.rm=TRUE)/n(), 'ECOG34'=100*sum (b.ecogvalue>2, na.rm=TRUE)/n(), 'NMiss'=sum (is.na (b.ecogvalue)), 'N0'=sum (b.ecogvalue==0, na.rm=TRUE), 'N1'=sum (b.ecogvalue==1, na.rm=TRUE), 'N2'=sum (b.ecogvalue==2, na.rm=TRUE), 'N34'=sum (b.ecogvalue>2, na.rm=TRUE))

ecogCat$Missing = round (ecogCat$Missing, 1)
ecogCat$ECOG0 = round (ecogCat$ECOG0, 1)
ecogCat$ECOG1 = round (ecogCat$ECOG1, 1)
ecogCat$ECOG2 = round (ecogCat$ECOG2, 1)
ecogCat$ECOG34 = round (ecogCat$ECOG34, 1)

pctEcog = data.frame('Variable'=rep ("ECOG PS", 5),'subgroup'=c('Missing', 'ECOG PS = 0', 'ECOG PS = 1', 'ECOG PS = 2', 'ECOG PS > 2'), 'Chemo'=c (paste0 (ecogCat$NMiss[1], ' (', ecogCat$Missing[1], ')'), paste0 (ecogCat$N0[1], ' (', ecogCat$ECOG0[1], ')'), paste0 (ecogCat$N1[1], ' (', ecogCat$ECOG1[1], ')'), paste0 (ecogCat$N2[1], ' (', ecogCat$ECOG2[1], ')'), paste0 (ecogCat$N34[1], ' (', ecogCat$ECOG34[1], ')')), 'Immuno'=c (paste0 (ecogCat$NMiss[2], ' (', ecogCat$Missing[2], ')'), paste0 (ecogCat$N0[2], ' (', ecogCat$ECOG0[2], ')'), paste0 (ecogCat$N1[2], ' (', ecogCat$ECOG1[2], ')'), paste0 (ecogCat$N2[2], ' (', ecogCat$ECOG2[2], ')'), paste0 (ecogCat$N34[2], ' (', ecogCat$ECOG34[2], ')')))

if (completeCaseVersion){
    pctEcog = pctEcog[-1,]
}

sex = group_by(sdat.c, treat) %>% summarize (., 'Missing'=100*mean (is.na (genderf)), 'Female'=100*sum (genderf=='Female', na.rm=TRUE)/n(), 'Male'=100*sum (genderf=='Male', na.rm=TRUE)/n(), 'NMiss'=sum (is.na (genderf)), 'N0'=sum (genderf=='Female', na.rm=TRUE), 'N1'=sum (genderf=='Male', na.rm=TRUE))
sex$Missing = round (sex$Missing, 1)
sex$Female = round (sex$Female, 1)
sex$Male = round (sex$Male, 1)

pctSex = data.frame('Variable'=rep ("Sex", 3),'subgroup'=c('Missing', 'Female', 'Male'), 'Chemo'=c (paste0 (sex$NMiss[1], ' (', sex$Missing[1], ')'), paste0 (sex$N0[1], ' (', sex$Female[1], ')'), paste0 (sex$N1[1], ' (', sex$Male[1], ')')), 'Immuno'=c (paste0 (sex$NMiss[2], ' (', sex$Missing[2], ')'), paste0 (sex$N0[2], ' (', sex$Female[2], ')'), paste0 (sex$N1[2], ' (', sex$Male[2], ')')))

if (completeCaseVersion){
    pctSex = pctSex[-1,]
}

practypec = group_by(sdat.c, treat) %>% summarize (., 'Community'=100*mean (practypec=='Community'), 'Non'=100*mean (practypec!='Community'), 'N0'=sum (practypec=='Community'), 'N1'=sum (practypec!='Community'))
practypec$Community = round (practypec$Community, 1)
practypec$Non = round (practypec$Non, 1)

pctPrac = data.frame('Variable'=rep ("Practice Type", 2),'subgroup'=c('Community', 'Non-community'), 'Chemo'=c (paste0 (practypec$N0[1], ' (', practypec$Community[1], ')'), paste0 (practypec$N1[1], ' (', practypec$Non[1], ')')), 'Immuno'=c (paste0 (practypec$N0[2], ' (', practypec$Community[2], ')'), paste0 (practypec$N1[2], ' (', practypec$Non[2], ')')))


smokey = group_by(sdat.c, treat) %>% summarize (., 'Smoking'=100*mean (smokey=='History of smoking'), 'Non'=100*mean (smokey!='History of smoking'), 'N0'=sum (smokey=='History of smoking'), 'N1'=sum (smokey!='History of smoking'))
smokey$Smoking = round (smokey$Smoking, 1)
smokey$Non = round (smokey$Non, 1)

pctSmoke = data.frame('Variable'=rep ("Smoking Status", 2), 'subgroup'=c('History of smoking', 'No history of smoking'), 'Chemo'=c (paste0 (smokey$N0[1], ' (', smokey$Smoking[1], ')'), paste0 (smokey$N1[1], ' (', smokey$Non[1], ')')), 'Immuno'=c (paste0 (smokey$N0[2], ' (', smokey$Smoking[2], ')'), paste0 (smokey$N1[2], ' (', smokey$Non[2], ')')))


grade = group_by(sdat.c, treat) %>% summarize (., 'High'=100*mean (dgradeh=='High grade (G2/G3/G4)'), 'Low'=100*mean (dgradeh!='High grade (G2/G3/G4)'), 'N0'=sum (dgradeh=='High grade (G2/G3/G4)'), 'N1'=sum (dgradeh!='High grade (G2/G3/G4)'))
grade$High = round (grade$High, 1)
grade$Low = round (grade$Low, 1)

pctGrade = data.frame('Variable'=rep ("Grade", 2),'subgroup'=c('High (G2/G3/G4)', 'Low '), 'Chemo'=c (paste0 (grade$N0[1], ' (', grade$High[1], ')'), paste0 (grade$N1[1], ' (', grade$Low[1], ')')), 'Immuno'=c (paste0 (grade$N0[2], ' (', grade$High[2], ')'), paste0 (grade$N1[2], ' (', grade$Low[2], ')')))


surgery = group_by(sdat.c, treat) %>% summarize (., 'Surgery'=100*mean (surgery=='Surgery'), 'Non'=100*mean (surgery!='Surgery'), 'N0'=sum(surgery=='Surgery'), 'N1'=sum (surgery!='Surgery'))
surgery$Surgery = round (surgery$Surgery, 1)
surgery$Non = round (surgery$Non, 1)

pctSurg = data.frame('Variable'=rep ("Surgery", 2),'subgroup'=c('Surgery', 'No Surgery'), 'Chemo'=c (paste0 (surgery$N0[1], ' (', surgery$Surgery[1], ')'), paste0 (surgery$N1[1], ' (', surgery$Non[1], ')')), 'Immuno'=c (paste0 (surgery$N0[2], ' (', surgery$Surgery[2], ')'), paste0 (surgery$N1[2], ' (', surgery$Non[2], ')')))

race = group_by(sdat.c, treat) %>% summarize (., 'White'=100*mean (race.ethnicity=='White'), 'Black'=100*mean (race.ethnicity=='Black'), 'Hisp'=100*mean (race.ethnicity=='Hispanic or Latino'), 'Other'=100*mean (race.ethnicity=='Other'), 'N0'=sum (race.ethnicity=='White'), 'N1'=sum (race.ethnicity=='Black'), 'N2'=sum (race.ethnicity=='Hispanic or Latino'), 'N3'=sum (race.ethnicity=='Other'))
race$White = round (race$White, 1)
race$Black = round (race$Black, 1)
race$Hisp = round (race$Hisp, 1)
race$Other = round (race$Other, 1)

pctRace = data.frame('Variable'=rep ("Race", 4),'subgroup'=c('White', 'Black', 'Hispanic or Latino', 'Other'), 'Chemo'=c (paste0(race$N0[1], ' (', race$White[1], ')'), paste0(race$N1[1], ' (', race$Black[1], ')'), paste0(race$N2[1], ' (', race$Hisp[1], ')'), paste0(race$N3[1], ' (', race$Other[1], ')')), 'Immuno'=c (paste0(race$N0[2], ' (', race$White[2], ')'), paste0(race$N1[2], ' (', race$Black[2], ')'), paste0(race$N2[2], ' (', race$Hisp[2], ')'), paste0(race$N3[2], ' (', race$Other[2], ')')))

site = group_by(sdat.c, treat) %>% summarize (., 'Bladder'=100*mean (primarysite=='Bladder'), 'Ureter'=100*mean (primarysite=='Ureter'), 'Renal'=100*mean (primarysite=='Renal'), 'Urethra'=100*mean (primarysite=='Urethra'), 'N0'=sum (primarysite=='Bladder'), 'N1'=sum (primarysite=='Ureter'), 'N2'=sum (primarysite=='Renal'), 'N3'=sum (primarysite=='Urethra'))
site$Bladder = round (site$Bladder, 1)
site$Ureter = round (site$Ureter, 1)
site$Renal = round (site$Renal, 1)
site$Urethra = round (site$Urethra, 1)

pctSite = data.frame('Variable'=rep ("Primary Site", 4),'subgroup'=c('Bladder', 'Ureter', 'Renal', 'Urethra'), 'Chemo'=c (paste0(site$N0[1], ' (', site$Bladder[1], ')'), paste0(site$N1[1], ' (', site$Ureter[1], ')'), paste0(site$N2[1], ' (', site$Renal[1], ')'), paste0(site$N3[1], ' (', site$Urethra[1], ')')), 'Immuno'=c (paste0(site$N0[2], ' (', site$Bladder[2], ')'), paste0(site$N1[2], ' (', site$Ureter[2], ')'), paste0(site$N2[2], ' (', site$Renal[2], ')'), paste0(site$N3[2], ' (', site$Urethra[2], ')')))


binDf = rbind (pctSex, pctRace, pctSite, pctPrac, pctSmoke, pctGrade, pctSurg, pctEcog)
binDf$Chemo = as.character(binDf$Chemo)
binDf$Immuno = as.character(binDf$Immuno)
colnames (binDf)[3:4] = c (name1, name2)

#dfCont = data.frame ('Variable'=c ('N', 'Age', 'ECOG'), 'subgroup'=c ('N', 'Age (mean, SD)', 'ECOG PS (mean, SD)'), 'Chemo'=c (N[1], contDf[1,1], contDf[2,1]), 'Immuno'=c (N[2], contDf[2,1], contDf[2,2]))
dfCont = data.frame ('Variable'=c ('N', 'Age'), 'subgroup'=c ('N', 'Age (mean, SD)'), 'Chemo'=c (N[1], contDf[1,1]), 'Immuno'=c (N[2], contDf[1,2]))
colnames (dfCont)[2:4] = c ('subgroup', name1, name2)

allDf = rbind (dfCont, binDf)
vars = allDf$Variable
uvars = unique (vars)
allDf$Variable = NULL
colnames (allDf)[1] = " "

if (!completeCaseVersion){
kable(allDf, row.names=FALSE, align=c("l", "c", "c"), format = "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"), full_width = F) %>%
  pack_rows("Sex N (%)", 3, 5) %>% 
  pack_rows("Race N (%)", 6, 9) %>% 
  pack_rows("Primary Site N (%)", 10, 13) %>% 
  pack_rows("Practice Type N (%)", 14, 15) %>% 
  pack_rows("Smoking Status N (%)", 16, 17) %>% 
  pack_rows("Grade N (%)", 18, 19) %>% 
  pack_rows("Surgery N (%)", 20, 21) %>%
  pack_rows("ECOG N (%)", 22, 26)
} else{
kable(allDf, row.names=FALSE, align=c("l", "c", "c"), format = "latex", booktabs = T) %>%
  kable_styling(latex_options = c("striped", "hold_position"), full_width = F) %>%
  pack_rows("Sex N (%)", 3, 5) %>% 
  pack_rows("Race N (%)", 5, 8) %>% 
  pack_rows("Primary Site N (%)", 9, 12) %>% 
  pack_rows("Practice Type N (%)", 13, 14) %>% 
  pack_rows("Smoking Status N (%)", 15, 16) %>% 
  pack_rows("Grade N (%)", 17, 18) %>% 
  pack_rows("Surgery N (%)", 19, 20) %>%
  pack_rows("ECOG N (%)", 21, 24)

}