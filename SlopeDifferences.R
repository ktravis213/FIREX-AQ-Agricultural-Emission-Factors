require(lmodel2)
# ------------
cContent = 500 # Assuming 50 % C, 500 g/kg
cContentAG = 415.1 
cContentGrassShrub = 462.7
cContentPresc = 511.1

# ---- USE varsALL from MCEfig.R -------
runMCE = 0
if (runMCE == 1){source("MCEFig.R") }

# ------ Correct all species to 0.92 MCE -------
allBOTH.filter$FinalEF_MCE92 = NaN
ind = which(is.nan(varsALL$Rval_ag)); varsALL$Rval_ag[ind] = 0
ind = which(is.nan(varsALL$Rval_presc)); varsALL$Rval_presc[ind] = 0
ind = which(is.nan(varsALL$Rval_grass)); varsALL$Rval_grass[ind] = 0

# ----- Choose best MCE dependence ---------
# USE ISAF for Formaldehyde
indA = which(varsALL$vars == 'CH2O_ISAF_HANISCO')
ind = which(varsALL$names == 'Formaldehyde')
varsALL$Slope_ag[ind] = varsALL$Slope_ag[indA];varsALL$Slope_presc[ind] = varsALL$Slope_presc[indA]; varsALL$Slope_grass[ind] = varsALL$Slope_grass[indA]
varsALL$Intercept_ag[ind] = varsALL$Intercept_ag[indA];varsALL$Intercept_presc[ind] = varsALL$Intercept_presc[indA]; varsALL$Intercept_grass[ind] = varsALL$Intercept_grass[indA]
varsALL$Rval_ag[ind] = varsALL$Rval_ag[indA];varsALL$Rval_presc[ind] = varsALL$Rval_presc[indA]; varsALL$Rval_grass[ind] = varsALL$Rval_grass[indA]

# Using PTRMS MVK/MACR/2butenals
indA = which(varsALL$vars == 'MVKMAC_NOAAPTR_ppbv_WARNEKE')
ind = which(varsALL$names == 'Methyl vinyl ketone' | varsALL$names == '2-Butenals' | varsALL$names == 'Methacrolein')
varsALL$Slope_ag[ind] = varsALL$Slope_ag[indA];varsALL$Slope_presc[ind] = varsALL$Slope_presc[indA]; varsALL$Slope_grass[ind] = varsALL$Slope_grass[indA]
varsALL$Intercept_ag[ind] = varsALL$Intercept_ag[indA];varsALL$Intercept_presc[ind] = varsALL$Intercept_presc[indA]; varsALL$Intercept_grass[ind] = varsALL$Intercept_grass[indA]
varsALL$Rval_ag[ind] = varsALL$Rval_ag[indA];varsALL$Rval_presc[ind] = varsALL$Rval_presc[indA]; varsALL$Rval_grass[ind] = varsALL$Rval_grass[indA]

# Using PTRMS C4Carbonyls
indA = which(varsALL$vars == 'C4Carbonyls_NOAAPTR_ppbv_WARNEKE')
ind = which(varsALL$names == 'Methyl Ethyl Ketone' | varsALL$names == 'Butanal' | varsALL$names == 'Isobutanal')
varsALL$Slope_ag[ind] = varsALL$Slope_ag[indA];varsALL$Slope_presc[ind] = varsALL$Slope_presc[indA]; varsALL$Slope_grass[ind] = varsALL$Slope_grass[indA]
varsALL$Intercept_ag[ind] = varsALL$Intercept_ag[indA];varsALL$Intercept_presc[ind] = varsALL$Intercept_presc[indA]; varsALL$Intercept_grass[ind] = varsALL$Intercept_grass[indA]
varsALL$Rval_ag[ind] = varsALL$Rval_ag[indA];varsALL$Rval_presc[ind] = varsALL$Rval_presc[indA]; varsALL$Rval_grass[ind] = varsALL$Rval_grass[indA]

# Using iWAS
indA = which(varsALL$vars == 'x2Me1Butene_NOAAiWAS_GILMAN')
ind = which(varsALL$names == '2-Methyl-1-butene')
varsALL$Slope_ag[ind] = varsALL$Slope_ag[indA];varsALL$Slope_presc[ind] = varsALL$Slope_presc[indA]; varsALL$Slope_grass[ind] = varsALL$Slope_grass[indA]
varsALL$Intercept_ag[ind] = varsALL$Intercept_ag[indA];varsALL$Intercept_presc[ind] = varsALL$Intercept_presc[indA]; varsALL$Intercept_grass[ind] = varsALL$Intercept_grass[indA]
varsALL$Rval_ag[ind] = varsALL$Rval_ag[indA];varsALL$Rval_presc[ind] = varsALL$Rval_presc[indA]; varsALL$Rval_grass[ind] = varsALL$Rval_grass[indA]

# Using iWAS
indA = which(varsALL$vars == 'MeFormate_NOAAiWAS_GILMAN')
ind = which(varsALL$names == 'Methyl formate')
varsALL$Slope_ag[ind] = varsALL$Slope_ag[indA];varsALL$Slope_presc[ind] = varsALL$Slope_presc[indA]; varsALL$Slope_grass[ind] = varsALL$Slope_grass[indA]
varsALL$Intercept_ag[ind] = varsALL$Intercept_ag[indA];varsALL$Intercept_presc[ind] = varsALL$Intercept_presc[indA]; varsALL$Intercept_grass[ind] = varsALL$Intercept_grass[indA]
varsALL$Rval_ag[ind] = varsALL$Rval_ag[indA];varsALL$Rval_presc[ind] = varsALL$Rval_presc[indA]; varsALL$Rval_grass[ind] = varsALL$Rval_grass[indA]

# Using iWAS
indA = which(varsALL$vars == 'x3Me1Butene_NOAAiWAS_GILMAN')
ind = which(varsALL$names == '3-Methyl-1-butene')
varsALL$Slope_ag[ind] = varsALL$Slope_ag[indA];varsALL$Slope_presc[ind] = varsALL$Slope_presc[indA]; varsALL$Slope_grass[ind] = varsALL$Slope_grass[indA]
varsALL$Intercept_ag[ind] = varsALL$Intercept_ag[indA];varsALL$Intercept_presc[ind] = varsALL$Intercept_presc[indA]; varsALL$Intercept_grass[ind] = varsALL$Intercept_grass[indA]
varsALL$Rval_ag[ind] = varsALL$Rval_ag[indA];varsALL$Rval_presc[ind] = varsALL$Rval_presc[indA]; varsALL$Rval_grass[ind] = varsALL$Rval_grass[indA]

# Using WAS
indA = which(varsALL$vars == 'MeCycHexane_WAS_BLAKE')
ind = which(varsALL$names == 'Methylcyclohexane')
varsALL$Slope_ag[ind] = varsALL$Slope_ag[indA];varsALL$Slope_presc[ind] = varsALL$Slope_presc[indA]; varsALL$Slope_grass[ind] = varsALL$Slope_grass[indA]
varsALL$Intercept_ag[ind] = varsALL$Intercept_ag[indA];varsALL$Intercept_presc[ind] = varsALL$Intercept_presc[indA]; varsALL$Intercept_grass[ind] = varsALL$Intercept_grass[indA]
varsALL$Rval_ag[ind] = varsALL$Rval_ag[indA];varsALL$Rval_presc[ind] = varsALL$Rval_presc[indA]; varsALL$Rval_grass[ind] = varsALL$Rval_grass[indA]

# Using PTRMS C8 aromatics
indA = which(varsALL$vars == 'C8Aromatics_NOAAPTR_ppbv_WARNEKE')
ind = which(varsALL$names == 'Ethylbenzene' | varsALL$names == 'm,p-Xylene' | varsALL$names == 'o-Xylene')
varsALL$Slope_ag[ind] = varsALL$Slope_ag[indA];varsALL$Slope_presc[ind] = varsALL$Slope_presc[indA]; varsALL$Slope_grass[ind] = varsALL$Slope_grass[indA]
varsALL$Intercept_ag[ind] = varsALL$Intercept_ag[indA];varsALL$Intercept_presc[ind] = varsALL$Intercept_presc[indA]; varsALL$Intercept_grass[ind] = varsALL$Intercept_grass[indA]
varsALL$Rval_ag[ind] = varsALL$Rval_ag[indA];varsALL$Rval_presc[ind] = varsALL$Rval_presc[indA]; varsALL$Rval_grass[ind] = varsALL$Rval_grass[indA]

# Using PTRMS C9 aromatics
indA = which(varsALL$vars == 'C9Aromatics_NOAAPTR_ppbv_WARNEKE')
ind = which(varsALL$names == '1,2,4-Trimethylbenzene' | varsALL$names == '3-Ethyltoluene' | varsALL$names == '4-Ethyltoluene' |
              varsALL$names == 'i-Propylbenzene' | varsALL$names == 'n-Propylbenzene' | varsALL$names == '1,3,5-trimethylbenzene' |
              varsALL$names == '2-Ethyltoluene')
varsALL$Slope_ag[ind] = varsALL$Slope_ag[indA];varsALL$Slope_presc[ind] = varsALL$Slope_presc[indA]; varsALL$Slope_grass[ind] = varsALL$Slope_grass[indA]
varsALL$Intercept_ag[ind] = varsALL$Intercept_ag[indA];varsALL$Intercept_presc[ind] = varsALL$Intercept_presc[indA]; varsALL$Intercept_grass[ind] = varsALL$Intercept_grass[indA]
varsALL$Rval_ag[ind] = varsALL$Rval_ag[indA];varsALL$Rval_presc[ind] = varsALL$Rval_presc[indA]; varsALL$Rval_grass[ind] = varsALL$Rval_grass[indA]

# Using PTRMS Monoterpenes
indA = which(varsALL$vars == 'Monoterpenes_NOAAPTR_ppbv_WARNEKE')
ind = which(varsALL$names == 'beta-Pinene' | varsALL$names == 'Camphene'| varsALL$names == 'alpha-Pinene'
            | varsALL$names == 'beta-Pinene/Myrcene' | varsALL$names == 'Myrcene' | varsALL$names == 'Tricyclene')
varsALL$Slope_ag[ind] = varsALL$Slope_ag[indA];varsALL$Slope_presc[ind] = varsALL$Slope_presc[indA]; varsALL$Slope_grass[ind] = varsALL$Slope_grass[indA]
varsALL$Intercept_ag[ind] = varsALL$Intercept_ag[indA];varsALL$Intercept_presc[ind] = varsALL$Intercept_presc[indA]; varsALL$Intercept_grass[ind] = varsALL$Intercept_grass[indA]
varsALL$Rval_ag[ind] = varsALL$Rval_ag[indA];varsALL$Rval_presc[ind] = varsALL$Rval_presc[indA]; varsALL$Rval_grass[ind] = varsALL$Rval_grass[indA]

# Using CAMS Ethane
indA = which(varsALL$vars == 'C2H6_CAMS_pptv_FRIED' )
ind = which(varsALL$names == 'Ethane')
varsALL$Slope_ag[ind] = varsALL$Slope_ag[indA];varsALL$Slope_presc[ind] = varsALL$Slope_presc[indA]; varsALL$Slope_grass[ind] = varsALL$Slope_grass[indA]
varsALL$Intercept_ag[ind] = varsALL$Intercept_ag[indA];varsALL$Intercept_presc[ind] = varsALL$Intercept_presc[indA]; varsALL$Intercept_grass[ind] = varsALL$Intercept_grass[indA]
varsALL$Rval_ag[ind] = varsALL$Rval_ag[indA];varsALL$Rval_presc[ind] = varsALL$Rval_presc[indA]; varsALL$Rval_grass[ind] = varsALL$Rval_grass[indA]

# Using CIT-CIMS HCN
indA = which(varsALL$vars == 'HCN_WENNBERG')
ind = which(varsALL$names == 'Hydrogen cyanide')
varsALL$Slope_ag[ind] = varsALL$Slope_ag[indA];varsALL$Slope_presc[ind] = varsALL$Slope_presc[indA]; varsALL$Slope_grass[ind] = varsALL$Slope_grass[indA]
varsALL$Intercept_ag[ind] = varsALL$Intercept_ag[indA];varsALL$Intercept_presc[ind] = varsALL$Intercept_presc[indA]; varsALL$Intercept_grass[ind] = varsALL$Intercept_grass[indA]
varsALL$Rval_ag[ind] = varsALL$Rval_ag[indA];varsALL$Rval_presc[ind] = varsALL$Rval_presc[indA]; varsALL$Rval_grass[ind] = varsALL$Rval_grass[indA]

#Using HNCO from PTRMS
indA = which(varsALL$vars == 'HNCO_NOAAPTR_ppbv_WARNEKE')
ind = which(varsALL$names == 'Isocyanic acid')
varsALL$Slope_ag[ind] = varsALL$Slope_ag[indA];varsALL$Slope_presc[ind] = varsALL$Slope_presc[indA]; varsALL$Slope_grass[ind] = varsALL$Slope_grass[indA]
varsALL$Intercept_ag[ind] = varsALL$Intercept_ag[indA];varsALL$Intercept_presc[ind] = varsALL$Intercept_presc[indA]; varsALL$Intercept_grass[ind] = varsALL$Intercept_grass[indA]
varsALL$Rval_ag[ind] = varsALL$Rval_ag[indA];varsALL$Rval_presc[ind] = varsALL$Rval_presc[indA]; varsALL$Rval_grass[ind] = varsALL$Rval_grass[indA]

#Using HONO from ACES
indA = which(varsALL$vars == 'HNO2_ACES_WOMACK')
ind = which(varsALL$names == 'Nitrous acid')
varsALL$Slope_ag[ind] = varsALL$Slope_ag[indA];varsALL$Slope_presc[ind] = varsALL$Slope_presc[indA]; varsALL$Slope_grass[ind] = varsALL$Slope_grass[indA]
varsALL$Intercept_ag[ind] = varsALL$Intercept_ag[indA];varsALL$Intercept_presc[ind] = varsALL$Intercept_presc[indA]; varsALL$Intercept_grass[ind] = varsALL$Intercept_grass[indA]
varsALL$Rval_ag[ind] = varsALL$Rval_ag[indA];varsALL$Rval_presc[ind] = varsALL$Rval_presc[indA]; varsALL$Rval_grass[ind] = varsALL$Rval_grass[indA]

# Using PTRMS AcetonePropanal_NOAAPTR_ppbv_WARNEKE
indA = which(varsALL$vars == 'AcetonePropanal_NOAAPTR_ppbv_WARNEKE' )
ind = which(varsALL$names == 'Acetone' | varsALL$names == 'Propanal')
varsALL$Slope_ag[ind] = varsALL$Slope_ag[indA];varsALL$Slope_presc[ind] = varsALL$Slope_presc[indA]; varsALL$Slope_grass[ind] = varsALL$Slope_grass[indA]
varsALL$Intercept_ag[ind] = varsALL$Intercept_ag[indA];varsALL$Intercept_presc[ind] = varsALL$Intercept_presc[indA]; varsALL$Intercept_grass[ind] = varsALL$Intercept_grass[indA]
varsALL$Rval_ag[ind] = varsALL$Rval_ag[indA];varsALL$Rval_presc[ind] = varsALL$Rval_presc[indA]; varsALL$Rval_grass[ind] = varsALL$Rval_grass[indA]

corVAL = 0.5
for (i in 1:length(varsALL$vars)){
  # only correct if at least 50% of variability is driven by MCE
  if (abs( varsALL$Rval_ag[i]^2) >= corVAL ){
    ind = which(allBOTH.filter$variable == varsALL$vars[i] & allBOTH.filter$fuel2 == 'agriculture')
    EF1= allBOTH.filter$MCE[ind]*varsALL$Slope_ag[i] +varsALL$Intercept_ag[i]
    EF2= (0.92*varsALL$Slope_ag[i] +varsALL$Intercept_ag[i])#*cContent/cContentAG
    allBOTH.filter$FinalEF_MCE92[ind] = allBOTH.filter$FinalEF[ind]*(EF2/EF1)
  }

  if (abs( varsALL$Rval_presc[i]^2) >= corVAL ){ # if there IS a correlation
    ind = which(allBOTH.filter$variable == varsALL$vars[i] & allBOTH.filter$fuel2 == 'prescribed')
    EF1= allBOTH.filter$MCE[ind]*varsALL$Slope_presc[i] +varsALL$Intercept_presc[i]
      EF2= (0.92*varsALL$Slope_presc[i] +varsALL$Intercept_presc[i])#*cContent/cContentPresc
      allBOTH.filter$FinalEF_MCE92[ind] = allBOTH.filter$FinalEF[ind]*(EF2/EF1)
  }
  if (abs(varsALL$Rval_grass[i]^2) >= corVAL){
    ind = which(allBOTH.filter$variable == varsALL$vars[i] & allBOTH.filter$fuel2 == 'grass')
    EF1= allBOTH.filter$MCE[ind]*varsALL$Slope_grass[i] +varsALL$Intercept_grass[i]
    EF2= (0.92*varsALL$Slope_grass[i] +varsALL$Intercept_grass[i])#*cContent/cContentGrassShrub
    allBOTH.filter$FinalEF_MCE92[ind] = allBOTH.filter$FinalEF[ind]*(EF2/EF1)
  }
}
# If correction went below zero
ind = which(allBOTH.filter$FinalEF_MCE92 < 0)
allBOTH.filter$FinalEF_MCE92[ind] = allBOTH.filter$FinalEF[ind]
ind = which(varsALL$Rval_ag == 0);varsALL$Rval_ag[ind] = NaN
ind = which(varsALL$Rval_presc == 0);varsALL$Rval_presc[ind] = NaN
ind = which(varsALL$Rval_grass == 0);varsALL$Rval_grass[ind] = NaN

write.csv(varsALL, file='varsALL.csv')

ind = which(allBOTH.filter$variable == 'CH4_DACOM_DISKIN' & allBOTH.filter$fuel2 == 'agriculture')
mean(allBOTH.filter$FinalEF[ind], na.rm=TRUE)
mean(allBOTH.filter$MCE[ind], na.rm=TRUE)
mean(allBOTH.filter$FinalEF_MCE92[ind], na.rm=TRUE)
ind = which(allBOTH.filter$variable == 'CH4_DACOM_DISKIN' & allBOTH.filter$fuel2 == 'prescribed' & is.finite(allBOTH.filter$FinalEF))
mean(allBOTH.filter$FinalEF[ind], na.rm=TRUE)
mean(allBOTH.filter$MCE[ind], na.rm=TRUE)
mean(allBOTH.filter$FinalEF_MCE92[ind], na.rm=TRUE)
