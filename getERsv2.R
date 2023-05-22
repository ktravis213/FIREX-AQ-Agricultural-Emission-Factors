require(OrgMassSpecR)
TEMP = 298
source("SpeciesProperties.R")
# ------- ============ FUNCTIONS ================= ----    
# -- Main EF/ER calculation function ----
ERsEFsALLhzv2 = function(pass, background,xspecies, fire, fuel, res, SLOW, blake, blakeBG, 
                         doBlake,gilman, gilmanBG, doGilman,becky, beckyBG, doBecky){
  debugKT = TRUE
  useFASTBG = TRUE
  OAtoOC = NaN # for 1hz data
  HtoC = NaN
  OtoC = NaN
  f43= NaN
  f44 = NaN
  f57 = NaN
  f60 =  NaN
  f82 = NaN
  f91 = NaN
  
  # CO/CO2/CH4
  # ------------------------------------------ DISKIN --------------
  if (debugKT == 1){print('DISKIN')}
  nn = colnames(pass)
  
  pass$CO2_ppb = pass$CO2_7000_ppm_DISKIN*1E3 ; background$CO2_ppb = background$CO2_7000_ppm_DISKIN*1E3
  #UNCERTAINTY: CO2<500 ppm: 0.25 ppm accuracy, 0.1 ppm precision @ 1 Hz; CO2>500 ppm: 2% total error @ 1 Hz
  pass$CO2_error = 250
  pass$CO_error = pass$CO_Unc_DACOM_DISKIN
  ind = which(pass$CO2_7000_ppm_DISKIN > 500)
  if (length(ind) > 0){pass$CO2_error[ind] = pass$CO2_7000_ppm_DISKIN[ind]*0.02*1E3}
  
  pass$CH4_error = pass$CH4_DACOM_DISKIN*0.01 #$CH4 - 1% or 20 ppbv
  # Are we doing dCO2 or dCO? 
  if (debugKT == 1){print(xspecies)}
  if (xspecies == 'CO2_ppb'){ xerror = 'CO2_error'}
  if (xspecies == 'CO_DACOM_DISKIN'){ xerror = 'CO_error'}
  
  pP = smallfun(pass,background,xspecies,'CO2_ppb', xerror, 'CO2_error',SLOW)
  pPCO = smallfun(pass,background,'CO_DACOM_DISKIN','CO2_ppb', 'CO_error', 'CO2_error',SLOW) #CO2/CO
  
  p1 = smallfun(pass,background,xspecies,'CO_DACOM_DISKIN', xerror, 'CO_error',SLOW) # CO/CO2
  p1CO = smallfun(pass,background,'CO_DACOM_DISKIN','CO_DACOM_DISKIN', 'CO_error', 'CO_error',SLOW) #CO/CO
  # BAD PASS IF max CO < 200 or r2 < 0.9? or CH4 E/R is less than zero
  p2 = smallfun(pass,background,xspecies,'CH4_DACOM_DISKIN', xerror, 'CH4_error',SLOW)
  p2CO = smallfun(pass,background,'CO_DACOM_DISKIN', 'CH4_DACOM_DISKIN', 'CO_error','CH4_error',SLOW)
  
  if (debugKT == 1){print(c(max(pass$CO_DACOM_DISKIN, na.rm=TRUE) ,pPCO$r_sq))}
  # don't calculate a pass if there is no CO or CO2 or CH4 data
  if (is.nan(p2$slope) |  is.nan(pP$slope) |  is.nan(p1$slope)){
    mult = NaN
  } else{
    mult = 1
    # don't plot bad passes
    #plotpass5hzSAVE(pass, fire)
  }
  
  # --------- (MCE) 
  mce = 1/(1+p1$slope)
  mce_uncert = abs(mce - 1/(1+( p1$uncertainty_slope)))
  mce_int = p1$ExcessX/(p1$ExcessX+p1$ExcessS)
  if (length(mce_int) == 0){mce_int = NaN}
  
  passERs   = c( pP$slope, p1$slope, p2$slope) # CO2, CO, CH4
  passERsCO   = c( pPCO$slope, p1CO$slope, p2CO$slope) # CO2, CO, CH4
  passErrors = c(mean(pass$CO2_error, na.rm=TRUE),mean(pass$CH4_error,na.rm=TRUE), mean(pass$CO_error, na.rm=TRUE))
  passERsint = c( pP$ERint, p1$ERint, p2$ERint)
  passERsintfill = c( pP$ERintfill, p1$ERintfill, p2$ERintfill)
  passERsCOint = c( pPCO$ERint, p1CO$ERint, p2CO$ERint)
  passERsCOintfill = c( pPCO$ERintfill, p1CO$ERintfill, p2CO$ERintfill)
  # maybe these would be useful?
  passERszero = c( pP$intercept, p1$intercept, p2$intercept)
  passERsCOzero = c( pPCO$intercept, p1CO$intercept, p2CO$intercept)
  
  passBGX      = c( pP$BGX, p1$BGX, p2$BGX) # mce, CO2 /CO2, CO/CO2
  passBGCO     = c( pPCO$BGX, p1CO$BGX, p2CO$BGX) # mce, CO2 /CO2, CO/CO2
  passBGS      = c( pP$BGS,   p1$BGS,   p2$BGS) # mce, CO2 /CO2, CO/CO2
  
  passRsqsCO   = c(pPCO$r_sq, p1CO$r_sq,  p2CO$r_sq) # CO2/CO, CO/CO, CH4/CO 
  passRsqsCO2  = c(pPCO$r_sq, p1$r_sq,  p2$r_sq) # CO2/CO, CO/CO2, CH4/CO 
  
  passExcessX  = c(pP$ExcessX, p1$ExcessX,  p2$ExcessX) #
  passExcessS  = c(pP$ExcessS, p1$ExcessS,  p2$ExcessS) # 
  
  passsigma = c( pP$uncertainty_slope, p1$uncertainty_slope, p2$uncertainty_slope)
  passMAX   = c(maxfun(pass, 'CO2_ppb'),maxfun(pass, 'CO_DACOM_DISKIN'),maxfun(pass, 'CH4_DACOM_DISKIN') )
  names1    = c("CO2_ppb","CO_DACOM_DISKIN","CH4_DACOM_DISKIN")
  formulas = c("CO2","CO","CH4")
  usenames = c("Carbon Dioxide","Carbon Monoxide","Methane")
  mWs       = c( mWCO2, mWCO, mWCH4)             
  nCs       = c(1, 1, 1) 
  M =(6.022E23*1.01325E5)/(8.31*TEMP)/1E6
  passOHrate = c(0,1.44E-13*(1+(M/4.2E19)),1.85E-12*exp(-1690/TEMP)) # at STP, 298K, 101325Pa,f rom the MCM
  isNMVOC = c(0,0,0) # 0 is inorganic, , 1 is VOC, 2 is nitrogen containing, 3 is halogen, 4 aerosol
  kind = c("CO2","CO","CH4")
  
  if (res != 5){ # not present in 5 Hz
    # ------------------------------------------ RYERSON ---------------
    if (debugKT == 1){print("RYERSON")}
    pass$NO2_error = 0.04 ;  pass$NO_error = 0.08
    pass$O3_error = 2.0 # it is really a column variable
    pass$NOy_error = pass$NOy_CL_RYERSON *0.13
    
    p3  = smallfun(pass,background,xspecies,'NO2_CL_RYERSON', xerror, 'NO2_error',SLOW)
    p4  = smallfun(pass,background,xspecies,'NO_CL_RYERSON', xerror, 'NO_error',SLOW)
    p4b = smallfun(pass,background,xspecies,'O3_CL_RYERSON', xerror, 'O3_error',SLOW)
    p4c = smallfun(pass,background,xspecies,'NOy_CL_RYERSON', xerror, 'NOy_error',SLOW)
    
    passERs   = c(passERs,p3$slope, p4$slope, p4b$slope, p4c$slope)
    passERsint   = c(passERsint,p3$ERint, p4$ERint, p4b$ERint, p4c$ERint)
    passERsintfill   = c(passERsintfill,p3$ERintfill, p4$ERintfill, p4b$ERintfill, p4c$ERintfill)
    
    passBGX    = c( passBGX,p3$BGX, p4$BGX, p4b$BGX, p4c$BGX) # 
    passBGS      = c( passBGS,  p3$BGS,   p4$BGS,   p4b$BGS,   p4c$BGS) # 
    passERszero = c( passERszero,p3$intercept, p4$intercept, p4b$intercept, p4c$intercept) 
    
    passRsqsCO2  = c(passRsqsCO2,p3$r_sq,p4$r_sq,p4b$r_sq,p4c$r_sq)
    p3  = smallfun(pass,background,'CO_DACOM_DISKIN','NO2_CL_RYERSON', xerror, 'NO2_error',SLOW)
    p4  = smallfun(pass,background,'CO_DACOM_DISKIN','NO_CL_RYERSON', xerror, 'NO_error',SLOW)
    p4b = smallfun(pass,background,'CO_DACOM_DISKIN','O3_CL_RYERSON', xerror, 'O3_error',SLOW)
    p4c = smallfun(pass,background,'CO_DACOM_DISKIN','NOy_CL_RYERSON', xerror, 'NOy_error',SLOW)
    passBGCO    = c( passBGCO,p3$BGX, p4$BGX, p4b$BGX, p4c$BGX) # 
    passERsCO  = c(passERsCO,p3$slope,p4$slope,p4b$slope,p4c$slope)
    passErrors = c(passErrors, pass$NO2_error, pass$NO_error, pass$O3_error, pass$NOy_error)
    passRsqsCO  = c(passRsqsCO,p3$r_sq,p4$r_sq,p4b$r_sq,p4c$r_sq)
    passERsCOint   = c(passERsCOint,p3$ERint, p4$ERint, p4b$ERint, p4c$ERint)
    passERsCOintfill   = c(passERsCOintfill,p3$ERintfill, p4$ERintfill, p4b$ERintfill, p4c$ERintfill)  
    passExcessX  = c(passExcessX,p3$ExcessX,p4$ExcessX,p4b$ExcessX,p4c$ExcessX)
    passExcessS  = c(passExcessS,p3$ExcessS,p4$ExcessS,p4b$ExcessS,p4c$ExcessS)
    passERsCOzero = c( passERsCOzero,p3$intercept, p4$intercept, p4b$intercept, p4c$intercept) 
    
    passsigma = c(passsigma, p3$uncertainty_slope, p4$uncertainty_slope, p4b$uncertainty_slope, p4c$uncertainty_slope)
    passMAX   = c(passMAX, maxfun(pass, 'NO2_CL_RYERSON'),maxfun(pass, 'NO_CL_RYERSON'),maxfun(pass, 'O3_CL_RYERSON'),maxfun(pass, 'NOy_CL_RYERSON') )
    passOHrate = c(passOHrate, NaN, NaN, NaN, NaN)
    isNMVOC = c(isNMVOC,2,2,5,5)
    kind = c(kind,"NOy","NOy","O3","NOy")
    names1    = c(names1,"NO2_RYERSON","NO_RYERSON","O3_RYERSON","NOy_RYERSON")
    formulas = c(formulas,"NO2","NO","O3","NOy")
    usenames = c(usenames,"Nitrogen dioxide","Nitrogen oxide","Ozone","NOy")
    mWs       = c(mWs, mWNO2, mWNO, mWO3, NaN)
    nCs       = c(nCs,  0,      0,    0,   NaN)
    # ------------------------------------------ JIMENEZ EESI ----------------
    if (debugKT == 1){print("JIMENEZ EESI")}
    convOC = (1.29/28.97)*12E-3  ; 
    # if (debugKT == 1){('WARNING, need to fix OAtoOC')
    convNO3 = (1.29/28.97)*62E-3  ; convCl = (1.29/28.97)*35E-3 ;convSSA = (1.29/28.97)*31.4E-3
    convClO4 = (1.29/28.97)*99.45E-3 ; convMSA = (1.29/28.97)*96.1E-3; convNH4 = (1.29/28.97)*18E-3; convSO4 = (1.29/28.97)*96E-3 
    convBr = (1.29/28.97)*79.9E-3 ; convI = (1.29/28.97)*126.9E-3 ; convC6H10O5 = (1.29/28.97)*162.141E-3; convC6H5NO4= (1.29/28.97)*155.1081E-3
    convK = (1.29/28.97)*39.1E-3 
    nair_STP = (6.022E23*1.013E5)/(273*8.31)
    passT = mean(pass$Static_Air_Temp, na.rm=TRUE)+273.15
    passP = mean(pass$Static_Pressure.x, na.rm=TRUE)
    passPatm = passP * 0.000986923
    pass$nair = (6.022E23*passP*100)/(passT*8.31)
    RHS = pass$nair * (1/6.022E23)*298/passT * passPatm/1/1E3 
    # BG
    backgroundT = mean(background$Static_Air_Temp, na.rm=TRUE)+273.15
    backgroundP = mean(background$Static_Pressure.x, na.rm=TRUE)
    backgroundPatm = backgroundP * 0.000986923
    background$nair = (6.022E23*backgroundP*100)/(backgroundT*8.31)
    BG_RHS = background$nair * (1/6.022E23)*298/backgroundT * backgroundPatm/1/1E3 
    
    pass$C6H10O5_JIMENEZ_ppb = pass$C6H10O5_JIMENEZ/pass$StdToVol_EESI_JIMENEZ/RHS/mWC6H10O5
    background$C6H10O5_JIMENEZ_ppb = background$C6H10O5_JIMENEZ/background$StdToVol_EESI_JIMENEZ/BG_RHS/mWC6H10O5
    pass$C6H10O5_JIMENEZ_error = pass$C6H10O5_JIMENEZ_ppb *0.5
    mmJ = maxfun(pass,'C6H10O5_JIMENEZ')
    p6k = smallfun(pass,background,xspecies,'C6H10O5_JIMENEZ_ppb', xerror, 'C6H10O5_JIMENEZ_error',SLOW)
    p6kCO = smallfun(pass,background,'CO_DACOM_DISKIN','C6H10O5_JIMENEZ_ppb', xerror, 'C6H10O5_JIMENEZ_error',SLOW)
    
    pass$C6H5NO4_JIMENEZ_ppb = pass$C6H5NO4_JIMENEZ/pass$StdToVol_EESI_JIMENEZ/RHS/mWC6H5NO4
    background$C6H5NO4_JIMENEZ_ppb = background$C6H5NO4_JIMENEZ/background$StdToVol_EESI_JIMENEZ/BG_RHS/mWC6H5NO4
    pass$C6H5NO4_JIMENEZ_error = pass$C6H5NO4_JIMENEZ_ppb *0.5
    mmJ2 = maxfun(pass,'C6H5NO4_JIMENEZ')
    p6l = smallfun(pass,background,xspecies,'C6H5NO4_JIMENEZ_ppb', xerror, 'C6H5NO4_JIMENEZ_error',SLOW)
    p6lCO = smallfun(pass,background,'CO_DACOM_DISKIN','C6H5NO4_JIMENEZ_ppb', xerror, 'C6H5NO4_JIMENEZ_error',SLOW)
    
    passERs = c(passERs,  p6k$slope,p6l$slope)
    passERsCO = c(passERsCO, p6kCO$slope,p6lCO$slope)
    passERsint = c(passERsint,  p6k$ERint, p6l$ERint)
    passERsintfill = c(passERsintfill, p6k$ERintfill, p6l$ERintfill)
    passERsCOint = c(passERsCOint, p6kCO$ERint,p6lCO$ERint)
    passERsCOintfill = c(passERsCOintfill,p6kCO$ERintfill, p6lCO$ERintfill)
    
    passBGX = c(passBGX,  p6k$BGX, p6l$BGX)
    passBGCO = c(passBGCO,p6kCO$BGX, p6lCO$BGX)
    
    passBGS = c(passBGS, p6k$BGS, p6l$BGS)
    
    passERszero = c( passERszero, p6k$intercept, p6l$intercept)
    passERsCOzero = c( passERsCOzero,p6kCO$intercept, p6lCO$intercept)
    passRsqsCO2 = c(passRsqsCO2, p6k$r_sq, p6l$r_sq)
    passRsqsCO = c(passRsqsCO, p6kCO$r_sq, p6lCO$r_sq)
    passExcessX  = c(passExcessX, p6k$ExcessX, p6l$ExcessX)
    passExcessS  = c(passExcessS,p6k$ExcessS,p6l$ExcessS)
    
    passsigma = c(passsigma, p6k$uncertainty_slope, p6l$uncertainty_slope)
    passMAX   = c( passMAX,  mmJ, mmJ2)
    passErrors = c(passErrors, pass$C6H10O5_JIMENEZ_error, pass$C6H5NO4_JIMENEZ_error)
    passOHrate = c(passOHrate, NaN, NaN)
    isNMVOC = c(isNMVOC,4,4)
    kind = c(kind,"aerosol","aerosol")
    
    names1=c(names1,"C6H10O5_JIMENEZ","C6H5NO4_JIMENEZ")
    formulas = c(formulas,"C6H10O5","C6H5NO4")
    usenames = c(usenames,"Levoglucosan","4-nitrocatechol")
    mWs = c(mWs,  mWC6H10O5, mWC6H5NO4)
    nCs  = c(nCs,6, 6)
    
    if (debugKT == 1){print(c(length(passERs), length(passERsint), length(passERsintfill),length(passERsCO),length(passsigma), length(passRsqsCO2), 
                              length(passRsqsCO),length(passMAX),length(passOHrate),length(passErrors),
                              length(names1), length(mWs), length(nCs), length(passBGS), length(passBGX), 
                              length(passExcessX), length(passExcessS), length(passERsCOint), length(passERsCOintfill),
                              length(passERszero), length(passERsCOzero), length(isNMVOC), length(kind)))}
     
    # ------------------------------------------ SCHWARZ ---------
    if (debugKT == 1){print("SCHWARZ")}
    # random zeros in BC
    ind = which(pass$BC_mass_90_550_nm_SCHWARZ == 0)
    if (length(ind) > 0){ pass$BC_mass_90_550_nm_SCHWARZ[ind] = NaN}
    pass$BC_mass_90_550_nm_SCHWARZ_ppb = pass$BC_mass_90_550_nm_SCHWARZ/pass$StdToVol_EESI_JIMENEZ/RHS/mWOC/1E3
    background$BC_mass_90_550_nm_SCHWARZ_ppb = background$BC_mass_90_550_nm_SCHWARZ/background$StdToVol_EESI_JIMENEZ/BG_RHS/mWOC/1E3
    pass$BC_mass_90_550_nm_SCHWARZ_error = pass$BC_mass_90_550_nm_SCHWARZ_ppb*0.2
    p7 = smallfun(pass,background,xspecies,'BC_mass_90_550_nm_SCHWARZ_ppb', xerror, 'BC_mass_90_550_nm_SCHWARZ_error',SLOW)
    p7CO = smallfun(pass,background,'CO_DACOM_DISKIN','BC_mass_90_550_nm_SCHWARZ_ppb', xerror, 'BC_mass_90_550_nm_SCHWARZ_error',SLOW)
    passERs = c(passERs, p7$slope)
    passERsCO = c(passERsCO, p7CO$slope)
    
    passERsint = c(passERsint, p7$ERint)
    passERsintfill = c(passERsintfill, p7$ERintfill)
    passERsCOint = c(passERsCOint, p7CO$ERint)
    passERsCOintfill = c(passERsCOintfill, p7CO$ERintfill)
    passBGX = c(passBGX,  p7$BGX)
    passBGCO = c(passBGCO,  p7CO$BGX)
    
    passBGS = c(passBGS,  p7$BGS)
    passERszero = c( passERszero,p7$intercept)
    passERsCOzero = c( passERsCOzero,p7CO$intercept)
    
    passRsqsCO2 = c(passRsqsCO2, p7$r_sq)
    passRsqsCO = c(passRsqsCO, p7CO$r_sq)
    passExcessX  = c(passExcessX,p7$ExcessX)
    passExcessS  = c(passExcessS,p7$ExcessS)
    
    passsigma = c(passsigma, p7$uncertainty_slope)
    passMAX   = c(passMAX, maxfun(pass, 'BC_mass_90_550_nm_SCHWARZ'))
    passOHrate = c(passOHrate, NaN)
    passErrors=c(passErrors,pass$BC_mass_90_550_nm_SCHWARZ_error )
    isNMVOC = c(isNMVOC,4)
    kind = c(kind,'aerosol')
    names1=c(names1,"BC_SCHWARZ")
    formulas = c(formulas,"BC")
    usenames = c(usenames,"Black carbon")
    mWs = c(mWs, mWOC)
    nCs  = c(nCs,   1)
    # ------------------------------------------ FRIED ------------------------------------------ 
    if (debugKT == 1){print('FRIED')}
    pass$CH2O_CAMS_pptv_FRIED_ppb = pass$CH2O_CAMS_pptv_FRIED/1E3
    background$CH2O_CAMS_pptv_FRIED_ppb = background$CH2O_CAMS_pptv_FRIED/1E3
    pass$tmperror = (pass$CH2O_CAMS_pptv_FRIED*0.02)/1E3
    passErrors=c(passErrors,pass$tmperror )
    
    pT = smallfun(pass,background,xspecies,'CH2O_CAMS_pptv_FRIED_ppb', xerror, 'tmperror',SLOW)
    pTCO = smallfun(pass,background,'CO_DACOM_DISKIN','CH2O_CAMS_pptv_FRIED_ppb', xerror, 'tmperror',SLOW)
    passERs = c(passERs, pT$slope) 
    passERsCO = c(passERsCO, pTCO$slope) 
    
    passERsint = c(passERsint, pT$ERint)
    passERsintfill = c(passERsintfill, pT$ERintfill)
    passERsCOint = c(passERsCOint, pTCO$ERint)
    passERsCOintfill = c(passERsCOintfill, pTCO$ERintfill)
    passBGX = c(passBGX,  pT$BGX)
    passBGCO = c(passBGCO,  pTCO$BGX)
    passBGS = c(passBGS,  pT$BGS)
    
    passERszero = c( passERszero,pT$intercept)
    passERsCOzero = c( passERsCOzero,pTCO$intercept)
    
    passRsqsCO2 = c(passRsqsCO2,pT$r_sq)
    passRsqsCO = c(passRsqsCO,pTCO$r_sq)
    passExcessX  = c(passExcessX,pT$ExcessX)
    passExcessS  = c(passExcessS,pT$ExcessS)
    
    passsigma = c(passsigma,pT$uncertainty_slope)
    passMAX = c(passMAX, maxfun(pass,'CH2O_CAMS_pptv_FRIED'))
    passOHrate = c(passOHrate, OHch2o)
    isNMVOC = c(isNMVOC,1)
    kind = c(kind, 'CH2O')
    
    names1=c(names1,'CH2O_CAMS_pptv_FRIED')
    formulas = c(formulas, "CH2O")
    usenames = c(usenames, "Formaldehyde")
    mWs       = c(mWs, mWCH2O)
    nCs = c(nCs,1)
    
    pass$C2H6_CAMS_pptv_FRIED_ppb = pass$C2H6_CAMS_pptv_FRIED/1E3
    background$C2H6_CAMS_pptv_FRIED_ppb = background$C2H6_CAMS_pptv_FRIED/1E3
    pass$tmperror = (pass$C2H6_CAMS_pptv_FRIED/1E3)*0.02
    passErrors=c(passErrors,pass$tmperror )
    pT = smallfun(pass,background,xspecies,'C2H6_CAMS_pptv_FRIED_ppb', xerror, 'tmperror',SLOW)
    pTCO = smallfun(pass,background,'CO_DACOM_DISKIN','C2H6_CAMS_pptv_FRIED_ppb', xerror, 'tmperror',SLOW)
    passERs = c(passERs, pT$slope)    
    passERsCO = c(passERsCO, pTCO$slope) 
    
    passERsint = c(passERsint, pT$ERint)
    passERsintfill = c(passERsintfill, pT$ERintfill)
    passERsCOint = c(passERsCOint, pTCO$ERint)
    passERsCOintfill = c(passERsCOintfill, pTCO$ERintfill)
    passBGX = c(passBGX,  pT$BGX)
    passBGCO = c(passBGCO,  pTCO$BGX)
    passBGS = c(passBGS,  pT$BGS)
    
    passERszero = c( passERszero,pT$intercept)
    passERsCOzero = c( passERsCOzero,pTCO$intercept)
    
    passRsqsCO2 = c(passRsqsCO2,pT$r_sq)
    passRsqsCO = c(passRsqsCO,pTCO$r_sq)
    passExcessX  = c(passExcessX,pT$ExcessX)
    passExcessS  = c(passExcessS,pT$ExcessS)
    passsigma = c(passsigma,pT$uncertainty_slope)
    passMAX = c(passMAX,maxfun(pass,'C2H6_CAMS_pptv_FRIED' ) )
    passOHrate = c(passOHrate, OHc2h6)
    names1=c(names1,'C2H6_CAMS_pptv_FRIED')
    formulas = c(formulas, "C2H6")
    usenames = c(usenames, "Ethane")
    isNMVOC = c(isNMVOC,1)
    kind = c(kind, 'alkane')
    
    mWs       = c(mWs, mWC2H6)
    nCs = c(nCs,2)

    # ------------------------------------------ WOMACK ------------------------------------------ 
    if (debugKT == 1){print('WOMACK')}
    #NO2, HNO2, CH3COCHO, CHOCHO
    pass$no2_error = pass$NO2_ACES_WOMACK * 0.04 + 0.26
    pass$hno2_error = pass$HNO2_ACES_WOMACK*0.04 + 0.60 
    pass$ch3cocho_error = pass$CH3COCHO_ACES_WOMACK*0.04 + 0.52 #MGLY
    pass$chocho_error = pass$CHOCHO_ACES_WOMACK * .04 + 0.14 #GLY
    pW1 = smallfun(pass,background,xspecies,'NO2_ACES_WOMACK',      xerror, 'no2_error',SLOW) 
    pW2 = smallfun(pass,background,xspecies,'HNO2_ACES_WOMACK',     xerror, 'no2_error',SLOW) 
    pW3 = smallfun(pass,background,xspecies,'CH3COCHO_ACES_WOMACK', xerror, 'ch3cocho_error',SLOW) 
    pW4 = smallfun(pass,background,xspecies,'CHOCHO_ACES_WOMACK',   xerror, 'chocho_error',SLOW) 
    passERs = c(passERs, pW1$slope,pW2$slope,pW3$slope,pW4$slope)    
    passERsint = c(passERsint, pW1$ERint,pW2$ERint,pW3$ERint,pW4$ERint)    
    passERsintfill = c(passERsintfill, pW1$ERintfill,pW2$ERintfill,pW3$ERintfill,pW4$ERintfill)    
    passErrors=c(passErrors,pass$no2_error,pass$hno2_error,pass$ch3cocho_error,pass$chocho_error)
    
    passBGX = c(passBGX,  pW1$BGX,pW2$BGX,pW3$BGX,pW4$BGX)   
    passBGS = c(passBGS,  pW1$BGS,pW2$BGS,pW3$BGS,pW4$BGS)   
    passERszero = c( passERszero,pW1$intercept,pW2$intercept,pW3$intercept,pW4$intercept)
    
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq,pW2$r_sq,pW3$r_sq,pW4$r_sq) 
    pW1 = smallfun(pass,background,'CO_DACOM_DISKIN','NO2_ACES_WOMACK',      xerror, 'no2_error',SLOW) 
    pW2 = smallfun(pass,background,'CO_DACOM_DISKIN','HNO2_ACES_WOMACK',     xerror, 'no2_error',SLOW) 
    pW3 = smallfun(pass,background,'CO_DACOM_DISKIN','CH3COCHO_ACES_WOMACK', xerror, 'ch3cocho_error',SLOW) 
    pW4 = smallfun(pass,background,'CO_DACOM_DISKIN','CHOCHO_ACES_WOMACK',   xerror, 'chocho_error',SLOW) 
    passERsCO = c(passERsCO, pW1$slope,pW2$slope,pW3$slope,pW4$slope) 
    passERsCOint = c(passERsCOint, pW1$ERint,pW2$ERint,pW3$ERint,pW4$ERint)    
    passERsCOintfill = c(passERsCOintfill, pW1$ERintfill,pW2$ERintfill,pW3$ERintfill,pW4$ERintfill)  

    passBGCO = c(passBGCO,  pW1$BGX,pW2$BGX,pW3$BGX,pW4$BGX)   
    passRsqsCO = c(passRsqsCO, pW1$r_sq,pW2$r_sq,pW3$r_sq,pW4$r_sq) 
    passExcessX  = c(passExcessX,pW1$ExcessX,pW2$ExcessX,pW3$ExcessX,pW4$ExcessX) 
    passExcessS  = c(passExcessS,pW1$ExcessS,pW2$ExcessS,pW3$ExcessS,pW4$ExcessS) 
    passERsCOzero = c( passERsCOzero,pW1$intercept,pW2$intercept,pW3$intercept,pW4$intercept)
    
    passsigma = c(passsigma, pW1$uncertainty_slope,pW2$uncertainty_slope,pW3$uncertainty_slope,pW4$uncertainty_slope)    
    passMAX = c(passMAX, maxfun(pass,'NO2_ACES_WOMACK'),maxfun(pass,'HNO2_ACES_WOMACK'),maxfun(pass,'CH3COCHO_ACES_WOMACK'),maxfun(pass,'CHOCHO_ACES_WOMACK'))    
    passOHrate = c(passOHrate, NaN, NaN, OHch3chocho, OHchocho )
    
    isNMVOC = c(isNMVOC,2, 2,1,1)
    kind = c(kind, 'nitrogen','nitrogen','oVOC','oVOC')
    
    names1=c(names1, 'NO2_ACES_WOMACK','HNO2_ACES_WOMACK','CH3COCHO_ACES_WOMACK','CHOCHO_ACES_WOMACK')
    formulas = c(formulas, "NO2","HNO2","CH3COCHO","CHOCHO")
    usenames = c(usenames, "Nitrogen dioxide","Nitrous acid","Methylglyoxal","Glyoxal")
    
    mWs    = c(mWs, mWNO2, mWHNO2, mWCH3COCHO, mWCHOCHO)
    nCs    = c(nCs,0, 0, 3, 2)
    
    if (debugKT == 1){print(c(length(passERs), length(passERsint), length(passERsintfill),length(passERsCO),length(passsigma), length(passRsqsCO2), length(passRsqsCO),length(passMAX),length(passOHrate),
                              length(names1), length(mWs), length(nCs), length(passBGS), length(passBGX), 
                              length(passExcessX), length(passExcessS), length(passERsCOint), length(passERsCOintfill),
                              length(passERszero), length(passERsCOzero), length(isNMVOC), length(kind)))}
    
    # ------------------------------------------ STCLAIR ------------------------------------------ 
    # if (debugKT == 1){print('STCLAIr',SLOW)}
    # #  Calibration uncertainty is +/-10 percent. Offset uncertainty is +/- 100 pptv.
    # ind = which(colnames(pass) == 'NO2_CANOE')
    # if (length(ind) > 0){pass$NO2_CANOE_STCLAIR = pass$NO2_CANOE; background$NO2_CANOE_STCLAIR = background$NO2_CANOE}
    # pass$no2_error2 = pass$NO2_CANOE_STCLAIR * 0.1 + 100
    # pass$no2_error2 = pass$no2_error2/1E3
    # pass$NO2_CANOE_STCLAIR_ppb = pass$NO2_CANOE_STCLAIR/1E3
    # background$NO2_CANOE_STCLAIR_ppb = background$NO2_CANOE_STCLAIR/1E3
    # pW1 = smallfun(pass,background,xspecies,'NO2_CANOE_STCLAIR_ppb',      xerror, 'no2_error2', SLOW) 
    # pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','NO2_CANOE_STCLAIR_ppb',      xerror, 'no2_error2', SLOW) 
    # passERs = c(passERs, pW1$slope)
    # passERsCO = c(passERsCO, pW1CO$slope)
    # passERsint = c(passERsint, pW1$ERint)
    # passERsintfill = c(passERsintfill, pW1$ERintfill)
    # passERsCOint = c(passERsCOint, pW1CO$ERint)
    # passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)
    # passBGX = c(passBGX,  pW1$BGX)
    # passBGCO = c(passBGCO,  pW1CO$BGX)
    # passBGS = c(passBGS,  pW1$BGS)
    # passERszero = c( passERszero,pW1$intercept)
    # passERsCOzero = c( passERsCOzero,pW1CO$intercept)
    # 
    # passRsqsCO2 = c(passRsqsCO2, pW1$r_sq) 
    # passRsqsCO = c(passRsqsCO, pW1CO$r_sq) 
    # passExcessX  = c(passExcessX,pW1$ExcessX)
    # passExcessS  = c(passExcessS,pW1$ExcessS)
    # 
    # passsigma = c(passsigma, pW1$uncertainty_slope)
    # passMAX = c(passMAX, maxfun(pass,'NO2_CANOE_STCLAIR')) 
    # passOHrate = c(passOHrate, NaN)
    # isNMVOC = c(isNMVOC,2)
    # kind = c(kind,'NOy')
    # names1=c(names1, 'NO2_CANOE_STCLAIR')
    # formulas = c(formulas,"NO2")
    # usenames = c(usenames,"Nitrogen dioxide")
    # mWs       = c(mWs, mWNO2)
    # nCs = c(nCs,0)
    
    # ------------------------------------------ VERES -------------------------
    if (debugKT == 1){print('VERES')}
    #HNO2
    pass$tmpy = (pass$HNO2_NOAACIMS_VERES )/1E3 # need to convert to ppb
    background$tmpy = (background$HNO2_NOAACIMS_VERES )/1E3
    pass$tmperror = (pass$HNO2_NOAACIMS_VERES * .15 + 3)/1E3
    passErrors=c(passErrors,pass$tmperror)
    
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)
    passERsintfill = c(passERsintfill, pW1$ERintfill)
    passBGX = c(passBGX,  pW1$BGX)
    passBGCO = c(passBGCO,  pW1CO$BGX)
    passBGS = c(passBGS,  pW1$BGS)
    passExcessX  = c(passExcessX,pW1$ExcessX)
    passExcessS  = c(passExcessS,pW1$ExcessS)
    passERszero = c( passERszero,pW1$intercept)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passERsCO = c(passERsCO, pW1CO$slope)
    passERsCOint = c(passERsCOint, pW1CO$ERint)
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)
    passERsCOzero = c( passERsCOzero,pW1CO$intercept)
    passRsqsCO = c(passRsqsCO, pW1CO$r_sq)   
    passsigma = c(passsigma , pW1$uncertainty_slope)   
    passMAX = c(passMAX, maxfun(pass,'tmpy'))
    passOHrate = c(passOHrate, NaN)
    isNMVOC = c(isNMVOC,2)
    kind = c(kind,'NOy')
    mWs       = c(mWs, mWHNO2)
    nCs = c(nCs,0)
    if (debugKT == 1){print('N2O5')}
    
    #N2O5_NOAACIMS_VERES# (15% + 2 pptv)
    pass$tmpy = (pass$N2O5_NOAACIMS_VERES )/1E3 # need to convert to ppb
    background$tmpy = (background$N2O5_NOAACIMS_VERES )/1E3 # need to convert to ppb
    pass$tmperror = (pass$N2O5_NOAACIMS_VERES * .15 + 2)/1E3
    passErrors=c(passErrors,pass$tmperror)
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)
    passERsintfill = c(passERsintfill, pW1$ERintfill)
    passBGX = c(passBGX,  pW1$BGX)
    passBGCO = c(passBGCO,  pW1CO$BGX)
    passBGS = c(passBGS,  pW1$BGS)
    passExcessX  = c(passExcessX,pW1$ExcessX)
    passExcessS  = c(passExcessS,pW1$ExcessS)
    passERszero = c( passERszero,pW1$intercept)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passERsCO = c(passERsCO, pW1CO$slope)
    passERsCOint = c(passERsCOint, pW1CO$ERint)
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)
    passERsCOzero = c( passERsCOzero,pW1CO$intercept)
    passRsqsCO = c(passRsqsCO, pW1CO$r_sq)   
    passsigma = c(passsigma , pW1$uncertainty_slope)   
    passMAX = c(passMAX, maxfun(pass,'tmpy'))
    passOHrate = c(passOHrate, NaN)
    isNMVOC = c(isNMVOC,2)
    kind = c(kind,'NOy')
    
    mWs       = c(mWs, mWN2O5)
    nCs = c(nCs,0)
    if (debugKT == 1){print('BrCN')}
    
    #BrCN_NOAACIMS_VERES # (25% + 0.4 ppt)
    pass$tmpy = (pass$BrCN_NOAACIMS_VERES )/1E3 # need to convert to ppb
    background$tmpy = (background$BrCN_NOAACIMS_VERES )/1E3 # need to convert to ppb
    pass$tmperror = (pass$BrCN_NOAACIMS_VERES * .25 + .4)/1E3
    passErrors=c(passErrors,pass$tmperror)
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)
    passERsintfill = c(passERsintfill, pW1$ERintfill)
    passBGX = c(passBGX,  pW1$BGX)
    passBGCO = c(passBGCO,  pW1CO$BGX)
    passBGS = c(passBGS,  pW1$BGS)
    passExcessX  = c(passExcessX,pW1$ExcessX)
    passExcessS  = c(passExcessS,pW1$ExcessS)
    passERszero = c( passERszero,pW1$intercept)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passERsCO = c(passERsCO, pW1CO$slope)
    passERsCOint = c(passERsCOint, pW1CO$ERint)
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)
    passERsCOzero = c( passERsCOzero,pW1CO$intercept)
    passRsqsCO = c(passRsqsCO, pW1CO$r_sq)   
    passsigma = c(passsigma , pW1$uncertainty_slope)   
    passMAX = c(passMAX, maxfun(pass,'tmpy'))
    passOHrate = c(passOHrate, NaN)
    isNMVOC = c(isNMVOC,3)
    mWs = c(mWs, mWBrCN)
    kind = c(kind,'nitrogen')
    
    nCs  = c(nCs,1)
    if (debugKT == 1){print('BrCl')}
    
    #BrCl_NOAACIMS_VERES #(25% + 0.4 ppt)
    pass$tmpy = (pass$BrCl_NOAACIMS_VERES )/1E3 # need to convert to ppb
    background$tmpy = (background$BrCl_NOAACIMS_VERES )/1E3 # need to convert to ppb
    pass$tmperror = (pass$BrCl_NOAACIMS_VERES * .25 + .4)/1E3
    passErrors=c(passErrors,pass$tmperror)
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)
    passERsintfill = c(passERsintfill, pW1$ERintfill)
    passBGX = c(passBGX,  pW1$BGX)
    passBGCO = c(passBGCO,  pW1CO$BGX)
    passBGS = c(passBGS,  pW1$BGS)
    passExcessX  = c(passExcessX,pW1$ExcessX)
    passExcessS  = c(passExcessS,pW1$ExcessS)
    passERszero = c( passERszero,pW1$intercept)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passERsCO = c(passERsCO, pW1CO$slope)
    passERsCOint = c(passERsCOint, pW1CO$ERint)
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)
    passERsCOzero = c( passERsCOzero,pW1CO$intercept)
    passRsqsCO = c(passRsqsCO, pW1CO$r_sq)   
    passsigma = c(passsigma , pW1$uncertainty_slope)   
    passMAX = c(passMAX, maxfun(pass,'tmpy'))
    passOHrate = c(passOHrate, NaN)
    isNMVOC = c(isNMVOC,3)
    kind = c(kind,'halogen')
    mWs = c(mWs, mWBrCl)
    
    nCs  = c(nCs,0)
    if (debugKT == 1){print('BrO')}
    
    #BrO_NOAACIMS_VERES # (15% + 0.4 ppt)
    pass$tmpy = (pass$BrO_NOAACIMS_VERES )/1E3 # need to convert to ppb
    background$tmpy = (background$BrO_NOAACIMS_VERES )/1E3 # need to convert to ppb
    pass$tmperror = (pass$BrO_NOAACIMS_VERES * .15 + .4)/1E3
    passErrors=c(passErrors,pass$tmperror)
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW)  
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)
    passERsintfill = c(passERsintfill, pW1$ERintfill)
    passBGX = c(passBGX,  pW1$BGX)
    passBGCO = c(passBGCO,  pW1CO$BGX)
    passBGS = c(passBGS,  pW1$BGS)
    passExcessX  = c(passExcessX,pW1$ExcessX)
    passExcessS  = c(passExcessS,pW1$ExcessS)
    passERszero = c( passERszero,pW1$intercept)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passERsCO = c(passERsCO, pW1CO$slope)
    passERsCOint = c(passERsCOint, pW1CO$ERint)
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)
    passERsCOzero = c( passERsCOzero,pW1CO$intercept)
    passRsqsCO = c(passRsqsCO, pW1CO$r_sq)   
    passsigma = c(passsigma , pW1$uncertainty_slope)   
    passMAX = c(passMAX, maxfun(pass,'tmpy'))
    passOHrate = c(passOHrate, NaN)
    isNMVOC = c(isNMVOC,3)
    kind = c(kind,'halogen')
    
    mWs = c(mWs, mWBrO)
    nCs  = c(nCs,0)
    if (debugKT == 1){print('CH3COOCl')}
    
    #CH3COOCl_NOAACIMS_VERES # (15% + 0.4 ppt)
    pass$tmpy = (pass$CH3COOCl_NOAACIMS_VERES )/1E3 # need to convert to ppb
    background$tmpy = (background$CH3COOCl_NOAACIMS_VERES )/1E3 # need to convert to ppb
    pass$tmperror = (pass$CH3COOCl_NOAACIMS_VERES * .15 + .4)/1E3
    passErrors=c(passErrors,pass$tmperror)
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)
    passERsintfill = c(passERsintfill, pW1$ERintfill)
    passBGX = c(passBGX,  pW1$BGX)
    passBGCO = c(passBGCO,  pW1CO$BGX)
    passBGS = c(passBGS,  pW1$BGS)
    passExcessX  = c(passExcessX,pW1$ExcessX)
    passExcessS  = c(passExcessS,pW1$ExcessS)
    passERszero = c( passERszero,pW1$intercept)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passERsCO = c(passERsCO, pW1CO$slope)
    passERsCOint = c(passERsCOint, pW1CO$ERint)
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)
    passERsCOzero = c( passERsCOzero,pW1CO$intercept)
    passRsqsCO = c(passRsqsCO, pW1CO$r_sq)   
    passsigma = c(passsigma , pW1$uncertainty_slope)   
    passMAX = c(passMAX, maxfun(pass,'tmpy'))
    passOHrate = c(passOHrate, NaN)
    isNMVOC = c(isNMVOC,3)
    kind = c(kind,'halogen')
    
    mWs = c(mWs, mWCH3COOCl)
    nCs  = c(nCs,2)
    if (debugKT == 1){print('Cl2')}
    
    #Cl2_NOAACIMS_VERES #  (15% + 0.4 ppt)
    pass$tmpy = (pass$Cl2_NOAACIMS_VERES )/1E3 # need to convert to ppb
    background$tmpy = (background$Cl2_NOAACIMS_VERES )/1E3 # need to convert to ppb
    pass$tmperror = (pass$Cl2_NOAACIMS_VERES * .15 + .4)/1E3
    passErrors=c(passErrors,pass$tmperror)
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW)  
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)
    passERsintfill = c(passERsintfill, pW1$ERintfill)
    passBGX = c(passBGX,  pW1$BGX)
    passBGCO = c(passBGCO,  pW1CO$BGX)
    passBGS = c(passBGS,  pW1$BGS)
    passExcessX  = c(passExcessX,pW1$ExcessX)
    passExcessS  = c(passExcessS,pW1$ExcessS)
    passERszero = c( passERszero,pW1$intercept)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passERsCO = c(passERsCO, pW1CO$slope)
    passERsCOint = c(passERsCOint, pW1CO$ERint)
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)
    passERsCOzero = c( passERsCOzero,pW1CO$intercept)
    passRsqsCO = c(passRsqsCO, pW1CO$r_sq)   
    passsigma = c(passsigma , pW1$uncertainty_slope)   
    passOHrate = c(passOHrate, NaN)
    isNMVOC = c(isNMVOC,3)
    kind = c(kind,'halogen')
    
    passMAX = c(passMAX, maxfun(pass,'tmpy'))
    mWs = c(mWs, mWCl2)
    nCs  = c(nCs,0)
    if (debugKT == 1){print('ClNO2')}
    
    #ClNO2_NOAACIMS_VERES #(15% + 0.05 pptv)
    pass$tmpy = (pass$ClNO2_NOAACIMS_VERES )/1E3 # need to convert to ppb
    background$tmpy = (background$ClNO2_NOAACIMS_VERES )/1E3 # need to convert to ppb
    pass$tmperror = (pass$ClNO2_NOAACIMS_VERES * .15 + .05)/1E3
    passErrors=c(passErrors,pass$tmperror)
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW)  
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)
    passERsintfill = c(passERsintfill, pW1$ERintfill)
    passBGX = c(passBGX,  pW1$BGX)
    passBGCO = c(passBGCO,  pW1CO$BGX)
    passBGS = c(passBGS,  pW1$BGS)
    passExcessX  = c(passExcessX,pW1$ExcessX)
    passExcessS  = c(passExcessS,pW1$ExcessS)
    passERszero = c( passERszero,pW1$intercept)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passERsCO = c(passERsCO, pW1CO$slope)
    passERsCOint = c(passERsCOint, pW1CO$ERint)
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)
    passERsCOzero = c( passERsCOzero,pW1CO$intercept)
    passRsqsCO = c(passRsqsCO, pW1CO$r_sq)   
    passsigma = c(passsigma , pW1$uncertainty_slope)   
    passOHrate = c(passOHrate, NaN)
    isNMVOC = c(isNMVOC,3)
    kind = c(kind,'halogen')
    
    passMAX = c(passMAX, maxfun(pass,'tmpy'))
    mWs = c(mWs, mWClNO2)
    nCs  = c(nCs,0)
    
    if (debugKT == 1){print('HCN')}
    #HCN_NOAACIMS_VERES#(15% + 50 ppt)
    pass$tmpy = (pass$HCN_NOAACIMS_VERES )/1E3 # need to convert to ppb
    background$tmpy = (background$HCN_NOAACIMS_VERES )/1E3 # need to convert to ppb
    pass$tmperror = (pass$HCN_NOAACIMS_VERES * .15 + 50)/1E3
    passErrors=c(passErrors,pass$tmperror)
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW)  
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)
    passERsintfill = c(passERsintfill, pW1$ERintfill)
    passBGX = c(passBGX,  pW1$BGX)
    passBGCO = c(passBGCO,  pW1CO$BGX)
    passBGS = c(passBGS,  pW1$BGS)
    passExcessX  = c(passExcessX,pW1$ExcessX)
    passExcessS  = c(passExcessS,pW1$ExcessS)
    passERszero = c( passERszero,pW1$intercept)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passERsCO = c(passERsCO, pW1CO$slope)
    passERsCOint = c(passERsCOint, pW1CO$ERint)
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)
    passERsCOzero = c( passERsCOzero,pW1CO$intercept)
    passRsqsCO = c(passRsqsCO, pW1CO$r_sq)   
    passsigma = c(passsigma , pW1$uncertainty_slope)   
    passOHrate = c(passOHrate, NaN)
    isNMVOC = c(isNMVOC,2)
    kind = c(kind,'nitrogen')
    
    passMAX = c(passMAX, maxfun(pass,'tmpy'))
    mWs = c(mWs, mWHCN)
    nCs  = c(nCs,1)
    
    if (debugKT == 1){print('HCOOH')}
    #HCOOH_NOAACIMS_VERES #(15% + 30 ppt)
    pass$tmpy = (pass$HCOOH_NOAACIMS_VERES )/1E3 # need to convert to ppb
    background$tmpy = (background$HCOOH_NOAACIMS_VERES )/1E3 # need to convert to ppb
    pass$tmperror = (pass$HCOOH_NOAACIMS_VERES * .15 + 30)/1E3
    passErrors=c(passErrors,pass$tmperror)
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW)  
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)
    passERsintfill = c(passERsintfill, pW1$ERintfill)
    passBGX = c(passBGX,  pW1$BGX)
    passBGCO = c(passBGCO,  pW1CO$BGX)
    passBGS = c(passBGS,  pW1$BGS)
    passExcessX  = c(passExcessX,pW1$ExcessX)
    passExcessS  = c(passExcessS,pW1$ExcessS)
    passERszero = c( passERszero,pW1$intercept)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passERsCO = c(passERsCO, pW1CO$slope)
    passERsCOint = c(passERsCOint, pW1CO$ERint)
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)
    passERsCOzero = c( passERsCOzero,pW1CO$intercept)
    passRsqsCO = c(passRsqsCO, pW1CO$r_sq)   
    passsigma = c(passsigma , pW1$uncertainty_slope)   
    passOHrate = c(passOHrate,OHhcooh)
    isNMVOC = c(isNMVOC,1)
    kind = c(kind,'oVOC')
    
    passMAX = c(passMAX, maxfun(pass,'tmpy'))
    mWs = c(mWs, mWHCOOH)
    nCs  = c(nCs,1)
    
    if (debugKT ==1){ 'HNCO'}
    #HNCO_NOAACIMS_VERES #(30% + 16 ppt)
    pass$tmpy = (pass$HNCO_NOAACIMS_VERES )/1E3 # need to convert to ppb
    background$tmpy = (background$HNCO_NOAACIMS_VERES )/1E3 # need to convert to ppb
    pass$tmperror = (pass$HNCO_NOAACIMS_VERES * .3 + 16)/1E3
    passErrors=c(passErrors,pass$tmperror)
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW)  
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)
    passERsintfill = c(passERsintfill, pW1$ERintfill)
    passBGX = c(passBGX,  pW1$BGX)
    passBGCO = c(passBGCO,  pW1CO$BGX)
    passBGS = c(passBGS,  pW1$BGS)
    passExcessX  = c(passExcessX,pW1$ExcessX)
    passExcessS  = c(passExcessS,pW1$ExcessS)
    passERszero = c( passERszero,pW1$intercept)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passERsCO = c(passERsCO, pW1CO$slope)
    passERsCOint = c(passERsCOint, pW1CO$ERint)
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)
    passERsCOzero = c( passERsCOzero,pW1CO$intercept)
    passRsqsCO = c(passRsqsCO, pW1CO$r_sq)   
    passsigma = c(passsigma , pW1$uncertainty_slope)   
    passOHrate = c(passOHrate,OHhnco)
    isNMVOC = c(isNMVOC,2)
    kind = c(kind,'nitrogen')
    
    passMAX = c(passMAX, maxfun(pass,'tmpy'))
    mWs = c(mWs, mWHNCO)
    nCs  = c(nCs,1)
    
    if (debugKT ==1){ 'HPMTF'}
    #HPMTF_NOAACIMS_VERES#(55% + 3 pptv)
    pass$tmpy = (pass$HPMTF_NOAACIMS_VERES )/1E3 # need to convert to ppb
    background$tmpy = (background$HPMTF_NOAACIMS_VERES )/1E3 # need to convert to ppb
    pass$tmperror = (pass$HPMTF_NOAACIMS_VERES * .55 + 3)/1E3
    passErrors=c(passErrors,pass$tmperror)
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW)  
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)
    passERsintfill = c(passERsintfill, pW1$ERintfill)
    passBGX = c(passBGX,  pW1$BGX)
    passBGCO = c(passBGCO,  pW1CO$BGX)
    passBGS = c(passBGS,  pW1$BGS)
    passExcessX  = c(passExcessX,pW1$ExcessX)
    passExcessS  = c(passExcessS,pW1$ExcessS)
    passERszero = c( passERszero,pW1$intercept)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passERsCO = c(passERsCO, pW1CO$slope)
    passERsCOint = c(passERsCOint, pW1CO$ERint)
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)
    passERsCOzero = c( passERsCOzero,pW1CO$intercept)
    passRsqsCO = c(passRsqsCO, pW1CO$r_sq)   
    passsigma = c(passsigma , pW1$uncertainty_slope)   
    passMAX = c(passMAX, maxfun(pass,'tmpy'))
    passOHrate = c(passOHrate,OHHPMTF) 
    isNMVOC = c(isNMVOC,6)
    kind = c(kind,'sulfur')
    
    mWs = c(mWs, mWHPMTF)
    nCs  = c(nCs, 2)
    names1=c(names1,   'HNO2_NOAACIMS_VERES','N2O5_NOAACIMS_VERES','BrCN_NOAACIMS_VERES',
             'BrCl_NOAACIMS_VERES', 'BrO_NOAACIMS_VERES', 'CH3COOCl_NOAACIMS_VERES',
             'Cl2_NOAACIMS_VERES',  'ClNO2_NOAACIMS_VERES', 'HCN_NOAACIMS_VERES',
             'HCOOH_NOAACIMS_VERES', 'HNCO_NOAACIMS_VERES','HPMTF_NOAACIMS_VERES')
    formulas = c(formulas,"HNO2","N2O5","BrCN","BrCl", 'BrO', 'CH3COOCl', 'Cl2',  'ClNO2', 'HCN','HCOOH', 'HNCO','C2H4O3S')
    usenames = c(usenames,"Nitrous acid","Dinitrogen pentoxide","Cyanogen bromide","Bromine chloride",
                 "Bromine monoxide","Chloroacetic acid","Chlorine","Nitryl chloride","Hydrogen cyanide","Formic acid",
                 "Isocyanic acid","Hydroperoxymethyl thioformate")
    
    # # ------------------------------------------ WISTHALER -----------
    # if (debugKT == 1){print('WISTHALER')}
    # pass$tmperror =0.5  # actually a column value
    # passErrors=c(passErrors,pass$tmperror)
    # pass$tmpy = pass$NH3_UIOPTR_ppbV_WISTHALER
    # background$tmpy = background$NH3_UIOPTR_ppbV_WISTHALER
    # pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW)  
    # pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    # passERs = c(passERs, pW1$slope)
    # passERsint = c(passERsint, pW1$ERint)
    # passERsintfill = c(passERsintfill, pW1$ERintfill)
    # passBGX = c(passBGX,  pW1$BGX)
    # passBGCO = c(passBGCO,  pW1CO$BGX)
    # passBGS = c(passBGS,  pW1$BGS)
    # passExcessX  = c(passExcessX,pW1$ExcessX)
    # passExcessS  = c(passExcessS,pW1$ExcessS)
    # passERszero = c( passERszero,pW1$intercept)
    # passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    # passERsCO = c(passERsCO, pW1CO$slope)
    # passERsCOint = c(passERsCOint, pW1CO$ERint)
    # passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)
    # passERsCOzero = c( passERsCOzero,pW1CO$intercept)
    # passRsqsCO = c(passRsqsCO, pW1CO$r_sq)   
    # passsigma = c(passsigma , pW1$uncertainty_slope)   
    # passMAX = c(passMAX, maxfun(pass,'tmpy'))
    # passOHrate = c(passOHrate,NaN)
    # isNMVOC = c(isNMVOC,2)
    # kind = c(kind,'NH3')
    # 
    # names1=c( names1,'NH3_WISTHALER')
    # formulas=c(formulas,"NH3")
    # usenames = c(usenames, "Ammonia")
    # mWs   = c(mWs, mWNH3)
    # nCs  = c(nCs, 0)
    # 
    # # nn = colnames(pass)
    # # ind1A = which(nn == 'CH2O_UIOPTR_ppbV')
    # # ind1B = which(nn == 'CH2O_UIOPTR_ppbV_WISTHALEr',SLOW)
    # # ind1 = max(c(ind1A,ind1B))
    # # if (debugKT == 1){(ind1)}
    # # pass$tmperror = (pass[,ind1] * .25) # Volume mixing ratios accuracy +/- 25 %
    # # pW1 = smallfun(pass,background,xspecies,nn[ind1],      xerror, 'tmperror',SLOW) 
    # # passERs = c(passERs, pW1$slope)
    # # passERsCO = c(passERsCO, pW1CO$slope)
    # # passERsint = c(passERsint, pW1$ERint)
    # # passERsintfill = c(passERsintfill, pW1$ERintfill)
    # #passERsCOint = c(passERsCOint, pW1CO$ERint)
    # #passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)
    # # passBGX = c(passBGX,  pW1$BGX)
    # #passBGCO = c(passBGCO,  pW1CO$BGX)
    # # passBGS = c(passBGS,  pW1$BGS)
    # # passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    # # pW1 = smallfun(pass,background,'CO_DACOM_DISKIN',nn[ind1],      xerror, 'tmperror',SLOW) 
    # # passRsqsCO = c(passRsqsCO, pW1$r_sq)   
    # # passExcessX  = c(passExcessX,pW1$ExcessX)
    # # passExcessS  = c(passExcessS,pW1$ExcessS)
    # # passsigma = c(passsigma, pW1$uncertainty_slope)   
    # # passMAX = c(passMAX, maxfun(pass,nn[ind1]))
    # # names1=c( names1, 'CH2O_WISTHALEr',SLOW)
    # # mWs   = c(mWs, mWCH2O)
    # # nCs  = c(nCs, 1)
    # 
    # if (debugKT == 1){print(c(length(passERs), length(passERsint), length(passERsintfill),length(passERsCO),length(passsigma), length(passRsqsCO2), length(passRsqsCO),length(passMAX),length(passOHrate),
    #                           length(names1), length(mWs), length(nCs), length(passBGS), length(passBGX), 
    #                           length(passExcessX), length(passExcessS), length(passERsCOint), length(passERsCOintfill),
    #                           length(passERszero), length(passERsCOzero), length(isNMVOC), length(kind)))}
    # # ------------------------------------------ BLAKE ------------------------------------------
    if (debugKT == 1){print('BLAKE')}
    bnames=c('OCS_WAS_BLAKE', 'DMS_WAS_BLAKE', 'CFC12_WAS_BLAKE', 'CFC11_WAS_BLAKE', 'CFC113_WAS_BLAKE', 'CFC114_WAS_BLAKE', 
             'HFC152a_WAS_BLAKE', 'HFC134a_WAS_BLAKE', 'HFC365mfc_WAS_BLAKE', 'HCFC22_WAS_BLAKE', 'HCFC142b_WAS_BLAKE', 'HCFC141b_WAS_BLAKE', 
             'H1301_WAS_BLAKE', 'H2402_WAS_BLAKE', 'H1211_WAS_BLAKE', 'CH3CCl3_WAS_BLAKE', 'CCl4_WAS_BLAKE', 'CHCl3_WAS_BLAKE', 'CH2Cl2_WAS_BLAKE', 
             'C2HCl3_WAS_BLAKE', 'C2Cl4_WAS_BLAKE', 'CH3Cl_WAS_BLAKE', 'CH3Br_WAS_BLAKE', 'CH3I_WAS_BLAKE', 'CH2Br2_WAS_BLAKE', 'CHBrCl2_WAS_BLAKE', 
             'CHBr2Cl_WAS_BLAKE', 'CHBr3_WAS_BLAKE', 'CH2ClCH2Cl_WAS_BLAKE', 'C2H5Cl_WAS_BLAKE', 'MeONO2_WAS_BLAKE', 'EthONO2_WAS_BLAKE', 
             'iPropONO2_WAS_BLAKE', 'nPropONO2_WAS_BLAKE', 'x2ButONO2_WAS_BLAKE', 'x3PentONO2_WAS_BLAKE', 'x2PentONO2_WAS_BLAKE', 
             'x3Me2ButONO2_WAS_BLAKE', 'Ethane_WAS_BLAKE', 'Ethene_WAS_BLAKE', 'Ethyne_WAS_BLAKE', 'Propene_WAS_BLAKE', 'Propane_WAS_BLAKE', 
             'Propadiene_WAS_BLAKE', 'Propyne_WAS_BLAKE', 'iButane_WAS_BLAKE', 'nButane_WAS_BLAKE', 'x1Butene_WAS_BLAKE', 'iButene_WAS_BLAKE', 
             't2Butene_WAS_BLAKE', 'c2Butene_WAS_BLAKE', 'x13Butadiene_WAS_BLAKE', 'x12Butadiene_WAS_BLAKE', 'x1Buten3yne_WAS_BLAKE', 
             'x13Butadyine_WAS_BLAKE', 'x1Butyne_WAS_BLAKE', 'x2Butyne_WAS_BLAKE', 'iPentane_WAS_BLAKE', 'nPentane_WAS_BLAKE', 
             'Isoprene_WAS_BLAKE', 'x1Pentene_WAS_BLAKE', 't2Pentene_WAS_BLAKE', 'c2Pentene_WAS_BLAKE', 'X3Me1Butene_WAS_BLAKE', 
             'X2Me1Butene_WAS_BLAKE', 'X2Me2Butene_WAS_BLAKE', 'x13Pentadienes_WAS_BLAKE', 'x3Me1PenteneAnd4Me1Pentene_WAS_BLAKE', 
             'x1Hexene_WAS_BLAKE', 'x1Heptene_WAS_BLAKE', 'x1Octene_WAS_BLAKE', 'x1Nonene_WAS_BLAKE', 'x1Decene_WAS_BLAKE', 'nHexane_WAS_BLAKE', 
             'nHeptane_WAS_BLAKE', 'nOctane_WAS_BLAKE', 'nNonane_WAS_BLAKE', 'nDecane_WAS_BLAKE', 'nUndecane_WAS_BLAKE', 'x22Dimebutane_WAS_BLAKE', 
             'x23Dimebutane_WAS_BLAKE', 'x2MePentane_WAS_BLAKE', 'x3MePentane_WAS_BLAKE', 'x2MeHexane_WAS_BLAKE', 'x3MeHexane_WAS_BLAKE', 
             'x23DimePentane_BLAKE', 'x224TrimePentane_WAS_BLAKE', 'x234TrimePentane_WAS_BLAKE', 'CycPentane_WAS_BLAKE', 
             'MeCycPentane_WAS_BLAKE', 'CycHexane_WAS_BLAKE', 'MeCycHexane_WAS_BLAKE', 'CycPentene_WAS_BLAKE', 'Benzene_WAS_BLAKE', 
             'Toluene_WAS_BLAKE', 'EthBenzene_WAS_BLAKE', 'mpXylene_WAS_BLAKE', 'oXylene_WAS_BLAKE', 'Styrene_WAS_BLAKE', 
             'EthynylBenzene_WAS_BLAKE', 'iPropBenzene_WAS_BLAKE', 'nPropBenzene_WAS_BLAKE', 'x3EthToluene_WAS_BLAKE', 'x4EthToluene_WAS_BLAKE', 
             'x2EthToluene_WAS_BLAKE', 'x135rimeBenzene_WAS_BLAKE', 'x124rimeBenzene_WAS_BLAKE', 'ClBenzene_WAS_BLAKE', 'aPinene_WAS_BLAKE', 
             'bPinene_WAS_BLAKE', 'Tricyclene_WAS_BLAKE', 'Camphene_WAS_BLAKE', 'Myrcene_WAS_BLAKE', 'Limonene_WAS_BLAKE', 'Furan_WAS_BLAKE', 
             'x2MeFuran_WAS_BLAKE', 'x3MeFuran_WAS_BLAKE', 'BenzFuran_WAS_BLAKE', 'iButanal_WAS_BLAKE', 'Butanal_WAS_BLAKE', 
             'AcetonePropanal_WAS_BLAKE', 'MEK_WAS_BLAKE', 'MAC_WAS_BLAKE', 'MVK_WAS_BLAKE', 'Acrolein_WAS_BLAKE', 'iPropanol_WAS_BLAKE', 
             'Nitromethane_WAS_BLAKE', 'Acrylonitrile_WAS_BLAKE', 'PropNitrile_WAS_BLAKE', 'MeAcetate_WAS_BLAKE')
    kind=c(kind,'sulfur', 'sulfur', 'halocarbon', 'halocarbon', 'halocarbon', 'halocarbon', 
             'halocarbon', 'halocarbon', 'halocarbon', 'halocarbon', 'halocarbon', 'halocarbon', 
             'halocarbon', 'halocarbon', 'halocarbon', 'halocarbon', 'halocarbon', 'halocarbon', 'halocarbon', 
             'halocarbon', 'halocarbon', 'halocarbon', 'halocarbon', 'halocarbon', 'halocarbon', 'halocarbon', 
             'halocarbon', 'halocarbon', 'halocarbon', 'halocarbon', 'alkyl nitrate', 'alkyl nitrate', 
             'alkyl nitrate', 'alkyl nitrate', 'alkyl nitrate', 'alkyl nitrate', 'alkyl nitrate', 
             'alkyl nitrate', 'alkane', 'alkene', 'alkyne', 'alkene', 'alkane', 
             'alkene', 'alkyne', 'alkane', 'alkane', 'alkene', 'alkene', 
             'alkene', 'alkene', 'alkene', 'alkene', 'alkyne', 
             'alkyne', 'alkyne', 'alkyne', 'alkane', 'alkane', 
             'alkene', 'alkene', 'alkene', 'alkene', 'alkene', 
             'alkene', 'alkene', 'alkene', 'alkene', 
             'alkene', 'alkene', 'alkene', 'alkene', 'alkene', 'alkane', 
             'alkane', 'alkane', 'alkane', 'alkane', 'alkane', 'alkane', 
             'alkane', 'alkane', 'alkane', 'alkane', 'alkane', 
             'alkane', 'alkane', 'alkane', 'alkane', 
             'alkane', 'aromatic', 'aromatic', 'alkene', 'aromatic', 
             'aromatic', 'aromatic', 'aromatic', 'aromatic', 'aromatic', 
             'aromatic', 'aromaticE', 'aromatic', 'aromatic', 'aromatic', 
             'aromatic', 'aromatic', 'aromatic', 'aromatic', 'alkene', 
             'alkene', 'alkene', 'alkene', 'alkene', 'alkene', 'oVOC', 
             'oVOC', 'oVOC', 'oVOC', 'oVOC', 'oVOC', 
             'oVOC', 'oVOC', 'oVOC', 'oVOC', 'oVOC', 'oVOC', 
             'nitrogen', 'nitrogen', 'nitrogen', 'oVOC')
    berrors=c(0.10, 0.20, 0.2, 0.2, 0.2, 0.2, 0.10, 0.5, 0.15, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.10, 0.10, 0.20, 0.10, 0.5, 0.10, 0.10, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.25, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.25, 0.25, 0.25, 0.50, 0.40, 0.40, 0.30, 0.30, 0.30, 0.30, 0.30, 0.20, 0.30, 0.40, 0.40, 0.40)
    passErrors=c(passErrors,berrors)
    bMW=c(mWOCS,  mWDMS, mWCFC12, mWCFC11, mWCFC113, mWCFC114, 
          mWHFC152a, mWHFC134a, mWHFC365mfc, mWHCFC22, mWHCFC142b, mWHCFC141b,
          mWH1301, mWH2402, mWH1211, mWCH3CCl3, mWCCl4, mWCHCl3, mWCH2Cl2,
          mWC2HCl3, mWC2Cl4, mWCH3Cl, mWCH3Br, mWCH3I, mWCH2Br2, mWCHBrCl2,
          mWCHBr2Cl, mWCHBr3, mWCH2ClCH2Cl, mWC2H5Cl, mWMeONO2, mWEthONO2,
          mWiPropONO2, mWnPropONO2, mWx2ButONO2, mWx3PentONO2, mWx2PentONO2,
          mWx3Me2ButONO2, mWEthane, mWEthene, mWEthyne, mWPropene, mWPropane,
          mWPropadiene, mWPropyne, mWiButane, mWnButane, mWx1Butene, mWiButene,
          mWt2Butene, mWc2Butene, mWx13Butadiene, mWx12Butadiene, mWx1Buten3yne,
          mWx13Butadyine, mWx1Butyne, mWx2Butyne, mWiPentane, mWnPentane,
          mWIsoprene, mWx1Pentene, mWt2Pentene, mWc2Pentene, mWX3Me1Butene,
          mWX2Me1Butene, mWX2Me2Butene, mWx13Pentadienes, mWx3Me1PenteneAnd4Me1Pentene,
          mWx1Hexene, mWx1Heptene, mWx1Octene, mWx1Nonene, mWx1Decene, mWnHexane,
          mWnHeptane, mWnOctane, mWnNonane, mWnDecane, mWnUndecane, mWx22Dimebutane,
          mWx23Dimebutane, mWx2MePentane, mWx3MePentane, mWx2MeHexane, mWx3MeHexane,
          mWx23DimePentane, mWx224TrimePentane, mWx234TrimePentane, mWCycPentane,
          mWMeCycPentane, mWCycHexane, mWMeCycHexane, mWCycPentene, mWBenzene,
          mWToluene, mWEthBenzene, mWmpXylene, mWoXylene, mWStyrene,
          mWEthynylBenzene, mWiPropBenzene, mWnPropBenzene, mWx3EthToluene, mWx4EthToluene,
          mWx2EthToluene, mWx135rimeBenzene, mWx124rimeBenzene, mWClBenzene, mWaPinene,
          mWbPinene, mWTricyclene, mWCamphene, mWMyrcene, mWLimonene, mWFuran,
          mWx2MeFuran, mWx3MeFuran, mWBenzFuran, mWiButanal, mWButanal,
          mWAcetonePropanal, mWMEK, mWMAC, mWMVK, mWAcrolein, mWiPropanol,
          mWNitromethane, mWAcrylonitrile, mWPropNitrile, mWMeAcetate)
    bNC=c(nCOCS,  nCDMS, nCCFC12, nCCFC11, nCCFC113, nCCFC114, 
          nCHFC152a, nCHFC134a, nCHFC365mfc, nCHCFC22, nCHCFC142b, nCHCFC141b,
          nCH1301, nCH2402, nCH1211, nCCH3CCl3, nCCCl4, nCCHCl3, nCCH2Cl2,
          nCC2HCl3, nCC2Cl4, nCCH3Cl, nCCH3Br, nCCH3I, nCCH2Br2, nCCHBrCl2,
          nCCHBr2Cl, nCCHBr3, nCCH2ClCH2Cl, nCC2H5Cl, nCMeONO2, nCEthONO2,
          nCiPropONO2, nCnPropONO2, nCx2ButONO2, nCx3PentONO2, nCx2PentONO2,
          nCx3Me2ButONO2, nCEthane, nCEthene, nCEthyne, nCPropene, nCPropane,
          nCPropadiene, nCPropyne, nCiButane, nCnButane, nCx1Butene, nCiButene,
          nCt2Butene, nCc2Butene, nCx13Butadiene, nCx12Butadiene, nCx1Buten3yne,
          nCx13Butadyine, nCx1Butyne, nCx2Butyne, nCiPentane, nCnPentane,
          nCIsoprene, nCx1Pentene, nCt2Pentene, nCc2Pentene, nCX3Me1Butene,
          nCX2Me1Butene, nCX2Me2Butene, nCx13Pentadienes, nCx3Me1PenteneAnd4Me1Pentene,
          nCx1Hexene, nCx1Heptene, nCx1Octene, nCx1Nonene, nCx1Decene, nCnHexane,
          nCnHeptane, nCnOctane, nCnNonane, nCnDecane, nCnUndecane, nCx22Dimebutane,
          nCx23Dimebutane, nCx2MePentane, nCx3MePentane, nCx2MeHexane, nCx3MeHexane,
          nCx23DimePentane, nCx224TrimePentane, nCx234TrimePentane, nCCycPentane,
          nCMeCycPentane, nCCycHexane, nCMeCycHexane, nCCycPentene, nCBenzene,
          nCToluene, nCEthBenzene, nCmpXylene, nCoXylene, nCStyrene,
          nCEthynylBenzene, nCiPropBenzene, nCnPropBenzene, nCx3EthToluene, nCx4EthToluene,
          nCx2EthToluene, nCx135rimeBenzene, nCx124rimeBenzene, nCClBenzene, nCaPinene,
          nCbPinene, nCTricyclene, nCCamphene, nCMyrcene, nCLimonene, nCFuran,
          nCx2MeFuran, nCx3MeFuran, nCBenzFuran, nCiButanal, nCButanal,
          nCAcetonePropanal, nCMEK, nCMAC, nCMVK, nCAcrolein, nCiPropanol,
          nCNitromethane, nCAcrylonitrile, nCPropNitrile, nCMeAcetate)
    #START
    passOHrate = c(passOHrate,OHOCS,  OHDimeFurans, OHCFC12, OHCFC11, OHCFC113, OHCFC114, 
                   OHHFC152a, OHHFC134a, OHHFC365mfc, OHHCFC22, OHHCFC142b, OHHCFC141b,
                   OHH1301, OHH2402, OHH1211, OHCH3CCl3, OHCCl4, OHCHCl3, OHCH2Cl2,
                   OHC2HCl3, OHC2Cl4, OHCH3Cl, OHCH3Br, OHCH3I, OHCH2Br2, OHCHBrCl2,
                   OHCHBr2Cl, OHCHBr3, OHCH2ClCH2Cl, OHC2H5Cl, OHMeONO2, OHEthONO2,
                   OHiPropONO2, OHnPropONO2, OHx2ButONO2, OHx3PentONO2, OHx2PentONO2,
                   OHx3Me2ButONO2, OHc2h6, OHEthene, OHEthyne, OHPropene, OHPropane,
                   OHPropadiene, OHPropyne, OHiButane, OHnButane, OHx1Butene, OHiButene,
                   OHt2Butene, OHc2Butene, OHx13Butadiene, OHx12Butadiene, OHx1Buten3yne,
                   OHx13Butadyine, OHx1Butyne, OHx2Butyne, OHiPentane, OHnPentane,
                   OHIsoprene, OHx1Pentene, OHt2Pentene, OHc2Pentene, OHX3Me1Butene,
                   OHX2Me1Butene, OHX2Me2Butene, OHx13Pentadienes, OHx3Me1PenteneAnd4Me1Pentene,
                   OHx1Hexene, OHx1Heptene, OHx1Octene, OHx1Nonene, OHx1Decene, OHnHexane,
                   OHnHeptane, OHnOctane, OHnNonane, OHnDecane, OHnUndecane, OHx22Dimebutane,
                   OHx23Dimebutane, OHx2MePentane, OHx3MePentane, OHx2MeHexane, OHx3MeHexane,
                   OHx23DimePentane, OHx224TrimePentane, OHx234TrimePentane, OHCycPentane,
                   OHMeCycPentane, OHCycHexane, OHMeCycHexane, OHCycPentene, OHBenzene,
                   OHToluene, OHEthBenzene, OHmpXylene, OHoXylene, OHStyrene,
                   OHEthynylBenzene, OHiPropBenzene, OHnPropBenzene, OHx3EthToluene, OHx4EthToluene,
                   OHx2EthToluene, OHx135rimeBenzene, OHx124rimeBenzene, OHClBenzene, OHaPinene,
                   OHbPinene, OHTricyclene, OHCamphene, OHMyrcene, OHLimonene, OHFuran,
                   OHx2MeFuran, OHx3MeFuran, OHBenzFuran, OHiButanal, OHButanal,
                   OHAcetonePropanal, OHMEK, OHMAC, OHMVK,OHAcrolein, OHiPropanol,
                   OHNitromethane, OHAcrylonitrile, OHPropNitrile, OHMeAcetate)
    if (doBlake == 1){
      for (i in 1:length(bnames)){
        ccs = colnames(blake) ; ccs2=colnames(blakeBG)
        #print(bnames[i])
        # if (debugKT == 1){print(bnames[i])}
        ind = which(ccs == bnames[i]) ; ind2 = which(ccs2 == bnames[i])
        # everything in ppb
        t1 = blake[,ind]/1E3 ; t1B = blakeBG[,ind2]/1E3
        # Assume if a species background is NaN we can set it to zero
        if (is.nan(t1B)){t1B = 0}
        
        t1CO = blake$CO_DACOM_DISKIN_BLAKE; t1COB = blakeBG$CO_DACOM_DISKIN_BLAKE
        t1CO2 = blake$CO2_7000_ppm_DISKIN_BLAKE*1E3; t1CO2B = blakeBG$CO2_7000_ppm_DISKIN_BLAKE*1E3
        
        # Are we using the fast background?
        if (useFASTBG  == 1 ){
          t1COB = mean(background$CO_DACOM_DISKIN, na.rm=TRUE)
          t1CO2B = mean(background$CO2_7000_ppm_DISKIN*1E3)
        }
        
        if (i == 1){
          passERsB = sum(t1 - t1B)/sum(t1CO2 - t1CO2B)
          passERsBCO = sum(t1 - t1B)/sum(t1CO - t1COB)
          
          passERsBI  = NaN
          passERsBIfill  = NaN
          passERsCOBI   = NaN
          passERsCOBIfill = NaN
          passRBCO2 =t1CO2B
          passRBCOB  = t1COB
          passRBS =t1B
          passRB = 1# so they dont get thrown away by R2 filter
          passRBCO =1 # so they dont get thrown away by R2 filter
          passRBXX  =sum(t1CO2 - t1CO2B)
          passRBSS  =sum(t1 - t1B)
          passSB = NaN
          passSBQ = NaN
          passSBQCO =NaN
          passMB = max(t1)
        }
        if (i > 1){ 
          passERsB = c(passERsB,sum(t1 - t1B)/sum(t1CO2 - t1CO2B))
          passERsBCO = c(passERsBCO,sum(t1 - t1B)/sum(t1CO - t1COB))
          passERsBI = c(passERsBI,NaN)
          passERsBIfill = c(passERsBIfill,NaN)
          passERsCOBI = c(passERsCOBI, NaN)
          passERsCOBIfill = c(passERsCOBIfill,NaN)
          passRBCO2 = c(passRBCO2,t1CO2B)
          passRBCOB  = c(passRBCOB,t1COB)
          passRBS = c(passRBS,t1B)
          
          passRB = c(passRB,1)# so they dont get thrown away by R2 filter
          passRBCO =c(passRBCO,1) # so they dont get thrown away by R2 filter
          passRBXX  = c(passRBXX,sum(t1CO2 - t1CO2B))
          passRBSS  = c(passRBSS,sum(t1 - t1B))
          passSBQ = c( passSBQ,NaN)
          passSBQCO = c(passSBQCO,NaN)
          passSB = c(passSB,NaN)
          passMB = c(passMB, max(t1))
        }
      }
      passERs = c(passERs, passERsB)  
      passERsCO = c(passERsCO,passERsBCO) 
      
      passERsint = c(passERsint, passERsBI)    
      passERsintfill = c(passERsintfill, passERsBIfill)   
      
      passERszero   = c(passERszero, passSBQ)
      passERsCOzero = c(passERsCOzero, passSBQCO)
      
      passERsCOint = c(passERsCOint, passERsCOBI)
      passERsCOintfill = c(passERsCOintfill, passERsCOBIfill)
      
      passBGX = c(passBGX,  passRBCO2)
      passBGCO = c(passBGCO,passRBCOB)
      passBGS = c(passBGS,  passRBS)
      passRsqsCO2 = c(passRsqsCO2, passRB)   
      passRsqsCO  = c(passRsqsCO, passRBCO)   
      passExcessX  = c(passExcessX,passRBXX)
      passExcessS  = c(passExcessS,passRBSS)
      
      passsigma = c(passsigma, passSB)
      passMAX = c(passMAX, passMB)
    } else{
      passERs = c(passERs, rep(NaN, length(bnames)))
      passERsCO = c(passERsCO,rep(NaN, length(bnames)))
      
      passERsint = c(passERsint,rep(NaN, length(bnames)))
      passERsintfill = c(passERsintfill,rep(NaN, length(bnames)))
      
      passERszero   = c(passERszero,rep(NaN, length(bnames)))
      passERsCOzero = c(passERsCOzero,rep(NaN, length(bnames)))
      
      passERsCOint = c(passERsCOint,rep(NaN, length(bnames)))
      passERsCOintfill = c(passERsCOintfill,rep(NaN, length(bnames)))
      
      passBGX = c(passBGX,  rep(NaN, length(bnames)))
      passBGCO = c(passBGCO,rep(NaN, length(bnames)))
      passBGS = c(passBGS,  rep(NaN, length(bnames)))
      passRsqsCO2 = c(passRsqsCO2,rep(NaN, length(bnames))) 
      passRsqsCO  = c(passRsqsCO,rep(NaN, length(bnames)))
      passExcessX  = c(passExcessX,rep(NaN, length(bnames)))
      passExcessS  = c(passExcessS,rep(NaN, length(bnames)))
      passsigma = c(passsigma,rep(NaN, length(bnames)))
      passMAX = c(passMAX, rep(NaN, length(bnames)))
    }
    isNMVOC=c(isNMVOC,6, 6, 3, 3, 3,3, 3,3, 3, 3, 3, 3, 
              3, 3, 3, 3, 3,3, 3,
              3, 3,3, 3,  3, 3,3, 3,3,3, 3,
              2,2, 2,2, 2,2,  2,2,
              1, 1,  1,  1, 1, 
              1, 1,  1,  1, 1, 1,  1,  1,1, 1,  1,  1, 1, 1,  1,  1,
              1, 1,  1,  1, 1, 1,  1,  1,1, 1,  1,  1, 1, 1,  1,  1,1, 1,  1,  1,
              1, 1,  1,  1,1, 1,  1,  1,1, 1,  1,  1, 1, 1,  1,  1, 1, 1,  1,  1,
              1, 1,  1,  1,1, 1,  1,  1,1, 1,  1,  1,1, 1,  1,  1,
              1, 1,  1,  1,1, 1,  1,  1,1, 1,  1,  1, 1, 2,  1)
    names1=c(names1,bnames)
    formulas = c(formulas, 'OCS',	'C2H6S', 'CFC12', 'CFC11', 'CFC113', 'CFC114', 'HFC152a', 'HFC134a', 'HFC365mfc', 'HCFC22', 
    'HCFC142b', 'HCFC141b', 'H1301', 'H2402', 'H1211', 'CH3CCl3', 'CCl4', 'CHCl3', 'CH2Cl2', 'C2HCl3', 'C2Cl4', 'CH3Cl', 
    'CH3Br', 'CH3I', 'CH2Br2', 'CHBrCl2', 'CHBr2Cl', 'CHBr3', 'CH2ClCH2Cl', 'C2H5Cl', 'CH3NO3', 'C2H5NO3', 'C3H7NO3', 
    'C3H7NO3', 'C4H9NO3', 'C5H11NO3', 'C5H11NO3', 'C5H11NO3', 'C2H6', 'C2H4', 'C2H2', 'C3H6', 'C3H8', 'C3H4', 'C3H4',
    'C4H10', 'C4H10', 'C4H8', 'C4H8', 'C4H8', 'C4H8', 'C4H6', 'C4H6', 'C4H4', 'C4H2', 'C4H6', 'C4H6', 'C5H12', 'C5H12', 
    'C5H8', 'C5H10', 'C5H10', 'C5H10', 'C5H10', 'C5H10', 'C5H10', 'C5H8', 'C6H12', 'C6H12', 'C7H14', 'C8H16', 'C9H18', 
    'C10H20', 'C₆H₁₄', 'C7H16', 'C8H18', 'C9H20', 'C10H22', 'C11H24', 'C6H14', 'C6H14', 'C6H14', 'C6H14', 'C7H16', 'C7H16',
    'C7H16', 'C8H18', 'C8H18', 'C5H10', 'C6H12', 'C6H12', 'C7H14', 'C5H8', 'C6H6', 'C7H8', 'C8H10', 'C8H10', 'C8H10', 'C8H8',
    'C8H6', 'C9H12', 'C9H12', 'C9H12', 'C9H12', 'C9H12', 'C9H12', 'C9H12', 'C6H5Cl', 'C10H16', 'C10H16', 'C10H16', 'C10H16', 
    'C10H16', 'C10H16', 'C4H4O', 'C5H6O', 'C5H6O', 'C8H6O', 'C4H8O', 'C4H8O', 'C3H6O', 'C4H8O', 'C4H6O', 'C4H6O', 'C3H4O',
    'C3H8O', 'CH3NO2', 'C3H3N', 'CH3CH2CN','CH3COOCH3')
                 
    usenames = c(usenames,'Carbonyl sulfide', 'Dimethyl sulfide', 'CFC12', 'CFC11', 'CFC113', 'CFC114', 'HFC152a',
    'HFC134a', 'HFC365mfc', 'HCFC22', 'HCFC142b', 'HCFC141b', 'H1301', 'H2402', 'H1211', 'Methyl Chloroform', 
    'Carbon tetrachloride', 'Chloroform', 'Dichloromethane', 'Trichloroethylene', 'Tetrachloroethylene', 'Methyl chloride',
    'Methyl Bromide', 'Methyl Iodide', 'Dibromomethane', 'Bromodichloromethane', 'Chlorodibromomethane', 'Bromoform',
    '1,2-Dichloroethane', 'Chloroethane', 'Methyl nitrate', 'Ethyl nitrate', 'i-propyl nitrate', 'N-propyl nitrate',
    '2-Butyl nitrate', 'Pentan-3-yl nitrate', '2-Pentanyl nitrate', '3-Methyl-2-butanyl nitrate', 'Ethane', 'Ethene', 
    'Ethyne', 'Propene', 'Propane', 'Propadiene', 'Propyne', 'Isobutane', 'n-Butane', '1-Butene', 'i-Butene',
    'trans-2-Butene', 'cis-2-Butene', '1,3-Butadiene', '1,2-Butadiene', '1-Buten-3-yne', '1,3-Butadiyne', '1-Butyne',
    '2-Butyne', 'Isopentane', 'n-Pentane', 'Isoprene', '1-Pentene', 'trans-2-Pentene', 'cis-2-Pentene', '3-Methyl-1-butene',
    '2-Methyl-1-butene', '2-Methyl-2-butene', '1,3-pentadiene', '3,4-Methyl-1-pentene', '1-Hexene', '1-Heptene', '1-Octene',
    '1-Nonene', '1-Decene', 'n-Hexane', 'n-Heptane', 'n-Octane', 'n-Nonane', 'n-Decane', 'n-Undecane', '2,2-Dimethylbutane', 
    '2,3-Dimethylbutane', '2-Methylpentane', '3-Methylpentane', '2-Methylhexane', '3-Methylhexane', '2,3-Dimethylpentane',
    '2,2,4-Trimethylpentane', '2,3,4-Trimethylpentane', 'Cyclopentane', 'Methylcyclopentane', 'Cyclohexane', 
    'Methylcyclohexane', 'Cyclopentene', 'Benzene', 'Toluene', 'Ethylbenzene', 'm,p-Xylene', 'o-Xylene', 'Styrene', 
    'Ethynylbenzene', 'i-Propylbenzene', 'n-Propylbenzene', '3-Ethyltoluene', '4-Ethyltoluene', '2-Ethyltoluene', 
    '1,3,5-Trimethylbenzene', '1,2,4-Trimethylbenzene', 'Chlorobenzene', '⍺-Pinene', 'β-Pinene', 'Tricyclene',
    'Camphene', 'Myrcene', 'Limonene', 'Furan', '2-Methylfuran', '3-Methylfuran', 'Benzofuran', 'Isobutanal', 'Butanal',
    'Acetone/Propanal', 'Methyl ethyl ketone', 'Methacrolein', 'Methyl vinyl ketone', 'Acrolein', 'Isopropanol', 
    'Nitromethane', 'Acrylonitrile', 'Propionitrile','Methyl acetate')
    mWs = c(mWs, bMW) 
    nCs = c(nCs, bNC) 
    # # ------------------------------------------ APEL ----------------------
    # #if (debugKT == 1){print("APEL")}
    # anames=c('HFC134a_TOGA_APEL', 'HCFC141b_TOGA_APEL', 'HCFC142b_TOGA_APEL', 'HCFC22_TOGA_APEL', 'CH2Cl2_TOGA_APEL',
    #          'CHCl3_TOGA_APEL', 'CH2ClCH2Cl_TOGA_APEL', 'CH3CCl3_TOGA_APEL', 'C2Cl4_TOGA_APEL', 'ClBenzene_TOGA_APEL',
    #          'CHBrCl2_TOGA_APEL', 'CHBr2Cl_TOGA_APEL', 'CH3Br_TOGA_APEL', 'CH2Br2_TOGA_APEL', 'CHBr3_TOGA_APEL',
    #          'CH2ClI_TOGA_APEL', 'CH3I_TOGA_APEL', 'CS2_TOGA_APEL', 'CH3SH_TOGA_APEL', 'DMS_TOGA_APEL', 'Propane_TOGA_APEL',
    #          'iButane_TOGA_APEL', 'nButane_TOGA_APEL', 'iPentane_TOGA_APEL', 'nPentane_TOGA_APEL', 'x2MePentane_TOGA_APEL',
    #          'x3MePentane_TOGA_APEL', 'nHexane_TOGA_APEL', 'x224TrimePentane_TOGA_APEL', 'nHeptane_TOGA_APEL',
    #          'nOctane_TOGA_APEL', 'Propene_TOGA_APEL', 'iButene1Butene_TOGA_APEL', 'Isoprene_TOGA_APEL',
    #          'Tricyclene_TOGA_APEL', 'aPinene_TOGA_APEL', 'Camphene_TOGA_APEL', 'bPineneMyrcene_TOGA_APEL',
    #          'LimoneneD3Carene_TOGA_APEL', 'Benzene_TOGA_APEL', 'Toluene_TOGA_APEL', 'EthBenzene_TOGA_APEL',
    #          'mpXylene_TOGA_APEL', 'oXylene_TOGA_APEL', 'Styrene_TOGA_APEL', 'EthynylBenzene_TOGA_APEL', 'CH2O_TOGA_APEL',
    #          'CH3CHO_TOGA_APEL', 'Propanal_TOGA_APEL', 'Butanal_TOGA_APEL', 'iButanal_TOGA_APEL', 'Acrolein_TOGA_APEL',
    #          'x2Butenals_TOGA_APEL', 'Acetone_TOGA_APEL', 'MEK_TOGA_APEL', 'CH3OH_TOGA_APEL', 'C2H5OH_TOGA_APEL',
    #          'iPropanol_TOGA_APEL', 'MBO_TOGA_APEL', 'MAC_TOGA_APEL', 'MVK_TOGA_APEL', 'MeFormate_TOGA_APEL',
    #          'MeAcetate_TOGA_APEL', 'Furan_TOGA_APEL', 'x2MeFuran_TOGA_APEL', 'x3MeFuran_TOGA_APEL', 'Furfural_TOGA_APEL',
    #          'HCN_TOGA_APEL', 'CH3CN_TOGA_APEL', 'PropNitrile_TOGA_APEL', 'Acrylonitrile_TOGA_APEL',
    #          'MeAcrylonitrile_TOGA_APEL', 'Pyrrole_TOGA_APEL', 'Nitromethane_TOGA_APEL', 'MeONO2_TOGA_APEL',
    #          'EthONO2_TOGA_APEL', 'iPropONO2_TOGA_APEL', 'x2ButONO2iButONO2_TOGA_APEL')
    # kind = c(kind, rep("NaN", length(anames)))
    # aMW = c(  mWHFC134a, mWHCFC141b, mWHCFC142b, mWHCFC22, mWCH2Cl2,
    #           mWCHCl3, mWCH2ClCH2Cl, mWCH3CCl3, mWC2Cl4, mWClBenzene,
    #           mWCHBrCl2, mWCHBr2Cl, mWCH3Br, mWCH2Br2, mWCHBr3,
    #           mWCH2ClI, mWCH3I, mWCS2, mWCH3SH, mWDMS, mWPropane,
    #           mWiButane, mWnButane, mWiPentane, mWnPentane, mWx2MePentane,
    #           mWx3MePentane, mWnHexane, mWx224TrimePentane, mWnHeptane,
    #           mWnOctane, mWPropene, mWiButene, mWIsoprene,
    #           mWTricyclene, mWaPinene, mWCamphene, mWbPinene,
    #           mWLimonene, mWBenzene, mWToluene, mWEthBenzene,
    #           mWmpXylene, mWoXylene, mWStyrene, mWEthynylBenzene, mWCH2O,
    #           mWAcetaldehyde, mWpropanal, mWButanal, mWiButanal, mWAcrolein,
    #           mW2Butenals, mWAcetonePropanal, mWMEK, mWCH3OH, mWC2H5OH,
    #           mWiPropanol, mWMBO, mWMAC, mWMVK, mWMeFormate,
    #           mWMeAcetate, mWFuran, mWx2MeFuran, mWx3MeFuran, mWfurfural,
    #           mWHCN, mWch3cn, mWPropNitrile, mWAcrylonitrile,
    #           mWMeAcrylonitrile, mwPyrrole, mWNitromethane, mWMeONO2,
    #           mWEthONO2, mWiPropONO2, MWx2ButONO2iButONO2)
    # aNC =c(  nCHFC134a, nCHCFC141b, nCHCFC142b, nCHCFC22, nCCH2Cl2,
    #          nCCHCl3, nCCH2ClCH2Cl, nCCH3CCl3, nCC2Cl4, nCClBenzene,
    #          nCCHBrCl2, nCCHBr2Cl, nCCH3Br, nCCH2Br2, nCCHBr3,
    #          nCCH2ClI, nCCH3I, nCCS2, nCCH3SH, nCDMS, nCPropane,
    #          nCiButane, nCnButane, nCiPentane, nCnPentane, nCx2MePentane,
    #          nCx3MePentane, nCnHexane, nCx224TrimePentane, nCnHeptane,
    #          nCnOctane, nCPropene, nCiButene, nCIsoprene,
    #          nCTricyclene, nCaPinene, nCCamphene, nCbPinene,
    #          nCLimonene, nCBenzene, nCToluene, nCEthBenzene,
    #          nCmpXylene, nCoXylene, nCStyrene, nCEthynylBenzene, nCCH2O,
    #          nCAcetaldehyde, nCpropanal, nCButanal, nCButanal, nCAcrolein,
    #          nC2Butenals, nCAcetonePropanal, nCMEK, nCCH3OH, nCC2H5OH,
    #          nCiPropanol, nCMBO, nCMAC, nCMVK, nCMeFormate,
    #          nCMeAcetate, nCFuran, nCx2MeFuran, nCx3MeFuran, nCfurfural,
    #          nCHCN, nCch3cn, nCPropNitrile, nCAcrylonitrile,
    #          nCMeAcrylonitrile, nCPyrrole, nCNitromethane, nCMeONO2,
    #          nCEthONO2, nCiPropONO2, nCx2ButONO2iButONO2)
    # isNMVOC=c(isNMVOC,3, 3, 3, 3,3, 3, 3, 3,3, 3,
    #           3, 3, 3, 3,3, 3, 3, 6, 6, 6, 1,
    #           1, 1, 1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    #           1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    #           2, 2, 2, 2,
    #           2, 2, 2, 2,2, 2, 2)
    # usenames=c(usenames,'HFC134a', 'HCFC141b', 'HCFC142b', 'HCFC22', 'Dichloromethane', 'Chloroform', '1,2-Dichloroethane',
    #            'Methyl Chloroform', 'Tetrachloroethylene', 'Chlorobenzene', 'Bromodichloromethane', 'Chlorodibromomethane', 
    #            'Methyl bromide', 'Dibromomethane', 'Bromoform', 'Chloroiodomethane', 'Iodomethane', 'Carbon disulfide', 'Methanethiol',
    #            'Dimethyl sulfate', 'Propane', 'Isobutane', 'n-Butane', 'Isopentane', 'n-Pentane', '2-Methylpentane', '3-Methylpentane',
    #            'n-Hexane', '2,2,4-Trimethylpentane', 'n-Heptane', 'n-Octane', 'Propene', 'i-Butene/1Butene', 'Isoprene', 'Tricyclene', 
    #            'alpha-Pinene', 'Camphene', 'beta-Pinene/Myrcene', 'Limonene/D3Carene', 'Benzene', 'Toluene', 'Ethylbenzene', 'm,p-Xylene', 
    #            'o-Xylene', 'Styrene', 'Ethynylbenzene', 'Formaldehyde', 'Acetaldehyde', 'Propanal', 'Butanal', 'Isobutanol', 'Acrolein', 
    #            '2-Butenals', 'Acetone/Propanal', 'Methyl Ethyl Ketone', 'Methanol', 'Ethanol', 'Isopropanol', '2-methyl-3-butene-2-ol',
    #            'Methacrolein', 'Methyl vinyl ketone', 'Methyl formate', 'Methyl acetate', 'Furan', '2-Methylfuran', '3-Methylfuran', 
    #            'Furfural', 'Hydrogen cyanide', 'Acetonitrile', 'Propionitrile', 'Acrylonitrile', 'Methacrylonitrile', 'Pyrrole',
    #            'Nitromethane', 'Methyl nitrate', 'Ethyl nitrate', 'i-propyl nitrate','2/i-Butyl nitrate')
    # formulas=c(formulas,'HFC134a', 'HCFC141b', 'HCFC142b', 'HCFC22', 'CH2Cl2', 'CHCl3', 'CH2ClCH2Cl', 'CH3CCl3', 
    #            'C2Cl4', 'C6H5Cl', 'CHBrCl2', 'CHBr2Cl', 'CH3Br', 'CH2Br2', 'CHBr3', 'CH2ClI', 'CH3I', 'CS2', 'CH3SH',
    #            'DMS', 'C3H8O', 'C4H10', 'C4H10', 'C5H12', 'C5H12', 'C6H14', 'C6H14', 'C6H14', 'C8H18', 'C7H16', 'C8H18',
    #            'C3H6', 'C4H8', 'C5H8', 'C10H16', 'C10H16', 'C10H16', 'C10H16', 'C10H16', 'C6H6', 'C7H8', 'C8H10', 'C8H10',
    #            'C8H10', 'C8H8', 'C8H6', 'CH2O', 'C2H4O', 'C3H6O', 'C4H8O', 'C4H8O', 'C3H4O', 'C4H6O', 'C3H6O', 'C4H8O', 
    #            'CH3OH', 'C2H5OH', 'C3H8O', 'C5H10O', 'C4H6O', 'C4H6O', 'C2H4O2', 'CH3COOCH3', 'C4H4O', 'C5H6O', 'C5H6O',
    #            'C5H4O2', 'HCN', 'CH3CN', 'CH3CH2CN', 'C3H3N', 'C4H5N', 'C4H5N', 'CH3NO2', 'CH3NO3', 'C2H5NO3', 'C3H7NO3',
    #            'C4H9NO3')
    # passOHrate = c(passOHrate,OHHFC134a, OHHCFC141b, OHHCFC142b, OHHCFC22, OHCH2Cl2,
    #                OHCHCl3, OHCH2ClCH2Cl, OHCH3CCl3, OHC2Cl4, OHClBenzene,
    #                OHCHBrCl2, OHCHBr2Cl, OHCH3Br, OHCH2Br2, OHCHBr3,
    #                OHCH2ClI, OHCH3I, OHCS2, OHCH3SH, OHDMS, OHPropane,
    #                OHiButane, OHnButane, OHiPentane, OHnPentane, OHx2MePentane,
    #                OHx3MePentane, OHnHexane, OHx224TrimePentane, OHnHeptane,
    #                OHnOctane, OHPropene, OHiButene, OHIsoprene,
    #                OHTricyclene, OHaPinene, OHCamphene, OHbPinene,
    #                OHLimonene, OHBenzene, OHToluene, OHEthBenzene,
    #                OHmpXylene, OHoXylene, OHStyrene, OHEthynylBenzene, OHch2o,
    #                OHCH3CHO,  OHpropanal, OHButanal, OHiButanal, OHAcrolein,
    #                OHx2Butenals, OHAcetonePropanal, OHMEK, OHCH3OH, OHC2H5OH,
    #                OHiPropanol, OHMBO, OHMAC, OHMVK, OHMeFormate,
    #                OHMeAcetate, OHFuran, OHx2MeFuran, OHx3MeFuran, OHFurfural,
    #                OHHCN, OHCH3CN, OHPropNitrile, OHAcrylonitrile,
    #                OHMeAcrylonitrile, OHpyrrole, OHNitromethane, OHMeONO2,
    #                OHEthONO2, OHiPropONO2, OHx2ButONO2)
    # if (debugKT == 1){print("APEL")}
    # 
    # #TOADD "UnknownC2H4O_TOGA"  ,        "C3O2_TOGA"       ,           "THF_TOGA"       ,
    # #"x23Butanedione_TOGA"  ,"EthAcetate_TOGA" ,"MePropionate_TOGA","x2EthFuran_TOGA" ,"DimeFurans_TOGA" ,"UnknownC6H8O_TOGA"   ,
    # #"VinylFuranA_TOGA"  ,         "VinylFuranB_TOGA"     ,      "x3Furaldehyde_TOGA"
    # 
    # aerror1 = c(2,2,2,2,2,0.1,1,1,0.2,0.2,0.1,0.1,2,0.04,0.04,0.2,0.1,0.4,2,2,10,2,2,2,2,1,1,1,1,2,1,10,4,2,0.2,0.2,
    #             0.2,0.2,0.2,0.6,0.6,0.4,0.4,0.4,0.2,2,40,10,2,1,1,2,2,10,1,10,4,10,1,2,1,2,2,2,1,1,2,10,2,2,2,4,2,0.4,
    #             2,2,1,1)
    # aerror2 = c(0.19,0.30,0.15,0.17,0.22,.20,0.17,0.15,0.15,0.24,0.15,0.15,0.17,0.22,.40,0.40,.30,0.15,0.30,0.21,0.35,0.25,0.21,0.24,0.38,0.24,0.15,0.27,0.18,0.20,0.15,0.40,0.30,0.28,0.38,0.24,0.31,0.35,0.26,0.17,0.35,0.38,0.40,0.40,0.40,0.40,0.35,0.28,0.21,0.30,0.30,0.20,0.44,0.22,0.20,0.40,0.29,0.30,0.25,0.20,0.20,0.30,0.30,0.22,0.20,0.20,0.40,0.30,0.29,0.24,0.19,0.22,0.40,0.40,0.40,0.40,0.16,0.25)
    # if (doAPEL == 1){
    #   for (i in 1:length(anames)){
    #     ccs = colnames(apel) ; ccs2=colnames(apelBG)
    #     #print(anames[i])
    #     ind = which(ccs == anames[i]) ; ind2 = which(ccs2 == anames[i])
    #     # everything in ppb
    #     t1 = apel[,ind]/1E3 ; t1B = apelBG[,ind2]/1E3
    #     #print(t1B)
    #     # Assume if a species background is NaN we can set it to zero
    #     if (is.nan(t1B)){t1B = 0}
    # 
    #     t1CO = apel$CO_DACOM_DISKIN_APEL; t1COB = apelBG$CO_DACOM_DISKIN_APEL
    #     t1CO2 = apel$CO2_7000_ppm_DISKIN_APEL*1E3; t1CO2B = apelBG$CO2_7000_ppm_DISKIN_APEL*1E3
    #     if (i == 1){
    #       passERsB = sum(t1 - t1B)/sum(t1CO2 - t1CO2B)
    #       passERsBCO = sum(t1 - t1B)/sum(t1CO - t1COB)
    # 
    #       passERsBI  = NaN
    #       passERsBIfill  = NaN
    #       passERsCOBI   = NaN
    #       passERsCOBIfill = NaN
    #       passRBCO2 =t1CO2B
    #       passRBCOB = t1COB
    #       passRBS =t1B
    #       passRB = 1# so they dont get thrown away by R2 filter
    #       passRBCO =1 # so they dont get thrown away by R2 filter
    #       passRBXX  =sum(t1CO2 - t1CO2B)
    #       passRBSS  =sum(t1 - t1B)
    #       passSB = NaN
    #       passSBQ = NaN
    #       passSBQCO =NaN
    #       passMB = max(t1)
    #     }
    #     if (i > 1){
    #       passERsB = c(passERsB,sum(t1 - t1B)/sum(t1CO2 - t1CO2B))
    #       passERsBCO = c(passERsBCO,sum(t1 - t1B)/sum(t1CO - t1COB))
    #       passERsBI = c(passERsBI,NaN)
    #       passERsBIfill = c(passERsBIfill,NaN)
    #       passERsCOBI = c(passERsCOBI, NaN)
    #       passERsCOBIfill = c(passERsCOBIfill,NaN)
    #       passRBCO2 = c(passRBCO2,t1CO2B)
    #       passRBCOB = c(passRBCOB,t1COB)
    #       passRBS = c(passRBS,t1B)
    #       passRB = c(passRB,1)# so they dont get thrown away by R2 filter
    #       passRBCO =c(passRBCO,1) # so they dont get thrown away by R2 filter
    #       passRBXX  = c(passRBXX,sum(t1CO2 - t1CO2B))
    #       passRBSS  = c(passRBSS,sum(t1 - t1B))
    #       passSBQ = c( passSBQ,NaN)
    #       passSBQCO = c(passSBQCO,NaN)
    #       passSB = c(passSB,NaN)
    #       passMB = c(passMB, max(t1))
    #     }
    #   }
    # 
    #   passERs = c(passERs, passERsB)
    #   passERsCO = c(passERsCO,passERsBCO)
    # 
    #   passERsint = c(passERsint, passERsBI)
    #   passERsintfill = c(passERsintfill, passERsBIfill)
    # 
    #   passERszero   = c(passERszero, passSBQ)
    #   passERsCOzero = c(passERsCOzero, passSBQCO)
    # 
    #   passERsCOint = c(passERsCOint, passERsCOBI)
    #   passERsCOintfill = c(passERsCOintfill, passERsCOBIfill)
    # 
    #   passBGX = c(passBGX,  passRBCO2)
    #   passBGCO = c(passBGCO, passRBCOB)
    #   passBGS = c(passBGS,  passRBS)
    #   passRsqsCO2 = c(passRsqsCO2, passRB)
    #   passRsqsCO  = c(passRsqsCO, passRBCO)
    #   passExcessX  = c(passExcessX,passRBXX)
    #   passExcessS  = c(passExcessS,passRBSS)
    # 
    #   passsigma = c(passsigma, passSB)
    #   passMAX = c(passMAX, passMB)
    # }   else{
    #   passERs = c(passERs,rep(NaN,length(anames)))  
    #   passERsCO = c(passERsCO,rep(NaN,length(anames)))  
    #   passERsint = c(passERsint, rep(NaN,length(anames)))  
    #   passERsintfill = c(passERsintfill, rep(NaN,length(anames)))  
    #   passERszero   = c(passERszero, rep(NaN,length(anames)))  
    #   passERsCOzero = c(passERsCOzero,rep(NaN,length(anames)))  
    #   passERsCOint = c(passERsCOint, rep(NaN,length(anames)))  
    #   passERsCOintfill = c(passERsCOintfill, rep(NaN,length(anames)))  
    #   passBGX = c(passBGX, rep(NaN,length(anames)))  
    #   passBGCO = c(passBGCO, rep(NaN,length(anames)))  
    #   passBGS = c(passBGS, rep(NaN,length(anames)))  
    #   passRsqsCO2 = c(passRsqsCO2, rep(NaN,length(anames)))  
    #   passRsqsCO  = c(passRsqsCO,rep(NaN,length(anames)))  
    #   passExcessX  = c(passExcessX,rep(NaN,length(anames)))  
    #   passExcessS  = c(passExcessS,rep(NaN,length(anames)))  
    #   passsigma = c(passsigma, rep(NaN,length(anames)))  
    #   passMAX = c(passMAX, rep(NaN,length(anames)))  
    # }   
    # names1=c(names1,anames)
    # mWs = c(mWs, aMW)
    # nCs = c(nCs, aNC)
    # 
    # 
    # if (debugKT == 1){print(c(length(passERs), length(passERsint), length(passERsintfill),length(passERsCO),length(passsigma), length(passRsqsCO2), length(passRsqsCO),length(passMAX),length(passOHrate),
    #                           length(names1), length(mWs), length(nCs), length(passBGS), length(passBGX), 
    #                           length(passExcessX), length(passExcessS), length(passERsCOint), length(passERsCOintfill),
    #                           length(passERszero), length(passERsCOzero), length(isNMVOC), length(kind)))}
    # 
    # ------------------------------------------ BECKY ----------------------
    if (debugKT == 1){print("BECKY")}
    anames=c('HFC134a_ppt', 'HCFC141b_ppt', 'HCFC142b_ppt', 'HCFC22_ppt', 'CH2Cl2_ppt', 
             'CHCl3_ppt', 'CH2ClCH2Cl_ppt', 'CH3CCl3_ppt', 'C2Cl4_ppt', 'ClBenzene_ppt',
             'CHBrCl2_ppt', 'CHBr2Cl_ppt', 'CH3Br_ppt', 'CH2Br2_ppt', 'CHBr3_ppt',
             'CH2ClI_ppt', 'CH3I_ppt', 'CS2_ppt', 'CH3SH_ppt', 'DMS_ppt', 'Propane_ppt',
             'iButane_ppt', 'nButane_ppt', 'iPentane_ppt', 'nPentane_ppt', 'x2MePentane_ppt',
             'x3MePentane_ppt', 'nHexane_ppt', 'x224TrimePentane_ppt', 'nHeptane_ppt',
             'nOctane_ppt', 'Propene_ppt', 'iButene1Butene_ppt', 'Isoprene_ppt', 
             'Tricyclene_ppt', 'aPinene_ppt', 'Camphene_ppt', 'bPineneMyrcene_ppt', 
             'LimoneneD3Carene_ppt', 'Benzene_ppt', 'Toluene_ppt', 'EthBenzene_ppt', 
             'mpXylene_ppt', 'oXylene_ppt', 'Styrene_ppt', 'EthynylBenzene_ppt', 'CH2O_ppt', 
             'CH3CHO_ppt', 'Propanal_ppt', 'Butanal_ppt', 'iButanal_ppt', 'Acrolein_ppt', 
             'x2Butenals_ppt', 'Acetone_ppt', 'MEK_ppt', 'CH3OH_ppt', 'C2H5OH_ppt', 
             'iPropanol_ppt', 'MBO_ppt', 'MAC_ppt', 'MVK_ppt', 'MeFormate_ppt', 
             'MeAcetate_ppt', 'Furan_ppt', 'x2MeFuran_ppt', 'x3MeFuran_ppt', 'Furfural_ppt', 
             'HCN_ppt', 'CH3CN_ppt', 'PropNitrile_ppt', 'Acrylonitrile_ppt',
             'MeAcrylonitrile_ppt', 'Pyrrole_ppt', 'Nitromethane_ppt', 'MeONO2_ppt', 
             'EthONO2_ppt', 'iPropONO2_ppt', 'iButONO2and2ButONO2_ppt')
    aMW = c(  mWHFC134a, mWHCFC141b, mWHCFC142b, mWHCFC22, mWCH2Cl2, 
              mWCHCl3, mWCH2ClCH2Cl, mWCH3CCl3, mWC2Cl4, mWClBenzene,
              mWCHBrCl2, mWCHBr2Cl, mWCH3Br, mWCH2Br2, mWCHBr3,
              mWCH2ClI, mWCH3I, mWCS2, mWCH3SH, mWDMS, mWPropane,
              mWiButane, mWnButane, mWiPentane, mWnPentane, mWx2MePentane,
              mWx3MePentane, mWnHexane, mWx224TrimePentane, mWnHeptane,
              mWnOctane, mWPropene, mWiButene, mWIsoprene, 
              mWTricyclene, mWaPinene, mWCamphene, mWbPinene, 
              mWLimonene, mWBenzene, mWToluene, mWEthBenzene, 
              mWmpXylene, mWoXylene, mWStyrene, mWEthynylBenzene, mWCH2O, 
              mWAcetaldehyde, mWpropanal, mWButanal, mWiButanal, mWAcrolein, 
              mW2Butenals, mWAcetonePropanal, mWMEK, mWCH3OH, mWC2H5OH, 
              mWiPropanol, mWMBO, mWMAC, mWMVK, mWMeFormate, 
              mWMeAcetate, mWFuran, mWx2MeFuran, mWx3MeFuran, mWfurfural, 
              mWHCN, mWch3cn, mWPropNitrile, mWAcrylonitrile,
              mWMeAcrylonitrile, mwPyrrole, mWNitromethane, mWMeONO2, 
              mWEthONO2, mWiPropONO2, MWx2ButONO2iButONO2)
    aNC =c(  nCHFC134a, nCHCFC141b, nCHCFC142b, nCHCFC22, nCCH2Cl2, 
             nCCHCl3, nCCH2ClCH2Cl, nCCH3CCl3, nCC2Cl4, nCClBenzene,
             nCCHBrCl2, nCCHBr2Cl, nCCH3Br, nCCH2Br2, nCCHBr3,
             nCCH2ClI, nCCH3I, nCCS2, nCCH3SH, nCDMS, nCPropane,
             nCiButane, nCnButane, nCiPentane, nCnPentane, nCx2MePentane,
             nCx3MePentane, nCnHexane, nCx224TrimePentane, nCnHeptane,
             nCnOctane, nCPropene, nCiButene, nCIsoprene, 
             nCTricyclene, nCaPinene, nCCamphene, nCbPinene, 
             nCLimonene, nCBenzene, nCToluene, nCEthBenzene, 
             nCmpXylene, nCoXylene, nCStyrene, nCEthynylBenzene, nCCH2O, 
             nCAcetaldehyde, nCpropanal, nCButanal, nCButanal, nCAcrolein, 
             nC2Butenals, nCAcetonePropanal, nCMEK, nCCH3OH, nCC2H5OH, 
             nCiPropanol, nCMBO, nCMAC, nCMVK, nCMeFormate, 
             nCMeAcetate, nCFuran, nCx2MeFuran, nCx3MeFuran, nCfurfural, 
             nCHCN, nCch3cn, nCPropNitrile, nCAcrylonitrile,
             nCMeAcrylonitrile, nCPyrrole, nCNitromethane, nCMeONO2, 
             nCEthONO2, nCiPropONO2, nCx2ButONO2iButONO2)
    kind = c(kind,'halocarbon', 'halocarbon', 'halocarbon', 'halocarbon', 'halocarbon', 
                    'halocarbon', 'CH2ClCH2Cl_ppt', 'halocarbon', 'halocarbon', 'ClBenzene_ppt',
                    'CHBrCl2_ppt', 'CHBr2Cl_ppt', 'CH3Br_ppt', 'CH2Br2_ppt', 'CHBr3_ppt',
                    'CH2ClI_ppt', 'CH3I_ppt', 'sulfur', 'sulfur', 'sulfur', 'alkane',
                    'alkane', 'alkane', 'alkane', 'alkane', 'alkane',
                    'alkane', 'alkane', 'alkane', 'alkane',
                    'alkane', 'alkene', 'alkene', 'alkene', 
                    'alkene', 'alkene', 'alkene', 'alkene', 
             'alkene', 'aromatic', 'aromatic', 'aromatic', 
             'aromatic', 'aromatic', 'aromatic', 'aromatic', 'CH2O', 
                    'oVOC', 'oVOC', 'oVOC', 'oVOC', 'oVOC', 
                    'oVOC', 'oVOC', 'oVOC', 'oVOC', 'oVOC', 
                    'oVOC', 'oVOC', 'oVOC', 'oVOC', 'oVOC', 
                    'oVOC', 'oVOC', 'x2MeFuran_ppt', 'oVOC', 'oVOC', 
                    'nitrogen', 'nitrogen', 'nitrogen', 'nitrogen',
                    'nitrogen', 'nitrogen', 'nitro', 'alkyl nitrate', 
                    'alkyl nitrate', 'alkyl nitrate', 'alkyl nitrate')
    isNMVOC=c(isNMVOC,3, 3, 3, 3,3, 3, 3, 3,3, 3,
              3, 3, 3, 3,3, 3, 3, 6, 6, 6, 1,
              1, 1, 1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
              2, 2, 2, 2,
              2, 2, 2, 2,2, 2, 2)
    passOHrate = c(passOHrate,OHHFC134a, OHHCFC141b, OHHCFC142b, OHHCFC22, OHCH2Cl2, 
                   OHCHCl3, OHCH2ClCH2Cl, OHCH3CCl3, OHC2Cl4, OHClBenzene,
                   OHCHBrCl2, OHCHBr2Cl, OHCH3Br, OHCH2Br2, OHCHBr3,
                   OHCH2ClI, OHCH3I, OHCS2, OHCH3SH, OHDMS, OHPropane,
                   OHiButane, OHnButane, OHiPentane, OHnPentane, OHx2MePentane,
                   OHx3MePentane, OHnHexane, OHx224TrimePentane, OHnHeptane,
                   OHnOctane, OHPropene, OHiButene, OHIsoprene, 
                   OHTricyclene, OHaPinene, OHCamphene, OHbPinene, 
                   OHLimonene, OHBenzene, OHToluene, OHEthBenzene, 
                   OHmpXylene, OHoXylene, OHStyrene, OHEthynylBenzene, OHch2o, 
                   OHCH3CHO,  OHpropanal, OHButanal, OHiButanal, OHAcrolein, 
                   OHx2Butenals, OHAcetonePropanal, OHMEK, OHCH3OH, OHC2H5OH, 
                   OHiPropanol, OHMBO, OHMAC, OHMVK, OHMeFormate, 
                   OHMeAcetate, OHFuran, OHx2MeFuran, OHx3MeFuran, OHFurfural, 
                   OHHCN, OHCH3CN, OHPropNitrile, OHAcrylonitrile,
                   OHMeAcrylonitrile, OHpyrrole, OHNitromethane, OHMeONO2, 
                   OHEthONO2, OHiPropONO2, OHx2ButONO2)
    if (debugKT == 1){print("Becky")}
    
    #TOADD "UnknownC2H4O_TOGA"  ,        "C3O2_TOGA"       ,           "THF_TOGA"       ,   
    #"x23Butanedione_TOGA"  ,"EthAcetate_TOGA" ,"MePropionate_TOGA","x2EthFuran_TOGA" ,"DimeFurans_TOGA" ,"UnknownC6H8O_TOGA"   ,      
    #"VinylFuranA_TOGA"  ,         "VinylFuranB_TOGA"     ,      "x3Furaldehyde_TOGA"    
    
    aerror1 = c(2,2,2,2,2,0.1,1,1,0.2,0.2,0.1,0.1,2,0.04,0.04,0.2,0.1,0.4,2,2,10,2,2,2,2,1,1,1,1,2,1,10,4,2,0.2,0.2,
                0.2,0.2,0.2,0.6,0.6,0.4,0.4,0.4,0.2,2,40,10,2,1,1,2,2,10,1,10,4,10,1,2,1,2,2,2,1,1,2,10,2,2,2,4,2,0.4,
                2,2,1,1)
    aerror2 = c(0.19,0.30,0.15,0.17,0.22,.20,0.17,0.15,0.15,0.24,0.15,0.15,0.17,0.22,.40,0.40,.30,0.15,0.30,0.21,0.35,0.25,0.21,0.24,0.38,0.24,0.15,0.27,0.18,0.20,0.15,0.40,0.30,0.28,0.38,0.24,0.31,0.35,0.26,0.17,0.35,0.38,0.40,0.40,0.40,0.40,0.35,0.28,0.21,0.30,0.30,0.20,0.44,0.22,0.20,0.40,0.29,0.30,0.25,0.20,0.20,0.30,0.30,0.22,0.20,0.20,0.40,0.30,0.29,0.24,0.19,0.22,0.40,0.40,0.40,0.40,0.16,0.25)
    passErrors=c(passErrors,aerror1)
    if (doBecky == 1){
      for (i in 1:length(anames)){
        ccs = colnames(becky) ; ccs2=colnames(beckyBG)
        #print(anames[i])
        ind = which(ccs == anames[i]) ; ind2 = which(ccs2 == anames[i])
        # everything in ppb
        t1 = as.numeric(becky[,ind]/1E3) ; t1B = as.numeric(beckyBG[,ind2]/1E3)
        #print(t1B)
        # Assume if a species background is NaN we can set it to zero
        if (is.nan(t1B)){t1B = 0}
        
        t1CO = as.numeric(becky$CO_DACOM_DISKIN_BECKY); t1COB = as.numeric(beckyBG$CO_DACOM_DISKIN_BECKY)
        t1CO2 = as.numeric(becky$CO2_7000_ppm_DISKIN_BECKY*1E3); t1CO2B = as.numeric(beckyBG$CO2_7000_ppm_DISKIN_BECKY*1E3)
        
        # Are we using the fast background?
        if (useFASTBG  == 1 ){
          t1COB = mean(background$CO_DACOM_DISKIN, na.rm=TRUE)
          t1CO2B = mean(background$CO2_7000_ppm_DISKIN*1E3)
        }
        
        if (i == 1){
          passERsB = sum(t1 - t1B)/sum(t1CO2 - t1CO2B)
          passERsBCO = sum(t1 - t1B)/sum(t1CO - t1COB)
          
          passERsBI  = NaN
          passERsBIfill  = NaN
          passERsCOBI   = NaN
          passERsCOBIfill = NaN
          passRBCO2 =t1CO2B
          passRBCOB = t1COB
          passRBS =t1B
          passRB = 1# so they dont get thrown away by R2 filter
          passRBCO =1 # so they dont get thrown away by R2 filter
          passRBXX  =sum(t1CO2 - t1CO2B)
          passRBSS  =sum(t1 - t1B)
          passSB = NaN
          passSBQ = NaN
          passSBQCO =NaN
          passMB = max(t1)
        }
        if (i > 1){ 
          passERsB = c(passERsB,sum(t1 - t1B)/sum(t1CO2 - t1CO2B))
          passERsBCO = c(passERsBCO,sum(t1 - t1B)/sum(t1CO - t1COB))
          passERsBI = c(passERsBI,NaN)
          passERsBIfill = c(passERsBIfill,NaN)
          passERsCOBI = c(passERsCOBI, NaN)
          passERsCOBIfill = c(passERsCOBIfill,NaN)
          passRBCO2 = c(passRBCO2,t1CO2B)
          passRBS = c(passRBS,t1B)
          passRBCOB = c(passRBCOB,t1COB)
          
          passRB = c(passRB,1)# so they dont get thrown away by R2 filter
          passRBCO =c(passRBCO,1) # so they dont get thrown away by R2 filter
          passRBXX  = c(passRBXX,sum(t1CO2 - t1CO2B))
          passRBSS  = c(passRBSS,sum(t1 - t1B))
          passSBQ = c( passSBQ,NaN)
          passSBQCO = c(passSBQCO,NaN)
          passSB = c(passSB,NaN)
          passMB = c(passMB, max(t1))
        }
      }
      
      passERs = c(passERs, passERsB)  
      passERsCO = c(passERsCO,passERsBCO) 
      
      passERsint = c(passERsint, passERsBI)    
      passERsintfill = c(passERsintfill, passERsBIfill)   
      
      passERszero   = c(passERszero, passSBQ)
      passERsCOzero = c(passERsCOzero, passSBQCO)
      
      passERsCOint = c(passERsCOint, passERsCOBI)
      passERsCOintfill = c(passERsCOintfill, passERsCOBIfill)
      
      passBGX = c(passBGX,  passRBCO2)
      passBGCO = c(passBGCO, passRBCOB)
      passBGS = c(passBGS,  passRBS)
      passRsqsCO2 = c(passRsqsCO2, passRB)   
      passRsqsCO  = c(passRsqsCO, passRBCO)   
      passExcessX  = c(passExcessX,passRBXX)
      passExcessS  = c(passExcessS,passRBSS)
      
      passsigma = c(passsigma, passSB)
      passMAX = c(passMAX, passMB)
    } else{
      passERs = c(passERs,rep(NaN,length(anames)))  
      passERsCO = c(passERsCO,rep(NaN,length(anames)))  
      passERsint = c(passERsint, rep(NaN,length(anames)))  
      passERsintfill = c(passERsintfill, rep(NaN,length(anames)))  
      passERszero   = c(passERszero, rep(NaN,length(anames)))  
      passERsCOzero = c(passERsCOzero,rep(NaN,length(anames)))  
      passERsCOint = c(passERsCOint, rep(NaN,length(anames)))  
      passERsCOintfill = c(passERsCOintfill, rep(NaN,length(anames)))  
      passBGX = c(passBGX, rep(NaN,length(anames)))  
      passBGCO = c(passBGCO, rep(NaN,length(anames)))  
      passBGS = c(passBGS, rep(NaN,length(anames)))  
      passRsqsCO2 = c(passRsqsCO2, rep(NaN,length(anames)))  
      passRsqsCO  = c(passRsqsCO,rep(NaN,length(anames)))  
      passExcessX  = c(passExcessX,rep(NaN,length(anames)))  
      passExcessS  = c(passExcessS,rep(NaN,length(anames)))  
      passsigma = c(passsigma, rep(NaN,length(anames)))  
      passMAX = c(passMAX, rep(NaN,length(anames)))  
    }    
    names1=c(names1,anames)
    usenames=c(usenames,'HFC134a', 'HCFC141b', 'HCFC142b', 'HCFC22', 'Dichloromethane', 'Chloroform', '1,2-Dichloroethane',
               'Methyl Chloroform', 'Tetrachloroethylene', 'Chlorobenzene', 'Bromodichloromethane', 'Chlorodibromomethane', 
               'Methyl bromide', 'Dibromomethane', 'Bromoform', 'Chloroiodomethane', 'Iodomethane', 'Carbon disulfide', 'Methanethiol',
               'Dimethyl sulfide', 'Propane', 'Isobutane', 'n-Butane', 'Isopentane', 'n-Pentane', '2-Methylpentane', '3-Methylpentane',
               'n-Hexane', '2,2,4-Trimethylpentane', 'n-Heptane', 'n-Octane', 'Propene', 'i-Butene/1Butene', 'Isoprene', 'Tricyclene', 
               '⍺-Pinene', 'Camphene', 'β-Pinene/Myrcene', 'Limonene/D3Carene', 'Benzene', 'Toluene', 'Ethylbenzene', 'm,p-Xylene', 
               'o-Xylene', 'Styrene', 'Ethynylbenzene', 'Formaldehyde', 'Acetaldehyde', 'Propanal', 'Butanal', 'Isobutanal', 'Acrolein', 
               '2-Butenals', 'Acetone/Propanal', 'Methyl ethyl ketone', 'Methanol', 'Ethanol', 'Isopropanol', '2-methyl-3-butene-2-ol',
               'Methacrolein', 'Methyl vinyl ketone', 'Methyl formate', 'Methyl acetate', 'Furan', '2-Methylfuran', '3-Methylfuran', 
               'Furfural', 'Hydrogen cyanide', 'Acetonitrile', 'Propionitrile', 'Acrylonitrile', 'Methacrylonitrile', 'Pyrrole',
               'Nitromethane', 'Methyl nitrate', 'Ethyl nitrate', 'i-propyl nitrate','2/i-Butyl nitrate')
    formulas=c(formulas,'HFC134a', 'HCFC141b', 'HCFC142b', 'HCFC22', 'CH2Cl2', 'CHCl3', 'CH2ClCH2Cl', 'CH3CCl3', 
               'C2Cl4', 'C6H5Cl', 'CHBrCl2', 'CHBr2Cl', 'CH3Br', 'CH2Br2', 'CHBr3', 'CH2ClI', 'CH3I', 'CS2', 'CH3SH',
               'C2H6S', 'C3H8', 'C4H10', 'C4H10', 'C5H12', 'C5H12', 'C6H14', 'C6H14', 'C6H14', 'C8H18', 'C7H16', 'C8H18',
               'C3H6', 'C4H8', 'C5H8', 'C10H16', 'C10H16', 'C10H16', 'C10H16', 'C10H16', 'C6H6', 'C7H8', 'C8H10', 'C8H10',
               'C8H10', 'C8H8', 'C8H6', 'CH2O', 'C2H4O', 'C3H6O', 'C4H8O', 'C4H8O', 'C3H4O', 'C4H6O', 'C3H6O', 'C4H8O', 
               'CH3OH', 'C2H5OH', 'C3H8O', 'C5H10O', 'C4H6O', 'C4H6O', 'C2H4O2', 'CH3COOCH3', 'C4H4O', 'C5H6O', 'C5H6O',
               'C5H4O2', 'HCN', 'CH3CN', 'CH3CH2CN', 'C3H3N', 'C4H5N', 'C4H5N', 'CH3NO2', 'CH3NO3', 'C2H5NO3', 'C3H7NO3',
               'C4H9NO3')
    mWs = c(mWs, aMW) 
    nCs = c(nCs, aNC)
    if (debugKT == 1){print(c(length(passERs), length(passERsint), length(passERsintfill),length(passERsCO),length(passsigma), length(passRsqsCO2), length(passRsqsCO),length(passMAX),length(passOHrate),
                              length(names1), length(mWs), length(nCs), length(passBGS), length(passBGX), 
                              length(passExcessX), length(passExcessS), length(passERsCOint), length(passERsCOintfill),
                              length(passERszero), length(passERsCOzero), length(isNMVOC), length(kind), length(formulas), length(usenames)))}
    # ------------------------------------------ GILMAN ----------------------
    if (debugKT == 1){print("GILMAN")}
    gnames = c('C2Cl4_NOAAiWAS_GILMAN', 'CHCl3_NOAAiWAS_GILMAN', 'Ethane_NOAAiWAS_GILMAN', 'Propane_NOAAiWAS_GILMAN', 
               'nButane_NOAAiWAS_GILMAN', 'iButane_NOAAiWAS_GILMAN', 'nPentane_NOAAiWAS_GILMAN', 'iPentane_NOAAiWAS_GILMAN',
               'nHexane_NOAAiWAS_GILMAN', 'x2MePentane_NOAAiWAS_GILMAN', 'x3MePentane_NOAAiWAS_GILMAN',
               'x22DiMeButane_NOAAiWAS_GILMAN', 'x24DiMePentane_NOAAiWAS_GILMAN', 'nOctane_NOAAiWAS_GILMAN',
               'x224TriMePentane_NOAAiWAS_GILMAN', 'nNonane_NOAAiWAS_GILMAN', 'nDecane_NOAAiWAS_GILMAN', 
               'MeCycPentane_NOAAiWAS_GILMAN', 'CycHexane_NOAAiWAS_GILMAN', 'MeCycHexane_NOAAiWAS_GILMAN', 
               'Ethyne_NOAAiWAS_GILMAN', 'Ethene_NOAAiWAS_GILMAN', 'Propene_NOAAiWAS_GILMAN', 'x1Butene_NOAAiWAS_GILMAN',
               'c2Butene_NOAAiWAS_GILMAN', 't2butene_NOAAiWAS_GILMAN', 'iButene_NOAAiWAS_GILMAN', 'x1Pentene_NOAAiWAS_GILMAN', 
               'c2Pentene_NOAAiWAS_GILMAN', 't2Pentene_NOAAiWAS_GILMAN', 'x2Me1Butene_NOAAiWAS_GILMAN', 'x3Me1Butene_NOAAiWAS_GILMAN', 
               't13Pentadiene_NOAAiWAS_GILMAN', 'Isoprene_NOAAiWAS_GILMAN', 'aPinene_NOAAiWAS_GILMAN', 'Benzene_NOAAiWAS_GILMAN', 
               'Toluene_NOAAiWAS_GILMAN', 'EthBenzene_NOAAiWAS_GILMAN', 'oXylene_NOAAiWAS_GILMAN', 'mpXylene_NOAAiWAS_GILMAN', 
               'Acetone_NOAAiWAS_GILMAN', 'MEK_NOAAiWAS_GILMAN', 'MeFormate_NOAAiWAS_GILMAN', 'Furan_NOAAiWAS_GILMAN', 
               'CH3CN_NOAAiWAS_GILMAN', 'Acrylonitrile_NOAAiWAS_GILMAN')
    formulas=c(formulas,'C2Cl4', 'CHCl3', 'C2H6', 'C3H8', 'C4H10', 'C4H10', 'C5H12', 'C5H12', 'C6H14', 'C6H14', 'C6H14', 
    'C6H14', 'C7H16', 'C8H18', 'C8H18', 'C9H20', 'C10H22', 'C6H12', 'C6H12', 'C7H14', 'C2H2', 'C2H4', 'C3H6', 'C4H8',
    'C4H8', 'C4H8', 'C4H8', 'C5H10', 'C5H10', 'C5H10', 'C5H10', 'C5H10', 'C5H8', 'C5H8', 'C10H16', 'C6H6', 'C7H8', 'C8H10',
    'C8H10', 'C8H10', 'C3H6O', 'C4H8O', 'C2H4O2', 'C4H4O', 'CH3CN','C3H3N')
    usenames=c(usenames,'Tetrachloroethylene', 'Chloroform', 'Ethane', 'Propane', 'n-Butane', 'Isobutane', 'n-Pentane',
               'Isopentane', 'n-Hexane', '2-Methylpentane', '3-Methylpentane', '2,2-Dimethylbutane', '2,4-Dimethylpentane',
               'n-Octane', '2,2,4-Trimethylpentane', 'n-Nonane', 'n-Decane', 'Methylcyclopentane', 'Cyclohexane', 
               'Methylcyclohexane', 'Ethyne', 'Ethene', 'Propene', '1-Butene', 'cis-2-Butene', 'trans-2-Butene', 
               'i-Butene', '1-Pentene', 'cis-2-Pentene', 'trans-2-Pentene', '2-Methyl-1-butene', '3-Methyl-1-butene',
               'trans-1,3-Pentadiene', 'Isoprene', '⍺-Pinene', 'Benzene', 'Toluene', 'Ethylbenzene', 'o-Xylene', 'm,p-Xylene',
               'Acetone', 'Methyl ethyl ketone', 'Methyl formate', 'Furan', 'Acetonitrile',	'Acrylonitrile')
    kind = c(kind,'halogen', 'halogen', 'alkane', 'alkane', 
               'alkane', 'alkane', 'alkane', 'alkane',
               'alkane', 'alkane', 'alkane',
               'alkane', 'alkane', 'alkane',
               'alkane', 'alkane', 'alkane', 
               'alkane', 'alkane', 'alkane', 
               'alkyne', 'alkene', 'alkene', 'alkene',
               'alkene', 'alkene', 'alkene', 'alkene', 
               'alkene', 'alkene', 'alkene', 'alkene', 
               'alkene', 'alkene', 'alkene', 'aromatic', 
               'aromatic', 'aromatic', 'aromatic', 'aromatic', 
               'oVOC', 'oVOC', 'oVOC', 'oVOC', 
               'nitrogen', 'nitrogen')
    gMW = c(mWC2Cl4, mWCHCl3, mWEthane, mWPropane, mWnButane, mWiButane,mWnPentane, mWiPentane,
               mWnHexane, mWx2MePentane, mWx3MePentane,
               mWx22Dimebutane, mWdimepentane, mWnOctane,
               mWx224TrimePentane, mWnNonane, mWnDecane, 
               mWMeCycPentane, mWCycHexane, mWMeCycHexane, 
               mWEthyne, mWEthene, mWPropene, mWx1Butene,
               mWc2Butene, mWt2Butene, mWiButene, mWx1Pentene, 
               mWc2Pentene, mWt2Pentene, mWX2Me1Butene, mWX3Me1Butene, 
               mW13Pentadiene, mWIsoprene, mWaPinene, mWBenzene, 
               mWToluene, mWEthBenzene, mWoXylene, mWmpXylene, 
               mWAcetonePropanal, mWMEK, mWMeFormate, mWFuran, 
               mWch3cn, mWAcrylonitrile)
    gNC =c(nCC2Cl4, nCCHCl3, nCEthane, nCPropane, nCnButane, nCiButane,nCnPentane, nCiPentane,
                 nCnHexane, nCx2MePentane, nCx3MePentane,
                 nCx22Dimebutane, nCx23DimePentane, nCnOctane,
                 nCx224TrimePentane, nCnNonane, nCnDecane, 
                 nCMeCycPentane, nCCycHexane, nCMeCycHexane, 
                 nCEthyne, nCEthene, nCPropene, nCx1Butene,
                 nCc2Butene, nCt2Butene, nCiButene, nCx1Pentene, 
                 nCc2Pentene, nCt2Pentene, nCX2Me1Butene, nCX3Me1Butene, 
                 nCx13Pentadienes, nCIsoprene, nCaPinene, nCBenzene, 
                 nCToluene, nCEthBenzene, nCoXylene, nCmpXylene, 
                 nCAcetonePropanal, nCMEK, nCMeFormate, nCFuran, 
                 2, nCAcrylonitrile)
    isNMVOC=c(isNMVOC,3, 3, 1, 1, 
              1, 1,1, 1, 
              1, 1, 1, 1,  1, 1, 
              1, 1,  1, 1, 1, 1, 
              1, 1, 1, 1, 
              1, 1, 1, 1, 
              1, 1,1, 1, 
              1, 1, 1, 1, 
              1, 1,1, 1, 
              1, 1,  1, 1, 2,2)
    gerror = c(seq(46))*0 # error not listed in ICARTT file, need to get this eventually
    passErrors=c(passErrors,gerror)
    
    if (doGilman== 1){
      for (i in 1:length(gnames)){
        ccs = colnames(gilman) ; ccs2=colnames(gilmanBG)
        #print(gnames[i])
        ind = which(ccs == gnames[i]) ; ind2 = which(ccs2 == gnames[i])
        # everything in ppb
        t1 = gilman[,ind] ; t1B = gilmanBG[,ind2]
        # Assume if a species background is NaN we can set it to zero
        if (is.nan(t1B)){t1B = 0}
        
        t1CO = gilman$CO_DACOM_DISKIN_GILMAN; t1COB = gilmanBG$CO_DACOM_DISKIN_GILMAN
        t1CO2 = gilman$CO2_7000_ppm_DISKIN_GILMAN*1E3; t1CO2B = gilmanBG$CO2_7000_ppm_DISKIN_GILMAN*1E3
        # Are we using the fast background?
        if (useFASTBG  == 1 ){
          t1COB = mean(background$CO_DACOM_DISKIN, na.rm=TRUE)
          t1CO2B = mean(background$CO2_7000_ppm_DISKIN*1E3)
        }
        
        if (i == 1){
          passERsB = sum(t1 - t1B)/sum(t1CO2 - t1CO2B)
          passERsBCO = sum(t1 - t1B)/sum(t1CO - t1COB)
          
          passERsBI  = NaN
          passERsBIfill  = NaN
          passERsCOBI   = NaN
          passERsCOBIfill = NaN
          passRBCO2 =t1CO2B
          passRBCOB = t1COB
          passRBS =t1B
          passRB = 1# so they dont get thrown away by R2 filter
          passRBCO =1 # so they dont get thrown away by R2 filter
          passRBXX  =sum(t1CO2 - t1CO2B)
          passRBSS  =sum(t1 - t1B)
          passSB = NaN
          passSBQ = NaN
          passSBQCO =NaN
          passMB = max(t1)
        }
        if (i > 1){ 
          passERsB = c(passERsB,sum(t1 - t1B)/sum(t1CO2 - t1CO2B))
          passERsBCO = c(passERsBCO,sum(t1 - t1B)/sum(t1CO - t1COB))
          passERsBI = c(passERsBI,NaN)
          passERsBIfill = c(passERsBIfill,NaN)
          passERsCOBI = c(passERsCOBI, NaN)
          passERsCOBIfill = c(passERsCOBIfill,NaN)
          passRBCO2 = c(passRBCO2,t1CO2B)
          passRBCOB  = c(passRBCOB,t1COB)
          passRBS = c(passRBS,t1B)
          passRB = c(passRB,1)# so they dont get thrown away by R2 filter
          passRBCO =c(passRBCO,1) # so they dont get thrown away by R2 filter
          passRBXX  = c(passRBXX,sum(t1CO2 - t1CO2B))
          passRBSS  = c(passRBSS,sum(t1 - t1B))
          passSBQ = c( passSBQ,NaN)
          passSBQCO = c(passSBQCO,NaN)
          passSB = c(passSB,NaN)
          passMB = c(passMB, max(t1))
        }
      }
      
      passERs = c(passERs, passERsB)  
      passERsCO = c(passERsCO,passERsBCO) 
      
      passERsint = c(passERsint, passERsBI)    
      passERsintfill = c(passERsintfill, passERsBIfill)   
      
      passERszero   = c(passERszero, passSBQ)
      passERsCOzero = c(passERsCOzero, passSBQCO)
      
      passERsCOint = c(passERsCOint, passERsCOBI)
      passERsCOintfill = c(passERsCOintfill, passERsCOBIfill)
      
      passBGX = c(passBGX,  passRBCO2)
      passBGCO = c(passBGCO,  passRBCOB)
      passBGS = c(passBGS,  passRBS)
      passRsqsCO2 = c(passRsqsCO2, passRB)   
      passRsqsCO  = c(passRsqsCO, passRBCO)   
      passExcessX  = c(passExcessX,passRBXX)
      passExcessS  = c(passExcessS,passRBSS)
      
      passsigma = c(passsigma, passSB)
      passMAX = c(passMAX, passMB)
    } else{
      
      passERs = c(passERs, rep(NaN, length(gnames)))  
      passERsCO = c(passERsCO,rep(NaN, length(gnames)))  
      
      passERsint = c(passERsint, rep(NaN, length(gnames)))  
      passERsintfill = c(passERsintfill, rep(NaN, length(gnames)))  
      
      passERszero   = c(passERszero, rep(NaN, length(gnames)))  
      passERsCOzero = c(passERsCOzero,rep(NaN, length(gnames)))  
      
      passERsCOint = c(passERsCOint,rep(NaN, length(gnames)))  
      passERsCOintfill = c(passERsCOintfill,rep(NaN, length(gnames)))  
      
      passBGX = c(passBGX,  rep(NaN, length(gnames)))  
      passBGCO = c(passBGCO, rep(NaN, length(gnames)))  
      passBGS = c(passBGS,  rep(NaN, length(gnames)))  
      passRsqsCO2 = c(passRsqsCO2, rep(NaN, length(gnames)))  
      passRsqsCO  = c(passRsqsCO, rep(NaN, length(gnames)))  
      passExcessX  = c(passExcessX,rep(NaN, length(gnames)))  
      passExcessS  = c(passExcessS,rep(NaN, length(gnames)))  
      
      passsigma = c(passsigma,rep(NaN, length(gnames)))  
      passMAX = c(passMAX, rep(NaN, length(gnames)))  
    }      
    passOHrate = c(passOHrate, OHC2Cl4, OHCHCl3, OHc2h6, OHPropane, 
                   OHnButane, OHiButane,OHnPentane, OHiPentane,
                   OHnHexane, OHx2MePentane, OHx3MePentane,
                   OHx22Dimebutane, OHx24DiMePentane, OHnOctane,
                   OHx224TrimePentane, OHnNonane, OHnDecane, 
                   OHMeCycPentane, OHCycHexane, OHMeCycHexane,
                   OHEthyne, OHEthene, OHPropene, OHx1Butene ,
                   OHc2Butene, OHt2Butene, OHiButene, OHx1Pentene, 
                   OHc2Pentene, OHt2Pentene, OHX2Me1Butene, OHX3Me1Butene, 
                   OHt13Pentadiene, OHIsoprene, OHaPinene, OHBenzene, 
                   OHToluene, OHEthBenzene, OHoXylene, OHmpXylene, 
                   OHAcetonePropanal, OHMEK, OHMeFormate, OHFuran, 
                   OHCH3CN, OHAcrylonitrile)
    names1=c(names1,gnames)
    mWs = c(mWs, gMW) 
    nCs = c(nCs, gNC) 
    
  }
  # Present in fast data
  
  if (debugKT == 1){print(c(length(passERs), length(passERsint), length(passERsintfill),length(passERsCO),length(passsigma), length(passRsqsCO2), length(passRsqsCO),length(passMAX),length(passOHrate),
                            length(names1), length(mWs), length(nCs), length(passBGS), length(passBGX), 
                            length(passExcessX), length(passExcessS), length(passERsCOint), length(passERsCOintfill),
                            length(passERszero), length(passERsCOzero), length(isNMVOC), length(kind)))}
  # ------------------------------------------ JIMENEZ AMS ----------------
  if (res == 5) {
    if (debugKT == 1){print("JIMENEZ AMS")}
    convOC = (1.29/28.97)*12E-3  ; 
    # if (debugKT == 1){('WARNING, need to fix OAtoOC')
    convNO3 = (1.29/28.97)*62E-3  ; convCl = (1.29/28.97)*35E-3 ;convSSA = (1.29/28.97)*31.4E-3
    convClO4 = (1.29/28.97)*99.45E-3 ; convMSA = (1.29/28.97)*96.1E-3; convNH4 = (1.29/28.97)*18E-3; convSO4 = (1.29/28.97)*96E-3 
    convBr = (1.29/28.97)*79.9E-3 ; convI = (1.29/28.97)*126.9E-3 ; convC6H10O5 = (1.29/28.97)*162.141E-3; convC6H5NO4= (1.29/28.97)*155.1081E-3
    convK = (1.29/28.97)*39.1E-3 
    nair_STP = (6.022E23*1.013E5)/(273*8.31)
    passT = mean(pass$Static_Air_Temp.x, na.rm=TRUE)+273.15
    passP = mean(pass$Static_Pressure.x, na.rm=TRUE)
    passPatm = passP * 0.000986923
    pass$nair = (6.022E23*passP*100)/(passT*8.31)
    RHS = pass$nair * (1/6.022E23)*298/passT * passPatm/1/1E3 
    # BG
    backgroundT = mean(background$Static_Air_Temp.x, na.rm=TRUE)+273.15
    backgroundP = mean(background$Static_Pressure.x, na.rm=TRUE)
    backgroundPatm = backgroundP * 0.000986923
    background$nair = (6.022E23*backgroundP*100)/(backgroundT*8.31)
    BG_RHS = background$nair * (1/6.022E23)*298/backgroundT * backgroundPatm/1/1E3 
    
    #Inorganics +/-34%, Organics +/-38%,
    cc2 = colnames(pass)
    ind = which(cc2 == 'OAtoOC_PM1_AMS_JIMENEZ' ) # are we in the 5 hz data or 1hz merge?
    if (length(ind) > 0){
      OAtoOC = pass$OAtoOC_PM1_AMS_JIMENEZ
      OAtoOCB = background$OAtoOC_PM1_AMS_JIMENEZ
      OtoC = pass$OtoC_Ratio_PM1_AMS_JIMENEZ
      HtoC = pass$HtoC_Ratio_PM1_AMS_JIMENEZ
      f43 = pass$f43_PM1_AMS_JIMENEZ
      f44 = pass$f44_PM1_AMS_JIMENEZ
      f57 = pass$f57_PM1_AMS_JIMENEZ
      f60 = pass$f60_PM1_AMS_JIMENEZ
      f82 = pass$f82_PM1_AMS_JIMENEZ
      f91 = pass$f91_PM1_AMS_JIMENEZ
      
    } else{ 
      OAtoOC = 1.8; OAtoOCB = 1.8
    }
    
    ind = which(cc2 == 'OC_PM1_AMS_JIMENEZ') # did we convert already?
    if (length(ind) > 0){
      pass$OC_PM1_AMS_JIMENEZ_ppb = pass$OC_PM1_AMS_JIMENEZ/pass$StdtoVol_AMS_JIMENEZ/RHS/mWOC
      background$OC_PM1_AMS_JIMENEZ_ppb = background$OC_PM1_AMS_JIMENEZ/background$StdtoVol_AMS_JIMENEZ/BG_RHS/mWOC
    } else{
      pass$OC_PM1_AMS_JIMENEZ_ppb = pass$OA_PM1_AMS_JIMENEZ/OAtoOC/pass$StdtoVol_AMS_JIMENEZ/RHS/mWOC
      background$OC_PM1_AMS_JIMENEZ_ppb = background$OA_PM1_AMS_JIMENEZ/OAtoOCB/background$StdtoVol_AMS_JIMENEZ/BG_RHS/mWOC
    }
    
    pass$OC_PM1_AMS_JIMENEZ_error = pass$OC_PM1_AMS_JIMENEZ_ppb*0.38
    p6 = smallfun(pass,background,xspecies,'OC_PM1_AMS_JIMENEZ_ppb', xerror, 'OC_PM1_AMS_JIMENEZ_error',SLOW)
    p6CO = smallfun(pass,background,'CO_DACOM_DISKIN','OC_PM1_AMS_JIMENEZ_ppb', xerror, 'OC_PM1_AMS_JIMENEZ_error',SLOW)
    
    pass$Sulfate_PM1_AMS_JIMENEZ_ppb = pass$Sulfate_PM1_AMS_JIMENEZ/pass$StdtoVol_AMS_JIMENEZ/RHS/mWSO4
    #pass$Sulfate_PM1_AMS_JIMENEZ_ppb = pass$Sulfate_PM1_AMS_JIMENEZ/convSO4/1E3
    background$Sulfate_PM1_AMS_JIMENEZ_ppb = background$Sulfate_PM1_AMS_JIMENEZ/background$StdtoVol_AMS_JIMENEZ/BG_RHS/mWSO4
    pass$Sulfate_PM1_AMS_JIMENEZ_error = pass$Sulfate_PM1_AMS_JIMENEZ_ppb * 0.34
    p6b = smallfun(pass,background,xspecies,'Sulfate_PM1_AMS_JIMENEZ_ppb', xerror, 'Sulfate_PM1_AMS_JIMENEZ_error',SLOW)
    p6bCO = smallfun(pass,background,'CO_DACOM_DISKIN','Sulfate_PM1_AMS_JIMENEZ_ppb', xerror, 'Sulfate_PM1_AMS_JIMENEZ_error',SLOW)
    
    pass$Nitrate_PM1_AMS_JIMENEZ_ppb = pass$Nitrate_PM1_AMS_JIMENEZ/pass$StdtoVol_AMS_JIMENEZ/RHS/mWNO3
    background$Nitrate_PM1_AMS_JIMENEZ_ppb = background$Nitrate_PM1_AMS_JIMENEZ/background$StdtoVol_AMS_JIMENEZ/BG_RHS/mWNO3
    pass$Nitrate_PM1_AMS_JIMENEZ_error = pass$Nitrate_PM1_AMS_JIMENEZ_ppb*0.34
    
    p6c = smallfun(pass,background,xspecies,'Nitrate_PM1_AMS_JIMENEZ_ppb', xerror, 'Nitrate_PM1_AMS_JIMENEZ_error',SLOW)
    p6cCO = smallfun(pass,background,'CO_DACOM_DISKIN','Nitrate_PM1_AMS_JIMENEZ_ppb', xerror, 'Nitrate_PM1_AMS_JIMENEZ_error',SLOW)
    
    pass$HNO3_NO3 = pass$Nitrate_PM1_AMS_JIMENEZ_ppb + (pass$HNO3.10Hz_CIT_WENNBERG/1E3) * mWNO3/mWNO3
    background$HNO3_NO3 = background$Nitrate_PM1_AMS_JIMENEZ_ppb + (background$HNO3.10Hz_CIT_WENNBERG/1E3)* mWNO3/mWNO3
    p6c2 = smallfun(pass,background,xspecies,'HNO3_NO3', xerror, 'Nitrate_PM1_AMS_JIMENEZ_error',SLOW)
    p6cCO2 = smallfun(pass,background,'CO_DACOM_DISKIN','HNO3_NO3', xerror, 'Nitrate_PM1_AMS_JIMENEZ_error',SLOW)
    
    pass$Ammonium_PM1_AMS_JIMENEZ_ppb = pass$Ammonium_PM1_AMS_JIMENEZ/pass$StdtoVol_AMS_JIMENEZ/RHS/mWNH4
    background$Ammonium_PM1_AMS_JIMENEZ_ppb = background$Ammonium_PM1_AMS_JIMENEZ/background$StdtoVol_AMS_JIMENEZ/BG_RHS/mWNH4
    pass$Ammonium_PM1_AMS_JIMENEZ_error = pass$Ammonium_PM1_AMS_JIMENEZ_ppb * 0.34
    p6d = smallfun(pass,background,xspecies,'Ammonium_PM1_AMS_JIMENEZ_ppb', xerror, 'Ammonium_PM1_AMS_JIMENEZ_error',SLOW)
    p6dCO = smallfun(pass,background,'CO_DACOM_DISKIN','Ammonium_PM1_AMS_JIMENEZ_ppb', xerror, 'Ammonium_PM1_AMS_JIMENEZ_error',SLOW)
    
    pass$NR_Chloride_PM1_AMS_JIMENEZ_ppb = pass$NR_Chloride_PM1_AMS_JIMENEZ/pass$StdtoVol_AMS_JIMENEZ/RHS/mWCl
    background$NR_Chloride_PM1_AMS_JIMENEZ_ppb = background$NR_Chloride_PM1_AMS_JIMENEZ/background$StdtoVol_AMS_JIMENEZ/BG_RHS/mWCl
    pass$NR_Chloride_PM1_AMS_JIMENEZ_error = pass$NR_Chloride_PM1_AMS_JIMENEZ_ppb*0.34
    p6e = smallfun(pass,background,xspecies,'NR_Chloride_PM1_AMS_JIMENEZ_ppb', xerror, 'NR_Chloride_PM1_AMS_JIMENEZ_error',SLOW)
    p6eCO = smallfun(pass,background,'CO_DACOM_DISKIN','NR_Chloride_PM1_AMS_JIMENEZ_ppb', xerror, 'NR_Chloride_PM1_AMS_JIMENEZ_error',SLOW)
    
    pass$Seasalt_PM1_AMS_JIMENEZ_ppb = pass$Seasalt_PM1_AMS_JIMENEZ/pass$StdtoVol_AMS_JIMENEZ/RHS/mWSSA
    background$Seasalt_PM1_AMS_JIMENEZ_ppb = background$Seasalt_PM1_AMS_JIMENEZ/background$StdtoVol_AMS_JIMENEZ/BG_RHS/mWSSA
    pass$Seasalt_PM1_AMS_JIMENEZ_error = pass$Seasalt_PM1_AMS_JIMENEZ_ppb * 0.34
    p6f = smallfun(pass,background,xspecies,'Seasalt_PM1_AMS_JIMENEZ_ppb', xerror, 'Seasalt_PM1_AMS_JIMENEZ_error',SLOW)
    p6fCO = smallfun(pass,background,'CO_DACOM_DISKIN','Seasalt_PM1_AMS_JIMENEZ_ppb', xerror, 'Seasalt_PM1_AMS_JIMENEZ_error',SLOW)
    
    pass$MSA_PM1_AMS_JIMENEZ_ppb = pass$MSA_PM1_AMS_JIMENEZ/pass$StdtoVol_AMS_JIMENEZ/RHS/mWMSA
    background$MSA_PM1_AMS_JIMENEZ_ppb = background$MSA_PM1_AMS_JIMENEZ/background$StdtoVol_AMS_JIMENEZ/BG_RHS/mWMSA
    pass$MSA_PM1_AMS_JIMENEZ_error = pass$MSA_PM1_AMS_JIMENEZ_ppb*0.34
    p6g = smallfun(pass,background,xspecies,'MSA_PM1_AMS_JIMENEZ_ppb', xerror, 'MSA_PM1_AMS_JIMENEZ_error',SLOW)
    p6gCO = smallfun(pass,background,'CO_DACOM_DISKIN','MSA_PM1_AMS_JIMENEZ_ppb', xerror, 'MSA_PM1_AMS_JIMENEZ_error',SLOW)
    
    pass$ClO4_PM1_AMS_JIMENEZ_ppb = pass$ClO4_PM1_AMS_JIMENEZ/pass$StdtoVol_AMS_JIMENEZ/RHS/mWClO4
    background$ClO4_PM1_AMS_JIMENEZ_ppb = background$ClO4_PM1_AMS_JIMENEZ/background$StdtoVol_AMS_JIMENEZ/BG_RHS/mWClO4
    pass$ClO4_PM1_AMS_JIMENEZ_error = pass$ClO4_PM1_AMS_JIMENEZ_ppb*0.34
    p6h = smallfun(pass,background,xspecies,'ClO4_PM1_AMS_JIMENEZ_ppb', xerror, 'ClO4_PM1_AMS_JIMENEZ_error',SLOW)
    p6hCO = smallfun(pass,background,'CO_DACOM_DISKIN','ClO4_PM1_AMS_JIMENEZ_ppb', xerror, 'ClO4_PM1_AMS_JIMENEZ_error',SLOW)
    
    pass$Bromine_PM1_AMS_JIMENEZ_ppb = pass$Bromine_PM1_AMS_JIMENEZ/pass$StdtoVol_AMS_JIMENEZ/RHS/mWBr
    background$Bromine_PM1_AMS_JIMENEZ_ppb = background$Bromine_PM1_AMS_JIMENEZ/background$StdtoVol_AMS_JIMENEZ/BG_RHS/mWBr
    pass$Bromine_PM1_AMS_JIMENEZ_error = pass$Bromine_PM1_AMS_JIMENEZ_ppb * 0.34
    p6i = smallfun(pass,background,xspecies,'Bromine_PM1_AMS_JIMENEZ_ppb', xerror, 'Bromine_PM1_AMS_JIMENEZ_error',SLOW)
    p6iCO = smallfun(pass,background,'CO_DACOM_DISKIN','Bromine_PM1_AMS_JIMENEZ_ppb', xerror, 'Bromine_PM1_AMS_JIMENEZ_error',SLOW)
    
    pass$Potassium_PM1_AMS_JIMENEZ_ppb = pass$Potassium_PM1_AMS_JIMENEZ/pass$StdtoVol_AMS_JIMENEZ/RHS/mWK
    background$Potassium_PM1_AMS_JIMENEZ_ppb= background$Potassium_PM1_AMS_JIMENEZ/background$StdtoVol_AMS_JIMENEZ/BG_RHS/mWK
    pass$Potassium_PM1_AMS_JIMENEZ_error = pass$Potassium_PM1_AMS_JIMENEZ_ppb * 0.34
    p6z = smallfun(pass,background,xspecies,'Potassium_PM1_AMS_JIMENEZ_ppb', xerror, 'Potassium_PM1_AMS_JIMENEZ_error',SLOW)
    p6zCO = smallfun(pass,background,'CO_DACOM_DISKIN','Potassium_PM1_AMS_JIMENEZ_ppb', xerror, 'Potassium_PM1_AMS_JIMENEZ_error',SLOW)
    
    pass$Iodine_PM1_AMS_JIMENEZ_ppb = pass$Iodine_PM1_AMS_JIMENEZ/pass$StdtoVol_AMS_JIMENEZ/RHS/mWI
    background$Iodine_PM1_AMS_JIMENEZ_ppb = background$Iodine_PM1_AMS_JIMENEZ/background$StdtoVol_AMS_JIMENEZ/BG_RHS/mWI
    pass$Iodine_PM1_AMS_JIMENEZ_error = pass$Iodine_PM1_AMS_JIMENEZ_ppb*0.34
    p6j = smallfun(pass,background,xspecies,'Iodine_PM1_AMS_JIMENEZ_ppb', xerror, 'Iodine_PM1_AMS_JIMENEZ_error',SLOW)
    p6jCO = smallfun(pass,background,'CO_DACOM_DISKIN','Iodine_PM1_AMS_JIMENEZ_ppb', xerror, 'Iodine_PM1_AMS_JIMENEZ_error',SLOW)
    
    passERs = c(passERs, p6$slope, p6b$slope, p6c$slope,p6c2$slope,  p6d$slope, p6e$slope, p6f$slope, p6g$slope, p6h$slope,p6i$slope, p6j$slope,  p6z$slope)
    passERsCO = c(passERsCO,p6CO$slope, p6bCO$slope, p6cCO$slope,p6cCO$slope,  p6dCO$slope, p6eCO$slope, p6fCO$slope, p6gCO$slope, p6hCO$slope, p6iCO$slope, p6jCO$slope, p6zCO$slope)
    passERsint = c(passERsint, p6$ERint, p6b$ERint, p6c$ERint, p6c2$ERint, p6d$ERint, p6e$ERint, p6f$ERint, p6g$ERint, p6h$ERint,  p6i$ERint, p6j$ERint,  p6z$ERint)
    passERsintfill = c(passERsintfill, p6$ERintfill, p6b$ERintfill, p6c$ERintfill,p6c2$ERintfill, p6d$ERintfill, p6e$ERintfill, p6f$ERintfill, p6g$ERintfill, p6h$ERintfill,
                         p6i$ERintfill, p6j$ERintfill,p6z$ERintfill)
    passERsCOint = c(passERsCOint, p6CO$ERint, p6bCO$ERint, p6cCO$ERint, p6cCO2$ERint, p6dCO$ERint, p6eCO$ERint, p6fCO$ERint, p6gCO$ERint, p6hCO$ERint, p6iCO$ERint, p6jCO$ERint,p6z$ERintfill)
    passERsCOintfill = c(passERsCOintfill, p6CO$ERintfill, p6bCO$ERintfill, p6cCO$ERintfill,p6cCO2$ERintfill, p6dCO$ERintfill, p6eCO$ERintfill, p6fCO$ERintfill, p6gCO$ERintfill, p6hCO$ERintfill,
                         p6iCO$ERintfill, p6jCO$ERintfill,p6zCO$ERintfill)
    
    passBGX = c(passBGX, p6$BGX, p6b$BGX, p6c$BGX,p6c2$BGX, p6d$BGX, p6e$BGX, p6f$BGX, p6g$BGX, p6h$BGX, p6i$BGX, p6j$BGX,  p6z$BGX)
    passBGCO = c(passBGCO, p6CO$BGX, p6bCO$BGX, p6cCO$BGX, p6cCO2$BGX, p6dCO$BGX, p6eCO$BGX, p6fCO$BGX, p6gCO$BGX, p6hCO$BGX, p6iCO$BGX, p6jCO$BGX, p6zCO$BGX)
    
    passBGS = c(passBGS, p6$BGS, p6b$BGS, p6c$BGS, p6c2$BGS, p6d$BGS, p6e$BGS, p6f$BGS, p6g$BGS, p6h$BGS, p6i$BGS, p6j$BGS,p6z$BGS)
    
    passERszero = c( passERszero,p6$intercept, p6b$intercept, p6c$intercept, p6c2$intercept, p6d$intercept, p6e$intercept, p6f$intercept, p6g$intercept, p6h$intercept,p6i$intercept, p6j$intercept,  p6z$intercept)
    passERsCOzero = c( passERsCOzero,p6CO$intercept, p6bCO$intercept, p6cCO$intercept, p6cCO2$intercept,p6dCO$intercept, p6eCO$intercept, p6fCO$intercept, p6gCO$intercept, p6hCO$intercept,
                       p6iCO$intercept, p6jCO$intercept, p6zCO$intercept)
    passRsqsCO2 = c(passRsqsCO2,p6$r_sq, p6b$r_sq, p6c$r_sq, p6c2$r_sq, p6d$r_sq, p6e$r_sq, p6f$r_sq, p6g$r_sq, p6h$r_sq,p6i$r_sq, p6j$r_sq, p6z$r_sq)
    passRsqsCO = c(passRsqsCO,p6CO$r_sq, p6bCO$r_sq, p6cCO$r_sq,  p6cCO2$r_sq, p6dCO$r_sq, p6eCO$r_sq, p6fCO$r_sq, p6gCO$r_sq, p6hCO$r_sq,p6iCO$r_sq, 
                   p6jCO$r_sq, p6zCO$r_sq)
    passExcessX  = c(passExcessX,p6$ExcessX, p6b$ExcessX, p6c$ExcessX,p6c2$ExcessX, p6d$ExcessX, p6e$ExcessX, p6f$ExcessX, p6g$ExcessX, p6h$ExcessX,
                     p6i$ExcessX, p6j$ExcessX,p6z$ExcessX)
    passExcessS  = c(passExcessS,p6$ExcessS, p6b$ExcessS, p6c$ExcessS, p6c2$ExcessS, p6d$ExcessS, p6e$ExcessS, p6f$ExcessS, p6g$ExcessS, p6h$ExcessS,
                     p6i$ExcessS, p6j$ExcessS,p6z$ExcessS)
    
    passsigma = c(passsigma,p6$uncertainty_slope, p6b$uncertainty_slope, p6c$uncertainty_slope,p6c2$uncertainty_slope, p6d$uncertainty_slope, p6e$uncertainty_slope, p6f$uncertainty_slope, p6g$uncertainty_slope, p6h$uncertainty_slope,
                  p6i$uncertainty_slope, p6j$uncertainty_slope,  p6z$uncertainty_slope)
    passMAX   = c( passMAX, maxfun(pass, 'OA_PM1_AMS_JIMENEZ'),maxfun(pass, 'Sulfate_PM1_AMS_JIMENEZ'),maxfun(pass, 'Nitrate_PM1_AMS_JIMENEZ'),
                   maxfun(pass, 'HNO3_NO3'), maxfun(pass, 'Ammonium_PM1_AMS_JIMENEZ'), maxfun(pass, 'NR_Chloride_PM1_AMS_JIMENEZ'), maxfun(pass, 'Seasalt_PM1_AMS_JIMENEZ'),
                   maxfun(pass,'MSA_PM1_AMS_JIMENEZ'), maxfun(pass,'ClO4_PM1_AMS_JIMENEZ'), maxfun(pass,'Bromine_PM1_AMS_JIMENEZ'), maxfun(pass,'Iodine_PM1_AMS_JIMENEZ'),
                    maxfun(pass,'Potassium_PM1_AMS_JIMENEZ'))
    passErrors = c(passErrors,pass$OA_PM1_AMS_JIMENEZ_error, pass$Sulfate_PM1_AMS_JIMENEZ_error, pass$Nitrate_PM1_AMS_JIMENEZ_error,pass$Nitrate_PM1_AMS_JIMENEZ_error, pass$Ammonium_PM1_AMS_JIMENEZ_error, pass$Chloride_error, pass$Seasalt_error, 
                   pass$MSA_PM1_AMS_JIMENEZ_error, pass$ClO4_PM1_AMS_JIMENEZ_error, pass$Bromine_PM1_AMS_JIMENEZ_error,
                   pass$Iodine_PM1_AMS_JIMENEZ_error,pass$Potassium_PM1_AMS_JIMENEZ_error)
    passOHrate = c(passOHrate, NaN, NaN,NaN, NaN, NaN, NaN, NaN, NaN, NaN,NaN,NaN,NaN)
    isNMVOC = c(isNMVOC,4,4,4,4,4,4,4,4,4,4,4,4)
    kind = c(kind,"aerosol","aerosol","aerosol","aerosol","aerosol","aerosol","aerosol","aerosol","aerosol","aerosol","aerosol","aerosol")
    
    names1=c(names1,"OC_JIMENEZ","Sulfate_JIMENEZ","Nitrate_JIMENEZ","HNO3_NO3", "Ammonium_JIMENEZ","NR_Chloride_JIMENEZ",
             "Seasalt_JIMENEZ", "MSA_JIMENEZ","ClO4_JIMENEZ","Bromine_JIMENEZ", "Iodine_JIMENEZ","Potassium_JIMENEZ")
    formulas = c(formulas,"OC","SO4","NO3","NO3","NH4","Cl","Na","MSA","ClO4","Br","I","K")
    usenames = c(usenames,"Organic carbon","Sulfate","Nitrate","HNO3_NO3","Ammonium","Chloride","Seasalt","Methanesulfonic acid",
                 "Perchlorate","Bromine","Iodine","Potassium")
    mWs = c(mWs, mWOC, mWSO4, mWNO3,mWNO3, mWNH4, mWCl, mWSSA, mWMSA, mWClO4, mWBr, mWI, mWK)
    nCs  = c(nCs,  1,     0,     0,   0,  0,     0,     0,   1,       0,     0,   0,  0)
     
  }
  if (debugKT == 1){print(c(length(passERs), length(passERsint), length(passERsintfill),length(passERsCO),length(passsigma), length(passRsqsCO2), 
                            length(passRsqsCO),length(passMAX),length(passOHrate),length(passErrors),
                            length(names1), length(mWs), length(nCs), length(passBGS), length(passBGX), 
                            length(passExcessX), length(passExcessS), length(passERsCOint), length(passERsCOintfill),
                            length(passERszero), length(passERsCOzero), length(isNMVOC), length(kind)))}
  
  # ------------------------------------------ WARNEKE ------------------------------------------ 
  if (debugKT ==1){ print("WARNEKE")}
  wnames = c('HCN_NOAAPTR_ppbv_WARNEKE', 'CH2O_NOAAPTR_ppbv_WARNEKE','CH3OH_NOAAPTR_ppbv_WARNEKE','CH3CN_NOAAPTR_ppbv_WARNEKE','HNCO_NOAAPTR_ppbv_WARNEKE',
             'CH3CHO_NOAAPTR_ppbv_WARNEKE','C2H5OH_NOAAPTR_ppbv_WARNEKE','HCOOH_NOAAPTR_ppbv_WARNEKE','Acrylonitrile_NOAAPTR_ppbv_WARNEKE','Acrolein_NOAAPTR_ppbv_WARNEKE',
             'AcetonePropanal_NOAAPTR_ppbv_WARNEKE','GlycolaldehydeCH3COOH_NOAAPTR_ppbv_WARNEKE','CH3NO2_NOAAPTR_ppbv_WARNEKE','DMS_NOAAPTR_ppbv_WARNEKE','C4H5N_NOAAPTR_ppbv_WARNEKE',
             'Furan_NOAAPTR_ppbv_WARNEKE','MVKMAC_NOAAPTR_ppbv_WARNEKE','C4Carbonyls_NOAAPTR_ppbv_WARNEKE','C3H6O2_NOAAPTR_ppbv_WARNEKE','Benzene_NOAAPTR_ppbv_WARNEKE',
             'x2MeFuranx3MeFuran_NOAAPTR_ppbv_WARNEKE','x2Furanone_NOAAPTR_ppbv_WARNEKE','x23Butanedione_NOAAPTR_ppbv_WARNEKE','Toluene_NOAAPTR_ppbv_WARNEKE','Phenol_NOAAPTR_ppbv_WARNEKE',
             'Furfural_NOAAPTR_ppbv_WARNEKE','DimeFurans_NOAAPTR_ppbv_WARNEKE','MaleicAnhyd_NOAAPTR_ppbv_WARNEKE','BenzNitrile_NOAAPTR_ppbv_WARNEKE','Styrene_NOAAPTR_ppbv_WARNEKE',
             'Benzaldehyde_NOAAPTR_ppbv_WARNEKE','C8Aromatics_NOAAPTR_ppbv_WARNEKE','C7H8O_NOAAPTR_ppbv_WARNEKE','Catecholx5MeFurfural_NOAAPTR_ppbv_WARNEKE','BenzFuran_NOAAPTR_ppbv_WARNEKE',
             'C9Aromatics_NOAAPTR_ppbv_WARNEKE','C6H4O3_NOAAPTR_ppbv_WARNEKE','Guaiacol_NOAAPTR_ppbv_WARNEKE','Naphthalene_NOAAPTR_ppbv_WARNEKE','Monoterpenes_NOAAPTR_ppbv_WARNEKE',
             'Creosols_NOAAPTR_ppbv_WARNEKE','Syringol_NOAAPTR_ppbv','Isoprene_NOAAPTR_ppbv_WARNEKE')
  isNMVOC = c(isNMVOC, 2, 1,1,2,2, 1,1,1,2,1,1,1,2,6,2,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1, 1,1,1,1,1,1,1,1,1,1,1,1,1)
  kind = c(kind,'nitrogen', 'CH2O','oVOC','nitrogen','nitrogen',
           'oVOC','oVOC','oVOC','nitrogen','oVOC',
           'oVOC','oVOC','nitrate','sulfur','nitrogen',
           'oVOC','oVOC','oVOC','oVOC','aromatic',
           'oVOC','oVOC','oVOC','aromatic','oVOC',
           'oVOC','oVOC','oVOC','nitrogen','aromatic',
           'oVOC','aromatic','oVOC','oVOC','oVOC',
           'aromatic','oVOC','oVOC','PAH','alkene',
           'oVOC','oVOC','alkene')
  # get ALL_WARNEKE
  #ind = which(colnames(pass)  %in% wnames)
  #tmpall = pass[,ind]
  #pass$ALL_WARNEKE = rowSums(tmpall, na.rm=TRUE)
  
  werrors = c(0.30,0.20,0.20,0.20,0.30,0.20,0.30,0.30,0.30,0.50,0.20,0.50,0.30,0.30,0.50,0.30,0.20,0.30,0.50,0.20,0.30,0.50,0.50,0.20,0.30,0.30,0.30,0.50,0.30,0.30,0.30,0.20,
              0.50,0.50,0.50,0.20,0.50,0.50,0.30,0.30,0.50,0.50,0.30)
  passErrors=c(passErrors,werrors)
  
  wform = c('HCN',	'CH2O',	'CH3OH',	'CH3CN',	'HNCO',	'C2H4O',	'C2H5OH','HCOOH',	'C3H3N','C3H4O','C3H6O','C2H4O2','CH3NO2','C2H6S','C4H5N','C4H4O',
            'C4H6O',	'C4H8O','C3H6O2',	'C6H6',	'C5H6O', 'C4H4O2',	'C4H6O2','C7H8','C6H6O','C5H4O2',	'C6H8O','C4H2O3','C7H5N','C8H8','C7H6O','C8H10',
            'C7H8O',	'C6H6O2','C8H6O','C9H12','C6H4O3','C7H8O2','C10H8',	'C10H16','C8H10O2','C8H10O3','C5H8')
  formulas=c(formulas,wform)
  usenames=c(usenames, 'Hydrogen cyanide', 'Formaldehyde','Methanol','Acetonitrile','Isocyanic acid','Acetaldehyde', 'Ethanol','Formic acid', 'Acrylonitrile',
             'Acrolein','Acetone/Propanal','Acetic acid/Glycolaldehyde',
             'Nitromethane','Dimethyl sulfide','Pyrrole/Butenenitrile',' Furan and fragments','MVK/MACR',
             'MEK/ 2-methyl propanal', 'Methyl acetate/Ethyl formate/Hydroxyacetone','Benzene',
             'sum of 2-methylfuran 3-methylfuran and fragments',
             '2(3H)-Furanone',
             '2,3-Butanedione/2-Oxobutanal/1,4-Butanedial','Toluene','Phenol', 'sum of 2-furfural (dominant isomer) 3-furfural and fragments',
             '2,5-Dimethylfuran/2-Ethylfuran/Other unidentified organic compounds',
             'Maleic anhydride','Benzonitrile','Styrene','Benzaldehyde', 
             'sum of m-xylene p-xylene o-xylene and ethyl benzene','2-Methylphenol (=o-cresol)/Anisol',
             '5-Methylfurfural + Benzene diols (=Catechol, Resorcinol)','Benzofuran','C9 aromatics','Hydroxybenzoquinone','Guaiacol','Naphthalene','Monoterpenes','Creosols','Syringol','Isoprene')
  
  wOH = c(OHHCN,    OHch2o,    OHCH3OH,    OHCH3CN,    OHhnco,    OHCH3CHO,    OHC2H5OH,    OHhcooh,
          OHAcrylonitrile,    OHAcrolein,    OHAcetonePropanal,    OHGlycolaldehydeCH3COOH,
          OHCH3NO2,    OHDMS,    OHC4H5N,    OHFuran,    OHMVKMAC,    OHMEK,#C4Carbonyls,
          OHC3H6O2,    OHBenzene,    OHx2MeFuranx3MeFuran,    OHx2Furanone,    OHx23Butanedione,
          OHToluene,    OHPhenol,    OHFurfural,    OHDimeFurans,    OHMaleicAnhyd,    OHBenzNitrile,    OHStyrene,  
          OHBenzaldehyde,    OHC8Aromatics,    OHC7H8O,
          OHCatecholx5MeFurfural,    OHBenzFuran,    OHC9Aromatics,
          OHC6H4O3,    OHGuaiacol,    OHNaphthalene,    OHMonoterpenes,    OHCreosols,    OHSyringol,OHIsoprene)
  
  wtau = 1/(5E6*wOH)/60/60 # hours
  # All from Koss et al., 2018 supplement except CH3OH from MCM at 298 K,
  # maleic anhydride, from https://doi.org/10.1002/kin.21387
  #wMWs = c(27.0,30.0,32.0,41.1,43.0,44.0,46.1,46.0,55.1,56.1,58.1,60.0,61.0,62.1,67.1,68.1,70.1,74.1,74.1,78.1,82.1,
  #         84.1,86.1,92.1,94.1,96.1,96.1,98.1,103.1,104.1,106.1,106.2,108.1,110.1,118.1,120.2,124.1,124.1,128.2,136.2,138.1,154.1)
  wMWs = c(mWHCN, mWCH2O,mWCH3OH,mWch3cn, mWHNCO,
             mWAcetaldehyde, mWC2H5OH, mWHCOOH, mWAcrylonitrile, mWAcrolein,
             mWAcetonePropanal, mWGlycolaldehydeCH3COOH, mWNitromethane, mWDMS, mwPyrrole,
             mWFuran, mWMVK, mWMEK, mWHAC, mWBenzene,
             mWx2MeFuran, mW2Furanone, mW23butanedione, mWToluene, mWPHENOL,
             mWfurfural, mWDimeFuran, mWMaleicAnhyd, mWBenzNitrile, mWStyrene,
             mWBenzaldehyde, mWC8Aromatics, mWC7H8O, mWCatecholx5MeFurfural, mWBenzFuran,
            mWC9Aromatics, mWC6H4O3, mWGuaiacol, mWNaphthalene,mWMonoterpenes,
             mWCresol,mWSyringol, mWIsoprene)
  
  ccs = colnames(pass)
  nn = length(wnames) 
  for (i in 1:nn){ # don't do ALL_WARNEKE
    #if (debugKT == 1){print(wnames[i])}
    ind = which(ccs == wnames[i])
    if (length(ind) == 0){print("ERROR")}
    pass$tmperror = pass[,ind]*werrors[i]
    pTMP = smallfun(pass,background,xspecies,wnames[i], xerror, 'tmperror',SLOW)
    pTMPCO = smallfun(pass,background,'CO_DACOM_DISKIN', wnames[i], 'CO_error','tmperror',SLOW)
    
    if (i == 1){
      passERsW = pTMP$slope
      passERsWCO = pTMPCO$slope
      passERsWI = pTMP$ERint
      passERsWIfill = pTMP$ERintfill
      passERsCOWI = pTMPCO$ERint
      passERsCOWIfill = pTMPCO$ERintfill
      passRWCO2 = pTMP$BGX
      passRWCOB = pTMPCO$BGX
      
      passRWS   = pTMP$BGS
      passRW   = pTMP$r_sq
      passRWCO   = pTMPCO$r_sq
      passRWQ = pTMP$intercept
      passRWQCO = pTMPCO$intercept
      passRWXX = pTMP$ExcessX
      passRWSS = pTMP$ExcessS
      passUWK   = pTMP$uncertainty_slope
      passUWM   = maxfun(pass, wnames[i])
    }
    if (i > 1){ 
      passERsW = c(passERsW,pTMP$slope)
      passERsWCO =c(passERsWCO, pTMPCO$slope)
      passERsWI =c(passERsWI,pTMP$ERint)
      passERsWIfill =c(passERsWIfill,pTMP$ERintfill)
      passERsCOWI = c(passERsCOWI,pTMPCO$ERint)
      passERsCOWIfill = c(passERsCOWIfill,pTMPCO$ERintfill)
      passRWCO2 = c(passRWCO2,pTMP$BGX)
      passRWS   = c(passRWS,pTMP$BGS)
      passRWCOB = c(passRWCOB,pTMPCO$BGX)
      
      passRW   = c(passRW,pTMP$r_sq)
      passRWCO   = c(passRWCO,pTMPCO$r_sq)
      passRWXX = c(passRWXX,pTMP$ExcessX)
      passRWSS = c(passRWSS,pTMP$ExcessS)
      passRWQ = c(passRWQ,pTMP$intercept)
      passRWQCO = c(passRWQCO,pTMPCO$intercept)
      passUWK   = c(passUWK,pTMP$uncertainty_slope)
      passUWM   = c(passUWM,maxfun(pass, wnames[i]))
    }
  }
  # make a total VOC ER/EF
  passERs   = c(passERs, passERsW)
  passERsCO = c(passERsCO,passERsWCO)
  passERsint       = c(passERsint, passERsWI)
  passERsintfill   = c(passERsintfill, passERsWIfill)    
  passERsCOint     = c(passERsCOint, passERsCOWI)
  passERsCOintfill = c(passERsCOintfill, passERsCOWIfill )
  passBGX          = c(passBGX, passRWCO2)
  passBGCO         = c(passBGCO,passRWCOB) 
  passBGS          = c(passBGS, passRWS)
  passERszero      = c(passERszero, passRWQ)
  passERsCOzero    = c(passERsCOzero, passRWQCO)
  passRsqsCO2      = c(passRsqsCO2,passRW)
  passRsqsCO       = c(passRsqsCO,passRWCO)
  passExcessS = c(passExcessS,passRWSS)
  passExcessX = c(passExcessX,passRWXX)
  passsigma = c(passsigma,passUWK)
  passMAX   = c(passMAX, passUWM)
  passOHrate = c(passOHrate, wOH)
  
  names1    = c(names1,wnames)
  mWs       = c(mWs, wMWs)
  nCs = c(nCs,1,1,1,2,1,2,2,1,3,3,3,2,1,2,4,4,4,4,3,6,5,4,4,7,6,5,6,4,7,8,7,8,7,6,8,9,6,7,10,10,8,8,5)
  
  
  if (debugKT == 1){print(c(length(passERs), length(passERsint), length(passERsintfill),length(passERsCO),length(passsigma), length(passRsqsCO2), length(passRsqsCO),length(passMAX),length(passOHrate),
                            length(names1), length(mWs), length(nCs), length(passBGS), length(passBGX), 
                            length(passExcessX), length(passExcessS), length(passERsCOint), length(passERsCOintfill),
                            length(passERszero), length(passERsCOzero), length(isNMVOC), length(kind), length(formulas),length(usenames)))}
  # ------------------------------------------ HANISCO ------------------------------------------ 
  if (debugKT == 1){print("HANISCO")}
  #10% + 10 pptv
  pass$CH2O_ISAF_HANISCO_ppb = pass$CH2O_ISAF_HANISCO/1E3
  background$CH2O_ISAF_HANISCO_ppb = background$CH2O_ISAF_HANISCO/1E3
  pass$tmperror = (pass$CH2O_ISAF_HANISCO*0.1+10)/1E3
  passErrors=c(passErrors,mean(pass$tmperror, na.rm=TRUE))
  
  # pass$tmperror = (pass$CH2O_ISAF_precision_HANISCO)/1E3
  
  pT = smallfun(pass,background,xspecies,'CH2O_ISAF_HANISCO_ppb', xerror, 'tmperror',SLOW)
  pTCO = smallfun(pass,background,'CO_DACOM_DISKIN','CH2O_ISAF_HANISCO_ppb', xerror, 'tmperror',SLOW)
  passERs = c(passERs, pT$slope)    
  passERsCO = c(passERsCO, pTCO$slope)    
  passERsint = c(passERsint, pT$ERint)    
  passERsintfill = c(passERsintfill, pT$ERintfill)    
  passERsCOint = c(passERsCOint, pTCO$ERint)    
  passERsCOintfill = c(passERsCOintfill, pTCO$ERintfill)  
  passBGX  = c(passBGX, pT$BGX)
  passBGCO  = c(passBGCO, pTCO$BGX)
  passBGS  = c(passBGS, pT$BGS)
  passRsqsCO2 = c(passRsqsCO2,pT$r_sq)
  passRsqsCO = c(passRsqsCO,pTCO$r_sq)
  passExcessX = c(passExcessX,pT$ExcessX)
  passExcessS = c(passExcessS,pT$ExcessS)
  passERszero= c(passERszero,pT$intercept)
  passERsCOzero = c(passERsCOzero,pTCO$intercept)
  
  passsigma = c(passsigma, pT$uncertainty_slope)
  passMAX   = c(passMAX, maxfun(pass,'CH2O_ISAF_HANISCO_ppb' ))
  passOHrate = c(passOHrate, 5.4E-12*exp(135/TEMP))
  isNMVOC = c(isNMVOC,1)
  kind = c(kind,'CH2O')
  names1=c(names1, 'CH2O_ISAF_HANISCO')
  formulas=c(formulas,'CH2O')
  usenames=c(usenames,'Formaldehyde')
  mWs       = c(mWs, mWCH2O)
  nCs = c(nCs,1)
  # ------------------------------------------ ROLLINS ------------------------------------------ 
  if (debugKT == 1){print("ROLLINS")}
  pass$SO2_LIF_ROLLINS_ppb = pass$SO2_LIF_ROLLINS/1E3
  background$SO2_LIF_ROLLINS_ppb = background$SO2_LIF_ROLLINS/1E3
  pass$tmperror = (pass$SO2_LIF_ROLLINS*0.09 + 2)/1E3
  passErrors=c(passErrors,pass$tmperror)
  
  pT = smallfun(pass,background,xspecies,'SO2_LIF_ROLLINS_ppb', xerror, 'tmperror',SLOW)
  pTCO = smallfun(pass,background,'CO_DACOM_DISKIN','SO2_LIF_ROLLINS_ppb', xerror, 'tmperror',SLOW)
  passERs = c(passERs, pT$slope)    
  passERsCO = c(passERsCO, pTCO$slope)    
  passERsint = c(passERsint, pT$ERint)    
  passERsintfill = c(passERsintfill, pT$ERintfill)    
  passERsCOint = c(passERsCOint, pTCO$ERint)    
  passERsCOintfill = c(passERsCOintfill, pTCO$ERintfill)  
  passBGX  = c(passBGX, pT$BGX)
  passBGCO  = c(passBGCO, pTCO$BGX)
  passBGS  = c(passBGS, pT$BGS)
  passRsqsCO2 = c(passRsqsCO2,pT$r_sq)
  passRsqsCO = c(passRsqsCO,pTCO$r_sq)
  passExcessX = c(passExcessX,pT$ExcessX)
  passExcessS = c(passExcessS,pT$ExcessS)
  passERszero= c(passERszero,pT$intercept)
  passERsCOzero = c(passERsCOzero,pTCO$intercept)
  passsigma = c(passsigma, pT$uncertainty_slope)
  passMAX   = c(passMAX, maxfun(pass,'SO2_LIF_ROLLINS_ppb' ))
  names1=c(names1, 'SO2_LIF_ROLLINS')
  formulas=c(formulas,'SO2')
  usenames=c(usenames,'Sulfur dioxide')
  passOHrate = c(passOHrate, NaN)
  mWs       = c(mWs, mWSO2)
  nCs = c(nCs,0)
  isNMVOC = c(isNMVOC,6)
  kind = c(kind,'sulfur')
  
  pass$NO_LIF_ROLLINS_ppb = pass$NO_LIF_ROLLINS/1E3
  background$NO_LIF_ROLLINS_ppb = background$NO_LIF_ROLLINS/1E3
  pass$tmperror = (pass$NO_LIF_ROLLINS*0.08)/1E3
  passErrors=c(passErrors,pass$tmperror)
  
  pT = smallfun(pass,background,xspecies,'NO_LIF_ROLLINS_ppb', xerror, 'tmperror',SLOW)
  pTCO = smallfun(pass,background,'CO_DACOM_DISKIN','NO_LIF_ROLLINS_ppb', xerror, 'tmperror',SLOW)
  passERs = c(passERs, pT$slope)    
  passERsCO = c(passERsCO, pTCO$slope)    
  passERsint = c(passERsint, pT$ERint)    
  passERsintfill = c(passERsintfill, pT$ERintfill)   
  passERsCOint = c(passERsCOint, pTCO$ERint)    
  passERsCOintfill = c(passERsCOintfill, pTCO$ERintfill)  
  passBGX  = c(passBGX, pT$BGX)
  passBGCO  = c(passBGCO, pTCO$BGX)
  passBGS  = c(passBGS, pT$BGS)
  passRsqsCO2 = c(passRsqsCO2,pT$r_sq)
  passRsqsCO = c(passRsqsCO,pTCO$r_sq)
  passERszero= c(passERszero,pT$intercept)
  passERsCOzero = c(passERsCOzero,pTCO$intercept)
  passExcessX = c(passExcessX,pT$ExcessX)
  passExcessS = c(passExcessS,pT$ExcessS)
  passsigma = c(passsigma, pT$uncertainty_slope)
  passMAX   = c(passMAX, maxfun(pass,'NO_LIF_ROLLINS_ppb' ))
  passOHrate = c(passOHrate, NaN)
  names1=c(names1, 'NO_LIF_ROLLINS')
  formulas=c(formulas,'NO')
  usenames=c(usenames,'Nitrogen oxide')
  isNMVOC = c(isNMVOC,2)
  mWs       = c(mWs, mWNO)
  kind = c(kind,'NOy')
  
  nCs = c(nCs,0)
  
  # ------------------------------------------ WENNBERG ------------------------------------------ 
  if (debugKT == 1){print("WENNBERG")}
  print('BUTENEHN')
  #BUTENE-HN-1Hz_CIT_WENNBERG#C4O4H9N, +-(25% of measurement value + 3 pptv),  pptv
  if (res == 1){
    pass$tmperror = (pass$BUTENE.HN.1Hz_CIT_WENNBERG*0.25 + 3)/1E3
  } else{pass$tmperror = (pass$BUTENE.HN.10Hz_CIT_WENNBERG*0.25 + 3)/1E3}
  passErrors=c(passErrors,pass$tmperror)
  if (res == 1){
    pass$tmpy = pass$BUTENE.HN.1Hz_CIT_WENNBERG/1E3
    background$tmpy = background$BUTENE.HN.1Hz_CIT_WENNBERG/1E3
  }else {
    pass$tmpy = pass$BUTENE.HN.10Hz_CIT_WENNBERG/1E3
    background$tmpy = background$BUTENE.HN.10Hz_CIT_WENNBERG/1E3
  }
  pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
  pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
  passERs = c(passERs, pW1$slope)
  passERsint = c(passERsint, pW1$ERint)    
  passERsintfill = c(passERsintfill, pW1$ERintfill)    
  passBGX  = c(passBGX, pW1$BGX)
  passBGCO  = c(passBGCO, pW1CO$BGX)
  passBGS  = c(passBGS, pW1$BGS)
  passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
  passExcessX = c(passExcessX,pW1$ExcessX)
  passExcessS = c(passExcessS,pW1$ExcessS)
  passsigma = c(passsigma, pW1$uncertainty_slope)
  passERszero= c(passERszero,pW1$intercept)
  passERsCO = c(passERsCO, pW1CO$slope)
  passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
  passERsCOint = c(passERsCOint, pW1CO$ERint)    
  passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
  passERsCOzero = c(passERsCOzero,pW1CO$intercept)
  passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
  passOHrate = c(passOHrate, NaN)
  names1=c(names1,'BUTENE.HN_WENNBERG')
  kind = c(kind,'nitrate')
  formulas=c(formulas,'C4H9NO4')
  usenames=c(usenames,'Butene Hydroxynitrates')
  
  isNMVOC = c(isNMVOC,2)
  mWs   = c(mWs, mWBUTENEHN)
  nCs  = c(nCs, 4)
  
  print('BUTENEHP')
  #   #BUTENE-HP-1Hz_CIT_WENNBERG#C4O3H10, +-(30% of measurement value + 7 pptv),  pptv
  if (res == 1){
    pass$tmperror = (pass$BUTENE.HP.1Hz_CIT_WENNBERG*0.3 + 7)/1E3
  } else{ pass$tmperror = (pass$BUTENE.HP.10Hz_CIT_WENNBERG*0.3 + 7)/1E3}
  passErrors=c(passErrors,pass$tmperror)
  if (res == 1){
    pass$tmpy = pass$BUTENE.HP.1Hz_CIT_WENNBERG/1E3
    background$tmpy = background$BUTENE.HP.1Hz_CIT_WENNBERG/1E3
  } else{
    pass$tmpy = pass$BUTENE.HP.10Hz_CIT_WENNBERG/1E3
    background$tmpy = background$BUTENE.HP.10Hz_CIT_WENNBERG/1E3
  }
  pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
  pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
  passERs = c(passERs, pW1$slope)
  passERsint = c(passERsint, pW1$ERint)    
  passERsintfill = c(passERsintfill, pW1$ERintfill)    
  passBGX  = c(passBGX, pW1$BGX)
  passBGCO  = c(passBGCO, pW1CO$BGX)
  passBGS  = c(passBGS, pW1$BGS)
  passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
  passExcessX = c(passExcessX,pW1$ExcessX)
  passExcessS = c(passExcessS,pW1$ExcessS)
  passsigma = c(passsigma, pW1$uncertainty_slope)
  passERszero= c(passERszero,pW1$intercept)
  passERsCO = c(passERsCO, pW1CO$slope)
  passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
  passERsCOint = c(passERsCOint, pW1CO$ERint)    
  passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
  passERsCOzero = c(passERsCOzero,pW1CO$intercept)
  passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
  passOHrate = c(passOHrate, OHbuteneHP)
  names1=c(names1,'BUTENE.HP_WENNBERG')
  kind = c(kind,'oVOC')
  formulas=c(formulas,'C4H10O3')
  usenames=c(usenames,'C4 Hydroxyperoxide')
  isNMVOC = c(isNMVOC,1)
  mWs   = c(mWs, mWBUTENEHP)
  nCs  = c(nCs, 4)
  
  #ETHENE-HN-1Hz_CIT_WENNBERG# C2O4H5N,+-(25% of measurement value + 7 pptv)  pptv
  if (res == 1){
    pass$tmperror = (pass$ETHENE.HN.1Hz_CIT_WENNBERG*0.25 + 7)/1E3
  } else{    pass$tmperror = (pass$ETHENE.HN.10Hz_CIT_WENNBERG*0.25 + 7)/1E3}
  passErrors=c(passErrors,pass$tmperror)
  if (res == 1){
    pass$tmpy = pass$ETHENE.HN.1Hz_CIT_WENNBERG/1E3
    background$tmpy = background$ETHENE.HN.1Hz_CIT_WENNBERG/1E3
  } else{
    pass$tmpy = pass$ETHENE.HN.10Hz_CIT_WENNBERG/1E3
    background$tmpy = background$ETHENE.HN.10Hz_CIT_WENNBERG/1E3
  }
  pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
  pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
  passERs = c(passERs, pW1$slope)
  passERsint = c(passERsint, pW1$ERint)    
  passERsintfill = c(passERsintfill, pW1$ERintfill)    
  passBGX  = c(passBGX, pW1$BGX)
  passBGCO  = c(passBGCO, pW1CO$BGX)
  passBGS  = c(passBGS, pW1$BGS)
  passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
  passExcessX = c(passExcessX,pW1$ExcessX)
  passExcessS = c(passExcessS,pW1$ExcessS)
  passsigma = c(passsigma, pW1$uncertainty_slope)
  passERszero= c(passERszero,pW1$intercept)
  passERsCO = c(passERsCO, pW1CO$slope)
  passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
  passERsCOint = c(passERsCOint, pW1CO$ERint)    
  passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
  passERsCOzero = c(passERsCOzero,pW1CO$intercept)
  passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
  passOHrate = c(passOHrate, NaN)
  names1=c(names1,  'ETHENE.HN_WENNBERG')
  kind = c(kind,'nitrate')
  formulas=c(formulas,'C2O4H5N')
  usenames=c(usenames,'Ethene Hydroxynitrate')
  
  isNMVOC = c(isNMVOC,2)
  mWs   = c(mWs, mWETHENEHN)
  nCs  = c(nCs, 2)
  print('ETHENEHN')
  
  #ETHENE-HP-1Hz_CIT_WENNBERG#C2O3H6, +-(30% of measurement value + 9 pptv),  pptv
  if (res == 1){
    pass$tmperror = (pass$ETHENE.HP.1Hz_CIT_WENNBERG*0.3 + 9)/1E3
  } else{pass$tmperror = (pass$ETHENE.HP.10Hz_CIT_WENNBERG*0.3 + 9)/1E3}
  passErrors=c(passErrors,pass$tmperror)
  if (res == 1){
    pass$tmpy = pass$ETHENE.HP.1Hz_CIT_WENNBERG/1E3
    background$tmpy = background$ETHENE.HP.1Hz_CIT_WENNBERG/1E3
  } else{
    pass$tmpy = pass$ETHENE.HP.10Hz_CIT_WENNBERG/1E3
    background$tmpy = background$ETHENE.HP.10Hz_CIT_WENNBERG/1E3
  }
  pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
  pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
  passERs = c(passERs, pW1$slope)
  passERsint = c(passERsint, pW1$ERint)    
  passERsintfill = c(passERsintfill, pW1$ERintfill)    
  passBGX  = c(passBGX, pW1$BGX)
  passBGCO  = c(passBGCO, pW1CO$BGX)
  
  passBGS  = c(passBGS, pW1$BGS)
  passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
  passExcessX = c(passExcessX,pW1$ExcessX)
  passExcessS = c(passExcessS,pW1$ExcessS)
  passsigma = c(passsigma, pW1$uncertainty_slope)
  passERszero= c(passERszero,pW1$intercept)
  passERsCO = c(passERsCO, pW1CO$slope)
  passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
  passERsCOint = c(passERsCOint, pW1CO$ERint)    
  passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
  passERsCOzero = c(passERsCOzero,pW1CO$intercept)
  passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
  passOHrate = c(passOHrate, OHetheneHP)
  names1=c(names1,'ETHENE.HP_WENNBERG')
  kind = c(kind,'oVOC')
  formulas=c(formulas,'C2H6O3')
  usenames=c(usenames,'Ethene Hydroxyperoxide')
  isNMVOC = c(isNMVOC,1)
  mWs   = c(mWs, mWETHENEHP)
  nCs  = c(nCs, 2)
  print('ETHENEHP')
  
  #PHENOL-1Hz_CIT_WENNBERG#C6OH6,  +-(25% of measurement value + 20 pptv) pptv
  if (res == 1){
    pass$tmperror = (pass$PHENOL.1Hz_CIT_WENNBERG*0.25 + 20)/1E3
  } else{
    pass$tmperror = (pass$PHENOL.10Hz_CIT_WENNBERG*0.25 + 20)/1E3
  }
  passErrors=c(passErrors,pass$tmperror)
  if (res == 1){
    pass$tmpy = pass$PHENOL.1Hz_CIT_WENNBERG/1E3
    background$tmpy = background$PHENOL.1Hz_CIT_WENNBERG/1E3
  } else{
    pass$tmpy = pass$PHENOL.10Hz_CIT_WENNBERG/1E3
    background$tmpy = background$PHENOL.10Hz_CIT_WENNBERG/1E3
  }
  pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
  pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
  passERs = c(passERs, pW1$slope)
  passERsint = c(passERsint, pW1$ERint)    
  passERsintfill = c(passERsintfill, pW1$ERintfill)    
  passBGX  = c(passBGX, pW1$BGX)
  passBGCO  = c(passBGCO, pW1CO$BGX)
  
  passBGS  = c(passBGS, pW1$BGS)
  passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
  passExcessX = c(passExcessX,pW1$ExcessX)
  passExcessS = c(passExcessS,pW1$ExcessS)
  passsigma = c(passsigma, pW1$uncertainty_slope)
  passERszero= c(passERszero,pW1$intercept)
  passERsCO = c(passERsCO, pW1CO$slope)
  passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
  passERsCOint = c(passERsCOint, pW1CO$ERint)    
  passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
  passERsCOzero = c(passERsCOzero,pW1CO$intercept)
  passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
  passOHrate = c(passOHrate,4.7E-13*exp(1220/TEMP))
  names1=c(names1,'PHENOL_WENNBERG')
  isNMVOC = c(isNMVOC,1)
  kind = c(kind,'oVOC')
  formulas=c(formulas,'C6H6O')
  usenames=c(usenames,'Phenol')
  mWs   = c(mWs,mWPHENOL)
  nCs  = c(nCs, 6)
  print('PHENOL')
  
  #H2O2-1Hz_CIT_WENNBERG# +-(25% of measurement value + 67 pptv),  pptv
  if (res == 1){
    pass$tmperror = (pass$H2O2.1Hz_CIT_WENNBERG*0.25 + 67)/1E3
  } else{
    pass$tmperror = (pass$H2O2.10Hz_CIT_WENNBERG*0.25 + 67)/1E3
  }
  passErrors=c(passErrors,pass$tmperror)
  if (res == 1){
    pass$tmpy = pass$H2O2.1Hz_CIT_WENNBERG/1E3
    background$tmpy = background$H2O2.1Hz_CIT_WENNBERG/1E3
  } else{
    pass$tmpy = pass$H2O2.10Hz_CIT_WENNBERG/1E3
    background$tmpy = background$H2O2.10Hz_CIT_WENNBERG/1E3
  } 
  pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
  pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
  passERs = c(passERs, pW1$slope)
  passERsint = c(passERsint, pW1$ERint)    
  passERsintfill = c(passERsintfill, pW1$ERintfill)    
  passBGX  = c(passBGX, pW1$BGX)
  passBGCO  = c(passBGCO, pW1CO$BGX)
  
  passBGS  = c(passBGS, pW1$BGS)
  passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
  passExcessX = c(passExcessX,pW1$ExcessX)
  passExcessS = c(passExcessS,pW1$ExcessS)
  passsigma = c(passsigma, pW1$uncertainty_slope)
  passERszero= c(passERszero,pW1$intercept)
  passERsCO = c(passERsCO, pW1CO$slope)
  passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
  passERsCOint = c(passERsCOint, pW1CO$ERint)    
  passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
  passERsCOzero = c(passERsCOzero,pW1CO$intercept)
  passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
  passOHrate = c(passOHrate, NaN)
  names1=c(names1,'H2O2_WENNBERG')
  kind = c(kind,'H2O2')
  formulas=c(formulas,'H2O2')
  usenames=c(usenames,'Hydrogen Peroxide')
  
  isNMVOC = c(isNMVOC,5)
  mWs   = c(mWs,mWH2O2)
  nCs  = c(nCs, 0)
  print('H2O2')
  
  #ISOPN-1Hz_CIT_WENNBERG#C5O4H9N, +-(25% of measurement value + 3 pptv),  pptv
  if (res == 1){
    pass$tmperror = (pass$ISOPN.1Hz_CIT_WENNBERG*0.25 + 3)/1E3
  } else{
    pass$tmperror = (pass$ISOPN.10Hz_CIT_WENNBERG*0.25 + 3)/1E3
  }
  passErrors=c(passErrors,pass$tmperror)
  if (res == 1){
    pass$tmpy = pass$ISOPN.1Hz_CIT_WENNBERG/1E3
    background$tmpy = background$ISOPN.1Hz_CIT_WENNBERG/1E3
  } else{
    pass$tmpy = pass$ISOPN.10Hz_CIT_WENNBERG/1E3
    background$tmpy = background$ISOPN.10Hz_CIT_WENNBERG/1E3
  }
  pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
  pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
  passERs = c(passERs, pW1$slope)
  passERsint = c(passERsint, pW1$ERint)    
  passERsintfill = c(passERsintfill, pW1$ERintfill)    
  passBGX  = c(passBGX, pW1$BGX)
  passBGCO  = c(passBGCO, pW1CO$BGX)
  passBGS  = c(passBGS, pW1$BGS)
  passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
  passExcessX = c(passExcessX,pW1$ExcessX)
  passExcessS = c(passExcessS,pW1$ExcessS)
  passsigma = c(passsigma, pW1$uncertainty_slope)
  passERszero= c(passERszero,pW1$intercept)
  passERsCO = c(passERsCO, pW1CO$slope)
  passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
  passERsCOint = c(passERsCOint, pW1CO$ERint)    
  passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
  passERsCOzero = c(passERsCOzero,pW1CO$intercept)
  passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
  passOHrate = c(passOHrate, 1.12E-10)   # ISOPANO3 = 1.12D-10, ISOPBNO3 = 2.17D-11, ISOPCNO3 = 1.12D-10
  isNMVOC = c(isNMVOC,2)
  kind = c(kind,'alkyl nitrate')
  formulas=c(formulas,'C5O4H9N')
  usenames=c(usenames,'Isoprene Hydroxynitrates ')
  
  names1=c(names1,'ISOPN_WENNBERG')
  mWs   = c(mWs, mWISOPN)
  nCs  = c(nCs,5)
  print('ISOPN')
  
  
  #PROPENE-HN-1Hz_CIT_WENNBERG#C3O4H7N: +-(25% of measurement value + 2.5 pptv)
  if (res == 1){
    pass$tmperror = (pass$PROPENE.HN.1Hz_CIT_WENNBERG*0.25 + 2.5)/1E3
  } else{
    pass$tmperror = (pass$PROPENE.HN.10Hz_CIT_WENNBERG*0.25 + 2.5)/1E3
  }
  passErrors=c(passErrors,pass$tmperror)
  if (res == 1){
    pass$tmpy = pass$PROPENE.HN.1Hz_CIT_WENNBERG/1E3
    background$tmpy = background$PROPENE.HN.1Hz_CIT_WENNBERG/1E3
  } else{
    pass$tmpy = pass$PROPENE.HN.10Hz_CIT_WENNBERG/1E3
    background$tmpy = background$PROPENE.HN.10Hz_CIT_WENNBERG/1E3
  }
  pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
  pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
  passERs = c(passERs, pW1$slope)
  passERsint = c(passERsint, pW1$ERint)    
  passERsintfill = c(passERsintfill, pW1$ERintfill)    
  passBGX  = c(passBGX, pW1$BGX)
  passBGCO  = c(passBGCO, pW1CO$BGX)
  passBGS  = c(passBGS, pW1$BGS)
  passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
  passExcessX = c(passExcessX,pW1$ExcessX)
  passExcessS = c(passExcessS,pW1$ExcessS)
  passsigma = c(passsigma, pW1$uncertainty_slope)
  passERszero= c(passERszero,pW1$intercept)
  passERsCO = c(passERsCO, pW1CO$slope)
  passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
  passERsCOint = c(passERsCOint, pW1CO$ERint)    
  passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
  passERsCOzero = c(passERsCOzero,pW1CO$intercept)
  passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
  passOHrate = c(passOHrate, NaN)
  names1=c(names1,'PROPENE.HN_WENNBERG')
  isNMVOC = c(isNMVOC,2)
  kind = c(kind,'nitrate')
  
  formulas=c(formulas,'C3O4H7N')
  usenames=c(usenames,'Propene Hydroxynitrate')
  mWs   = c(mWs,  mWPROPENEHN)
  nCs  = c(nCs, 3)
  print('PROPENEHN')
  
  #PROPENE-HP-1Hz_CIT_WENNBERG#C3H8O3 +-(30% of measurement value + 15 pptv) ,  pptv
  if (res == 1){
    pass$tmperror = (pass$PROPENE.HP.1Hz_CIT_WENNBERG*0.3 + 15)/1E3
  } else{
    pass$tmperror = (pass$PROPENE.HP.10Hz_CIT_WENNBERG*0.3 + 15)/1E3
  }
  passErrors=c(passErrors,pass$tmperror)
  if (res == 1){
    pass$tmpy = pass$PROPENE.HP.1Hz_CIT_WENNBERG/1E3
    background$tmpy = background$PROPENE.HP.1Hz_CIT_WENNBERG/1E3
  } else{
    pass$tmpy = pass$PROPENE.HP.10Hz_CIT_WENNBERG/1E3
    background$tmpy = background$PROPENE.HP.10Hz_CIT_WENNBERG/1E3
  }
  pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
  pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
  passERs = c(passERs, pW1$slope)
  passERsint = c(passERsint, pW1$ERint)    
  passERsintfill = c(passERsintfill, pW1$ERintfill)    
  passBGX  = c(passBGX, pW1$BGX)
  passBGCO  = c(passBGCO, pW1CO$BGX)
  passBGS  = c(passBGS, pW1$BGS)
  passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
  passExcessX = c(passExcessX,pW1$ExcessX)
  passExcessS = c(passExcessS,pW1$ExcessS)
  passsigma = c(passsigma, pW1$uncertainty_slope)
  passERszero= c(passERszero,pW1$intercept)
  passERsCO = c(passERsCO, pW1CO$slope)
  passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
  passERsCOint = c(passERsCOint, pW1CO$ERint)    
  passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
  passERsCOzero = c(passERsCOzero,pW1CO$intercept)
  passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
  passOHrate = c(passOHrate,OHpropeneHP)

  names1=c(names1,'PROPENE.HP_WENNBERG')
  formulas=c(formulas,'C3O3H8')
  usenames=c(usenames,'Propene Hydroxyperoxide')
  isNMVOC = c(isNMVOC,1)
  kind = c(kind,'oVOC')
  
  mWs   = c(mWs,mWPROPENEHP)
  nCs  = c(nCs, 3)
  print('PROPENEHP')
  
  #HCN-1Hz_CIT_WENNBERG#+-(25% of measurement value + 61 pptv),  pptv
  if (res == 1){
    pass$tmperror = (pass$HCN.1Hz_CIT_WENNBERG*0.25 + 61)/1E3
  }else{
    pass$tmperror = (pass$HCN.10Hz_CIT_WENNBERG*0.25 + 61)/1E3
  }
  passErrors=c(passErrors,pass$tmperror)
  if (res == 1){
    pass$tmpy = pass$HCN.1Hz_CIT_WENNBERG/1E3
    background$tmpy = background$HCN.1Hz_CIT_WENNBERG/1E3
  } else{
    pass$tmpy = pass$HCN.10Hz_CIT_WENNBERG/1E3
    background$tmpy = background$HCN.10Hz_CIT_WENNBERG/1E3
  }
  pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
  pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
  passERs = c(passERs, pW1$slope)
  passERsint = c(passERsint, pW1$ERint)    
  passERsintfill = c(passERsintfill, pW1$ERintfill)    
  passBGX  = c(passBGX, pW1$BGX)
  passBGCO  = c(passBGCO, pW1CO$BGX)
  passBGS  = c(passBGS, pW1$BGS)
  passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
  passExcessX = c(passExcessX,pW1$ExcessX)
  passExcessS = c(passExcessS,pW1$ExcessS)
  passsigma = c(passsigma, pW1$uncertainty_slope)
  passERszero= c(passERszero,pW1$intercept)
  passERsCO = c(passERsCO, pW1CO$slope)
  passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
  passERsCOint = c(passERsCOint, pW1CO$ERint)    
  passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
  passERsCOzero = c(passERsCOzero,pW1CO$intercept)
  passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
  passOHrate = c(passOHrate, NaN)
  names1=c(names1, 'HCN_WENNBERG')
  formulas=c(formulas,'HCN')
  usenames=c(usenames,'Hydrogen cyanide')
  isNMVOC = c(isNMVOC,2)
  kind = c(kind,'nitrogen')
  
  mWs   = c(mWs, mWHCN)
  nCs  = c(nCs,1)
  print('HCN - Wennberg')
  
  # #HNO3-1Hz_CIT_WENNBERG#+-(25% of measurement value + 48 pptv),  pptv
  # if (res ==1){ pass$HNO3.1Hz_CIT_WENNBERG = pass$HNO3_CIT_WENNBERG ; background$HNO3.1Hz_CIT_WENNBERG = background$HNO3_CIT_WENNBERG}
  if (res == 1){
    pass$tmperror = (pass$HNO3.1Hz_CIT_WENNBERG*0.25 + 48)/1E3
  } else{
    pass$tmperror = (pass$HNO3.10Hz_CIT_WENNBERG*0.25 + 48)/1E3
  }
  passErrors=c(passErrors,pass$tmperror)
  if (res == 1){
    pass$tmpy = pass$HNO3.1Hz_CIT_WENNBERG/1E3
    background$tmpy = background$HNO3.1Hz_CIT_WENNBERG/1E3
  } else{
    pass$tmpy = pass$HNO3.10Hz_CIT_WENNBERG/1E3
    background$tmpy = background$HNO3.10Hz_CIT_WENNBERG/1E3
  }
  pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
  pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
  passERs = c(passERs, pW1$slope)
  passERsint = c(passERsint, pW1$ERint)    
  passERsintfill = c(passERsintfill, pW1$ERintfill)    
  passBGX  = c(passBGX, pW1$BGX)
  passBGCO  = c(passBGCO, pW1CO$BGX)
  passBGS  = c(passBGS, pW1$BGS)
  passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
  passExcessX = c(passExcessX,pW1$ExcessX)
  passExcessS = c(passExcessS,pW1$ExcessS)
  passsigma = c(passsigma, pW1$uncertainty_slope)
  passERszero= c(passERszero,pW1$intercept)
  passERsCO = c(passERsCO, pW1CO$slope)
  passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
  passERsCOint = c(passERsCOint, pW1CO$ERint)    
  passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
  passERsCOzero = c(passERsCOzero,pW1CO$intercept)
  passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
  passOHrate = c(passOHrate, NaN)
  names1=c(names1, 'HNO3_WENNBERG')
  isNMVOC = c(isNMVOC,5)
  kind = c(kind,'nitrogen')
  formulas=c(formulas,'HNO3')
  usenames=c(usenames,'Nitric acid')
  mWs   = c(mWs, mWHNO3)
  nCs  = c(nCs,0)
  print('HNO3 - Wennberg')
  
  # ------------------------------------------ HUEY ------------------------------------------ 
  # PAN
  if (debugKT == 1){ print("HUEY")}
  pass$tmperror = (pass$PAN_GTCIMS_HUEY*0.2)/1E3
  passErrors=c(passErrors,pass$tmperror)
  pass$tmpy = pass$PAN_GTCIMS_HUEY/1E3
  background$tmpy = background$PAN_GTCIMS_HUEY/1E3
  pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
  pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
  passERs = c(passERs, pW1$slope)
  passERsint = c(passERsint, pW1$ERint)    
  passERsintfill = c(passERsintfill, pW1$ERintfill)    
  passBGX  = c(passBGX, pW1$BGX)
  passBGCO  = c(passBGCO, pW1CO$BGX)
  passBGS  = c(passBGS, pW1$BGS)
  passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
  passExcessX = c(passExcessX,pW1$ExcessX)
  passExcessS = c(passExcessS,pW1$ExcessS)
  passsigma = c(passsigma, pW1$uncertainty_slope)
  passERszero= c(passERszero,pW1$intercept)
  passERsCO = c(passERsCO, pW1CO$slope)
  passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
  passERsCOint = c(passERsCOint, pW1CO$ERint)    
  passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
  passERsCOzero = c(passERsCOzero,pW1CO$intercept)
  passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
  passOHrate = c(passOHrate, NaN)
  
  if (debugKT == 1){ print("APAN")}
  # APAN_GTCIMS_HUEY,  pptv - no uncertainty given, use 30%
  pass$tmperror = (pass$APAN_GTCIMS_HUEY*0.3)/1E3
  passErrors=c(passErrors,pass$tmperror)
  pass$tmpy = pass$APAN_GTCIMS_HUEY/1E3
  background$tmpy = background$APAN_GTCIMS_HUEY/1E3
  pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
  pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
  passERs = c(passERs, pW1$slope)
  passERsint = c(passERsint, pW1$ERint)    
  passERsintfill = c(passERsintfill, pW1$ERintfill)    
  passBGX  = c(passBGX, pW1$BGX)
  passBGCO  = c(passBGCO, pW1CO$BGX)
  passBGS  = c(passBGS, pW1$BGS)
  passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
  passExcessX = c(passExcessX,pW1$ExcessX)
  passExcessS = c(passExcessS,pW1$ExcessS)
  passsigma = c(passsigma, pW1$uncertainty_slope)
  passERszero= c(passERszero,pW1$intercept)
  passERsCO = c(passERsCO, pW1CO$slope)
  passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
  passERsCOint = c(passERsCOint, pW1CO$ERint)    
  passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
  passERsCOzero = c(passERsCOzero,pW1CO$intercept)
  passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
  passOHrate = c(passOHrate, NaN)
  names1 = c(names1,'PAN_HUEY','APAN_HUEY')
  
  formulas=c(formulas,'PAN','APAN')
  usenames=c(usenames,'Peroxyacetyl nitrate','Peroxyacryloyl nitrate')
  
  isNMVOC = c(isNMVOC,5, 5)
  kind = c(kind,'NOy','NOy')
  
  mWs    = c(mWs,mWPAN, mWAPAN)
  nCs     = c(nCs, 2, 3)
  
  if (res != 5){ # not present in 5 Hz
    
    # PPN GTCIMS_HUEY,  pptv - no uncertainty given, use 30%
    pass$tmperror = (pass$PPN_GTCIMS_HUEY*0.3)/1E3
    passErrors=c(passErrors,pass$tmperror)
    pass$tmpy = pass$PPN_GTCIMS_HUEY/1E3
    background$tmpy = background$PPN_GTCIMS_HUEY/1E3
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)    
    passERsintfill = c(passERsintfill, pW1$ERintfill)    
    passBGX  = c(passBGX, pW1$BGX)
    passBGCO  = c(passBGCO, pW1CO$BGX)
    passBGS  = c(passBGS, pW1$BGS)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passExcessX = c(passExcessX,pW1$ExcessX)
    passExcessS = c(passExcessS,pW1$ExcessS)
    passsigma = c(passsigma, pW1$uncertainty_slope)
    passERszero= c(passERszero,pW1$intercept)
    passERsCO = c(passERsCO, pW1CO$slope)
    passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
    passERsCOint = c(passERsCOint, pW1CO$ERint)    
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
    passERsCOzero = c(passERsCOzero,pW1CO$intercept)
    passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
    passOHrate = c(passOHrate, NaN)
    
    # PBN_GTCIMS_HUEY,  pptv - no uncertainty given, use 30%
    pass$tmperror = (pass$PBN_GTCIMS_HUEY*0.3)/1E3
    passErrors=c(passErrors,pass$tmperror)
    pass$tmpy = pass$PBN_GTCIMS_HUEY/1E3
    background$tmpy = background$PBN_GTCIMS_HUEY/1E3
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)    
    passERsintfill = c(passERsintfill, pW1$ERintfill)    
    passBGX  = c(passBGX, pW1$BGX)
    passBGCO  = c(passBGCO, pW1CO$BGX)
    passBGS  = c(passBGS, pW1$BGS)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passExcessX = c(passExcessX,pW1$ExcessX)
    passExcessS = c(passExcessS,pW1$ExcessS)
    passsigma = c(passsigma, pW1$uncertainty_slope)
    passERszero= c(passERszero,pW1$intercept)
    passERsCO = c(passERsCO, pW1CO$slope)
    passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
    passERsCOint = c(passERsCOint, pW1CO$ERint)    
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
    passERsCOzero = c(passERsCOzero,pW1CO$intercept)
    passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
    passOHrate = c(passOHrate, NaN)
    names1 = c(names1,'PPN_HUEY','PBN_HUEY')
    
    formulas=c(formulas,'PPN','PBN')
    usenames=c(usenames,'Peroxyacetyl nitrate','Peroxybenzoyl nitrate')
    isNMVOC = c(isNMVOC,5,5)
    kind = c(kind,'NOy','NOy')
    
    mWs    = c(mWs, mWPPN, mWPBN)
    nCs     = c(nCs,  3, 7)
  }
  
  # ------------------------------------------ MOORE------------------------------------------ 
  if (res ==5){
    if (debugKT == 1){print("MOORE")}
    print('Fast CNgt6')
    #20% - stdPT,cm-3. Number concentration of aerosol particles with nominal diameters greater than 6nm
    pass$tmperror = (pass$FastCNgt6nm_stdPT_MOORE*0.2)
    passErrors=c(passErrors,pass$tmperror)
    pass$tmpy = pass$FastCNgt6nm_stdPT_MOORE
    background$tmpy = background$FastCNgt6nm_stdPT_MOORE
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)    
    passERsintfill = c(passERsintfill, pW1$ERintfill)    
    passBGX  = c(passBGX, pW1$BGX)
    passBGCO  = c(passBGCO, pW1CO$BGX)
    passBGS  = c(passBGS, pW1$BGS)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passExcessX = c(passExcessX,pW1$ExcessX)
    passExcessS = c(passExcessS,pW1$ExcessS)
    passsigma = c(passsigma, pW1$uncertainty_slope)
    passERszero= c(passERszero,pW1$intercept)
    passERsCO = c(passERsCO, pW1CO$slope)
    passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
    passERsCOint = c(passERsCOint, pW1CO$ERint)    
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
    passERsCOzero = c(passERsCOzero,pW1CO$intercept)
    passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
    passOHrate = c(passOHrate, NaN)
    names1=c(names1,'CNgt6nm')
    kind = c(kind,'aerosol')
    formulas=c(formulas,'N/A')
    usenames=c(usenames,'CN > 6nm')
    
    isNMVOC = c(isNMVOC,4)
    mWs   = c(mWs,NA)
    nCs  = c(nCs, NA)
  } else{
    print('CNgt6')
    #20% - stdPT,cm-3. Number concentration of aerosol particles with nominal diameters greater than 6nm
    pass$tmperror = (pass$CNgt6nm_stdPT*0.2)
    passErrors=c(passErrors,pass$tmperror)
    pass$tmpy = pass$CNgt6nm_stdPT
    background$tmpy = background$CNgt6nm_stdPT
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)    
    passERsintfill = c(passERsintfill, pW1$ERintfill)    
    passBGX  = c(passBGX, pW1$BGX)
    passBGCO  = c(passBGCO, pW1CO$BGX)
    passBGS  = c(passBGS, pW1$BGS)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passExcessX = c(passExcessX,pW1$ExcessX)
    passExcessS = c(passExcessS,pW1$ExcessS)
    passsigma = c(passsigma, pW1$uncertainty_slope)
    passERszero= c(passERszero,pW1$intercept)
    passERsCO = c(passERsCO, pW1CO$slope)
    passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
    passERsCOint = c(passERsCOint, pW1CO$ERint)    
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
    passERsCOzero = c(passERsCOzero,pW1CO$intercept)
    passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
    passOHrate = c(passOHrate, NaN)
    names1=c(names1,'CNgt6nm')
    kind = c(kind,'aerosol')
    formulas=c(formulas,'N/A')
    usenames=c(usenames,'CN > 6nm')
    
    isNMVOC = c(isNMVOC,4)
    mWs   = c(mWs,NA)
    nCs  = c(nCs, NA)
    
    print('CNgt3')
    #20% - stdPT,cm-3. Number concentration of aerosol particles with nominal diameters greater than 3nm
    pass$tmperror = (pass$CNgt3nm_stdPT*0.2)
    passErrors=c(passErrors,pass$tmperror)
    pass$tmpy = pass$CNgt3nm_stdPT
    background$tmpy = background$CNgt3nm_stdPT
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)    
    passERsintfill = c(passERsintfill, pW1$ERintfill)    
    passBGX  = c(passBGX, pW1$BGX)
    passBGCO  = c(passBGCO, pW1CO$BGX)
    passBGS  = c(passBGS, pW1$BGS)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passExcessX = c(passExcessX,pW1$ExcessX)
    passExcessS = c(passExcessS,pW1$ExcessS)
    passsigma = c(passsigma, pW1$uncertainty_slope)
    passERszero= c(passERszero,pW1$intercept)
    passERsCO = c(passERsCO, pW1CO$slope)
    passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
    passERsCOint = c(passERsCOint, pW1CO$ERint)    
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
    passERsCOzero = c(passERsCOzero,pW1CO$intercept)
    passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
    passOHrate = c(passOHrate, NaN)
    names1=c(names1,'CNgt3nm')
    kind = c(kind,'aerosol')
    formulas=c(formulas,'N/A')
    usenames=c(usenames,'CN > 3nm')
    
    isNMVOC = c(isNMVOC,4)
    mWs   = c(mWs,NA)
    nCs  = c(nCs, NA)
    
    print('CNgt20')
    #20% - stdPT,cm-3. Number concentration of aerosol particles with nominal diameters greater than 3nm
    pass$tmperror = (pass$CNgt20nm_stdPT*0.2)
    passErrors=c(passErrors,pass$tmperror)
    pass$tmpy = pass$CNgt20nm_stdPT
    background$tmpy = background$CNgt20nm_stdPT
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)    
    passERsintfill = c(passERsintfill, pW1$ERintfill)    
    passBGX  = c(passBGX, pW1$BGX)
    passBGCO  = c(passBGCO, pW1CO$BGX)
    passBGS  = c(passBGS, pW1$BGS)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passExcessX = c(passExcessX,pW1$ExcessX)
    passExcessS = c(passExcessS,pW1$ExcessS)
    passsigma = c(passsigma, pW1$uncertainty_slope)
    passERszero= c(passERszero,pW1$intercept)
    passERsCO = c(passERsCO, pW1CO$slope)
    passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
    passERsCOint = c(passERsCOint, pW1CO$ERint)    
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
    passERsCOzero = c(passERsCOzero,pW1CO$intercept)
    passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
    passOHrate = c(passOHrate, NaN)
    names1=c(names1,'CNgt20nm')
    kind = c(kind,'aerosol')
    formulas=c(formulas,'N/A')
    usenames=c(usenames,'CN > 20nm')
    
    isNMVOC = c(isNMVOC,4)
    mWs   = c(mWs,NA)
    nCs  = c(nCs, NA)
    
    print('CCN0.34')
    #20% - stdPT,cm-3. Number concentration of aerosol particles with nominal diameters greater than 3nm
    pass$tmperror = (pass$CCN_034_stdPT*0.2)
    passErrors=c(passErrors,pass$tmperror)
    pass$tmpy = pass$CCN_034_stdPT
    background$tmpy = background$CCN_034_stdPT
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)    
    passERsintfill = c(passERsintfill, pW1$ERintfill)    
    passBGX  = c(passBGX, pW1$BGX)
    passBGCO  = c(passBGCO, pW1CO$BGX)
    passBGS  = c(passBGS, pW1$BGS)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passExcessX = c(passExcessX,pW1$ExcessX)
    passExcessS = c(passExcessS,pW1$ExcessS)
    passsigma = c(passsigma, pW1$uncertainty_slope)
    passERszero= c(passERszero,pW1$intercept)
    passERsCO = c(passERsCO, pW1CO$slope)
    passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
    passERsCOint = c(passERsCOint, pW1CO$ERint)    
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
    passERsCOzero = c(passERsCOzero,pW1CO$intercept)
    passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
    passOHrate = c(passOHrate, NaN)
    names1=c(names1,'CCN_034_stdPT')
    kind = c(kind,'aerosol')
    formulas=c(formulas,'N/A')
    usenames=c(usenames,'CCN_034')
    
    isNMVOC = c(isNMVOC,4)
    mWs   = c(mWs,NA)
    nCs  = c(nCs, NA)
    
    print('Ngt100nm_LAS_stdPT')
    pass$tmperror = (pass$Ngt100nm_LAS_stdPT*0.2)
    passErrors=c(passErrors,pass$tmperror)
    pass$tmpy = pass$Ngt100nm_LAS_stdPT
    background$tmpy = background$Ngt100nm_LAS_stdPT
    pW1 = smallfun(pass,background,xspecies,'tmpy',      xerror, 'tmperror',SLOW) 
    pW1CO = smallfun(pass,background,'CO_DACOM_DISKIN','tmpy',      xerror, 'tmperror',SLOW) 
    passERs = c(passERs, pW1$slope)
    passERsint = c(passERsint, pW1$ERint)    
    passERsintfill = c(passERsintfill, pW1$ERintfill)    
    passBGX  = c(passBGX, pW1$BGX)
    passBGCO  = c(passBGCO, pW1CO$BGX)
    passBGS  = c(passBGS, pW1$BGS)
    passRsqsCO2 = c(passRsqsCO2, pW1$r_sq)   
    passExcessX = c(passExcessX,pW1$ExcessX)
    passExcessS = c(passExcessS,pW1$ExcessS)
    passsigma = c(passsigma, pW1$uncertainty_slope)
    passERszero= c(passERszero,pW1$intercept)
    passERsCO = c(passERsCO, pW1CO$slope)
    passRsqsCO  = c(passRsqsCO, pW1CO$r_sq)   
    passERsCOint = c(passERsCOint, pW1CO$ERint)    
    passERsCOintfill = c(passERsCOintfill, pW1CO$ERintfill)  
    passERsCOzero = c(passERsCOzero,pW1CO$intercept)
    passMAX   = c(passMAX, maxfun(pass,'tmpy' ))
    passOHrate = c(passOHrate, NaN)
    names1=c(names1,'Ngt100nm_LAS_stdPT')
    kind = c(kind,'aerosol')
    formulas=c(formulas,'N/A')
    usenames=c(usenames,'Ngt100nm_LAS_stdPT')
    
    isNMVOC = c(isNMVOC,4)
    mWs   = c(mWs,NA)
    nCs  = c(nCs, NA)
  }  
  # --------- PHOTOCHEMICAL CLOCK(s) -----------------
  # maleic anhydryde to furan
  pass$tmperrorx = pass$Furan_NOAAPTR_ppbv_WARNEKE * 0.3
  pass$tmperrory = pass$MaleicAnhyd_NOAAPTR_ppbv_WARNEKE* 0.5
  passErrors=c(passErrors,pass$tmperrory)
  
  pA = smallfun(pass,background,'Furan_NOAAPTR_ppbv_WARNEKE','MaleicAnhyd_NOAAPTR_ppbv_WARNEKE', 'tmperrorx', 'tmperrory', SLOW) 
  matoF = pA$slope
  matoF_rq = pA$r_sq
  
  passERs = c(passERs, NaN)    
  passERsCO = c(passERsCO,  NaN)
  passERsint = c(passERsint, NaN)    
  passERsintfill = c(passERsintfill,NaN)    
  passERsCOint = c(passERsCOint, NaN)     
  passERsCOintfill = c(passERsCOintfill, NaN)   
  passERszero = c(passERszero, NaN)     
  passERsCOzero = c(passERsCOzero, NaN)   
  
  passBGX  = c(passBGX,  NaN)
  passBGCO = c(passBGCO,NaN)
  passBGS  = c(passBGS,  NaN)
  passRsqsCO2 = c(passRsqsCO2,NaN)
  passRsqsCO = c(passRsqsCO,NaN) 
  passExcessX = c(passExcessX,NaN)
  passExcessS = c(passExcessS,NaN)
  passsigma = c(passsigma, pA$uncertainty_slope)
  passMAX   = c(passMAX, NaN)
  passOHrate = c(passOHrate, NaN)
  names1=c(names1,'MAtoF')
  formulas=c(formulas,'MAtoF')
  usenames=c(usenames,'MAtoF')
  isNMVOC = c(isNMVOC,NaN)
  kind = c(kind,NaN)
  mWs = c(mWs, NaN) # 
  nCs = c(nCs, NaN) #
  
  if (res != 5){ # not present in 5 Hz
    #  NOx to NOy
    pass$NOx_RYERSON=  (pass$NO_CL_RYERSON + pass$NO2_CL_RYERSON)
    pass$tmperrorx = pass$NOy_CL_RYERSON* 0.13
    pass$tmperrory = pass$NOx_RYERSON * 0.08
    passErrors=c(passErrors,pass$tmperrorx)
    pA = smallfun(pass,background, 'NOy_CL_RYERSON','NOx_RYERSON','tmperrorx', 'tmperrory', SLOW) 
    passERs = c(passERs, pA$slope)    
    passERsCO = c(passERsCO, NaN) 
    passERszero = c(passERszero, NaN)     
    passERsCOzero = c(passERsCOzero, NaN)   
    passERsint = c(passERsint, pA$ERint)   
    passERsintfill = c(passERsintfill, pA$ERintfill)    
    passERsCOint = c(passERsCOint, NaN)   
    passERsCOintfill = c(passERsCOintfill, NaN)
    passBGX  = c(passBGX, pA$BGX)
    passBGCO = c(passBGCO,NaN)
    
    passBGS  = c(passBGS, pA$BGS)
    passRsqsCO2 = c(passRsqsCO2,NaN)  
    passRsqsCO = c(passRsqsCO, NaN)
    
    passExcessX = c(passExcessX,NaN)
    passExcessS = c(passExcessS,NaN)
    passsigma = c(passsigma, pA$uncertainty_slope)
    passMAX = c(passMAX, NaN)
    passOHrate = c(passOHrate, NaN)
    names1=c(names1,'NOx/NOy')
    formulas=c(formulas,'NOx/NOy')
    usenames=c(usenames,'NOx/NOy')
    isNMVOC = c(isNMVOC,NaN)
    kind = c(kind,NaN)
    mWs = c(mWs, NaN) # 
    nCs = c(nCs, NaN) #
    
    #  PAN to NOy
    pass$tmperrory = (pass$PAN_GTCIMS_HUEY*0.2)/1E3
    passErrors=c(passErrors,pass$tmperror)
    pass$tmpy = pass$PAN_GTCIMS_HUEY/1E3
    background$tmpy = background$PAN_GTCIMS_HUEY/1E3
    pA = smallfun(pass,background, 'NOy_CL_RYERSON','tmpy','tmperrorx', 'tmperrory', SLOW) 
    passERs = c(passERs, pA$slope)    
    passERsCO = c(passERsCO, NaN) 
    passERszero = c(passERszero, NaN)     
    passERsCOzero = c(passERsCOzero, NaN)   
    passERsint = c(passERsint, pA$ERint)   
    passERsintfill = c(passERsintfill, pA$ERintfill)    
    passERsCOint = c(passERsCOint, NaN)   
    passERsCOintfill = c(passERsCOintfill, NaN)
    passBGX  = c(passBGX, pA$BGX)
    passBGCO = c(passBGCO,NaN)
    
    passBGS  = c(passBGS, pA$BGS)
    passRsqsCO2 = c(passRsqsCO2,NaN)  
    passRsqsCO = c(passRsqsCO, NaN)
    
    passExcessX = c(passExcessX,NaN)
    passExcessS = c(passExcessS,NaN)
    passsigma = c(passsigma, pA$uncertainty_slope)
    passMAX = c(passMAX, NaN)
    passOHrate = c(passOHrate, NaN)
    names1=c(names1,'PAN/NOy')
    
    formulas=c(formulas,'NOx/NOy')
    usenames=c(usenames,'NOx/NOy')
    isNMVOC = c(isNMVOC,NaN)
    kind = c(kind,NaN)
    mWs = c(mWs, NaN) # 
    nCs = c(nCs, NaN) #
    
    #  HNO3 to NOy
    # #HNO3-1Hz_CIT_WENNBERG#+-(25% of measurement value + 48 pptv),  pptv
    # if (res ==1){ pass$HNO3.1Hz_CIT_WENNBERG = pass$HNO3_CIT_WENNBERG ; background$HNO3.1Hz_CIT_WENNBERG = background$HNO3_CIT_WENNBERG}
    # pass$tmperror = (pass$HNO3.1Hz_CIT_WENNBERG*0.25 + 48)/1E3
    # pass$tmpy = pass$HNO3.1Hz_CIT_WENNBERG/1E3
    # background$tmpy = background$HNO3.1Hz_CIT_WENNBERG/1E3
    # 
    # pA = smallfun(pass,background, 'NOy_CL_RYERSON','tmpy','tmperrorx', 'tmperrory', SLOW) 
    # passERs = c(passERs, pA$slope)    
    # passERsCO = c(passERsCO, NaN) 
    # passERszero = c(passERszero, NaN)     
    # passERsCOzero = c(passERsCOzero, NaN)   
    # passERsint = c(passERsint, pA$ERint)   
    # passERsintfill = c(passERsintfill, pA$ERintfill)    
    # passERsCOint = c(passERsCOint, NaN)   
    # passERsCOintfill = c(passERsCOintfill, NaN)
    # passBGX  = c(passBGX, pA$BGX)
    #passBGCO = c(passBGCO,NaN)
    # passBGS  = c(passBGS, pA$BGS)
    # passRsqsCO2 = c(passRsqsCO2,NaN)  
    # passRsqsCO = c(passRsqsCO, NaN)
    # 
    # passExcessX = c(passExcessX,NaN)
    # passExcessS = c(passExcessS,NaN)
    # passsigma = c(passsigma, pA$uncertainty_slope)
    # passMAX = c(passMAX, NaN)
    # passOHrate = c(passOHrate, NaN)
    # names1=c(names1,'HNO3/NOy')
    # isNMVOC = c(isNMVOC,NaN)
    # mWs = c(mWs, NaN) # 
    # nCs = c(nCs, NaN) #
    # 
    
  }
  # ---- Combine all data ------
  print("data right length?")
  # only want intfill, int is not sueful
  print(c(length(passERs), length(passERsint), length(passERsintfill),length(passERsCO),
                            length(passsigma), length(passRsqsCO2), length(passRsqsCO),length(passMAX),length(passOHrate),
                            length(names1), length(mWs), length(nCs), length(passBGS), length(passBGX), length(passBGCO),
                            length(passExcessX), length(passExcessS), length(passERsCOint), length(passERsCOintfill),
                            length(passERszero), length(passERsCOzero), length(isNMVOC), length(kind),length(passErrors),
          length(usenames),length(formulas)))
  
  datareturn = as.data.frame(cbind(as.numeric(passERs),as.numeric(passERsCO), as.numeric(passERsintfill),as.numeric(passERsCOintfill),
                                   as.numeric(passsigma),
                                   round(passRsqsCO2, digits=2),round(passRsqsCO, digits=2),as.numeric(passMAX), as.numeric(passOHrate),
                                   as.numeric(passBGX), as.numeric(passBGCO),as.numeric(passBGS), as.numeric(passExcessX), 
                                   as.numeric(passExcessS), as.numeric(passERszero), as.numeric(passERsCOzero), as.numeric(isNMVOC), kind))
  colnames(datareturn)=c("ERtoX", "ERtoCO","ERtoXintfill","ERtoCOintfill",
                         "ERsigma","R2toX","R2toCO","maxval", "OHrate","BGX","BGCO","BGSpecies", "ExcessX", "ExcessS",
                         "intercept","COintercept", "Category","kind")
  datareturn$variable=names1
  datareturn$names = usenames
  datareturn$formula = formulas
  datareturn$mWs=mWs
  datareturn$nCs=nCs
  datareturn$Start = min(pass$Time_Start)
  datareturn$Stop  = max(pass$Time_Start)
  datareturn$npts = length(pass$CO_DACOM_DISKIN)
  datareturn$MAtoF = matoF
  datareturn$MAtoF_rq = matoF_rq
  datareturn$OAtoOC = mean(OAtoOC, na.rm=TRUE)
  datareturn$OtoC = mean(OtoC, na.rm=TRUE)
  datareturn$HtoC = mean(HtoC, na.rm=TRUE)
  datareturn$f43 = mean(f43, na.rm=TRUE)
  datareturn$f44 = mean(f44, na.rm=TRUE)
  datareturn$f57 = mean(f57, na.rm=TRUE)
  datareturn$f60 = mean(f60, na.rm=TRUE)
  datareturn$f82 = mean(f82, na.rm=TRUE)
  datareturn$f91 = mean(f91, na.rm=TRUE)
  
  # get lat/lon/alt of maximum smoke for this pass
  #ind = which(datareturn$CO_DACOM_DISKIN == max(datareturn$CO_DACOM_DISKIN, na.rm=TRUE)
  #             & is.finite(datareturn$Latitude_YANG))
  #if (length(ind) > 0 ){
  
  #} else{ datareturn$lon = NaN; datareturn$lat = NaN; datareturn$alt = NaN}
  # ---- ++++++++++++ Emission Factors ++++++++++++++++ ----
  require(imputeTS)
  # ------ Carbon content - for now just 50%
  datareturn$fire = fire
  datareturn$fuel = fuel
  cContent = 500 # Assuming 50 % C, 500 g/kg
  if (fuel == 'corn' | fuel == 'soybean' | fuel == 'winter wheat' | fuel == 'rice'){
   cContent = 415.1
  }
  if (fuel == 'grass' | fuel == 'shrub'){
   cContent = 462.7
  }
  if (fuel == 'pile' | fuel == 'slash') {
   cContent = 511.1
  }
  #00000 slope method 0000000
  datareturn$R2toCO  = as.numeric(datareturn$R2toCO) # R2 to CO
  #  ----- using CO2 as the xspecies
  datareturn$ERtoX = as.numeric(datareturn$ERtoX) # E/R to CO2
  datareturn$ERtoXint = as.numeric(datareturn$ERtoXint) # integrate over CO2
  datareturn$R2toX  = as.numeric(datareturn$R2toX) # R2 to CO2
  datareturn$ERsigma = as.numeric(datareturn$ERsigma) # slope uncertainty
  
  datareturn$TC1 = datareturn$ERtoX[1]*datareturn$nCs[1]+ datareturn$ERtoX[2]*datareturn$nCs[2]+  # just CO, CO2, CH4
    datareturn$ERtoX[3] *datareturn$nCs[3]
  datareturn$TC2 = sum(datareturn$ERtoX*datareturn$nCs, na.rm=TRUE) # all carbon containing species (that I have looked up mWs for at the moment)
  datareturn$PCTCarbon = (datareturn$ERtoX*datareturn$nCs)/datareturn$TC2 # all carbon containing species (that I have looked up mWs for at the moment)
  datareturn$EF1   = cContent * datareturn$mWs/12 * datareturn$ERtoX /datareturn$TC1 # EF slope to CO2
  datareturn$EF2   = cContent * datareturn$mWs/12 * datareturn$ERtoX /datareturn$TC2 #  ---- including all the carbon
  
  #  ----- using CO as the xspecies
 datareturn$ERtoCO = as.numeric(datareturn$ERtoCO)
 datareturn$ERtoXintfill = as.numeric(datareturn$ERtoXintfill)
 datareturn$ERtoCOint = as.numeric(datareturn$ERtoCOint)
 datareturn$ERtoCOintfill = as.numeric(datareturn$ERtoCOintfill)
 
  datareturn$TC1CO = datareturn$ERtoCO[1]*datareturn$nCs[1]+ datareturn$ERtoCO[2]*datareturn$nCs[2]+  # just CO, CO2, CH4
    datareturn$ERtoCO[3] *datareturn$nCs[3]
  datareturn$TC2CO = sum(datareturn$ERtoCO*datareturn$nCs, na.rm=TRUE) # all carbon containing species (that I have looked up mWs for at the moment)
  datareturn$EF1CO   = cContent * datareturn$mWs/12 * datareturn$ERtoCO /datareturn$TC1CO # EF slope to CO
  datareturn$EF2CO   = cContent * datareturn$mWs/12 * datareturn$ERtoCO /datareturn$TC2CO     # ---- including all the carbon 
  # For particle number
  ind = which(datareturn$variable == 'CNgt6nm')
  c1 = datareturn$ERtoCO[ind]/(datareturn$ERtoCO[1]*datareturn$nCs[1]+ datareturn$ERtoCO[2]*datareturn$nCs[2]+  # just CO, CO2, CH4
            datareturn$ERtoCO[3] *datareturn$nCs[3]) # particles/cm3/ppbC
  datareturn$EF1CO[ind] =  (c1*6.022E23/2.6E10/12)*cContent #conversion at 273K and 1atm
  # get error
  datareturn$TC1E = datareturn$ERsigma[1]*datareturn$nCs[1]+ datareturn$ERsigma[2]*datareturn$nCs[2]+  # just CO, CO2, CH4
    datareturn$ERsigma[3] *datareturn$nCs[3]
  datareturn$EF1sigmaMAX   = cContent * datareturn$mWs/12 * (datareturn$ERtoX + datareturn$ERsigma) /
    (datareturn$TC1 + datareturn$TC1E)
  datareturn$EF1sigmaMIN   = cContent * datareturn$mWs/12 * (datareturn$ERtoX - datareturn$ERsigma) /
    (datareturn$TC1 - datareturn$TC1E)
  # append useful information for interpretation
  datareturn$mce   = mce 
  datareturn$mce_int   = mce_int 
  
  # set EF & ER to zero if bad pass
  datareturn$EF1 = datareturn$EF1 * mult 
  datareturn$EF2 = datareturn$EF2 * mult 
  datareturn$ERtoX = datareturn$ERtoX * mult 
  datareturn$R2toCO  = datareturn$R2toCO * mult 
  # 00000 integration method 0000000
  datareturn$TC1int = datareturn$ERtoXint[1]*datareturn$nCs[1]+ datareturn$ERtoXint[2]*datareturn$nCs[2]+  # just CO, CO2, CH4
    datareturn$ERtoXint[3] *datareturn$nCs[3]
  datareturn$TC2int = sum(datareturn$ERtoXint*datareturn$nCs, na.rm=TRUE) # all carbon containing species (that I have looked up mWs for at the moment)
  datareturn$PCTCarbon = (datareturn$ERtoXint*datareturn$nCs)/datareturn$TC2int # all carbon containing species (that I have looked up mWs for at the moment)
  datareturn$EF1int   = cContent * datareturn$mWs/12 * datareturn$ERtoXint /datareturn$TC1int
  datareturn$EF2int   = cContent * datareturn$mWs/12 * datareturn$ERtoXint /datareturn$TC2int
  # 00000 integration method with datafilling0000000
  datareturn$TC1intfill = datareturn$ERtoXintfill[1]*datareturn$nCs[1]+ datareturn$ERtoXintfill[2]*datareturn$nCs[2]+  # just CO, CO2, CH4
    datareturn$ERtoXintfill[3] *datareturn$nCs[3]
  datareturn$TC2intfill = sum(datareturn$ERtoXintfill*datareturn$nCs, na.rm=TRUE) # all carbon containing species (that I have looked up mWs for at the moment)
  datareturn$EF1intfill   = cContent * datareturn$mWs/12 * datareturn$ERtoXintfill /datareturn$TC1intfill
  datareturn$EF2intfill   = cContent * datareturn$mWs/12 * datareturn$ERtoXintfill /datareturn$TC2intfill
  
  # 00000 integration method using CO 0000000
  datareturn$TC1COint = datareturn$ERtoCOint[1]*datareturn$nCs[1]+ datareturn$ERtoCOint[2]*datareturn$nCs[2]+  # just CO, CO2, CH4
    datareturn$ERtoCOint[3] *datareturn$nCs[3]
  datareturn$TC2COint = sum(datareturn$ERtoCOint*datareturn$nCs, na.rm=TRUE) # all carbon containing species (that I have looked up mWs for at the moment)
  datareturn$EF1COint   = cContent * datareturn$mWs/12 * datareturn$ERtoCOint /datareturn$TC1COint
  datareturn$EF2COint   = cContent * datareturn$mWs/12 * datareturn$ERtoCOint /datareturn$TC2COint
  # 00000 integration method using CO with datafilling0000000
  datareturn$TC1COintfill = datareturn$ERtoCOintfill[1]*datareturn$nCs[1]+ datareturn$ERtoCOintfill[2]*datareturn$nCs[2]+  # just CO, CO2, CH4
    datareturn$ERtoCOintfill[3] *datareturn$nCs[3]
  datareturn$TC2COintfill = sum(datareturn$ERtoCOintfill*datareturn$nCs, na.rm=TRUE) # all carbon containing species (that I have looked up mWs for at the moment)
  datareturn$EF1COintfill   = cContent * datareturn$mWs/12 * datareturn$ERtoCOintfill /datareturn$TC1COintfill
  datareturn$EF2COintfill   = cContent * datareturn$mWs/12 * datareturn$ERtoCOintfill /datareturn$TC2COintfill
  
  # set EF & ER to zero if bad pass
  datareturn$EF1 = datareturn$EF1 * mult 
  datareturn$EF2 = datareturn$EF2 * mult 
  datareturn$ERtoX = datareturn$ERtoX * mult 
  datareturn$R2toCO  = datareturn$R2toCO * mult 
  
  datareturn$EF1int = datareturn$EF1int * mult 
  datareturn$EF2int = datareturn$EF2int * mult 
  datareturn$ERtoXint = datareturn$ERtoXint * mult 
  datareturn$EF1intfill = datareturn$EF1intfill * mult 
  datareturn$EF2intfill = datareturn$EF2intfill * mult 
  datareturn$ERtoXintfill = datareturn$ERtoXintfill * mult 
  
  datareturn$EF1COint = datareturn$EF1COint * mult 
  datareturn$EF2COint = datareturn$EF2COint * mult 
  datareturn$ERtoCOint = datareturn$ERtoCOint * mult 
  datareturn$EF1COintfill = datareturn$EF1COintfill * mult 
  datareturn$EF2COintfill = datareturn$EF2COintfill * mult 
  datareturn$ERtoCOintfill = datareturn$ERtoCOintfill * mult 
  
  # useful variables
  
  datareturn$TEDR = mean(pass$TEDR, na.rm = TRUE) 
  datareturn$WS = mean(pass$Wind_Speed_YANG, na.rm=TRUE) #Wind_Speed, m/s, Met_WindSpeed_InSitu_None, limited to where Roll_Angle does not exceed 5 degrees
  
  datareturn$WD = mean(pass$Wind_Direction_YANG, na.rm=TRUE) #m/s   Wind_Direction, deg 0-360 cw from +y, Met_WindDirection_InSitu_None, limited to where Roll_Angle does not exceed 5 degrees
  datareturn$VS = mean(pass$Vertical_Speed_YANG, na.rm=TRUE)#m/s
  datareturn$lat = mean(pass$Latitude.x, na.rm=TRUE)
  datareturn$lon = mean(pass$Longitude.x, na.rm=TRUE)
  datareturn$RH  = median(pass$RHw_DLH_DISKIN, na.rm=TRUE) #  Relative_Humidity, percent, Met_RelativeHumidityWater_InSitu_None
  datareturn$H2O = median(pass$H2O_DLH_DISKIN, na.rm=TRUE) #  Mixing_Ratio, g/kg, Met_H2OMR_InSitu_None
  
  datareturn$Radar_Altitude =  median(pass$Radar_Altitude_YANG, na.rm=TRUE) #Radar_Altitude, ft, Platform_AltitudeAGL_InSitu_None
  datareturn$T= median(pass$Total_Air_Temp_YANG, na.rm=TRUE)  #Total_Air_Temp, degrees Celcius, None
  datareturn$Static_Air_Temp = median(pass$Static_Air_Temp_YANG)#, degrees Celcius, Met_StaticAirTemperature_InSitu_None
  #Potential_Temp, Kelvin, Met_PotentialTemperature_InSitu_None
  #Dew_Point, degrees Celcius, Met_DewPoint_InSitu_None
  
  datareturn$IR_Surf_Temp = median(pass$IR_Surf_Temp_YANG, na.rm=TRUE)#, degrees Celcius, Met_SurfaceTemperature_InSitu_None
  #Static_Pressure, hPa, Met_StaticPressure_InSitu_None
  
  
  #---- return ----
  return(datareturn)
}
smallfun = function(pass,background,xKT,sKT, xKTE, sKTE, SLOW){
  require(imputeTS)
  
  pass = as.data.frame(pass) ; background = as.data.frame(background)
  ccc = colnames(pass); cccB = colnames(background)
  tmpX = which(ccc == xKT); tmpS = which(ccc == sKT)
  tmpXB = which(cccB == xKT) ;tmpSB = which(cccB == sKT)
  indKTX = which(is.finite(pass[,tmpX])) # finite X
  indKTS = which(is.finite(pass[,tmpS])) # finite S
  indKTSX = which(is.finite(pass[,tmpS]) & is.finite(pass[,tmpX])) # finite both
  uval = length(unique(pass[,tmpS])) # seems like we need at least 4 unique points > NaN
  pctVAL = round(length(indKTS)/length(indKTX), digits=1) # number of finite points greater than zero
  # kludge for CO2
  if (pctVAL > 0.4 & sKT == 'CO2_ppb'){ pctVAL = 0.61}
  if (length(indKTSX) > 3 ){
    testcor = (cor(pass[,tmpS],pass[,tmpX], use='pairwise.complete.obs'))
    testcor = testcor^2
  } else{ testcor =0}
  if (is.na(testcor)){ testcor = 0}
  # need to have at least 4? points that are different (non zero standard deviation)
  # must have at least 50% of data points of CO?
  # do I probably need to fill in the missing data?
  # don't bother if R^2 is less than 0.75?
  # lets filter for R2 after
  
  useerror = 1 # consider observational uncertainty?
  ttS = pass[,tmpS] ; ttX = pass[,tmpX] ; 
  useYork = 0 # only for CO vs. CO2
  if (sKT == 'CO_DACOM_DISKIN'){ useYork == 1}
  # print(c('SMALLFUN',sKT, pctVAL,uval, testcor))
  if (SLOW == 0){
    if (pctVAL > 0.6 & sd(pass[,tmpS], na.rm=TRUE) != 0 & uval >= 4 & length(indKTSX) >= 4){ # 4 unique points, 4 points total, 60% of CO2 data 
      if (useerror == 1){
        if (useYork == 1){pKT = york_regression(pass, x=xKT,y=sKT, variance_x = xKTE, variance_y = sKTE )}
        if (useYork == 0){pKT = lm(pass[,tmpS]~pass[,tmpX])}
      } else{ 
        pass$xKTE = 0
        pass$sKTE = 0
        if (useYork == 1){pKT = york_regression(pass, x=xKT,y=sKT, variance_x = xKTE, variance_y = sKTE )}
        if (useYork == 0){pKT = lm(pass[,tmpS]~pass[,tmpX])}
        
      }
      # y-species
      
      # ---- Get background  -----
      pKT$BGS = mean(background[,tmpSB], na.rm=TRUE)
      pKT$BGX = mean(background[,tmpXB], na.rm=TRUE)
      if (is.nan(pKT$BGS) | is.na(pKT$BGS)){
        if (useYork == 1){
          if (pKT$intercept < 0){ pKT$intercept = 0}
          pKT$BGS = pKT$intercept
          
        }
        if (useYork == 0){
          if (pKT$coefficients[1] < 0){ pKT$coefficients[1] = 0}
          pKT$BGS = pKT$coefficients[1]
          pKT$intercept = pKT$coefficients[2]
          pKT$uncertainty_slope = NaN
          
        }
        # should I set to zero?
      }
      pKT$ExcessS = sum(ttS - pKT$BGS, na.rm=TRUE)
      pKT$ExcessX = sum(ttX - pKT$BGX, na.rm=TRUE)
      if (pKT$ExcessS < 0){pKT$ExcessS = NaN}
      
      pKT$ERint = pKT$ExcessS/pKT$ExcessX
      # see what happens if I datafill
      ind = which(is.nan(pass[,tmpS])) ; ind2 = which(is.nan(pass[,tmpX]))
      
      #if (length(ind) > 0){print(c("Datafilling", sKT))}
      #if (length(ind2) > 0){print(c("Datafilling", xKT))}
      ttS = na_interpolation(ttS) # fill NaNs
      ttX = na_interpolation(ttX) # fill NaNs
      
      pKT$ExcessSfill = sum(ttS - pKT$BGS, na.rm=TRUE)
      pKT$ExcessXfill = sum(ttX - pKT$BGX, na.rm=TRUE)
      pKT$ERintfill = pKT$ExcessSfill/pKT$ExcessXfill
      pKT$r_sq = (cor(pass[,tmpS],pass[,tmpX], use='pairwise.complete.obs'))^2
      # can't have a negative ER!
      if (useYork == 1){
        if (pKT$slope < 0){ pKT$slope = NaN}
      }
      if (useYork == 0){
        pKT$uncertainty_slope = NaN
        pKT$slope = pKT$coefficients[2]
        pKT$intercept = pKT$coefficients[1]
        if (pKT$slope < 0){ pKT$slope = NaN}
      }
    } else{
      pKT = data.frame('slope'=NaN,'intercept'=NaN,'uncertainty_slope'=NaN,'r_sq'=NaN,
                       'ERint' = NaN, 'ERintfill'=NaN,'ExcessS' = NaN, 'ExcessX'=NaN, 
                       'BGX' = NaN, 'BGS' = NaN) 
    }
  }
  if (SLOW == 1){
    if (sd(pass[,tmpS], na.rm=TRUE) != 0 & uval > 4 ){ 
      if (useerror == 1){
        if (useYork == 1){pKT = york_regression(pass, x=xKT,y=sKT, variance_x = xKTE, variance_y = sKTE )}
        if (useYork == 0){pKT = lm(pass[,tmpS]~pass[,tmpX])}
      } else{ 
        pass$xKTE = 0
        pass$sKTE = 0
        if (useYork == 1){pKT = york_regression(pass, x=xKT,y=sKT, variance_x = xKTE, variance_y = sKTE )}
        if (useYork == 0){pKT = lm(pass[,tmpS]~pass[,tmpX])}
        
      }
      # y-species
      
      # ---- Get background from first 5 points -----
      # old method
      #pKT$BGS = quantile(background[,tmpSB],probs=c(0.1,0.2,0.3,0.4,0.5,0.75,1), type = 8, na.rm=TRUE)[1] # 10% percentile
      #pKT$BGX = quantile(background[,tmpXB],probs=c(0.1,0.2,0.3,0.4,0.5,0.75,1),type = 8, na.rm=TRUE)[1]
      # new method
      pKT$BGS = mean(background[,tmpSB], na.rm=TRUE)
      pKT$BGX =mean(background[,tmpXB], na.rm=TRUE)
      
      pKT$ExcessS = sum(ttS - pKT$BGS, na.rm=TRUE)
      pKT$ExcessX = sum(ttX - pKT$BGX, na.rm=TRUE)
      pKT$ERint = pKT$ExcessS/pKT$ExcessX
      # see what happens if I datafill
      ind = which(is.nan(pass[,tmpS])) ; ind2 = which(is.nan(pass[,tmpX]))
      if (length(ind) > 0){print(c("Datafilling", sKT))}
      if (length(ind2) > 0){print(c("Datafilling", xKT))}
      ttS = na_interpolation(ttS) # fill NaNs
      ttX = na_interpolation(ttX) # fill NaNs
      
      pKT$ExcessSfill = sum(ttS - pKT$BGS, na.rm=TRUE)
      pKT$ExcessXfill = sum(ttX - pKT$BGX, na.rm=TRUE)
      pKT$ERintfill = pKT$ExcessSfill/pKT$ExcessXfill
      
      if (useYork == 1){
        if (pKT$slope < 0){ pKT$slope = NaN}
      }
      if (useYork == 0){
        pKT$slope = pKT$coefficients[2]
        pKT$intercept = pKT$coefficients[1]
        pKT$uncertainty_slope = NaN
        
        if (pKT$slope < 0){ pKT$slope = NaN}
      }
    } else{
      pKT = data.frame('slope'=NaN,'uncertainty_slope'=NaN,'r_sq'=NaN,
                       'ERint' = NaN, 'ERintfill'=NaN,'ExcessS' = NaN, 'ExcessX'=NaN, 
                       'BGX' = NaN, 'BGS' = NaN) 
    }
  }
  return(pKT)
}
# ----- ********* Functions for averaging instruments -----
mergelines = function(alldata, variable1, variable2){
  ind = which(alldata$variable == variable1)
  t1= alldata[ind,]
  ind = which(alldata$variable== variable2)
  t2=alldata[ind,]
  for (i in 1:length(t1$variable)){
    newline1 = t1[i,]
    newline2 = t2[i,]
    # Going to replace the relevant variables in newline1 with the avg of newline1 and newline2
    newline1$variable = paste0(t1$names[1],'_',t1$PI[1], '_',t2$PI[1])
    newline1$PI = paste(t1$PI[1],'+',t2$PI[1],sep='')
    
    # 1hz
    newline1$EF1.1hz = mean(c(newline1$EF1.1hz, newline2$EF1.1hz), na.rm=TRUE)
    newline1$EF1CO.1hz = mean(c(newline1$EF1CO.1hz, newline2$EF1CO.1hz), na.rm=TRUE)
    
    newline1$EF2.1hz = mean(c(newline1$EF2.1hz, newline2$EF2.1hz), na.rm=TRUE)
    newline1$EF2CO.1hz = mean(c(newline1$EF2CO.1hz, newline2$EF2CO.1hz), na.rm=TRUE)
    
    newline1$ERtoCO.1hz = mean(c(newline1$ERtoCO.1hz, newline2$ERtoCO.1hz), na.rm=TRUE)
    newline1$ERtoX.1hz = mean(c(newline1$ERtoX.1hz, newline2$ERtoX.1hz), na.rm=TRUE)
    
    newline1$ERtoCOintfill.1hz= mean(c(newline1$ERtoCOintfill.1hz, newline2$ERtoCOintfill.1hz), na.rm=TRUE)
    newline1$ERtoXintfill.1hz = mean(c(newline1$ERtoXintfill.1hz, newline2$ERtoXintfill.1hz), na.rm=TRUE)
    
    # 5hz
    newline1$EF1.5hz = mean(c(newline1$EF1.5hz, newline2$EF1.5hz), na.rm=TRUE)
    newline1$EF1CO.5hz = mean(c(newline1$EF1CO.5hz, newline2$EF1CO.5hz), na.rm=TRUE)
    
    newline1$EF2.5hz = mean(c(newline1$EF2.5hz, newline2$EF2.5hz), na.rm=TRUE)
    newline1$EF2CO.5hz = mean(c(newline1$EF2CO.5hz, newline2$EF2CO.5hz), na.rm=TRUE)
    
    newline1$ERtoCO.5hz = mean(c(newline1$ERtoCO.5hz, newline2$ERtoCO.5hz), na.rm=TRUE)
    newline1$ERtoX.5hz = mean(c(newline1$ERtoX.5hz, newline2$ERtoX.5hz), na.rm=TRUE)
    
    newline1$ERtoCOintfill.5hz= mean(c(newline1$ERtoCOintfill.5hz, newline2$ERtoCOintfill.5hz), na.rm=TRUE)
    newline1$ERtoXintfill.5hz = mean(c(newline1$ERtoXintfill.5hz, newline2$ERtoXintfill.5hz), na.rm=TRUE)
    
    # Final
    newline1$FinalEF= mean(c(newline1$FinalEF, newline2$FinalEF), na.rm=TRUE)
    newline1$FinalR2 = mean(c(newline1$FinalR2, newline2$FinalR2), na.rm=TRUE)
    newline1$FinalERtoCO = mean(c(newline1$FinalERtoCO, newline2$FinalERtoCO), na.rm=TRUE)
    
    alldata = rbind(alldata, newline1)
  }
  
  # set USEME for the indivudal species to 0
  # ---------------- Determine which EFs are averaged, used, or not used -------
  # ------- Taking the average of >1 measurement
  # Now using averaged HCHO, not using TOGA
  ind = which(alldata$variable == variable1 | alldata$variable == variable2)
  alldata$USEME[ind] = 0
  
  return(alldata)
}
mergelines3 = function(alldata, variable1, variable2,variable3){
  ind = which(alldata$variable == variable1)
  t1= alldata[ind,]
  ind = which(alldata$variable== variable2)
  t2=alldata[ind,]
  ind = which(alldata$variable== variable3)
  t3=alldata[ind,]
  for (i in 1:length(t1$variable)){
    newline1 = t1[i,]
    newline2 = t2[i,]
    newline3 = t3[i,]
    
    # Going to replace the relevant variables in newline1 with the avg of newline1 and newline2
    newline1$variable = paste0(t1$names[1],'_',t1$PI[1], '_',t2$PI[1],'_',t3$P1[1])
    
    newline1$PI = paste(t1$PI[1],'+',t2$PI[1],'+',t3$PI[1],sep='')
    
    # 1hz
    newline1$EF1.1hz = mean(c(newline1$EF1.1hz, newline2$EF1.1hz, newline3$EF1.1hz), na.rm=TRUE)
    newline1$EF1CO.1hz = mean(c(newline1$EF1CO.1hz, newline2$EF1CO.1hz, newline3$EF1CO.1hz), na.rm=TRUE)
    
    newline1$EF2.1hz = mean(c(newline1$EF2.1hz, newline2$EF2.1hz, newline3$EF2.1hz), na.rm=TRUE)
    newline1$EF2CO.1hz = mean(c(newline1$EF2CO.1hz, newline2$EF2CO.1hz, newline3$EF2CO.1hz), na.rm=TRUE)
    
    newline1$ERtoCO.1hz = mean(c(newline1$ERtoCO.1hz, newline2$ERtoCO.1hz, newline3$ERtoCO.1hz), na.rm=TRUE)
    newline1$ERtoX.1hz = mean(c(newline1$ERtoX.1hz, newline2$ERtoX.1hz, newline3$ERtoX.1hz), na.rm=TRUE)
    
    newline1$ERtoCOintfill.1hz= mean(c(newline1$ERtoCOintfill.1hz, newline2$ERtoCOintfill.1hz, newline3$ERtoCOintfill.1hz), na.rm=TRUE)
    newline1$ERtoXintfill.1hz = mean(c(newline1$ERtoXintfill.1hz, newline2$ERtoXintfill.1hz, newline3$ERtoXintfill.1hz), na.rm=TRUE)
    
    # 5hz
    newline1$EF1.5hz = mean(c(newline1$EF1.5hz, newline2$EF1.5hz, newline3$EF1.5hz), na.rm=TRUE)
    newline1$EF1CO.5hz = mean(c(newline1$EF1CO.5hz, newline2$EF1CO.5hz, newline3$EF1CO.5hz), na.rm=TRUE)
    
    newline1$EF2.5hz = mean(c(newline1$EF2.5hz, newline2$EF2.5hz, newline3$EF2.5hz), na.rm=TRUE)
    newline1$EF2CO.5hz = mean(c(newline1$EF2CO.5hz, newline2$EF2CO.5hz, newline3$EF2CO.5hz), na.rm=TRUE)
    
    newline1$ERtoCO.5hz = mean(c(newline1$ERtoCO.5hz, newline2$ERtoCO.5hz, newline3$ERtoCO.5hz), na.rm=TRUE)
    newline1$ERtoX.5hz = mean(c(newline1$ERtoX.5hz, newline2$ERtoX.5hz, newline3$ERtoX.5hz), na.rm=TRUE)
    
    newline1$ERtoCOintfill.5hz= mean(c(newline1$ERtoCOintfill.5hz, newline2$ERtoCOintfill.5hz, newline3$ERtoCOintfill.5hz), na.rm=TRUE)
    newline1$ERtoXintfill.5hz = mean(c(newline1$ERtoXintfill.5hz, newline2$ERtoXintfill.5hz,newline3$ERtoXintfill.5hz), na.rm=TRUE)
    # Chosen
    #newline1$ChosenEF.5hz= mean(c(newline1$ChosenEF.5hz, newline2$ChosenEF.5hz, newline3$ChosenEF.5hz), na.rm=TRUE)
    #newline1$ChosenEF.1hz= mean(c(newline1$ChosenEF.1hz, newline2$ChosenEF.1hz, newline3$ChosenEF.1hz), na.rm=TRUE)
    #newline1$ChosenEF.R2.5hz = mean(c(newline1$ChosenEF.R2.5hz, newline2$ChosenEF.R2.5hz, newline3$ChosenEF.R2.5hz), na.rm=TRUE)
    #newline1$ChosenEF.R2.1hz = mean(c(newline1$ChosenEF.R2.1hz, newline2$ChosenEF.R2.1hz, newline3$ChosenEF.R2.1hz), na.rm=TRUE)
    # Final
    newline1$FinalEF= mean(c(newline1$FinalEF, newline2$FinalEF,newline3$FinalEF), na.rm=TRUE)
    newline1$FinalR2 = mean(c(newline1$FinalR2, newline2$FinalR2, newline3$FinalR2), na.rm=TRUE)
    newline1$FinalERtoCO = mean(c(newline1$FinalERtoCO, newline2$FinalERtoCO,newline3$FinalERtoCO), na.rm=TRUE)
    
    alldata = rbind(alldata, newline1)
  }
  
  # set USEME for the indivudal species to 0
  # ---------------- Determine which EFs are averaged, used, or not used -------
  # ------- Taking the average of >1 measurement
  # Now using averaged HCHO, not using TOGA
  ind = which(alldata$variable == variable1 | alldata$variable == variable2 | alldata$variable == variable3)
  alldata$USEME[ind] = 0
  
  return(alldata)
}
mergelines4 = function(alldata, variable1, variable2,variable3, variable4){
  ind = which(alldata$variable == variable1)
  t1= alldata[ind,]
  ind = which(alldata$variable== variable2)
  t2=alldata[ind,]
  ind = which(alldata$variable== variable3)
  t3=alldata[ind,]
  ind = which(alldata$variable== variable4)
  t4=alldata[ind,]
  for (i in 1:length(t1$variable)){
    newline1 = t1[i,]
    newline2 = t2[i,]
    newline3 = t3[i,]
    newline4 = t4[i,]
    
    newline1$USEME = 1
    # Going to replace the relevant variables in newline1 with the avg of newline1 and newline2
    newline1$variable = paste0(t1$names[1],'_',t1$PI[1], '_',t2$PI[1],'_',t3$P1[1],'_',t4$P1[1])
    
    newline1$PI = paste(t1$PI[1],'+',t2$PI[1],'+',t3$PI[1],'+',t4$PI[1],sep='')
    
    # 1hz
    newline1$EF1.1hz = mean(c(newline1$EF1.1hz, newline2$EF1.1hz, newline3$EF1.1hz, newline4$EF1.1hz), na.rm=TRUE)
    newline1$EF1CO.1hz = mean(c(newline1$EF1CO.1hz, newline2$EF1CO.1hz, newline3$EF1CO.1hz, newline4$EF1CO.1hz), na.rm=TRUE)
    
    newline1$EF2.1hz = mean(c(newline1$EF2.1hz, newline2$EF2.1hz, newline3$EF2.1hz, newline4$EF2.1hz), na.rm=TRUE)
    newline1$EF2CO.1hz = mean(c(newline1$EF2CO.1hz, newline2$EF2CO.1hz, newline3$EF2CO.1hz, newline4$EF2CO.1hz), na.rm=TRUE)
    # fix commented lines eventually
    #newline1$ERtoCO.1hz = mean(c(newline1$ERtoCO.1hz, newline2$ERtoCO.1hz, newline3$ERtoCO.1hz), na.rm=TRUE)
    #newline1$ERtoX.1hz = mean(c(newline1$ERtoX.1hz, newline2$ERtoX.1hz, newline3$ERtoX.1hz), na.rm=TRUE)
    
    #newline1$ERtoCOintfill.1hz= mean(c(newline1$ERtoCOintfill.1hz, newline2$ERtoCOintfill.1hz, newline3$ERtoCOintfill.1hz), na.rm=TRUE)
    #newline1$ERtoXintfill.1hz = mean(c(newline1$ERtoXintfill.1hz, newline2$ERtoXintfill.1hz, newline3$ERtoXintfill.1hz), na.rm=TRUE)
    
    # 5hz
    #newline1$EF1.5hz = mean(c(newline1$EF1.5hz, newline2$EF1.5hz, newline3$EF1.5hz), na.rm=TRUE)
    #newline1$EF1CO.5hz = mean(c(newline1$EF1CO.5hz, newline2$EF1CO.5hz, newline3$EF1CO.5hz), na.rm=TRUE)
    
    #newline1$EF2.5hz = mean(c(newline1$EF2.5hz, newline2$EF2.5hz, newline3$EF2.5hz), na.rm=TRUE)
    #newline1$EF2CO.5hz = mean(c(newline1$EF2CO.5hz, newline2$EF2CO.5hz, newline3$EF2CO.5hz), na.rm=TRUE)
    
    #newline1$ERtoCO.5hz = mean(c(newline1$ERtoCO.5hz, newline2$ERtoCO.5hz, newline3$ERtoCO.5hz), na.rm=TRUE)
    #newline1$ERtoX.5hz = mean(c(newline1$ERtoX.5hz, newline2$ERtoX.5hz, newline3$ERtoX.5hz), na.rm=TRUE)
    
    #newline1$ERtoCOintfill.5hz= mean(c(newline1$ERtoCOintfill.5hz, newline2$ERtoCOintfill.5hz, newline3$ERtoCOintfill.5hz), na.rm=TRUE)
    #newline1$ERtoXintfill.5hz = mean(c(newline1$ERtoXintfill.5hz, newline2$ERtoXintfill.5hz,newline3$ERtoXintfill.5hz), na.rm=TRUE)
    # Chosen
    #newline1$ChosenEF.5hz= mean(c(newline1$ChosenEF.5hz, newline2$ChosenEF.5hz, newline3$ChosenEF.5hz, newline4$ChosenEF.5hz), na.rm=TRUE)
    #newline1$ChosenEF.1hz= mean(c(newline1$ChosenEF.1hz, newline2$ChosenEF.1hz, newline3$ChosenEF.1hz, newline4$ChosenEF.1hz), na.rm=TRUE)
    #newline1$ChosenEF.R2.5hz = mean(c(newline1$ChosenEF.R2.5hz, newline2$ChosenEF.R2.5hz, newline3$ChosenEF.R2.5hz, newline4$ChosenEF.R2.5hz), na.rm=TRUE)
    #newline1$ChosenEF.R2.1hz = mean(c(newline1$ChosenEF.R2.1hz, newline2$ChosenEF.R2.1hz, newline3$ChosenEF.R2.1hz, newline4$ChosenEF.R2.1hz), na.rm=TRUE)
    # Final
    newline1$FinalEF= mean(c(newline1$FinalEF, newline2$FinalEF,newline3$FinalEF,newline4$FinalEF), na.rm=TRUE)
    newline1$FinalR2 = mean(c(newline1$FinalR2, newline2$FinalR2, newline3$FinalR2, newline4$FinalR2), na.rm=TRUE)
    newline1$FinalERtoCO = mean(c(newline1$FinalERtoCO, newline2$FinalERtoCO,newline3$FinalERtoCO,newline4$FinalERtoCO), na.rm=TRUE)
    
    alldata = rbind(alldata, newline1)
  }
  
  # set USEME for the indivudal species to 0
  # ---------------- Determine which EFs are averaged, used, or not used -------
  # ------- Taking the average of >1 measurement
  # Now using averaged HCHO, not using TOGA
  ind = which(alldata$variable == variable1 | alldata$variable == variable2 | alldata$variable == variable3 | alldata$variable == variable4)
  alldata$USEME[ind] = 0
  
  return(alldata)
}

# -------*********  Plotting functions --------
plotpassSPECIES = function(pass,species, shortname){
  cc = colnames(pass)
  par(mfrow=c(1,2))
  ind = which(cc == species)
  yy=max(pass[,ind], na.rm=TRUE)
  plot(pass$Time_Start, pass[,ind], ylab = species,type='l',
       main = "", xlab = "Time",ylim=c(0, yy))
  lines(pass$Time_Start, pass$Furan_NOAAPTR_ppbv_WARNEKE/2, col='green')
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, xaxt = "n", yaxt = "n",type='l',
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,max(pass$CO_DACOM_DISKIN, na.rm=TRUE)))
  axis(side = 4, cex=cz)
  mtext("CO, ppb", side = 4, line = 3)
  legend("topright", c(shortname, "CO",'Furan'), 
         col = c("black", "red","green"), bty='n',lty = c(1,1, 2))
  
  plot(pass$CO_DACOM_DISKIN, pass$C2H5OH_NOAAPTR_ppbv_WARNEKE, pch=19, xlab='CO, ppb', ylab='C2H5OH, ppb')
}
plotpass5hz = function(pass){
  co_bg = quantile(pass$CO_DACOM_DISKIN,probs=c(0.1,0.2,0.3,0.4,0.5,0.75,1), type = 8, na.rm=TRUE)[1] # 10% percentile
  co2_bg = quantile(pass$CO2_7000_ppm_DISKIN,probs=c(0.1,0.2,0.3,0.4,0.5,0.75,1),type = 8, na.rm=TRUE)[1]
  furan_bg = quantile(pass$Furan_NOAAPTR_ppbv_WARNEKE,probs=c(0.1,0.2,0.3,0.4,0.5,0.75,1),type = 8, na.rm=TRUE)[1]
  so2_bg = quantile(pass$SO2_LIF_ROLLINS,probs=c(0.1,0.2,0.3,0.4,0.5,0.75,1),type = 8, na.rm=TRUE)[1]
  no_bg = quantile(pass$NO_LIF_ROLLINS,probs=c(0.1,0.2,0.3,0.4,0.5,0.75,1),type = 8, na.rm=TRUE)[1]
  ch2o_bg=quantile(pass$CH2O_ISAF_HANISCO,probs=c(0.1,0.2,0.3,0.4,0.5,0.75,1),type = 8, na.rm=TRUE)[1]
  phenol_bg=quantile(pass$Phenol_NOAAPTR_ppbv_WARNEKE,probs=c(0.1,0.2,0.3,0.4,0.5,0.75,1),type = 8, na.rm=TRUE)[1]
  ethanol_bg=quantile(pass$C2H5OH_NOAAPTR_ppbv_WARNEKE,probs=c(0.1,0.2,0.3,0.4,0.5,0.75,1),type = 8, na.rm=TRUE)[1]
  pan_bg=quantile(pass$PAN_GTCIMS_HUEY,probs=c(0.1,0.2,0.3,0.4,0.5,0.75,1),type = 8, na.rm=TRUE)[1]
  prpn_bg=quantile(pass$PROPENE.HN.1Hz_CIT_WENNBERG,probs=c(0.1,0.2,0.3,0.4,0.5,0.75,1),type = 8, na.rm=TRUE)[1]
  benzene_bg = quantile(pass$Acrolein_NOAAPTR_ppbv_WARNEKE,probs=c(0.1,0.2,0.3,0.4,0.5,0.75,1),type = 8, na.rm=TRUE)[1]
  
  cz = 1.2
  par(mfrow=c(4,2),mar = c(3, 4, 2, 4), cex=cz)
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  tt=3
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",type='o',lwd=tt,
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  abline(h=co2_bg)
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, xaxt = "n", yaxt = "n",type='o',lwd=tt,pch=16,cex=0.5,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,max(pass$CO_DACOM_DISKIN, na.rm=TRUE)))
  abline(h=co_bg, col='red')
  axis(side = 4, cex=cz)
  mtext("CO, ppb", side = 4, line = 2, cex=cz)
  legend("topright", c("CO2", "CO"), col = c("black", "red"), bty='n',lty = c(1, 2))
  
  
  # ---- correlations -----
  ind = which(is.finite(pass$CO2_7000_ppm_DISKIN))
  if (length(ind) > 3){
    cc1 = cor.test(pass$CO2_7000_ppm_DISKIN,pass$CO_DACOM_DISKIN)
    cc1 = cc1$estimate^2
  } else{
    cc1 = NaN
  }
  text(median(pass$Time_Start),max(pass$CO_DACOM_DISKIN)*0.9,paste("r2=",round(cc1, digits = 2) ))
  print(c("r^2, CO vs. CO2",cc1))
  
  ind = which(is.finite(pass$Acrolein_NOAAPTR_ppbv_WARNEKE))
  if (length(ind) > 3){
    ccFURAN = cor.test(pass$CO_DACOM_DISKIN, pass$Acrolein_NOAAPTR_ppbv_WARNEKE)
    ccFURAN_E=ccFURAN$estimate^2
    print(c("r^2, CO vs. Acrolein",ccFURAN_E))
  } else{  
    ccFURAN_E=NaN 
    print(c("r^2, CO vs. Acrolein",NaN))}
  
  ind = which(is.finite(pass$HCOOH_NOAAPTR_ppbv_WARNEKE))
  if (length(ind) > 3){
    ccHCOOH= cor.test(pass$CO_DACOM_DISKIN, pass$HCOOH_NOAAPTR_ppbv_WARNEKE)
    ccHCOOH_E=ccFURAN$estimate^2
    print(c("r^2, CO vs. HCOOH",ccHCOOH_E))
  } else{  
    ccHCOOH_E=NaN 
    print(c("r^2, CO vs. HCOOH",NaN))}
  
  ind = which(is.finite(pass$C2H5OH_NOAAPTR_ppbv_WARNEKE))
  if (length(ind) > 3){
    ccC2H5OH= cor.test(pass$CO_DACOM_DISKIN, pass$C2H5OH_NOAAPTR_ppbv_WARNEKE)
    ccC2H5OH_E=ccC2H5OH$estimate^2
    print(c("r^2, CO vs. C2H5OH",ccC2H5OH_E))
  } else{  
    ccHCOOH_E=NaN 
    print(c("r^2, CO vs. C2H5OH",NaN))}
  
  ind = which(is.finite(pass$CH3OH_NOAAPTR_ppbv_WARNEKE))
  if (length(ind) > 3){
    ccC2H5OH= cor.test(pass$CO_DACOM_DISKIN, pass$CH3OH_NOAAPTR_ppbv_WARNEKE)
    ccC2H5OH_E=ccC2H5OH$estimate^2
    print(c("r^2, CO vs. CH3OH",ccC2H5OH_E))
  } else{  
    ccHCOOH_E=NaN 
    print(c("r^2, CO vs. CH3OH",NaN))}
  
  
  ind = which(is.finite(pass$CH2O_ISAF_HANISCO))
  if (length(ind) > 3){
    ccCH2O = cor.test(pass$CO_DACOM_DISKIN, pass$CH2O_ISAF_HANISCO)
    ccCH2O_E = ccCH2O$estimate^2
    print(c("r^2, CO vs. CH2O",ccCH2O_E))
  } else{  
    ccCH2O_E = NaN
    print(c("r^2, CO vs. CH2O",NaN))
  }
  ind = which(is.finite(pass$NO_LIF_ROLLINS))
  if (length(ind) > 3){
    ccNO = cor.test(pass$CO_DACOM_DISKIN, pass$NO_LIF_ROLLINS)
    ccNO_E = ccNO$estimate^2
    print(c("r^2, CO vs. NO",ccNO_E))
  } else{   
    ccNO_E = NaN
    print(c("r^2, CO vs. NO",NaN))}
  
  ind = which(is.finite(pass$SO2_LIF_ROLLINS))
  if (length(ind) > 3){
    ccSO2 = cor.test(pass$CO_DACOM_DISKIN, pass$SO2_LIF_ROLLINS)
    ccSO2_E = ccSO2$estimate^2
    print(c("r^2, CO vs. SO2",ccSO2_E))
  } else{  
    ccSO2_E = NaN
    print(c("r^2, CO vs. SO2",ccSO2_E))}
  
  ind = which(is.finite(pass$PHENOL.10Hz_CIT_WENNBERG))
  if (length(ind) > 3){
    ccPHENOL = cor.test(pass$CO_DACOM_DISKIN, pass$PHENOL.10Hz_CIT_WENNBERG)
    ccPHENOL_E = ccPHENOL$estimate^2
    print(c("r^2, CO vs. PHENOL",ccPHENOL_E))
  } else{  
    ccPHENOL_E = NaN
    print(c("r^2, CO vs. PHENOL",NaN))}
  
  
  ind = which(is.finite(pass$PAN_GTCIMS_HUEY) & is.finite(pass$CO_DACOM_DISKIN))
  if (length(ind) > 3){
    ccPAN = cor.test(pass$CO_DACOM_DISKIN, pass$PAN_GTCIMS_HUEY)
    ccPAN_E = ccPAN$estimate^2
    print(c("r^2, CO vs. PAN",ccPAN_E))
  } else{ 
    ccPAN_E = NaN
    print(c("r^2, CO vs. PAN",NaN))}
  
  ind = which(is.finite(pass$PROPENE.HN.1Hz_CIT_WENNBERG) & is.finite(pass$CO_DACOM_DISKIN))
  if (length(ind) > 3){
    ccPROPENE = cor.test(pass$CO_DACOM_DISKIN, pass$PROPENE.HN.1Hz_CIT_WENNBERG)
    ccPROPENE_E=ccPROPENE$estimate^2
    print(c("r^2, CO vs. PROPENE-HN",ccPROPENE_E))
  } else{ 
    ccPROPENE_E = NaN
    print(c("r^2, CO vs. PROPENE-HN",NaN))}
  
  # --------- plot the rest -------
  plot(pass$CO2_7000_ppm, pass$CO_DACOM, xlab='CO2, ppm', ylab='CO, ppb')
  
  # add 
  doadd = 0
  if (doadd == 1){
    yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
    plot(pass$Time_Start, pass$CO_DACOM_DISKIN,ylab = "CO, ppb",type='o', lwd=tt,
         main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
    abline(h=co_bg)
    par(new = TRUE, cex=cz)
    yyz =max(pass$CH3OH_NOAAPTR_ppbv_WARNEKE, na.rm=TRUE)
    if (is.infinite(yyz)){ yyz = 10}
    plot(pass$Time_Start, pass$CH3OH_NOAAPTR_ppbv_WARNEKE,  xaxt = "n", yaxt = "n",type='o',
         cex=cz,
         ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yyz))
    axis(side = 4, cex=cz)
    
    mtext("CH3OH, ppb", side = 4, line = 2, cex=cz)
    legend("topright", c("CO", "CH3OH"), col = c("black", "red"), bty='n',lty = c(1, 2))
    yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
    
    plot(pass$Time_Start, pass$CO_DACOM_DISKIN,ylab = "CO, ppb",type='o', lwd=tt,
         main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
    abline(h=co_bg)
    par(new = TRUE, cex=cz)
    yyz =max(pass$HCOOH_NOAAPTR_ppbv_WARNEKE, na.rm=TRUE)
    if (is.infinite(yyz)){ yyz = 10}
    plot(pass$Time_Start, pass$HCOOH_NOAAPTR_ppbv_WARNEKE,  xaxt = "n", yaxt = "n",type='o',
         cex=cz,
         ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yyz))
    axis(side = 4, cex=cz)
    
    mtext("HCOOH, ppb", side = 4, line = 2, cex=cz)
    legend("topright", c("CO", "HCOOH"), col = c("black", "red"), bty='n',lty = c(1, 2))
  }
  # -----
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",type='o', lwd=tt,
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yyz = max(pass$SO2_LIF_ROLLINS, na.rm=TRUE)
  if (is.infinite(yyz)){ yyz = 10}
  plot(pass$Time_Start, pass$SO2_LIF_ROLLINS,  xaxt = "n", yaxt = "n",type='o', cex=cz,
       ylab = "", xlab = "", col = "red", lty = 2, lwd=tt,ylim=c(0,yyz))
  abline(h=so2_bg, col='red')
  axis(side = 4, cex=cz)
  mtext("SO2, ppt",side = 4, line = 2, cex=cz)
  legend("topright", c("CO", "SO2"), col = c("black", "red"), bty='n',lty = c(1, 2))
  text(median(pass$Time_Start),yyz*0.85,paste("r2=",round(ccSO2_E, digits = 2) ))
  
  yyz = max(pass$NO_LIF_ROLLINS, na.rm=TRUE)
  if (is.infinite(yyz)){ yyz = 10}
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",type='o',lwd=tt,
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$NO_LIF_ROLLINS,  xaxt = "n", yaxt = "n",type='o', cex=cz,
       ylab = "", xlab = "", col = "red", lty = 2, lwd=tt,ylim=c(0,yyz))
  axis(side = 4, cex=cz)
  abline(h=no_bg, col='red')
  mtext("NO, ppt", side = 4, line = 2, cex=cz)
  legend("topright", c("CO", "NO"), col = c("black", "red"), bty='n',lty = c(1, 2))
  text(median(pass$Time_Start),median(pass$NO_LIF_ROLLINS),paste("r2=",round(ccNO_E, digits = 2) ))
  
  yyz =max(pass$CH2O_ISAF_HANISCO, na.rm=TRUE)
  if (is.infinite(yyz)){ yyz = 10}
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",type='o',lwd=tt,
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$CH2O_ISAF_HANISCO,  xaxt = "n", yaxt = "n",type='o', cex=cz,
       ylab = "", xlab = "", col = "red", lwd=tt,lty = 2, ylim=c(0,yyz))
  abline(h=ch2o_bg, col='red')
  axis(side = 4, cex=cz)
  mtext("HCHO, ppt",side = 4, line = 2, cex=cz)
  legend("topright", c("CO", "HCHO"), col = c("black", "red"), bty='n',lty = c(1, 2))
  text(median(pass$Time_Start),median(pass$CH2O_ISAF_HANISCO),paste("r2=",round(ccCH2O_E, digits = 2) ))
  
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN,ylab = "CO, ppb",type='o', lwd=tt,
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  abline(h=co_bg)
  par(new = TRUE, cex=cz)
  yyz =max(pass$Acrolein_NOAAPTR_ppbv_WARNEKE , na.rm=TRUE)
  if (is.infinite(yyz)){ yyz = 10}
  plot(pass$Time_Start, pass$Acrolein_NOAAPTR_ppbv_WARNEKE,  xaxt = "n", yaxt = "n",type='o',
       cex=cz,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yyz))
  #abline(h=Acrolein_bg, col='red')
  axis(side = 4, cex=cz)
  
  mtext("Acrolein, ppt", side = 4, line = 2, cex=cz)
  legend("topright", c("CO", "Acrolein"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc1 = cor.test(pass$CO2_7000_ppm_DISKIN,pass$CO_DACOM_DISKIN)
  text(median(pass$Time_Start),yyz*0.85,paste("r2=",round(ccFURAN_E, digits = 2) ))
  
  yyz =max(pass$PHENOL.10Hz_CIT_WENNBERG, na.rm=TRUE)
  if (is.infinite(yyz)){ yyz = 10}
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",lwd=tt
       ,type='o', xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$PHENOL.10Hz_CIT_WENNBERG, xaxt = "n", yaxt = "n",type='o', cex=cz,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yyz))
  abline(h=phenol_bg, col='red')
  axis(side = 4, cex=cz)
  mtext("Phenol, ppt", side = 4, line = 2, cex=cz)
  legend("topright", c("CO", "Phenol"), col = c("black", "red"), bty='n',lty = c(1, 2))
  text(median(pass$Time_Start),yyz*0.991,paste("r2=",round(ccPHENOL_E, digits = 2) ))
  
  if (is.infinite(yyz) | is.nan(yyz)){ yyz = 10}
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",type='o',lwd=tt,
       xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yyz =max(pass$PAN_GTCIMS_HUEY, na.rm=TRUE)
  if (is.infinite(yyz) | is.nan(yyz)){ yyz = 10}
  plot(pass$Time_Start, pass$PAN_GTCIMS_HUEY,  xaxt = "n", yaxt = "n",type='o', cex=cz,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yyz))
  abline(h=pan_bg, col='red')
  abline(h=prpn_bg*100, col='blue')
  lines(pass$Time_Start, pass$PROPENE.HN.1Hz_CIT_WENNBERG*10, col='blue')
  axis(side = 4, cex=cz)
  mtext("PAN, ppt", side = 4, line = 2, cex=cz)
  legend("topright", c("CO", "PAN",'PROPENE-HN*10'), col = c("black", "red","blue"), bty='n',lty = c(1, 2,1))
  text(median(pass$Time_Start),yyz*0.991,paste("r2=",round(ccPAN_E, digits = 2) ))
  
}
plotpass5hzB = function(pass){
  
  cz = 1.2
  par(mfrow=c(4,2),mar = c(3, 4, 2, 4), cex=cz)
  yy=max(pass$CO2_7000_ppm, na.rm=TRUE)
  tt=3
  plot(pass$Time_Start, pass$CO2_7000_ppm_DISKIN, ylab = "CO2, ppm",type='o',lwd=tt,
       main = "", xlab = "Time", ylim=c(min(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE), yy))
  abline(h=co2_bg)
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, xaxt = "n", yaxt = "n",type='o',lwd=tt,pch=16,cex=0.5,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,max(pass$CO_DACOM_DISKIN, na.rm=TRUE)))
  abline(h=co_bg, col='red')
  axis(side = 4, cex=cz)
  mtext("CO, ppb", side = 4, line = 2, cex=cz)
  legend("topright", c("CO2", "CO"), col = c("black", "red"), bty='n',lty = c(1, 2))
  
  
  # ---- correlations -----
  ind = which(is.finite(pass$CO2_7000_ppm_DISKIN))
  if (length(ind) > 3){
    cc1 = cor.test(pass$CO2_7000_ppm_DISKIN,pass$CO_DACOM_DISKIN)
    cc1 = cc1$estimate^2
  } else{
    cc1 = NaN
  }
  text(median(pass$Time_Start),median(pass$CO_DACOM_DISKIN),paste("r2=",round(cc1, digits = 2) ))
  print(c("r^2, CO vs. CO2",cc1))
  
  ind = which(is.finite(pass$NO_CL_RYERSON))
  if (length(ind) > 3){
    ccNO_CL_RYERSON = cor.test(pass$CO_DACOM_DISKIN, pass$NO_CL_RYERSON)
    ccNO_CL_RYERSON=ccNO_CL_RYERSON$estimate^2
    print(c("r^2, CO vs. NO_CL_RYERSON",ccNO_CL_RYERSON))
  } else{  
    ccFURAN_E=NaN 
    print(c("r^2, CO vs. Furan",NaN))}
  
  ind = which(is.finite(pass$HCOOH_NOAAPTR_ppbv_WARNEKE))
  if (length(ind) > 3){
    ccHCOOH= cor.test(pass$CO_DACOM_DISKIN, pass$HCOOH_NOAAPTR_ppbv_WARNEKE)
    ccHCOOH_E=ccFURAN$estimate^2
    print(c("r^2, CO vs. HCOOH",ccHCOOH_E))
  } else{  
    ccHCOOH_E=NaN 
    print(c("r^2, CO vs. HCOOH",NaN))}
  
  ind = which(is.finite(pass$C2H5OH_NOAAPTR_ppbv_WARNEKE))
  if (length(ind) > 3){
    ccC2H5OH= cor.test(pass$CO_DACOM_DISKIN, pass$C2H5OH_NOAAPTR_ppbv_WARNEKE)
    ccC2H5OH_E=ccC2H5OH$estimate^2
    print(c("r^2, CO vs. C2H5OH",ccC2H5OH_E))
  } else{  
    ccHCOOH_E=NaN 
    print(c("r^2, CO vs. C2H5OH",NaN))}
  
  ind = which(is.finite(pass$CH2O_ISAF_HANISCO))
  if (length(ind) > 3){
    ccCH2O = cor.test(pass$CO_DACOM_DISKIN, pass$CH2O_ISAF_HANISCO)
    ccCH2O_E = ccCH2O$estimate^2
    print(c("r^2, CO vs. CH2O",ccCH2O_E))
  } else{  
    ccCH2O_E = NaN
    print(c("r^2, CO vs. CH2O",NaN))
  }
  ind = which(is.finite(pass$NO_LIF_ROLLINS))
  if (length(ind) > 3){
    ccNO = cor.test(pass$CO_DACOM_DISKIN, pass$NO_LIF_ROLLINS)
    ccNO_E = ccNO$estimate^2
    print(c("r^2, CO vs. NO",ccNO_E))
  } else{   
    ccNO_E = NaN
    print(c("r^2, CO vs. NO",NaN))}
  
  ind = which(is.finite(pass$SO2_LIF_ROLLINS))
  if (length(ind) > 3){
    ccSO2 = cor.test(pass$CO_DACOM_DISKIN, pass$SO2_LIF_ROLLINS)
    ccSO2_E = ccSO2$estimate^2
    print(c("r^2, CO vs. SO2",ccSO2_E))
  } else{  
    ccSO2_E = NaN
    print(c("r^2, CO vs. SO2",ccSO2_E))}
  
  ind = which(is.finite(pass$PHENOL.1Hz_CIT_WENNBERG))
  if (length(ind) > 3){
    ccPHENOL = cor.test(pass$CO_DACOM_DISKIN, pass$PHENOL.1Hz_CIT_WENNBERG)
    ccPHENOL_E = ccPHENOL$estimate^2
    print(c("r^2, CO vs. PHENOL",ccPHENOL_E))
  } else{  
    ccPHENOL_E = NaN
    print(c("r^2, CO vs. PHENOL",NaN))}
  
  ind = which(is.finite(pass$PAN_GTCIMS_HUEY) & is.finite(pass$CO_DACOM_DISKIN))
  if (length(ind) > 3){
    ccPAN = cor.test(pass$CO_DACOM_DISKIN, pass$PAN_GTCIMS_HUEY)
    ccPAN_E = ccPAN$estimate^2
    print(c("r^2, CO vs. PAN",ccPAN_E))
  } else{ 
    ccPAN_E = NaN
    print(c("r^2, CO vs. PAN",NaN))}
  
  ind = which(is.finite(pass$PROPENE.HN.1Hz_CIT_WENNBERG) & is.finite(pass$CO_DACOM_DISKIN))
  if (length(ind) > 3){
    ccPROPENE = cor.test(pass$CO_DACOM_DISKIN, pass$PROPENE.HN.1Hz_CIT_WENNBERG)
    ccPROPENE_E=ccPROPENE$estimate^2
    print(c("r^2, CO vs. PROPENE-HN",ccPROPENE_E))
  } else{ 
    ccPROPENE_E = NaN
    print(c("r^2, CO vs. PROPENE-HN",NaN))}
  
  # --------- plot the rest -------
  plot(pass$CO2_7000_ppm, pass$CO_DACOM, xlab='CO2, ppm', ylab='CO, ppb')
  
  # add 
  doadd = 0
  if (doadd == 1){
    yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
    plot(pass$Time_Start, pass$CO_DACOM_DISKIN,ylab = "CO, ppb",type='o', lwd=tt,
         main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
    abline(h=co_bg)
    par(new = TRUE, cex=cz)
    yyz =max(pass$C2H5OH_NOAAPTR_ppbv_WARNEKE, na.rm=TRUE)
    if (is.infinite(yyz)){ yyz = 10}
    plot(pass$Time_Start, pass$C2H5OH_NOAAPTR_ppbv_WARNEKE,  xaxt = "n", yaxt = "n",type='o',
         cex=cz,
         ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yyz))
    axis(side = 4, cex=cz)
    
    mtext("C2H5OH, ppb", side = 4, line = 2, cex=cz)
    legend("topright", c("CO", "C2H5OH"), col = c("black", "red"), bty='n',lty = c(1, 2))
    yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
    
    plot(pass$Time_Start, pass$CO_DACOM_DISKIN,ylab = "CO, ppb",type='o', lwd=tt,
         main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
    abline(h=co_bg)
    par(new = TRUE, cex=cz)
    yyz =max(pass$HCOOH_NOAAPTR_ppbv_WARNEKE, na.rm=TRUE)
    if (is.infinite(yyz)){ yyz = 10}
    plot(pass$Time_Start, pass$HCOOH_NOAAPTR_ppbv_WARNEKE,  xaxt = "n", yaxt = "n",type='o',
         cex=cz,
         ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yyz))
    axis(side = 4, cex=cz)
    
    mtext("HCOOH, ppb", side = 4, line = 2, cex=cz)
    legend("topright", c("CO", "HCOOH"), col = c("black", "red"), bty='n',lty = c(1, 2))
  }
  # -----
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",type='o', lwd=tt,
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yyz = max(pass$SO2_LIF_ROLLINS, na.rm=TRUE)
  if (is.infinite(yyz)){ yyz = 10}
  plot(pass$Time_Start, pass$SO2_LIF_ROLLINS,  xaxt = "n", yaxt = "n",type='o', cex=cz,
       ylab = "", xlab = "", col = "red", lty = 2, lwd=tt,ylim=c(0,yyz))
  abline(h=so2_bg, col='red')
  axis(side = 4, cex=cz)
  mtext("SO2, ppt",side = 4, line = 2, cex=cz)
  legend("topright", c("CO", "SO2"), col = c("black", "red"), bty='n',lty = c(1, 2))
  text(median(pass$Time_Start),yyz*0.85,paste("r2=",round(ccSO2_E, digits = 2) ))
  
  yyz = max(pass$NO_LIF_ROLLINS, na.rm=TRUE)
  if (is.infinite(yyz)){ yyz = 10}
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",type='o',lwd=tt,
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$NO_LIF_ROLLINS,  xaxt = "n", yaxt = "n",type='o', cex=cz,
       ylab = "", xlab = "", col = "red", lty = 2, lwd=tt,ylim=c(0,yyz))
  axis(side = 4, cex=cz)
  abline(h=no_bg, col='red')
  mtext("NO, ppt", side = 4, line = 2, cex=cz)
  legend("topright", c("CO", "NO"), col = c("black", "red"), bty='n',lty = c(1, 2))
  text(median(pass$Time_Start),median(pass$NO_LIF_ROLLINS),paste("r2=",round(ccNO_E, digits = 2) ))
  
  yyz =max(pass$CH2O_ISAF_HANISCO, na.rm=TRUE)
  if (is.infinite(yyz)){ yyz = 10}
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",type='o',lwd=tt,
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$CH2O_ISAF_HANISCO,  xaxt = "n", yaxt = "n",type='o', cex=cz,
       ylab = "", xlab = "", col = "red", lwd=tt,lty = 2, ylim=c(0,yyz))
  abline(h=ch2o_bg, col='red')
  axis(side = 4, cex=cz)
  mtext("HCHO, ppt",side = 4, line = 2, cex=cz)
  legend("topright", c("CO", "HCHO"), col = c("black", "red"), bty='n',lty = c(1, 2))
  text(median(pass$Time_Start),median(pass$CH2O_ISAF_HANISCO),paste("r2=",round(ccCH2O_E, digits = 2) ))
  
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN,ylab = "CO, ppb",type='o', lwd=tt,
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  abline(h=co_bg)
  par(new = TRUE, cex=cz)
  yyz =max(pass$Furan_NOAAPTR_ppbv_WARNEKE, na.rm=TRUE)
  if (is.infinite(yyz)){ yyz = 10}
  plot(pass$Time_Start, pass$Furan_NOAAPTR_ppbv_WARNEKE,  xaxt = "n", yaxt = "n",type='o',
       cex=cz,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yyz))
  abline(h=furan_bg, col='red')
  axis(side = 4, cex=cz)
  
  mtext("Furan, ppt", side = 4, line = 2, cex=cz)
  legend("topright", c("CO", "Furan"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc1 = cor.test(pass$CO2_7000_ppm_DISKIN,pass$CO_DACOM_DISKIN)
  text(max(pass$Time_Start)-2,yyz*0.85,paste("r2=",round(ccFURAN_E, digits = 2) ))
  
  yyz =max(pass$PHENOL.1Hz_CIT_WENNBERG, na.rm=TRUE)
  if (is.infinite(yyz)){ yyz = 10}
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",lwd=tt
       ,type='o', xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$PHENOL.1Hz_CIT_WENNBERG, xaxt = "n", yaxt = "n",type='o', cex=cz,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yyz))
  abline(h=phenol_bg, col='red')
  axis(side = 4, cex=cz)
  mtext("Phenol, ppt", side = 4, line = 2, cex=cz)
  legend("topright", c("CO", "Phenol"), col = c("black", "red"), bty='n',lty = c(1, 2))
  text(max(pass$Time_Start)-2,yyz*0.991,paste("r2=",round(ccPHENOL_E, digits = 2) ))
  
  yyz =max(pass$PAN_GTCIMS_HUEY, na.rm=TRUE)
  if (is.infinite(yyz)){ yyz = 10}
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",type='o',lwd=tt,
       xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$PAN_GTCIMS_HUEY,  xaxt = "n", yaxt = "n",type='o', cex=cz,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yyz))
  abline(h=pan_bg, col='red')
  abline(h=prpn_bg*100, col='blue')
  
  lines(pass$Time_Start, pass$PROPENE.HN.1Hz_CIT_WENNBERG*10, col='blue')
  axis(side = 4, cex=cz)
  mtext("PAN, ppt", side = 4, line = 2, cex=cz)
  legend("topright", c("CO", "PAN",'PROPENE-HN*10'), col = c("black", "red","blue"), bty='n',lty = c(1, 2,1))
  text(max(pass$Time_Start)-2,yyz*0.991,paste("r2=",round(ccPAN_E, digits = 2) ))
  
}
plotslowstuff = function(pass,startO, stopO, fire,passI, blake, blakeBG, GILMAN, GILMANBG, BECKY, BECKYBG){
  par(mfrow=c(2,1))
  ind = which(pass$Time_Start >= startO & pass$Time_Start <= stopO)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN,xlim=c(startO-5,stopO+10),main=paste(fire,passI),
       pch=19, type='o', ylab='CO_DACOM_DISKIN', ylim=c(0,max(pass$CO_DACOM_DISKIN[ind], na.rm=TRUE )))
  points(c(blake$Time_Start,blake$Time_Stop), c(blake$CO_DACOM_DISKIN_BLAKE,blake$CO_DACOM_DISKIN_BLAKE),type='o', col='red', pch=19)
  points(blake$Time_Start, blakeBG$CO_DACOM_DISKIN_BLAKE,type='o', col='red', pch=1, lwd=3)
  points(c(BECKY$Time_Start,BECKY$Time_Stop...5), c(BECKY$CO_DACOM_DISKIN_BECKY,BECKY$CO_DACOM_DISKIN_BECKY),type='o', col='green', pch=19)
  points(BECKY$Time_Start, BECKYBG$CO_DACOM_DISKIN_BECKY,type='o', col='green', pch=1, lwd=3)
  points(c(GILMAN$Time_Start,GILMAN$Time_Stop), c(GILMAN$CO_DACOM_DISKIN_GILMAN,GILMAN$CO_DACOM_DISKIN_GILMAN),type='o', col='purple', pch=19)
  points(GILMAN$Time_Start, GILMANBG$CO_DACOM_DISKIN_GILMAN,type='o', col='purple', pch=1, lwd=3)
  
  plot(pass$Time_Start, pass$Benzene_NOAAPTR_ppbv_WARNEKE*1E3,xlim=c(startO-5,stopO+10),
       pch=19, type='o', ylab='Benzene ppt', ylim=c(0,max(pass$Benzene_NOAAPTR_ppbv_WARNEKE[ind]*1E3, na.rm=TRUE )))
  points(c(blake$Time_Start,blake$Time_Stop), c(blake$Benzene_WAS_BLAKE,blake$Benzene_WAS_BLAKE),type='o', col='red', pch=19)
  points(blake$Time_Start, blakeBG$Benzene_WAS_BLAKE,type='o', col='red', pch=1, lwd=3)
  points(c(BECKY$Time_Start,BECKY$Time_Stop...5), c(BECKY$Benzene_ppt,BECKY$Benzene_ppt),type='o', col='green', pch=19)
  points(BECKY$Time_Start, BECKYBG$CO_DACOM_DISKIN_BECKY,type='o', col='green', pch=1, lwd=3)
  points(c(GILMAN$Time_Start,GILMAN$Time_Stop), c(GILMAN$Benzene_NOAAiWAS_GILMAN,GILMAN$Benzene_NOAAiWAS_GILMAN)*1E3,type='o', col='purple', pch=19)
  points(GILMAN$Time_Start, GILMANBG$Benzene_NOAAiWAS_GILMAN*1E3,type='o', col='purple', pch=1, lwd=3)
  legend("topright",c("PTRMS","WAS","TOGA","iWAS"), col=c("black","red","green","purple"), pch=19)
}
plotslopes5hz = function(pass){

  par(mfrow=c(4,2),mar = c(5, 5, 3, 5))
  plot(pass$CO2_7000_ppm, pass$CO_DACOM_DISKIN, xlab='CO2, ppm', ylab='CO, ppb')
  ll=lm(pass$CO_DACOM_DISKIN~pass$CO2_7000_ppm)
  text(x = median(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE), y=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)*0.5,labels = round(ll$coefficients[2])/1E3)
  plot(pass$CO2_7000_ppm, pass$CH4_DACOM_DISKIN, xlab='CO2, ppm', ylab='CH4, ppb')
  ind = which(is.finite(pass$Furan_NOAAPTR_ppbv_WARNEKE))
  if (length(ind) > 0){
    plot(pass$CO2_7000_ppm, pass$Furan_NOAAPTR_ppbv_WARNEKE, xlab='CO2, ppm', ylab='Furan, ppb')
  }
  plot(pass$CO2_7000_ppm, pass$SO2_LIF_ROLLINS/1E3, xlab='CO2, ppm', ylab='SO2, ppb')
  xx = pass$CO2_7000_ppm*1E3; yy= pass$SO2_LIF_ROLLINS/1E3
  ll=lm(yy~xx)
  text(x = median(xx/1E3, na.rm=TRUE), y=max(yy, na.rm=TRUE)*0.5,labels = round(ll$coefficients[2], digits=5))
  plot(pass$CO2_7000_ppm, pass$NO_LIF_ROLLINS,  xlab='CO2, ppm', ylab='NO, ppt')
  xx = pass$CO2_7000_ppm*1E3; yy= pass$NO_LIF_ROLLINS/1E3
  ll=lm(yy~xx)
  text(x = median(xx/1E3, na.rm=TRUE), y=max(yy, na.rm=TRUE)*0.5,labels = round(ll$coefficients[2], digits=5))
  plot(pass$CO2_7000_ppm, pass$CH2O_ISAF_HANISCO,  xlab='CO2, ppm', ylab='CH2O, ppt')
  plot(pass$CO2_7000_ppm, pass$PHENOL.1Hz_CIT_WENNBERG,  xlab='CO2, ppm', ylab='Phenol, ppt')
  plot(pass$CO2_7000_ppm, pass$PAN_GTCIMS_HUEY,  xlab='CO2, ppm', ylab='PAN, ppt')
 # plot(pass$CO2_7000_ppm, pass$PROPENE.HN.1Hz_CIT_WENNBERG,  xlab='CO2, ppm', ylab='PROPENE-HN, ppt')

}
plotpass1hzJUSTNH3 = function(pass){
  cz=1.2
  ind = which(is.finite(pass$CO2_7000_ppm_DISKIN))
  if (length(ind) > 3){
    cc1 = cor.test(pass$NH3_UIOPTR_ppbV_WISTHALER,pass$CO_DACOM_DISKIN)
    cc1 = cc1$estimate^2
  } else{
    cc1 = NaN
  }
  print(c("r^2, CO2 vs. Acrolein",cc1))
  
  par(mfrow=c(1,2),cex=1.5,mar = c(4, 4, 3, 4))
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",type='o',lwd=3,
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  text(mean(pass$Time_Start)-2,yy*0.996,paste("r_sq =",round(cc1, digits = 2) ))
  par(new = TRUE, cex=1.5, mar=c(4,4,3,4))
  plot(pass$Time_Start, pass$NH3_UIOPTR_ppbV_WISTHALER, xaxt = "n", yaxt = "n",type='o',lwd=3,pch=19,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(min(pass$NH3_UIOPTR_ppbV_WISTHALER, na.rm=TRUE),
                                                          max(pass$NH3_UIOPTR_ppbV_WISTHALER, na.rm=TRUE)))
  axis(side = 4, cex=cz)
  mtext("NH3, ppb",cex=1.5, side = 4, line = 3)
  legend("topright", c("CO", "NH3"),lwd=4,
         col = c("black", "red"), bty='n',lty = c(1, 2))
  
  plot(pass$CO_DACOM_DISKIN, pass$NH3_UIOPTR_ppbV_WISTHALER, 
       xlab="CO, ppb", ylab='Ammonia, ppb', pch=19)
  
}
plotpass5hzJUSTOA = function(pass){
  cz=1.2
  ind = which(is.finite(pass$CO_DACOM_DISKIN))
  if (length(ind) > 3){
    cc1 = cor.test(pass$OA_PM1_AMS_JIMENEZ,pass$CO_DACOM_DISKIN)
    cc1 = cc1$estimate^2
  } else{
    cc1 = NaN
  }
  print(c("r^2, CO vs. OA",cc1))
  
  par(mfrow=c(1,2),cex=1.5,mar = c(4, 4, 3, 4))
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",type='o',lwd=3,
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  text(mean(pass$Time_Start)-2,yy*0.996,paste("r_sq =",round(cc1, digits = 2) ))
  par(new = TRUE, cex=1.5, mar=c(4,4,3,4))
  plot(pass$Time_Start, pass$OA_PM1_AMS_JIMENEZ, xaxt = "n", yaxt = "n",type='o',lwd=3,pch=19,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(min(pass$OA_PM1_AMS_JIMENEZ, na.rm=TRUE),
                                                          max(pass$OA_PM1_AMS_JIMENEZ, na.rm=TRUE)))
  axis(side = 4, cex=cz)
  mtext("OA, ug m-3",cex=1.5, side = 4, line = 3)
  legend("topright", c("CO", "OA"),lwd=4,
         col = c("black", "red"), bty='n',lty = c(1, 2))
  
  plot(pass$CO_DACOM_DISKIN, pass$OA_PM1_AMS_JIMENEZ, 
       xlab="CO, ppb", ylab='Ammonia, ppb', pch=19)
  
}
plotpass5hzJUSTCO = function(pass){
  cz=1.2
  ind = which(is.finite(pass$CO2_7000_ppm_DISKIN))
  if (length(ind) > 3){
    cc1 = cor.test(pass$CO2_7000_ppm_DISKIN,pass$CO_DACOM_DISKIN)
    cc1 = cc1$estimate^2
  } else{
    cc1 = NaN
  }
  print(c("r^2, CO vs. CO2",cc1))
  
  par(mfrow=c(1,2),cex=1.5,mar = c(4, 4, 3, 4))
  yy=max(pass$CO2_7000_ppm, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO2_7000_ppm, ylab = expression(paste('CO'[2],',', ' ppm')),type='o',lwd=3,
       main = "", xlab = "Time", ylim=c(min(pass$CO2_7000_ppm, na.rm=TRUE), yy))
  text(mean(pass$Time_Start),yy*0.996,paste("r_sq =",round(cc1, digits = 2) ))
  par(new = TRUE, cex=1.5, mar=c(4,4,3,4))
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, xaxt = "n", yaxt = "n",type='o',lwd=3,pch=19,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,max(pass$CO_DACOM_DISKIN, na.rm=TRUE)))
  axis(side = 4, cex=cz)
  mtext("CO, ppb",cex=1.5, side = 4, line = 3)
  legend("topright", c(expression(paste("CO"[2])), "CO"),lwd=4,
         col = c("black", "red"), bty='n',lty = c(1, 2))

  ind = which(is.finite(pass$Furan_NOAAPTR_ppbv_WARNEKE))
  if (length(ind) > 3){
    cc = cor.test(pass$CO_DACOM_DISKIN, pass$Furan_NOAAPTR_ppbv_WARNEKE)
    print(c("r^2, CO vs. Furan",cc$estimate^2))
  } else{  print(c("r^2, CO vs. Furan",NaN))}
  ind = which(is.finite(pass$CH2O_ISAF_HANISCO))
  if (length(ind) > 3){
    cc = cor.test(pass$CO_DACOM_DISKIN, pass$CH2O_ISAF_HANISCO)
    print(c("r^2, CO vs. CH2O",cc$estimate^2))
  } else{  print(c("r^2, CO vs. CH2O",NaN))}
  
  
  ind = which(is.finite(pass$NO_LIF_ROLLINS))
  if (length(ind) > 3){
    cc = cor.test(pass$CO_DACOM_DISKIN, pass$NO_LIF_ROLLINS)
    print(c("r^2, CO vs. NO",cc$estimate^2))
  } else{  print(c("r^2, CO vs. NO",NaN))}
  
  ind = which(is.finite(pass$SO2_LIF_ROLLINS))
  if (length(ind) > 3){
    cc = cor.test(pass$CO_DACOM_DISKIN, pass$SO2_LIF_ROLLINS)
    print(c("r^2, CO vs. SO2",cc$estimate^2))
  } else{  print(c("r^2, CO vs. SO2",NaN))}
  
  ind = which(is.finite(pass$PHENOL.1Hz_CIT_WENNBERG))
  if (length(ind) > 3){
    cc = cor.test(pass$CO_DACOM_DISKIN, pass$PHENOL.1Hz_CIT_WENNBERG)
    print(c("r^2, CO vs. PHENOL",cc$estimate^2))
  } else{  print(c("r^2, CO vs. PHENOL",NaN))}
  
  ind = which(is.finite(pass$PAN_GTCIMS_HUEY) & is.finite(pass$CO_DACOM_DISKIN))
  if (length(ind) > 3){
    cc = cor.test(pass$CO_DACOM_DISKIN, pass$PAN_GTCIMS_HUEY)
    print(c("r^2, CO vs. PAN",cc$estimate^2))
  } else{  print(c("r^2, CO vs. PAN",NaN))}
  
  plot(pass$CO2_7000_ppm, pass$CO_DACOM, 
       xlab=expression(paste('CO'[2],',', ' ppm')), ylab='CO, ppb', pch=19)
  

}
plotpass5hzJUSTCN = function(pass){
  cz=1.2
  ind = which(is.finite(pass$CO_DACOM_DISKIN))
  if (length(ind) > 3){
    cc1 = cor.test(pass$FastCNgt6nm_stdPT_MOORE,pass$CO_DACOM_DISKIN)
    cc1 = cc1$estimate^2
  } else{
    cc1 = NaN
  }
  print(c("r^2, CN vs. CO",cc1))
  
  par(mfrow=c(1,2),cex=1.5,mar = c(4, 4, 3, 4))
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = expression(paste('CO, ppb')),type='o',lwd=3,
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  text(mean(pass$Time_Start),yy*0.66,paste("r_sq =",round(cc1, digits = 2) ))
  par(new = TRUE, cex=1.5, mar=c(4,4,3,4))
  plot(pass$Time_Start, pass$FastCNgt6nm_stdPT_MOORE, xaxt = "n", yaxt = "n",type='o',lwd=3,pch=19,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,max(pass$FastCNgt6nm_stdPT_MOORE, na.rm=TRUE)))
  axis(side = 4, cex=cz)
  mtext("CNgt6nm",cex=1.5, side = 4, line = 3)
  legend("topright", c("CO", "CNgt6nm"),lwd=4,
         col = c("black", "red"), bty='n',lty = c(1, 2))
  
  plot(pass$CO_DACOM_DISKIN, pass$FastCNgt6nm_stdPT_MOORE, 
       xlab='CO,ppb', ylab='CNgt6nm', pch=19)
  
  
}
#START MAKE BACKGROUND0.
plotpass5hzJUSTCOHNO2 = function(pass){
  cz=1.2
  ind = which(is.finite(pass$CO2_7000_ppm_DISKIN))
  if (length(ind) > 3){
    cc1 = cor.test(pass$CO2_7000_ppm_DISKIN,pass$CO_DACOM_DISKIN)
    cc1 = cc1$estimate^2
  } else{
    cc1 = NaN
  }
  print(c("r^2, CO vs. CO2",cc1))
  
  par(mfrow=c(1,2),cex=1.5,mar = c(4, 4, 3, 4))
  yy=max(pass$CO2_7000_ppm, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO2_7000_ppm, ylab = expression(paste('CO'[2],',', ' ppm')),type='o',lwd=3,
       main = "", xlab = "Time", ylim=c(min(pass$CO2_7000_ppm, na.rm=TRUE), yy))
  text(mean(pass$Time_Start),yy*0.56,paste("r_sq =",round(cc1, digits = 2) ))
  par(new = TRUE, cex=1.5, mar=c(4,4,3,4))
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, xaxt = "n", yaxt = "n",type='o',lwd=3,pch=19,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,max(pass$CO_DACOM_DISKIN, na.rm=TRUE)))
  axis(side = 4, cex=cz)
  mtext("CO, ppb",cex=1.5, side = 4, line = 3)
  legend("topright", c(expression(paste("CO"[2])), "CO"),lwd=4,
         col = c("black", "red"), bty='n',lty = c(1, 2))
  
 
  ind = which(is.finite(pass$PAN_GTCIMS_HUEY) & is.finite(pass$CO_DACOM_DISKIN))
  if (length(ind) > 3){
    cc = cor.test(pass$CO_DACOM_DISKIN, pass$PAN_GTCIMS_HUEY)
    print(c("r^2, CO vs. PAN",cc$estimate^2))
  } else{  print(c("r^2, CO vs. PAN",NaN))}
  
  plot( pass$Time_Start,pass$HNO2_NOAACIMS_VERES, 
       ylab=expression(paste('HNO'[2],',', ' ppt')), xlab='Time', pch=19)
  
  
}
plotpass5hzSAVE = function(pass, id){
  cz = 1.
  doplot1 = 0
  if (doplot1 == 1){
     file=paste(id,pass$Time_Start[1],".jpeg",sep='')
     jpeg(file=file)
  }
  par(mfrow=c(4,2),mar = c(5, 5, 3, 5))
  yy=max(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO2_7000_ppm_DISKIN, ylab = "CO2, ppm",type='l',
       main = "", xlab = "Time", ylim=c(min(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, xaxt = "n", yaxt = "n",type='l',
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,max(pass$CO_DACOM_DISKIN, na.rm=TRUE)))
  axis(side = 4, cex=cz)
  mtext("CO, ppb", side = 4, line = 3)
  legend("topright", c("CO2", "CO"), col = c("black", "red"), bty='n',lty = c(1, 2))
  ind = length(pass$CO2_7000_ppm_DISKIN)
  if (length(ind) > 2){
    cc1 = cor.test(pass$CO2_7000_ppm_DISKIN,pass$CO_DACOM_DISKIN)
    cc1 = cc1$estimate^2
  } else{ cc1 = NaN}
  print(c("r^2, CO vs. CO2",cc1))
  text(max(pass$Time_Start)-7,yy*0.991,paste("r2=",round(cc1, digits = 2) ))
  
  ind = which(is.finite(pass$Furan_NOAAPTR_ppbv_WARNEKE))
  if (length(ind) > 3){
    cc = cor.test(pass$CO_DACOM_DISKIN, pass$Furan_NOAAPTR_ppbv_WARNEKE)
    print(c("r^2, CO vs. Furan",cc$estimate^2))
  } else{  print(c("r^2, CO vs. Furan",NaN))}
  
  ind = which(is.finite(pass$CH2O_ISAF_HANISCO))
  if (length(ind) > 3){
    cc = cor.test(pass$CO_DACOM_DISKIN, pass$CH2O_ISAF_HANISCO)
    print(c("r^2, CO vs. CH2O",cc$estimate^2))
  } else{  print(c("r^2, CO vs. CH2O",NaN))}
  
  
  ind = which(is.finite(pass$NO_LIF_ROLLINS))
  if (length(ind) > 3){
    cc = cor.test(pass$CO_DACOM_DISKIN, pass$NO_LIF_ROLLINS)
    print(c("r^2, CO vs. NO",cc$estimate^2))
  } else{  print(c("r^2, CO vs. NO",NaN))}
  
  ind = which(is.finite(pass$SO2_LIF_ROLLINS))
  if (length(ind) > 3){
    cc = cor.test(pass$CO_DACOM_DISKIN, pass$SO2_LIF_ROLLINS)
    print(c("r^2, CO vs. SO2",cc$estimate^2))
  } else{  print(c("r^2, CO vs. SO2",NaN))}
  
  ind = which(is.finite(pass$PHENOL.1Hz_CIT_WENNBERG))
  if (length(ind) > 3){
    cc = cor.test(pass$CO_DACOM_DISKIN, pass$PHENOL.1Hz_CIT_WENNBERG)
    print(c("r^2, CO vs. PHENOL",cc$estimate^2))
  } else{  print(c("r^2, CO vs. PHENOL",NaN))}
  
  ind = which(is.finite(pass$PAN_GTCIMS_HUEY) & is.finite(pass$CO_DACOM_DISKIN))
  if (length(ind) > 3){
    cc = cor.test(pass$CO_DACOM_DISKIN, pass$PAN_GTCIMS_HUEY)
    print(c("r^2, CO vs. PAN",cc$estimate^2))
  } else{  print(c("r^2, CO vs. PAN",NaN))}
  
  ind = which(is.finite(pass$PROPENE.HN.1Hz_CIT_WENNBERG))
  if (length(ind) > 3){
    cc = cor.test(pass$CO_DACOM_DISKIN, pass$PROPENE.HN.1Hz_CIT_WENNBERG)
    print(c("r^2, CO vs. PROPENE-HN",cc$estimate^2))
  } else{  print(c("r^2, CO vs. PROPENE-HN",NaN))}
  
  plot(pass$CO2_7000_ppm, pass$CO_DACOM, xlab='CO2, ppm', ylab='CO, ppb')
  
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN,ylab = "CO, ppb",type='l', lwd=3,
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yyz =max(pass$Furan_NOAAPTR_ppbv_WARNEKE, na.rm=TRUE)
  if (is.infinite(yyz)){ yyz = 10}
  plot(pass$Time_Start, pass$Furan_NOAAPTR_ppbv_WARNEKE,  xaxt = "n", yaxt = "n",type='l',
       cex=cz,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yyz))
  axis(side = 4, cex=cz)
  mtext("Furan, ppt", side = 4, line = 3)
  legend("topright", c("CO", "Furan"), col = c("black", "red"), bty='n',lty = c(1, 2))
  
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",type='l', lwd=3,
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yyz = max(pass$SO2_LIF_ROLLINS, na.rm=TRUE)
  if (is.infinite(yyz)){ yyz = 10}
  plot(pass$Time_Start, pass$SO2_LIF_ROLLINS,  xaxt = "n", yaxt = "n",type='l', cex=cz,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yyz))
  axis(side = 4, cex=cz)
  mtext("SO2, ppt", side = 4, line = 3)
  legend("topright", c("CO", "SO2"), col = c("black", "red"), bty='n',lty = c(1, 2))
  
  yyz = max(pass$NO_LIF_ROLLINS, na.rm=TRUE)
  if (is.infinite(yyz)){ yyz = 10}
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",type='l',lwd=3,
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$NO_LIF_ROLLINS,  xaxt = "n", yaxt = "n",type='l', cex=cz,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yyz))
  axis(side = 4, cex=cz)
  mtext("NO, ppt", side = 4, line = 3)
  legend("topright", c("CO", "NO"), col = c("black", "red"), bty='n',lty = c(1, 2))
  
  yyz =max(pass$CH2O_ISAF_HANISCO, na.rm=TRUE)
  if (is.infinite(yyz)){ yyz = 10}
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",type='l',lwd=3,
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$CH2O_ISAF_HANISCO,  xaxt = "n", yaxt = "n",type='l', cex=cz,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yyz))
  axis(side = 4, cex=cz)
  mtext("HCHO, ppt", side = 4, line = 3)
  legend("topright", c("CO", "HCHO"), col = c("black", "red"), bty='n',lty = c(1, 2))
  
  yyz =max(pass$PHENOL.1Hz_CIT_WENNBERG, na.rm=TRUE)
  if (is.infinite(yyz)){ yyz = 10}
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",lwd=3
       ,type='l', xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$PHENOL.1Hz_CIT_WENNBERG, xaxt = "n", yaxt = "n",type='l', cex=cz,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yyz))
  axis(side = 4, cex=cz)
  mtext("Phenol, ppt", side = 4, line = 3)
  legend("topright", c("CO", "Phenol"), col = c("black", "red"), bty='n',lty = c(1, 2))
  
  yyz =max(pass$PAN_GTCIMS_HUEY, na.rm=TRUE)
  if (is.infinite(yyz)){ yyz = 10}
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab = "CO, ppb",type='l',lwd=3,
       xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$PAN_GTCIMS_HUEY,  xaxt = "n", yaxt = "n",type='l', cex=cz,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yyz))
  lines(pass$Time_Start, pass$PROPENE.HN.1Hz_CIT_WENNBERG*100, col='blue')
  axis(side = 4, cex=cz)
  mtext("PAN, ppt", side = 4, line = 3)
  legend("topright", c("CO", "PAN",'PROPENE-HN'), col = c("black", "red","blue"), bty='n',lty = c(1, 2,1))
 dev.off() 
}
plotpass1hz = function(pass){
  cz=1
  par(mfrow=c(4,2),mar = c(5, 5, 3, 5))

  yy=max(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO2_7000_ppm_DISKIN,type ="o", ylab = "CO2, ppm",
       main = "", xlab = "Time", ylim=c(min(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, type = "o", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,max(pass$CO_DACOM_DISKIN, na.rm=TRUE)))
  axis(side = 4, cex=cz)
  mtext("CO, ppb", side = 4, line = 3)
  legend("topright", c("CO2", "CO"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc = cor.test(pass$CO2_7000_ppm_DISKIN,pass$CO_DACOM_DISKIN, use='pairwise.complete.obs')
  print(cc$estimate^2)
  text(max(pass$Time_Start)-5,yy*0.9,paste("r2=",round(cc$estimate^2, digits = 2) ))
  
  
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  if (is.nan(yy)){yy=1}
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN,type ="o", ylab = "CO, ppb",
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yy=max(pass$C6H10O5_JIMENENZ, na.rm=TRUE);  if (is.infinite(yy)){yy=1}
  plot(pass$Time_Start, pass$C6H10O5_JIMENENZ, type = "o", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yy))
  axis(side = 4, cex=cz)
  mtext("C6H10O5, ug m-3", side = 4, line = 3)
  legend("topright", c("CO", "C6H10O5_JIMENENZ"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc = cor(pass$CO_DACOM_DISKIN,pass$C6H10O5_JIMENENZ, use='pairwise.complete.obs')
  print(cc^2)
  text(median(pass$Time_Start),yy*0.95,paste("r2=",round(cc^2, digits = 2) ))
  
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN,type ="o", ylab = "CO, ppb",
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yy=max(pass$BC_mass_90_550_nm_SCHWARZ, na.rm=TRUE)
  if (is.nan(yy) | is.infinite(yy)){yy=1}
  plot(pass$Time_Start, pass$BC_mass_90_550_nm_SCHWARZ, type = "o", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yy))
  axis(side = 4, cex=cz)
  mtext("BC, ng m-3", side = 4, line = 3)
  legend("topright", c("CO", "BC"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc = cor(pass$CO_DACOM_DISKIN,pass$BC_mass_90_550_nm_SCHWARZ, use='pairwise.complete.obs')
  print(cc^2)
  text(median(pass$Time_Start),yy*0.9,paste("r2=",round(cc^2, digits = 2) ))
  
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN,type ="o", ylab = "CO, ppb",
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yy=max(pass$HNO2_NOAACIMS_VERES, na.rm=TRUE);  if (is.infinite(yy)){yy=1}
  
  plot(pass$Time_Start, pass$HNO2_NOAACIMS_VERES, type = "o", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yy))
  axis(side = 4, cex=cz)
  mtext("HONO, ppt", side = 4, line = 3)
  legend("topright", c("CO", "HONO VERES"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc = cor(pass$CO_DACOM_DISKIN,pass$HNO2_NOAACIMS_VERES, use='pairwise.complete.obs')
  print(cc^2)
  text(median(pass$Time_Start),yy*0.9,paste("r2=",round(cc^2, digits = 2) ))
  
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN,type ="o", ylab = "CO, ppb",
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yy=max(pass$CH3COCHO_ACES_WOMACK, na.rm=TRUE)
  if (is.nan(yy) | is.infinite(yy)){ yy=1}
  plot(pass$Time_Start, pass$CH3COCHO_ACES_WOMACK, type = "o", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yy))
  axis(side = 4, cex=cz)
  mtext("CH3COCHO, ppt?", side = 4, line = 3)
  legend("topright", c("CO", "CH3COCHO WOMACK"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc = cor(pass$CO_DACOM_DISKIN,pass$CH3COCHO_ACES_WOMACK, use='pairwise.complete.obs')
  print(cc^2)
  text(median(pass$Time_Start),yy*0.9,paste("r2=",round(cc^2, digits = 2) ))
  
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN,type ="o", ylab = "CO, ppb",
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yy=max(pass$NO_CL_RYERSON, na.rm=TRUE)
  if (is.nan(yy) | is.infinite(yy)){yy=1}
  plot(pass$Time_Start, pass$NO_CL_RYERSON, type = "o", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yy))
  axis(side = 4, cex=cz)
  mtext("NO, ppt?", side = 4, line = 3)
  legend("topright", c("CO", "NO RYERSON"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc = cor(pass$CO_DACOM_DISKIN,pass$NO_CL_RYERSON)
  print(cc^2)
  text(median(pass$Time_Start),yy*0.9,paste("r2=",round(cc^2, digits = 2) ))
  
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN,type ="o", ylab = "CO, ppb",
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yy=max(pass$C2H6_CAMS_pptv_FRIED, na.rm=TRUE)
  if (is.nan(yy) | is.infinite(yy)){yy=1}
  plot(pass$Time_Start, pass$C2H6_CAMS_pptv_FRIED, type = "o", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yy))
  axis(side = 4, cex=cz)
  mtext("C2H6, ppt", side = 4, line = 3)
  legend("topright", c("CO", "CAMS C2H6"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc = cor(pass$CO_DACOM_DISKIN,pass$C2H6_CAMS_pptv_FRIED)
  print(cc^2)
  text(median(pass$Time_Start),yy*0.9,paste("r2=",round(cc^2, digits = 2) ))
  
  yy=max(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN,type ="o", ylab = "CO, ppb",
       main = "", xlab = "Time", ylim=c(min(pass$CO_DACOM_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yy=max(pass$NO2_CL_RYERSON, na.rm=TRUE)
  if (is.nan(yy) | is.infinite(yy)){yy=1}
  plot(pass$Time_Start, pass$NO2_CL_RYERSON, type = "o", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yy))
  axis(side = 4, cex=cz)
  mtext("NO2 Ryerson, ppb", side = 4, line = 3)
  legend("topright", c("CO", "NO2"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc = cor(pass$CO_DACOM_DISKIN,pass$NO2_CL_RYERSON)
  print(cc^2)
  text(median(pass$Time_Start),yy*0.9,paste("r2=",round(cc^2, digits = 2) ))
  
}
plotpass1hzBLAKE = function(pass,blake, blakeBG){
  cz=1
  par(mfrow=c(1,1),mar = c(5, 5, 3, 5))
  
  yy=max(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO2_7000_ppm_DISKIN,type ="o", ylab = "CO2, ppm",
       main = "", xlab = "Time", ylim=c(min(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE), yy))
  points(blake$Time_Start, blake$CO2_7000_ppm_DISKIN_BLAKE, col='red')
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, type = "o", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,max(pass$CO_DACOM_DISKIN, na.rm=TRUE)))
  axis(side = 4, cex=cz)
  mtext("CO, ppb", side = 4, line = 3)
  legend("topright", c("CO2", "CO"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc = cor.test(pass$CO2_7000_ppm_DISKIN,pass$CO_DACOM_DISKIN, use='pairwise.complete.obs')
  print(cc$estimate^2)
  text(max(pass$Time_Start)-5,yy*0.9,paste("r2=",round(cc$estimate^2, digits = 2) ))
  
  
  yy=max(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE)
  if (is.nan(yy)){yy=1}
  plot(pass$Time_Start, pass$CO2_7000_ppm_DISKIN,type ="o", ylab = "CO2, ppm",
       main = "", xlab = "Time", ylim=c(min(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yy=max(pass$OA_PM1_AMS_JIMENEZ, na.rm=TRUE);  if (is.infinite(yy)){yy=1}
  plot(pass$Time_Start, pass$OA_PM1_AMS_JIMENEZ, type = "o", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yy))
  axis(side = 4, cex=cz)
  mtext("OA, ug m-3", side = 4, line = 3)
  legend("topright", c("CO2", "OA"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc = cor(pass$CO2_7000_ppm_DISKIN,pass$OA_PM1_AMS_JIMENEZ, use='pairwise.complete.obs')
  print(cc^2)
  text(median(pass$Time_Start),yy*0.95,paste("r2=",round(cc^2, digits = 2) ))
  
  yy=max(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO2_7000_ppm_DISKIN,type ="o", ylab = "CO2, ppm",
       main = "", xlab = "Time", ylim=c(min(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yy=max(pass$BC_mass_90_550_nm_SCHWARZ, na.rm=TRUE)
  if (is.nan(yy) | is.infinite(yy)){yy=1}
  plot(pass$Time_Start, pass$BC_mass_90_550_nm_SCHWARZ, type = "o", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yy))
  axis(side = 4, cex=cz)
  mtext("BC, ng m-3", side = 4, line = 3)
  legend("topright", c("CO2", "BC"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc = cor(pass$CO2_7000_ppm_DISKIN,pass$BC_mass_90_550_nm_SCHWARZ, use='pairwise.complete.obs')
  print(cc^2)
  text(median(pass$Time_Start),yy*0.9,paste("r2=",round(cc^2, digits = 2) ))
  
  yy=max(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO2_7000_ppm_DISKIN,type ="o", ylab = "CO2, ppm",
       main = "", xlab = "Time", ylim=c(min(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yy=max(pass$HNO2_NOAACIMS_VERES, na.rm=TRUE);  if (is.infinite(yy)){yy=1}
  
  plot(pass$Time_Start, pass$HNO2_NOAACIMS_VERES, type = "o", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yy))
  axis(side = 4, cex=cz)
  mtext("HONO, ppt", side = 4, line = 3)
  legend("topright", c("CO2", "HONO VERES"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc = cor(pass$CO2_7000_ppm_DISKIN,pass$HNO2_NOAACIMS_VERES, use='pairwise.complete.obs')
  print(cc^2)
  text(median(pass$Time_Start),yy*0.9,paste("r2=",round(cc^2, digits = 2) ))
  
  yy=max(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO2_7000_ppm_DISKIN,type ="o", ylab = "CO2, ppm",
       main = "", xlab = "Time", ylim=c(min(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yy=max(pass$CH3COCHO_ACES_WOMACK, na.rm=TRUE)
  if (is.nan(yy) | is.infinite(yy)){ yy=1}
  plot(pass$Time_Start, pass$CH3COCHO_ACES_WOMACK, type = "o", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yy))
  axis(side = 4, cex=cz)
  mtext("CH3COCHO, ppt?", side = 4, line = 3)
  legend("topright", c("CO2", "CH3COCHO WOMACK"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc = cor(pass$CO2_7000_ppm_DISKIN,pass$CH3COCHO_ACES_WOMACK, use='pairwise.complete.obs')
  print(cc^2)
  text(median(pass$Time_Start),yy*0.9,paste("r2=",round(cc^2, digits = 2) ))
  
  yy=max(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO2_7000_ppm_DISKIN,type ="o", ylab = "CO2, ppm",
       main = "", xlab = "Time", ylim=c(min(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yy=max(pass$NO_CL_RYERSON, na.rm=TRUE)
  if (is.nan(yy) | is.infinite(yy)){yy=1}
  plot(pass$Time_Start, pass$NO_CL_RYERSON, type = "o", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yy))
  axis(side = 4, cex=cz)
  mtext("NO, ppt?", side = 4, line = 3)
  legend("topright", c("CO2", "NO RYERSON"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc = cor(pass$CO2_7000_ppm_DISKIN,pass$NO_CL_RYERSON)
  print(cc^2)
  text(median(pass$Time_Start),yy*0.9,paste("r2=",round(cc^2, digits = 2) ))
  
  yy=max(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO2_7000_ppm_DISKIN,type ="o", ylab = "CO2, ppm",
       main = "", xlab = "Time", ylim=c(min(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yy=max(pass$C2H6_CAMS_pptv_FRIED, na.rm=TRUE)
  if (is.nan(yy) | is.infinite(yy)){yy=1}
  plot(pass$Time_Start, pass$C2H6_CAMS_pptv_FRIED, type = "o", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yy))
  axis(side = 4, cex=cz)
  mtext("C2H6, ppt", side = 4, line = 3)
  legend("topright", c("CO2", "CAMS C2H6"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc = cor(pass$CO2_7000_ppm_DISKIN,pass$C2H6_CAMS_pptv_FRIED)
  print(cc^2)
  text(median(pass$Time_Start),yy*0.9,paste("r2=",round(cc^2, digits = 2) ))
  
  yy=max(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO2_7000_ppm_DISKIN,type ="o", ylab = "CO2, ppm",
       main = "", xlab = "Time", ylim=c(min(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yy=max(pass$NO2_CL_RYERSON, na.rm=TRUE)
  if (is.nan(yy) | is.infinite(yy)){yy=1}
  plot(pass$Time_Start, pass$NO2_CL_RYERSON, type = "o", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yy))
  axis(side = 4, cex=cz)
  mtext("NO2 Ryerson, ppb", side = 4, line = 3)
  legend("topright", c("CO2", "NO2"), col = c("black", "red"), bty='n',lty = c(1, 2))
  cc = cor(pass$CO2_7000_ppm_DISKIN,pass$NO2_CL_RYERSON)
  print(cc^2)
  text(median(pass$Time_Start),yy*0.9,paste("r2=",round(cc^2, digits = 2) ))
  
}
plotpass1hzHNO2 = function(pass){
  par(mfrow=c(1,1),mar = c(5, 5, 3, 5))
  cz=1
  yy=max(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE)
  plot(pass$Time_Start, pass$CO2_7000_ppm_DISKIN,type ="o", ylab = "CO2, ppm",
       main = "", xlab = "Time", ylim=c(min(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  yy2=max(pass$HNO2_NOAACIMS_VERES, na.rm=TRUE);  if (is.infinite(yy2)){yy2=1}
  yy1=max(pass$HNO2_ACES_WOMACK*1E3, na.rm=TRUE);  if (is.infinite(yy1)){yy1=1}
  yy = max(yy1, yy2, na.rm=TRUE)
  
  plot(pass$Time_Start, pass$HNO2_NOAACIMS_VERES, type = "o", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,yy))
  lines(pass$Time_Start, pass$HNO2_ACES_WOMACK*1E3, col='blue')
  axis(side = 4, cex=cz)
  mtext("HONO, ppt", side = 4, line = 3)
  legend("topright", c("CO2", "HONO VERES","HONO WOMACK"), col = c("black", "red","blue"), bty='n',lty = c(1, 2,2))
  cc = cor(pass$CO2_7000_ppm_DISKIN,pass$HNO2_NOAACIMS_VERES, use='pairwise.complete.obs')
  cc2 = cor(pass$CO2_7000_ppm_DISKIN,pass$HNO2_ACES_WOMACK*1E3, use='pairwise.complete.obs')
  print(cc^2)
  text(median(pass$Time_Start),yy*0.9,paste("r2=",round(cc^2, digits = 2) ), col='red')
  text(median(pass$Time_Start),yy*0.85,paste("r2=",round(cc2^2, digits = 2) ), col='blue')
  
}

corSLOW = function(pass){
  print(c('RYERSON NO',cor(pass$NO_CL_RYERSON, pass$CO_DACOM_DISKIN) ))
  print(c('JIMENEZ OA', cor(pass$OA_PM1_AMS_JIMENEZ, pass$CO_DACOM_DISKIN)))
  print(c('SCHWARZ BC', cor(pass$BC_mass_90_550_nm_SCHWARZ, pass$CO_DACOM_DISKIN)))
  print(c('FRIED CH2O', cor(pass$CH2O_CAMS_pptv_FRIED, pass$CO_DACOM_DISKIN)))
  print(c('WOMACK CHOCHO',cor(pass$CHOCHO_ACES_WOMACK, pass$CO_DACOM_DISKIN)))
  print(c('STCLAIR NO', cor(pass$NO2_CANOE_STCLAIR, pass$CO_DACOM_DISKIN)))
  print(c('VERES HONO', cor(pass$HNO2_NOAACIMS_VERES, pass$CO_DACOM_DISKIN)))
  print(c('WISTHALER NH3', cor(pass$NH3_UIOPTR_ppbV_WISTHALER, pass$CO_DACOM_DISKIN)))
  print(c('BLAKE FURAN', cor(pass$Furan_WAS_BLAKE, pass$CO_DACOM_DISKIN)))
  print(c('APEL CH3CHO', cor(pass$CH3CHO_TOGA_APEL, pass$CO_DACOM_DISKIN)))
  print(c('GILMAN FURAN', cor(pass$Furan_NOAAiWAS_GILMAN, pass$CO_DACOM_DISKIN)))
}
# ----- ********* Get number of plumes and correlation with MCE -------
getplumesANDmce = function(firedata, origdata){
  # what is the variability across a fire?
  firedata$COUNT_EFFINAL = NaN
  firedata$corMCE_FINAL = NaN
  firedata$corMCE_ERCO = NaN
  
  for (i in 1:length(firedata$Group.1)){
    
    ind = which(origdata$variable == firedata$Group.1[i] & is.finite(origdata$FinalERtoCO) &
                  is.finite(origdata$FinalEF) & origdata$USEME != 0)
    # remove the impact of outliers
    tmpDat = origdata[ind,]
    tmpMEAN = mean(tmpDat$FinalEF, na.rm=TRUE); tmpSD = sd(tmpDat$FinalEF, na.rm=TRUE)
    ind = which(tmpDat$FinalEF < (tmpMEAN+tmpSD*2)) 
    tmpDat = tmpDat[ind,]
    
    #firedata$COUNT_EF5hz[i] = length(ind)
    if (length(tmpDat$FinalEF) > 3){
      #tt = cor.test(tmpDat$FinalERtoCO, tmpDat$MCE)
      #if (tt$p.value < 0.05){ firedata$corMCE_ERCO[i] = tt$estimate }
      tt = cor.test(tmpDat$FinalEF, tmpDat$MCE)
      if (tt$p.value < 0.05){ firedata$corMCE_FINAL[i] = tt$estimate }
      firedata$COUNT_EFFINAL[i] = length(tmpDat$FinalEF)
    }
  }
  return(firedata)
}
getplumesANDmcebyfuel = function(firedatafuel, origdata){
  # what is the variability across a fire?
  firedatafuel$COUNT_EFFINAL = NaN
  firedatafuel$corMCE_FINAL = NaN
  firedatafuel$corMCE_ERCO = NaN
  firedatafuel$slopeMCE = NaN; firedatafuel$interceptMCE = NaN ; firedatafuel$slopeError = NaN
  for (i in 1:length(firedatafuel$Group.1)){
    ind = which(origdata$variable == firedatafuel$Group.2[i] &
                  origdata$fuelORIG == firedatafuel$Group.1[i] & is.finite(origdata$FinalEF) &
                  origdata$FinalEF != 0 & is.finite(origdata$FinalERtoCO))
    firedatafuel$COUNT_EFFINAL[i] = length(ind)
    if (length(ind) > 2){
      tt = cor.test(origdata$FinalEF[ind], origdata$MCE[ind])
      if (!is.na(tt$p.value) & !is.nan(tt$p.value)){ 
        if (tt$p.value < 0.05){
          firedatafuel$corMCE_FINAL[i] = tt$estimate
          tS = lm(origdata$FinalEF[ind] ~ origdata$MCE[ind])
          tSS = summary(tS)
          firedatafuel$slopeMCE[i] = tS$coefficients[1]; firedatafuel$interceptMCE[i] = tS$coefficients[2]
          firedatafuel$slopeError[i] = tSS$coefficients[2,2]
          
        }
      }
      
      tt = cor.test(origdata$FinalERtoCO[ind], origdata$MCE[ind])
      if (!is.na(tt$p.value) & !is.nan(tt$p.value)){ 
        if (tt$p.value < 0.05){ firedatafuel$corMCE_ERCO[i] = tt$estimate }
      }
    }
    
  }
  return(firedatafuel)
}
getplumesANDmcebyfuel1VAR = function(firedatafuel, origdata){
  # what is the variability across a fire?
  firedatafuel$COUNT_EFFINAL = NaN
  firedatafuel$corMCE_FINAL = NaN
  firedatafuel$corMCE_ERCO = NaN
  for (i in 1:length(firedatafuel$Group.1)){
    ind = which(origdata$fuelORIG == firedatafuel$Group.1[i] & is.finite(origdata$FinalEF) &
                  origdata$FinalEF != 0 & is.finite(origdata$FinalERtoCO))
    firedatafuel$COUNT_EFFINAL[i] = length(ind)
    if (length(ind) > 2){
      tt = cor.test(origdata$FinalEF[ind], origdata$MCE[ind])
      if (tt$p.value < 0.05){ firedatafuel$corMCE_FINAL[i] = tt$estimate }
      tt = cor.test(origdata$FinalERtoCO[ind], origdata$MCE[ind])
      if (is.na(tt$p.value) | is.nan(tt$p.value)){tt$p.value = 0}
      if (tt$p.value < 0.05){ firedatafuel$corMCE_ERCO[i] = tt$estimate }
      
    }
    
  }
  return(firedatafuel)
}

# ------- *********** Time alignment functions--------
time_alignBW = function(co2.5hz, shift, co.ch4.5hz, shift2,
                      warneke.5hz, shift3, isaf.5hz, shift4,  
                      rollins.5hz, shift5,rollinsS.5hz, shift6,
                      cit.5hz, shift7, gtcims.5hz, shift8,
                      met.5hz, mce.5hz){
  print(c('1'))
  co2.5hz.shift = co2.5hz
  co2.5hz.shift$Time_Start = co2.5hz.shift$Time_Start + shift
  print(c('2'))
  co.ch4.5hz.shift = co.ch4.5hz
  co.ch4.5hz$Time_Start = co.ch4.5hz$Time_Start + shift2
  print(c('3'))
  warneke.5hz.shift = warneke.5hz
  warneke.5hz.shift$Time_Start = warneke.5hz.shift$Time_Start+ shift3
  print(c('4'))
  isaf.5hz.shift = isaf.5hz
  isaf.5hz.shift$Time_Start = isaf.5hz.shift$Time_Start+ shift4
  print(c('5'))
  rollins.5hz.shift = rollins.5hz
  rollins.5hz.shift$Time_Start = rollins.5hz.shift$Time_Start+ shift5
  rollinsS.5hz.shift = rollinsS.5hz
  rollinsS.5hz.shift$Time_Start = rollinsS.5hz.shift$Time_Start+ shift6
  print(c('7'))
  cit.5hz.shift = cit.5hz
  cit.5hz.shift$Time_Start = cit.5hz.shift$Time_Start+ shift7
  print(c('6'))
  gtcims.5hz.shift = gtcims.5hz
  gtcims.5hz.shift$Time_Start = gtcims.5hz.shift$Time_Start+ shift8
  
  tmp = merge(co2.5hz.shift, co.ch4.5hz.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, warneke.5hz.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, isaf.5hz.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, rollins.5hz.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, rollinsS.5hz.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, cit.5hz.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, gtcims.5hz.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, met.5hz, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, mce.5hz, by='Time_Start', all.x=TRUE)
  
  ind = which(is.finite(tmp$Time_Start) )#& is.finite(tmp$CO2_7000_ppm))
  tmp= tmp[ind,]
  
  return(tmp)
}
time_align = function(start,stop,co2,  co.ch4,  warneke, isaf, 
                      rollins, rollinsS,cit, gtcims, 
                      moore, jimenez, met){
  
  start = start - 1; stop = stop + 1 # to allow space for shifts
  ind1 = which(co2$Time_Start >= start & co2$Time_Start <= stop )
  ind2 = which(co.ch4$Time_Start >= start & co.ch4$Time_Start <= stop )
  ind3 = which(warneke$Time_Start >= start & warneke$Time_Start <= stop )
  ind4 = which(isaf$Time_Start >= start & isaf$Time_Start <= stop )
  ind5 = which(rollins$Time_Start >= start & rollins$Time_Start <= stop )
  ind6 = which(rollinsS$Time_Start >= start & rollinsS$Time_Start <= stop )
  ind7 = which(cit$Time_Start >= start & cit$Time_Start <= stop )
  ind8 = which(gtcims$Time_Start >= start & gtcims$Time_Start <= stop )
  ind9 = which(moore$Time_mid >= start & moore$Time_mid <= stop )
  ind10 = which(jimenez$Time_Start >= start & jimenez$Time_Start <= stop )
  
  if (length(ind9) == 0){  ind9 = which(moore$Time_Start >= start & moore$Time_Start <= stop )}
  ind10 = which(met$Time_Start >= start & met$Time_Start <= stop )
  
  debug=0
  shifts = seq(-1,1,0.2)
  # Aligning everything to CO
  # ---- ++++++++ CO2 find best shift for pass - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    co2.shift  = co2[ind1,]
    co2.shift  <- co2.shift[, !duplicated(colnames(co2.shift))]
    co2.shift$Time_Start = co2.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],co2.shift, by='Time_Start',all.x=TRUE)
    ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$CO2_7000_ppm_DISKIN)
  #  if (debug == 1){print(c(ccM$estimate,shifts[i]))}
    corsM = c(corsM,ccM$estimate )
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("co2 shift",shifts[indM]))}
  co2.shift = co2
  co2.shift$Time_Start = co2.shift$Time_Start+ shifts[indM]
  tmp = merge( co.ch4,co2.shift, by='Time_Start', all.x=TRUE)


  #co2.shift  = co2[ind1,]
  #co2.shift  <- co2.shift[, !duplicated(colnames(co2.shift))]
  #ind1B = which(co2.shift$CO2_7000_ppm_DISKIN == max(co2.shift$CO2_7000_ppm_DISKIN, na.rm=TRUE))
  #ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  #hifts = co2.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  #f (debug==1){print(c("co2 shift",shifts))}
  #o2.shift  = co2
  #co2.shift$Time_Start = co2.shift$Time_Start - shifts
  #tmp = merge(co.ch4,co2.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Warneke find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    warneke.shift  = warneke[ind3,]
    warneke.shift  <- warneke.shift[, !duplicated(colnames(warneke.shift))]
    warneke.shift$Time_Start = warneke.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],warneke.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$HCN_NOAAPTR_ppbv_WARNEKE))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$HCN_NOAAPTR_ppbv_WARNEKE)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("Warneke shift",shifts[indM]))}
  warneke.shift  = warneke
  warneke.shift$Time_Start = warneke.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, warneke.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ ISAF find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    isaf.shift  = isaf[ind4,]
    isaf.shift  <- isaf.shift[, !duplicated(colnames(isaf.shift))]
    isaf.shift$Time_Start = isaf.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],isaf.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$CH2O_ISAF_HANISCO))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$CH2O_ISAF_HANISCO)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("isaf shift",shifts[indM]))}
  isaf.shift  = isaf
  isaf.shift$Time_Start = isaf.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, isaf.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Rollins N find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    rollins.shift  = rollins[ind5,]
    rollins.shift  <- rollins.shift[, !duplicated(colnames(rollins.shift))]
    rollins.shift$Time_Start = rollins.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],rollins.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$NO_LIF_ROLLINS))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$NO_LIF_ROLLINS)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("Rollins shift",shifts[indM]))}
  rollins.shift  = rollins
  rollins.shift$Time_Start = rollins.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, rollins.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Rollins S find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    rollinsS.shift  = rollinsS[ind6,]
    rollinsS.shift  <- rollinsS.shift[, !duplicated(colnames(rollinsS.shift))]
    rollinsS.shift$Time_Start = rollinsS.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],rollinsS.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$SO2_LIF_ROLLINS))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$SO2_LIF_ROLLINS)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("Rollins S shift",shifts[indM]))}
  rollinsS.shift  = rollinsS
  rollinsS.shift$Time_Start = rollinsS.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, rollinsS.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++CIT find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    cit.shift  = cit[ind7,]
    cit.shift  <- cit.shift[, !duplicated(colnames(cit.shift))]
    cit.shift$Time_Start = cit.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],cit.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$HCN.1Hz_CIT_WENNBERG))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$HCN.1Hz_CIT_WENNBERG)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
    
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("CIT shift",shifts[indM]))}
  cit.shift  = cit
  cit.shift$Time_Start = cit.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, cit.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++GTCIMS find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    gtcims.shift  = gtcims[ind8,]
    gtcims.shift  <- gtcims.shift[, !duplicated(colnames(gtcims.shift))]
    gtcims.shift$Time_Start = gtcims.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],gtcims.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$APAN_GTCIMS_HUEY))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$APAN_GTCIMS_HUEY)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("Gtcims shift",shifts[indM]))}
  gtcims.shift  = gtcims
  gtcims.shift$Time_Start = gtcims.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, gtcims.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ MOORE - find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    moore.shift     = moore[ind9,]; 
    moore.shift$Time_Start = moore.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],moore.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$FastCNgt6nm_stdPT_MOORE)  & is.finite(tmp2$CO_DACOM_DISKIN))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$FastCNgt6nm_stdPT_MOORE)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
    if (debug == 1){print(c(ccM$estimate,shifts[i]))}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM,na.rm=TRUE))){indM=6}
  if (debug==1){print(c("Moore shift",shifts[indM]))}
  moore.shift     = moore; 
  moore.shift$Time_Start = moore.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, moore.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ JIMENEZ - find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    jimenez.shift     = jimenez[ind10,]; 
    jimenez.shift$Time_Start = jimenez.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],jimenez.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$OA_PM1_AMS_JIMENEZ)  & is.finite(tmp2$CO_DACOM_DISKIN))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$OA_PM1_AMS_JIMENEZ)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
    if (debug == 1){print(c(ccM$estimate,shifts[i]))}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM,na.rm=TRUE))){indM=6}
  if (debug==1){print(c("Jimenez shift",shifts[indM]))}
  jimenez.shift     = jimenez; 
  jimenez.shift$Time_Start = jimenez.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, jimenez.shift, by='Time_Start', all.x=TRUE)
  
  # ++++++++ merge with met
  tmp = merge(tmp, met, by='Time_Start', all.x=TRUE)

  ind = which(is.finite(tmp$Time_Start) )#& is.finite(tmp$CO2_7000_ppm))
  tmp= tmp[ind,]
  return(tmp)
}
time_align23 = function(start,stop,co2,  co.ch4,  warneke, isaf, 
                      rollins, rollinsS,cit, gtcims, 
                      moore,fastBC, met){
  
  start = start - 1; stop = stop + 1 # to allow space for shifts
  ind1 = which(co2$Time_Start >= start & co2$Time_Start <= stop )
  ind2 = which(co.ch4$Time_Start >= start & co.ch4$Time_Start <= stop )
  ind3 = which(warneke$Time_Start >= start & warneke$Time_Start <= stop )
  ind4 = which(isaf$Time_Start >= start & isaf$Time_Start <= stop )
  ind5 = which(rollins$Time_Start >= start & rollins$Time_Start <= stop )
  ind6 = which(rollinsS$Time_Start >= start & rollinsS$Time_Start <= stop )
  ind7 = which(cit$Time_Start >= start & cit$Time_Start <= stop )
  ind8 = which(gtcims$Time_Start >= start & gtcims$Time_Start <= stop )
  ind9 = which(moore$Time_mid >= start & moore$Time_mid <= stop )
  if (length(ind9) == 0){  ind9 = which(moore$Time_Start >= start & moore$Time_Start <= stop )}
  ind10 = which(met$Time_Start >= start & met$Time_Start <= stop )
  ind11 = which(fastBC$Time_Start >= start & fastBC$Time_Start <= stop )
  
  debug=1
  shifts = seq(-1,1,0.2)
  # Aligning everything to CO
  # ---- ++++++++ CO2 find best shift for pass - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    co2.shift  = co2[ind1,]
    co2.shift  <- co2.shift[, !duplicated(colnames(co2.shift))]
    co2.shift$Time_Start = co2.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],co2.shift, by='Time_Start',all.x=TRUE)
    ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$CO2_7000_ppm_DISKIN)
    #  if (debug == 1){print(c(ccM$estimate,shifts[i]))}
    corsM = c(corsM,ccM$estimate )
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("co2 shift",shifts[indM]))}
  co2.shift = co2
  co2.shift$Time_Start = co2.shift$Time_Start+ shifts[indM]
  tmp = merge( co.ch4,co2.shift, by='Time_Start', all.x=TRUE)
  
  
  #co2.shift  = co2[ind1,]
  #co2.shift  <- co2.shift[, !duplicated(colnames(co2.shift))]
  #ind1B = which(co2.shift$CO2_7000_ppm_DISKIN == max(co2.shift$CO2_7000_ppm_DISKIN, na.rm=TRUE))
  #ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  #hifts = co2.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  #f (debug==1){print(c("co2 shift",shifts))}
  #o2.shift  = co2
  #co2.shift$Time_Start = co2.shift$Time_Start - shifts
  #tmp = merge(co.ch4,co2.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Warneke find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    warneke.shift  = warneke[ind3,]
    warneke.shift  <- warneke.shift[, !duplicated(colnames(warneke.shift))]
    warneke.shift$Time_Start = warneke.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],warneke.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$HCN_NOAAPTR_ppbv_WARNEKE))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$HCN_NOAAPTR_ppbv_WARNEKE)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("Warneke shift",shifts[indM]))}
  warneke.shift  = warneke
  warneke.shift$Time_Start = warneke.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, warneke.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ ISAF find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    isaf.shift  = isaf[ind4,]
    isaf.shift  <- isaf.shift[, !duplicated(colnames(isaf.shift))]
    isaf.shift$Time_Start = isaf.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],isaf.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$CH2O_ISAF_HANISCO))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$CH2O_ISAF_HANISCO)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("isaf shift",shifts[indM]))}
  isaf.shift  = isaf
  isaf.shift$Time_Start = isaf.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, isaf.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Rollins N find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    rollins.shift  = rollins[ind5,]
    rollins.shift  <- rollins.shift[, !duplicated(colnames(rollins.shift))]
    rollins.shift$Time_Start = rollins.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],rollins.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$NO_LIF_ROLLINS))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$NO_LIF_ROLLINS)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("Rollins shift",shifts[indM]))}
  rollins.shift  = rollins
  rollins.shift$Time_Start = rollins.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, rollins.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Rollins S find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    rollinsS.shift  = rollinsS[ind6,]
    rollinsS.shift  <- rollinsS.shift[, !duplicated(colnames(rollinsS.shift))]
    rollinsS.shift$Time_Start = rollinsS.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],rollinsS.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$SO2_LIF_ROLLINS))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$SO2_LIF_ROLLINS)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("Rollins S shift",shifts[indM]))}
  rollinsS.shift  = rollinsS
  rollinsS.shift$Time_Start = rollinsS.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, rollinsS.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++CIT find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    cit.shift  = cit[ind7,]
    cit.shift  <- cit.shift[, !duplicated(colnames(cit.shift))]
    cit.shift$Time_Start = cit.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],cit.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$HCN.1Hz_CIT_WENNBERG))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$HCN.1Hz_CIT_WENNBERG)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
    
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("CIT shift",shifts[indM]))}
  cit.shift  = cit
  cit.shift$Time_Start = cit.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, cit.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++GTCIMS find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    gtcims.shift  = gtcims[ind8,]
    gtcims.shift  <- gtcims.shift[, !duplicated(colnames(gtcims.shift))]
    gtcims.shift$Time_Start = gtcims.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],gtcims.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$APAN_GTCIMS_HUEY))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$APAN_GTCIMS_HUEY)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("Gtcims shift",shifts[indM]))}
  gtcims.shift  = gtcims
  gtcims.shift$Time_Start = gtcims.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, gtcims.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ MOORE - find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    moore.shift     = moore[ind9,]; 
    moore.shift$Time_Start = moore.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],moore.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$FastCNgt6nm_stdPT_MOORE)  & is.finite(tmp2$CO_DACOM_DISKIN))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$FastCNgt6nm_stdPT_MOORE)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
    if (debug == 1){print(c(ccM$estimate,shifts[i]))}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM,na.rm=TRUE))){indM=6}
  if (debug==1){print(c("Moore shift",shifts[indM]))}
  moore.shift     = moore; 
  moore.shift$Time_Start = moore.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, moore.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Fast BC - find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    fastBC.shift     = fastBC[ind11,]; 
    fastBC.shift$Time_Start = fastBC.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],fastBC.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$bc_ng_kg_5Hz)  & is.finite(tmp2$bc_ng_kg_5Hz))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$bc_ng_kg_5Hz)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
    if (debug == 1){print(c(ccM$estimate,shifts[i]))}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM,na.rm=TRUE))){indM=6}
  if (debug==1){print(c("fastBC shift",shifts[indM]))}
  fastBC.shift     = fastBC; 
  fastBC.shift$Time_Start = fastBC.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, fastBC.shift, by='Time_Start', all.x=TRUE)
  # ++++++++ merge with met
  tmp = merge(tmp, met, by='Time_Start', all.x=TRUE)
  
  ind = which(is.finite(tmp$Time_Start) )#& is.finite(tmp$CO2_7000_ppm))
  tmp= tmp[ind,]
  return(tmp)
}
time_alignPEAKS = function(start,stop,co2,  co.ch4,  warneke, isaf, 
                      rollins, rollinsS,cit, gtcims, 
                      moore, met){
  
  start = start - 1; stop = stop + 1 # to allow space for shifts
  ind1 = which(co2$Time_Start >= start & co2$Time_Start <= stop )
  ind2 = which(co.ch4$Time_Start >= start & co.ch4$Time_Start <= stop )
  ind3 = which(warneke$Time_Start >= start & warneke$Time_Start <= stop )
  ind4 = which(isaf$Time_Start >= start & isaf$Time_Start <= stop )
  ind5 = which(rollins$Time_Start >= start & rollins$Time_Start <= stop )
  ind6 = which(rollinsS$Time_Start >= start & rollinsS$Time_Start <= stop )
  ind7 = which(cit$Time_Start >= start & cit$Time_Start <= stop )
  ind8 = which(gtcims$Time_Start >= start & gtcims$Time_Start <= stop )
  ind9 = which(moore$Time_mid >= start & moore$Time_mid <= stop )
  ind10 = which(met$Time_Start >= start & met$Time_Start <= stop )
  
  debug=1
  # ---- ++++++++ CO2 find best shift - ++++++++ -----
  
  co2.shift  = co2[ind1,]
  co2.shift  <- co2.shift[, !duplicated(colnames(co2.shift))]
  ind1B = which(co2.shift$CO2_7000_ppm_DISKIN == max(co2.shift$CO2_7000_ppm_DISKIN, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = co2.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("co2 shift",round(shifts,digits = 2)))}
  co2.shift  = co2
  co2.shift$Time_Start = co2.shift$Time_Start - round(shifts,digits = 2)
  tmp = merge(co.ch4,co2.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Warneke find best shift - ++++++++ -----
  warneke.shift  = warneke[ind3,]
  warneke.shift  <- warneke.shift[, !duplicated(colnames(warneke.shift))]
  ind1B = which(warneke.shift$Benzene_NOAAPTR_ppbv_WARNEKE== max(warneke.shift$Benzene_NOAAPTR_ppbv_WARNEKE, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = warneke.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (length(shifts) == 0){shifts = 0}
  if (debug==1){print(c("warneke shift",round(shifts,digits = 2)))}
  warneke.shift     = warneke[ind3,]; warneke.shift$Time_Start = warneke.shift$Time_Start- round(shifts,digits = 2)
  tmp = merge(tmp,warneke.shift,  by='Time_Start', all.x=TRUE)
  # ---- ++++++++ ISAF find best shift - ++++++++ -----
  isaf.shift  = isaf[ind4,]
  isaf.shift  <- isaf.shift[, !duplicated(colnames(isaf.shift))]
  ind1B = which(isaf.shift$CH2O_ISAF_HANISCO== max(isaf.shift$CH2O_ISAF_HANISCO, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = isaf.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("isaf shift",round(shifts,digits = 2)))}
  isaf.shift     = isaf[ind4,]; isaf.shift$Time_Start = isaf.shift$Time_Start - round(shifts,digits = 2)
  tmp = merge(tmp,isaf.shift,  by='Time_Start', all.x=TRUE)
  # ---- ++++++++ Rollins N find best shift - ++++++++ -----
  rollins.shift  = rollins[ind5,]
  rollins.shift  <- rollins.shift[, !duplicated(colnames(rollins.shift))]
  ind1B = which(rollins.shift$NO_LIF_ROLLINS== max(rollins.shift$NO_LIF_ROLLINS, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = rollins.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("rollins shift",round(shifts,digits = 2)))}
  rollins.shift     = rollins[ind5,]; rollins.shift$Time_Start = rollins.shift$Time_Start - round(shifts,digits = 2)
  tmp = merge(tmp,rollins.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Rollins S find best shift - ++++++++ -----
  rollinsS.shift  = rollinsS[ind6,]
  rollinsS.shift  <- rollinsS.shift[, !duplicated(colnames(rollinsS.shift))]
  ind1B = which(rollinsS.shift$SO2_LIF_ROLLINS== max(rollinsS.shift$SO2_LIF_ROLLINS, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = rollinsS.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("rollinsS shift",round(shifts,digits = 2)))}
  rollinsS.shift     = rollinsS[ind6,]; rollinsS.shift$Time_Start = rollinsS.shift$Time_Start - round(shifts,digits = 2)
  tmp = merge(tmp,rollinsS.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++CIT find best shift - ++++++++ -----
  cit.shift  = cit[ind7,]
  cit.shift  <- cit.shift[, !duplicated(colnames(cit.shift))]
  ind1B = which(cit.shift$PHENOL.1Hz_CIT_WENNBERG== max(cit.shift$PHENOL.1Hz_CIT_WENNBERG, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = cit.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (length(shifts) == 0){shifts = 0}
  if (debug==1){print(c("CIT shift",round(shifts,digits = 2)))}
  cit.shift     = cit[ind7,]; cit.shift$Time_Start = cit.shift$Time_Start - round(shifts,digits = 2)
  tmp = merge(tmp,cit.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++GTCIMS find best shift - ++++++++ -----
  gtcims.shift  = gtcims[ind8,]
  gtcims.shift  <- gtcims.shift[, !duplicated(colnames(gtcims.shift))]
  ind1B = which(gtcims.shift$APAN_GTCIMS_HUEY== max(gtcims.shift$APAN_GTCIMS_HUEY, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = gtcims.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (length(shifts) == 0){shifts = 0}
  if (debug==1){print(c("gtcims shift",round(shifts,digits = 2)))}
  gtcims.shift     = gtcims[ind8,]; gtcims.shift$Time_Start = gtcims.shift$Time_Start - round(shifts,digits = 2)
  tmp = merge(tmp,gtcims.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ MOORE - find best shift - ++++++++ -----
  moore.shift     = moore[ind9,]; moore.shift$Time_Start = moore.shift$Time_mid #+ shift18
  moore.shift  <- moore.shift[, !duplicated(colnames(moore.shift))]
  ind1B = which(moore.shift$FastCNgt6nm_stdPT_MOORE == max(moore.shift$FastCNgt6nm_stdPT_MOORE, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = moore.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("Moore shift",round(shifts,digits = 2)))}
  moore.shift     = moore[ind9,]; moore.shift$Time_Start = moore.shift$Time_mid - round(shifts,digits = 2)
  tmp = merge(tmp,moore.shift,  by='Time_Start', all.x=TRUE)
  
    # ++++++++ merge with met
  tmp = merge(tmp, met, by='Time_Start', all.x=TRUE)
  
  ind = which(is.finite(tmp$Time_Start) )#& is.finite(tmp$CO2_7000_ppm))
  tmp= tmp[ind,]
  return(tmp)
}

time_alignORIG = function(co2.5hz, shift, co.ch4.5hz, shift2,
                      warneke.5hz, shift3, isaf.5hz, shift4,  
                      rollins.5hz, shift5,rollinsS.5hz, shift6,
                      cit.5hz, shift7, gtcims.5hz, shift8,
                      moore.5hz, shift9,
                      met.5hz){
  co2.5hz.shift = co2.5hz
  co2.5hz.shift$Time_Start = co2.5hz.shift$Time_Start + shift
  co.ch4.5hz.shift = co.ch4.5hz
  co.ch4.5hz$Time_Start = co.ch4.5hz$Time_Start + shift2
  warneke.5hz.shift = warneke.5hz
  warneke.5hz.shift$Time_Start = warneke.5hz.shift$Time_Start+ shift3
  isaf.5hz.shift = isaf.5hz
  isaf.5hz.shift$Time_Start = isaf.5hz.shift$Time_Start+ shift4
  rollins.5hz.shift = rollins.5hz
  rollins.5hz.shift$Time_Start = rollins.5hz.shift$Time_Start+ shift5
  rollinsS.5hz.shift = rollinsS.5hz
  rollinsS.5hz.shift$Time_Start = rollinsS.5hz.shift$Time_Start+ shift6
  cit.5hz.shift = cit.5hz
  cit.5hz.shift$Time_Start = cit.5hz.shift$Time_Start+ shift7
  gtcims.5hz.shift = gtcims.5hz
  gtcims.5hz.shift$Time_Start = gtcims.5hz.shift$Time_Start+ shift8
  moore.5hz.shift = moore.5hz
  moore.5hz.shift$Time_Start = moore.5hz.shift$Time_Start+ shift9
  
  tmp = merge(co2.5hz.shift, co.ch4.5hz.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, warneke.5hz.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, isaf.5hz.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, rollins.5hz.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, rollinsS.5hz.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, cit.5hz.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, gtcims.5hz.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, moore.5hz.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, met.5hz, by='Time_Start', all.x=TRUE)
  #tmp = merge(tmp, mce.5hz, by='Time_Start', all.x=TRUE)
  
  ind = which(is.finite(tmp$Time_Start) )#& is.finite(tmp$CO2_7000_ppm))
  tmp= tmp[ind,]
  
  return(tmp)
}

time_alignSLOW = function(co2,   shift,  co.ch4,   shift2,
                      warneke,   shift3, isaf,     shift4,  
                      rollins,   shift5, rollinsS, shift6,
                      cit,       shift7, gtcims,   shift8,
                      ryerson,   shift9, jimenez,   shift10,
                      schwarz,   shift11, freid.c2h6, shift12, freid.ch2o, shift13,
                      womack,    shift14, stclair, shift15, 
                      veres,     shift16, wisthaler, shift17,
                      blake,     shift18, apel, shift19, gilman, shift20,
                      met){

  co2.shift = co2 ;       co2.shift$Time_Start = co2.shift$Time_Start + shift
  co.ch4.shift = co.ch4;  co.ch4$Time_Start = co.ch4.shift$Time_Start + shift2
  warneke.shift = warneke ; warneke.shift$Time_Start = warneke.shift$Time_Start+ shift3
  isaf.shift = isaf;  isaf.shift$Time_Start = isaf.shift$Time_Start+ shift4
  rollins.shift = rollins; rollins.shift$Time_Start = rollins.shift$Time_Start+ shift5
  rollinsS.shift = rollinsS;  rollinsS.shift$Time_Start = rollinsS.shift$Time_Start+ shift6
  cit.shift = cit ;  cit.shift$Time_Start = cit.shift$Time_Start+ shift7
  gtcims.shift     = gtcims ;  gtcims.shift$Time_Start = gtcims.shift$Time_Start+ shift8
  ryerson.shift    = ryerson ; ryerson.shift$Time_Start = ryerson.shift$Time_Start+ shift9
  jimenez.shift    = jimenez ; jimenez.shift$Time_Start = jimenez.shift$Time_Start +   shift10
  schwarz.shift    = schwarz ;schwarz.shift$Time_Start = schwarz.shift$Time_Start + shift11
  freid.c2h6.shift = freid.c2h6; freid.c2h6.shift$Time_Start = freid.c2h6.shift$Time_Start+ shift12
  freid.ch2o.shift = freid.ch2o; freid.ch2o.shift$Time_Start = freid.ch2o.shift$Time_Start+ shift13
  womack.shift     = womack; womack.shift$Time_Start = womack.shift$Time_Start + shift14
  # remove duplicated columns
  womack.shift  <- womack.shift[, !duplicated(colnames(womack.shift))]
  stclair.shift    = stclair; stclair.shift$Time_Start = stclair.shift$Time_Start + shift15
  veres.shift      = veres; veres.shift$Time_Start = veres.shift$Time_Start + shift16
  veres.shift  <- veres.shift[, !duplicated(colnames(veres.shift))]
  wisthaler.shift  = wisthaler; wisthaler.shift$Time_Start = wisthaler.shift$Time_Start + shift17
  blake.shift      = blake; blake.shift$Time_Start = blake.shift$Time_Start + shift18
  apel.shift       = apel; apel.shift$Time_Start = apel.shift$Time_Start + shift19
  gilman.shift     = gilman; gilman.shift$Time_Start = gilman.shift$Time_Start + shift20
  
  tmp = merge(co2.shift, co.ch4.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, warneke.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, isaf.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, rollins.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, rollinsS.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, cit.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, gtcims.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, ryerson.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, jimenez.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, schwarz.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, freid.c2h6.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, freid.ch2o.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, womack.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, stclair.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, veres.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, wisthaler.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, blake.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, apel.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, gilman.shift, by='Time_Start', all.x=TRUE)
  tmp = merge(tmp, met, by='Time_Start', all.x=TRUE)
  tmp  <- tmp[, !duplicated(colnames(tmp))]
  
  ind = which(is.finite(tmp$Time_Start))# & is.finite(tmp$CO2_7000_ppm))
  tmp= tmp[ind,]
  
  return(tmp)
}
time_alignSLOWNOCANS = function(start,stop,co2,   co.ch4,  
                          warneke,   isaf,  rollins, rollinsS, 
                          cit,       gtcims, ryerson,   jimenez,  
                          schwarz,    freid.c2h6,  freid.ch2o, 
                          womack, stclair, veres,     moore, met){ # wisthaler,
  
  start = start - 5; stop = stop + 5 # to allow space for shifts
  ind1 = which(co2$Time_Start >= start & co2$Time_Start <= stop )
  ind2 = which(co.ch4$Time_Start >= start & co.ch4$Time_Start <= stop )
  ind3 = which(warneke$Time_Start >= start & warneke$Time_Start <= stop )
  ind4 = which(isaf$Time_Start >= start & isaf$Time_Start <= stop )
  ind5 = which(rollins$Time_Start >= start & rollins$Time_Start <= stop )
  ind6 = which(rollinsS$Time_Start >= start & rollinsS$Time_Start <= stop )
  ind7 = which(cit$Time_Start >= start & cit$Time_Start <= stop )
  ind8 = which(gtcims$Time_Start >= start & gtcims$Time_Start <= stop )
  ind9 = which(ryerson$Time_Start >= start & ryerson$Time_Start <= stop )
  ind10 = which(jimenez$Time_Start >= start & jimenez$Time_Start <= stop )
  ind11 = which(schwarz$Time_Start >= start & schwarz$Time_Start <= stop )
  ind12 = which(freid.c2h6$Time_Start >= start & freid.c2h6$Time_Start <= stop )
  ind13 = which(freid.ch2o$Time_Start >= start & freid.ch2o$Time_Start <= stop )
  ind14 = which(womack$Time_Start >= start & womack$Time_Start <= stop )
  ind15 = which(stclair$Time_Start >= start & stclair$Time_Start <= stop ) 
  ind16 = which(veres$Time_Start >= start & veres$Time_Start <= stop )
#  ind17 = which(wisthaler$Time_Start >= start & wisthaler$Time_Start <= stop )
  ind18 = which(moore$Time_mid >= start & moore$Time_mid <= stop )
  if (length(ind18) == 0){  ind18 = which(moore$Time_Start >= start & moore$Time_Start <= stop )}
  
  ind19 = which(met$Time_mid >= start & met$Time_mid <= stop )
  
  shifts = seq(-5,5,by=1)
  debug=1
  # Aligning everything to CO
  # ---- ++++++++ CO2 find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    co2.shift  = co2[ind1,]
    co2.shift  <- co2.shift[, !duplicated(colnames(co2.shift))]
    co2.shift$Time_Start = co2.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],co2.shift, by='Time_Start', all.x=TRUE)
    ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$CO2_7000_ppm_DISKIN)
    corsM = c(corsM,ccM$estimate )
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("co2 shift",shifts[indM]))}
  co2.shift  = co2
  co2.shift$Time_Start = co2.shift$Time_Start+ shifts[indM]
  tmp = merge(co.ch4,co2.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Warneke find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    warneke.shift  = warneke[ind3,]
    warneke.shift  <- warneke.shift[, !duplicated(colnames(warneke.shift))]
    warneke.shift$Time_Start = warneke.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],warneke.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$HCN_NOAAPTR_ppbv_WARNEKE))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$HCN_NOAAPTR_ppbv_WARNEKE)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("Warneke shift",shifts[indM]))}
  warneke.shift  = warneke
  warneke.shift$Time_Start = warneke.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, warneke.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ ISAF find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    isaf.shift  = isaf[ind4,]
    isaf.shift  <- isaf.shift[, !duplicated(colnames(isaf.shift))]
    isaf.shift$Time_Start = isaf.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],isaf.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$CH2O_ISAF_HANISCO))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$CH2O_ISAF_HANISCO)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }  
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("isaf shift",shifts[indM]))}
  isaf.shift  = isaf
  isaf.shift$Time_Start = isaf.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, isaf.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Rollins N find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    rollins.shift  = rollins[ind5,]
    rollins.shift  <- rollins.shift[, !duplicated(colnames(rollins.shift))]
    rollins.shift$Time_Start = rollins.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],rollins.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$NO_LIF_ROLLINS))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$NO_LIF_ROLLINS)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("Rollins shift",shifts[indM]))}
  rollins.shift  = rollins
  rollins.shift$Time_Start = rollins.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, rollins.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Rollins S find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    rollinsS.shift  = rollinsS[ind6,]
    rollinsS.shift  <- rollinsS.shift[, !duplicated(colnames(rollinsS.shift))]
    rollinsS.shift$Time_Start = rollinsS.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],rollinsS.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$SO2_LIF_ROLLINS))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$SO2_LIF_ROLLINS)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("Rollins S shift",shifts[indM]))}
  rollinsS.shift  = rollinsS
  rollinsS.shift$Time_Start = rollinsS.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, rollinsS.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++CIT find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    cit.shift  = cit[ind7,]
    cit.shift  <- cit.shift[, !duplicated(colnames(cit.shift))]
    cit.shift$Time_Start = cit.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],cit.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$HCN.1Hz_CIT_WENNBERG))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$HCN.1Hz_CIT_WENNBERG)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
    
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("CIT shift",shifts[indM]))}
  cit.shift  = cit
  cit.shift$Time_Start = cit.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, cit.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++GTCIMS find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    gtcims.shift  = gtcims[ind8,]
    gtcims.shift  <- gtcims.shift[, !duplicated(colnames(gtcims.shift))]
    gtcims.shift$Time_Start = gtcims.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],gtcims.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$APAN_GTCIMS_HUEY))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$APAN_GTCIMS_HUEY)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
    }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM, na.rm=TRUE))){indM=6} # shift of 0
  if (debug==1){print(c("Gtcims shift",shifts[indM]))}
  gtcims.shift  = gtcims
  gtcims.shift$Time_Start = gtcims.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, gtcims.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Ryerson find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    ryerson.shift  = ryerson[ind9,]
    ryerson.shift  <- ryerson.shift[, !duplicated(colnames(ryerson.shift))]
    ryerson.shift$Time_Start = ryerson.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],ryerson.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$NO_CL_RYERSON))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN,  tmp2$NO_CL_RYERSON)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}

  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM,na.rm=TRUE))){indM=6}
  if (debug==1){print(c("Ryerson shift",shifts[indM]))}
  ryerson.shift  = ryerson
  ryerson.shift$Time_Start = ryerson.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, ryerson.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Jimenez  find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    jimenez.shift  = jimenez[ind10,]
    jimenez.shift  <- jimenez.shift[, !duplicated(colnames(jimenez.shift))]
    jimenez.shift$Time_Start = jimenez.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],jimenez.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$C6H10O5_JIMENENZ))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$C6H10O5_JIMENENZ)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM,na.rm=TRUE))){indM=6}
  if (debug==1){print(c("Jimenez shift",shifts[indM]))}
  jimenez.shift  = jimenez
  jimenez.shift$Time_Start = jimenez.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, jimenez.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Schwarz-  find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    schwarz.shift  = schwarz[ind11,]
    schwarz.shift  <- schwarz.shift[, !duplicated(colnames(schwarz.shift))]
    schwarz.shift$Time_Start = schwarz.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],schwarz.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$BC_mass_90_550_nm_SCHWARZ)& is.finite(tmp2$CO_DACOM_DISKIN))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$BC_mass_90_550_nm_SCHWARZ)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM,na.rm=TRUE))){indM=6}
  if (debug==1){print(c("Schwarz shift",shifts[indM]))}
  schwarz.shift  = schwarz
  schwarz.shift$Time_Start = schwarz.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, schwarz.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Fried C2H6 -  find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    freid.c2h6.shift  = freid.c2h6[ind12,]
    freid.c2h6.shift  <- freid.c2h6.shift[, !duplicated(colnames(freid.c2h6.shift))]
    freid.c2h6.shift$Time_Start = freid.c2h6.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],freid.c2h6.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$C2H6_CAMS_pptv))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$C2H6_CAMS_pptv)
      corsM = c(corsM,ccM$estimate )
    } else{corsM = c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM,na.rm=TRUE))){indM=6}
  if (debug==1){print(c("Freid.c2h6 shift",shifts[indM]))}
  freid.c2h6.shift  = freid.c2h6
  freid.c2h6.shift$Time_Start = freid.c2h6.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, freid.c2h6.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ freid CH2O -  find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    freid.ch2o.shift  = freid.ch2o[ind13,]
    freid.ch2o.shift  <- freid.ch2o.shift[, !duplicated(colnames(freid.ch2o.shift))]
    freid.ch2o.shift$Time_Start = freid.ch2o.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],freid.ch2o.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$C2H6_CAMS_pptv))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$CH2O_CAMS_pptv)
      corsM = c(corsM,ccM$estimate )
    } else{corsM = c(corsM,NaN)}
    # if (debug==1){print(c(shifts[i],ccM$estimate))
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM,na.rm=TRUE))){indM=6}
  if (debug==1){print(c("freid.ch2o shift",shifts[indM]))}
  freid.ch2o.shift  = freid.ch2o
  freid.ch2o.shift$Time_Start = freid.ch2o.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, freid.ch2o.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Womack- find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    womack.shift  = womack[ind14,]
    womack.shift  <- womack.shift[, !duplicated(colnames(womack.shift))]
    womack.shift$Time_Start = womack.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],womack.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$CHOCHO_ACES_WOMACK))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$CHOCHO_ACES_WOMACK)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM,na.rm=TRUE))){indM=6}
  if (debug==1){print(c("Womack shift",shifts[indM]))}
  womack.shift  = womack
  womack.shift$Time_Start = womack.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, womack.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ St Clair- find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    stclair.shift  = stclair[ind15,]
    stclair.shift  <- stclair.shift[, !duplicated(colnames(stclair.shift))]
    stclair.shift$Time_Start = stclair.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],stclair.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$NO2_CANOE))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$NO2_CANOE)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM,na.rm=TRUE))){indM=6}
  if (debug==1){print(c("Stclair shift",shifts[indM]))}
  stclair.shift  = stclair
  stclair.shift$Time_Start = stclair.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, stclair.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Veres - find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    veres.shift  = veres[ind16,]
    veres.shift  <- veres.shift[, !duplicated(colnames(veres.shift))]
    veres.shift$Time_Start = veres.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],veres.shift, by='Time_Start', all.x=TRUE)
    if (length(which(is.finite(tmp2$HNO2_NOAACIMS_VERES))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$HNO2_NOAACIMS_VERES)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM,na.rm=TRUE))){indM=6}
  if (debug==1){print(c("Veres shift",shifts[indM]))}
  veres.shift  = veres
  veres.shift$Time_Start = veres.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, veres.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Wisthaler - find best shift - ++++++++ -----
 # corsM = c()
#  for (i in 1:length(shifts)){
#    wisthaler.shift  = wisthaler[ind17,]
#    wisthaler.shift$Time_Start = wisthaler.shift$Time_Start+ shifts[i]
#    tmp2 = merge(co.ch4[ind2,],wisthaler.shift, by='Time_Start', all.x=TRUE)
#    if (length(which(is.finite(tmp2$NH3_UIOPTR_ppbV_WISTHALER))) > 2){
#      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$NH3_UIOPTR_ppbV_WISTHALER)
#      corsM = c(corsM,ccM$estimate )
#    } else{corsM=c(corsM,NaN)}
#  }
#  indM = which(corsM == max(corsM, na.rm=TRUE))
##  if (!is.finite(max(corsM,na.rm=TRUE))){indM=6}
#  if (debug==1){print(c("Wisthaler shift",shifts[indM]))}
#  wisthaler.shift  = wisthaler
#  wisthaler.shift$Time_Start = wisthaler.shift$Time_Start+ shifts[indM]
#  tmp = merge(tmp, wisthaler.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ MOORE - find best shift - ++++++++ -----
  corsM = c()
  for (i in 1:length(shifts)){
    moore.shift     = moore[ind18,]; moore.shift$Time_Start = moore.shift$Time_mid #+ shift18
    moore.shift$Time_Start = moore.shift$Time_Start+ shifts[i]
    tmp2 = merge(co.ch4[ind2,],moore.shift, by='Time_Start', all.x=TRUE)
    ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$CNgt3nm_stdPT)
    if (length(which(is.finite(tmp2$CNgt3nm_stdPT))) > 2){
      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$CNgt3nm_stdPT)
      corsM = c(corsM,ccM$estimate )
    } else{corsM=c(corsM,NaN)}
  }
  indM = which(corsM == max(corsM, na.rm=TRUE))
  if (!is.finite(max(corsM,na.rm=TRUE))){indM=6}
  if (debug==1){print(c("Moore shift",shifts[indM]))}
  moore.shift     = moore; moore.shift$Time_Start = moore.shift$Time_mid #+ shift18
  moore.shift$Time_Start = moore.shift$Time_Start+ shifts[indM]
  tmp = merge(tmp, moore.shift, by='Time_Start', all.x=TRUE)
  
  # ++++++++ merge with met

  tmp = merge(tmp, met, by='Time_Start', all.x=TRUE)
  tmp  <- tmp[, !duplicated(colnames(tmp))]
  
  ind = which(is.finite(tmp$Time_Start))# & is.finite(tmp$CO2_7000_ppm))
  tmp= tmp[ind,]
  
  return(tmp)
}

time_alignSLOWNOCANSPEAKS = function(start,stop,co2,   co.ch4,  
                                warneke,   isaf,  rollins, rollinsS, 
                                cit,       gtcims, ryerson,   jimenez,  
                                schwarz,    freid.c2h6,  freid.ch2o, 
                                womack, stclair, veres,     moore, met){ # wisthaler,
  
  start = start - 5; stop = stop + 5 # to allow space for shifts
  ind1 = which(co2$Time_Start >= start & co2$Time_Start <= stop )
  ind2 = which(co.ch4$Time_Start >= start & co.ch4$Time_Start <= stop )
  ind3 = which(warneke$Time_Start >= start & warneke$Time_Start <= stop )
  ind4 = which(isaf$Time_Start >= start & isaf$Time_Start <= stop )
  ind5 = which(rollins$Time_Start >= start & rollins$Time_Start <= stop )
  ind6 = which(rollinsS$Time_Start >= start & rollinsS$Time_Start <= stop )
  ind7 = which(cit$Time_Start >= start & cit$Time_Start <= stop )
  ind8 = which(gtcims$Time_Start >= start & gtcims$Time_Start <= stop )
  ind9 = which(ryerson$Time_Start >= start & ryerson$Time_Start <= stop )
  ind10 = which(jimenez$Time_Start >= start & jimenez$Time_Start <= stop )
  ind11 = which(schwarz$Time_Start >= start & schwarz$Time_Start <= stop )
  ind12 = which(freid.c2h6$Time_Start >= start & freid.c2h6$Time_Start <= stop )
  ind13 = which(freid.ch2o$Time_Start >= start & freid.ch2o$Time_Start <= stop )
  ind14 = which(womack$Time_Start >= start & womack$Time_Start <= stop )
  ind15 = which(stclair$Time_Start >= start & stclair$Time_Start <= stop ) 
  ind16 = which(veres$Time_Start >= start & veres$Time_Start <= stop )
 # ind17 = which(wisthaler$Time_Start >= start & wisthaler$Time_Start <= stop )
  ind18 = which(moore$Time_mid >= start & moore$Time_mid <= stop )
  ind19 = which(met$Time_mid >= start & met$Time_mid <= stop )
  
  debug=1
  # Aligning everything to CO
  
  # ---- ++++++++ CO2 find best shift - ++++++++ -----
  
  co2.shift  = co2[ind1,]
  co2.shift  <- co2.shift[, !duplicated(colnames(co2.shift))]
  ind1B = which(co2.shift$CO2_7000_ppm_DISKIN == max(co2.shift$CO2_7000_ppm_DISKIN, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = co2.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("co2 shift",shifts))}
  co2.shift  = co2
  co2.shift$Time_Start = co2.shift$Time_Start - shifts
  tmp = merge(co.ch4,co2.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Warneke find best shift - ++++++++ -----
  warneke.shift  = warneke[ind3,]
  warneke.shift  <- warneke.shift[, !duplicated(colnames(warneke.shift))]
  ind1B = which(warneke.shift$Benzene_NOAAPTR_ppbv_WARNEKE== max(warneke.shift$Benzene_NOAAPTR_ppbv_WARNEKE, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = warneke.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("warneke shift",shifts))}
  warneke.shift     = warneke[ind3,]; warneke.shift$Time_Start = warneke.shift$Time_Start- shifts
  tmp = merge(tmp,warneke.shift,  by='Time_Start', all.x=TRUE)
  # ---- ++++++++ ISAF find best shift - ++++++++ -----
  isaf.shift  = isaf[ind4,]
  isaf.shift  <- isaf.shift[, !duplicated(colnames(isaf.shift))]
  ind1B = which(isaf.shift$CH2O_ISAF_HANISCO== max(isaf.shift$CH2O_ISAF_HANISCO, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = isaf.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("isaf shift",shifts))}
  isaf.shift     = isaf[ind4,]; isaf.shift$Time_Start = isaf.shift$Time_Start - shifts
  tmp = merge(tmp,isaf.shift,  by='Time_Start', all.x=TRUE)
  # ---- ++++++++ Rollins N find best shift - ++++++++ -----
  rollins.shift  = rollins[ind5,]
  rollins.shift  <- rollins.shift[, !duplicated(colnames(rollins.shift))]
  ind1B = which(rollins.shift$NO_LIF_ROLLINS== max(rollins.shift$NO_LIF_ROLLINS, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = rollins.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("rollins shift",shifts))}
  rollins.shift     = rollins[ind5,]; rollins.shift$Time_Start = rollins.shift$Time_Start - shifts
  tmp = merge(tmp,rollins.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Rollins S find best shift - ++++++++ -----
  rollinsS.shift  = rollinsS[ind6,]
  rollinsS.shift  <- rollinsS.shift[, !duplicated(colnames(rollinsS.shift))]
  ind1B = which(rollinsS.shift$SO2_LIF_ROLLINS== max(rollinsS.shift$SO2_LIF_ROLLINS, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = rollinsS.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("rollinsS shift",shifts))}
  rollinsS.shift     = rollinsS[ind6,]; rollinsS.shift$Time_Start = rollinsS.shift$Time_Start - shifts
  tmp = merge(tmp,rollinsS.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++CIT find best shift - ++++++++ -----
  cit.shift  = cit[ind7,]
  cit.shift  <- cit.shift[, !duplicated(colnames(cit.shift))]
  ind1B = which(cit.shift$PHENOL.1Hz_CIT_WENNBERG== max(cit.shift$PHENOL.1Hz_CIT_WENNBERG, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = cit.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("CIT shift",shifts))}
  cit.shift     = cit[ind7,]; cit.shift$Time_Start = cit.shift$Time_Start - shifts
  tmp = merge(tmp,cit.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++GTCIMS find best shift - ++++++++ -----
  gtcims.shift  = gtcims[ind8,]
  gtcims.shift  <- gtcims.shift[, !duplicated(colnames(gtcims.shift))]
  ind1B = which(gtcims.shift$APAN_GTCIMS_HUEY== max(gtcims.shift$APAN_GTCIMS_HUEY, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = gtcims.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("gtcims shift",shifts))}
  gtcims.shift     = gtcims[ind8,]; gtcims.shift$Time_Start = gtcims.shift$Time_Start - shifts
  tmp = merge(tmp,gtcims.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Ryerson find best shift - ++++++++ -----
  ryerson.shift  = ryerson[ind9,]
  ryerson.shift  <- ryerson.shift[, !duplicated(colnames(ryerson.shift))]
  ind1B = which(ryerson.shift$NO_CL_RYERSON== max(ryerson.shift$NO_CL_RYERSON, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = ryerson.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("Ryerson shift",shifts))}
  ryerson.shift     = ryerson[ind9,]; ryerson.shift$Time_Start = ryerson.shift$Time_Start - shifts
  tmp = merge(tmp,ryerson.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Jimenez  find best shift - ++++++++ -----
  jimenez.shift  = jimenez[ind10,]
  jimenez.shift  <- jimenez.shift[, !duplicated(colnames(jimenez.shift))]
  ind1B = which(jimenez.shift$OA_PM1_AMS_JIMENEZ == max(jimenez.shift$OA_PM1_AMS_JIMENEZ, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = jimenez.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("jimenez shift",shifts))}
  jimenez.shift     = jimenez[ind10,]; jimenez.shift$Time_Start = jimenez.shift$Time_Start - shifts
  tmp = merge(tmp,jimenez.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Schwarz-  find best shift - ++++++++ -----
  schwarz.shift  = schwarz[ind11,]
  schwarz.shift  <- schwarz.shift[, !duplicated(colnames(schwarz.shift))]
  ind1B = which(schwarz.shift$BC_mass_90_550_nm_SCHWARZ == max(schwarz.shift$BC_mass_90_550_nm_SCHWARZ, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = schwarz.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("schwarz shift",shifts))}
  schwarz.shift     = schwarz[ind11,]; schwarz.shift$Time_Start = schwarz.shift$Time_Start - shifts
  tmp = merge(tmp,schwarz.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Fried C2H6 -  find best shift - ++++++++ -----
  freid.c2h6.shift  = freid.c2h6[ind12,]
  freid.c2h6.shift  <- freid.c2h6.shift[, !duplicated(colnames(freid.c2h6.shift))]
  ind1B = which(freid.c2h6.shift$C2H6_CAMS_pptv == max(freid.c2h6.shift$C2H6_CAMS_pptv, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = freid.c2h6.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("freid.c2h6 shift",shifts))}
  freid.c2h6.shift     = freid.c2h6[ind12,]; freid.c2h6.shift$Time_Start = freid.c2h6.shift$Time_Start - shifts
  tmp = merge(tmp,freid.c2h6.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ freid CH2O -  find best shift - ++++++++ -----
  freid.ch2o.shift  = freid.ch2o[ind13,]
  freid.ch2o.shift  <- freid.ch2o.shift[, !duplicated(colnames(freid.ch2o.shift))]
  ind1B = which(freid.ch2o.shift$CH2O_CAMS_pptv == max(freid.ch2o.shift$CH2O_CAMS_pptv, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = freid.ch2o.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("freid.ch2o shift",shifts))}
  freid.ch2o.shift     = freid.ch2o[ind13,]; freid.ch2o.shift$Time_Start = freid.ch2o.shift$Time_Start - shifts
  tmp = merge(tmp,freid.ch2o.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Womack- find best shift - ++++++++ -----
  womack.shift  = womack[ind14,]
  womack.shift  <- womack.shift[, !duplicated(colnames(womack.shift))]
  ind1B = which(womack.shift$CHOCHO_ACES_WOMACK == max(womack.shift$CHOCHO_ACES_WOMACK, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = womack.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("womack shift",shifts))}
  womack.shift     = womack[ind14,]; womack.shift$Time_Start = womack.shift$Time_Start - shifts
  tmp = merge(tmp,womack.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ St Clair- find best shift - ++++++++ -----
  stclair.shift  = stclair[ind15,]
  stclair.shift  <- stclair.shift[, !duplicated(colnames(stclair.shift))]
  ind1B = which(stclair.shift$NO2_CANOE== max(stclair.shift$NO2_CANOE, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = stclair.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("stclair shift",shifts))}
  stclair.shift     = stclair[ind15,]; stclair.shift$Time_Start = stclair.shift$Time_Start - shifts
  tmp = merge(tmp,stclair.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Veres - find best shift - ++++++++ -----
  veres.shift  = veres[ind16,]
  veres.shift  <- veres.shift[, !duplicated(colnames(veres.shift))]
  ind1B = which(veres.shift$HNO2_NOAACIMS_VERES == max(veres.shift$HNO2_NOAACIMS_VERES, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = veres.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("veres shift",shifts))}
  veres.shift     = veres[ind16,]; veres.shift$Time_Start = veres.shift$Time_Start - shifts
  tmp = merge(tmp,veres.shift,  by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ Wisthaler - find best shift leave as correlation as peaks are bad - ++++++++ -----
 # corsM = c()
#  for (i in 1:length(shifts)){
#    wisthaler.shift  = wisthaler[ind17,]
##    wisthaler.shift$Time_Start = wisthaler.shift$Time_Start+ shifts[i]
#    tmp2 = merge(co.ch4[ind2,],wisthaler.shift, by='Time_Start', all.x=TRUE)
#    if (length(which(is.finite(tmp2$NH3_UIOPTR_ppbV_WISTHALER))) > 2){
#      ccM=cor.test(tmp2$CO_DACOM_DISKIN, tmp2$NH3_UIOPTR_ppbV_WISTHALER)
#      corsM = c(corsM,ccM$estimate )
##    } else{corsM=c(corsM,NaN)}
#  }
#  indM = which(corsM == max(corsM, na.rm=TRUE))
#  if (!is.finite(max(corsM,na.rm=TRUE))){indM=6}
#  if (debug==1){print(c("Wisthaler shift",shifts[indM]))}
#  wisthaler.shift  = wisthaler
#  wisthaler.shift$Time_Start = wisthaler.shift$Time_Start+ shifts[indM]
#  tmp = merge(tmp, wisthaler.shift, by='Time_Start', all.x=TRUE)
  
  # ---- ++++++++ MOORE - find best shift - ++++++++ -----
  
  moore.shift     = moore[ind18,]; moore.shift$Time_Start = moore.shift$Time_mid #+ shift18
  moore.shift  <- moore.shift[, !duplicated(colnames(moore.shift))]
  ind1B = which(moore.shift$CNgt3nm_stdPT == max(moore.shift$CNgt3nm_stdPT, na.rm=TRUE))
  ind2B = which(co.ch4$CO_DACOM_DISKIN[ind2] == max(co.ch4$CO_DACOM_DISKIN[ind2], na.rm=TRUE))
  shifts = moore.shift$Time_Start[ind1B] - co.ch4$Time_Start[ind2[ind2B]] 
  if (debug==1){print(c("Moore shift",shifts))}
  moore.shift     = moore[ind18,]; moore.shift$Time_Start = moore.shift$Time_mid - shifts
  tmp = merge(tmp,moore.shift,  by='Time_Start', all.x=TRUE)
  
  # ++++++++ merge with met
  
  tmp = merge(tmp, met, by='Time_Start', all.x=TRUE)
  tmp  <- tmp[, !duplicated(colnames(tmp))]
  
  ind = which(is.finite(tmp$Time_Start))# & is.finite(tmp$CO2_7000_ppm))
  tmp= tmp[ind,]
  
  return(tmp)
}


map_myfire = function(firedata){
  map.us <- get_map(c(left =min(firedata$Longitude_YANG, na.rm=TRUE)-.05, bottom =min(firedata$Latitude_YANG, na.rm=TRUE)-0.05,
                      right = max(firedata$Longitude_YANG, na.rm=TRUE)+.05, top = max(firedata$Latitude_YANG, na.rm=TRUE)+.05))
  ggmap(map.us) +
    geom_point(data = firedata, aes(x = Longitude_YANG, y = Latitude_YANG,
                                    colour =((CO_DACOM_DISKIN))),  size =1) +   ggtitle("FIREX-AQ")
}
plotpass = function(pass){
  par(mfrow=c(1,2))
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, ylab='CO, ppb', xlab='Time_Start', type='l', pch=19)
  tmp = cor.test(pass$CO_DACOM_DISKIN,  pass$CO2_7000_ppm_DISKIN)
  plot( pass$CO_DACOM_DISKIN,  pass$CO2_7000_ppm_DISKIN, xlab='CO, ppb', ylab='CO2, ppm', pch=19)
  print(round((tmp$estimate)^2, digits=2))
  text(390,415, paste("r2 = ",round((tmp$estimate)^2, digits=2)))
}
maxfun = function(pass,sKT){
  ccc = colnames(pass)
  tmp = which(ccc == sKT)
  indKT = which(is.finite(pass[,tmp]))
  pKT = max(pass[,tmp], na.rm=TRUE) 
  return(pKT)
}
plotfire=function(fire, flags, indA, name, xlim){
  plot(fire$Time_Start, fire$CO_DACOM_DISKIN, ylab='CO, ppb', type='l', xlim=xlim,main=paste(name,", DOY=",fire$Day_Of_Year[1], sep=''))
  for (i in 1:length(indA)){
    start = flags$`transect_start_time (UTC s from midnight)`[indA[i]]
    stop = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]
    ind = which(fire$Time_Start >= start & fire$Time_Start <= stop )
    lines(fire$Time_Start[ind], fire$CO_DACOM_DISKIN[ind], col='red')
  }
}

plotEFvsMCE = function(fire, xvar){
  ind = which(fire$variable == xvar)
  ggplot(fire[ind,]) + geom_point(aes(x=mce.5hz,y=ChosenEF.5hz,color=fuel, size=fuel)) + xlab('MCE')+ ylab(xvar) + theme_classic()+
    theme(text = element_text(size=20))+ scale_color_brewer(palette="Dark2")
}
plotEFvsMCE1hz = function(fire, xvar){
  ind = which(fire$variable == xvar)
  ggplot(fire[ind,]) + geom_point(aes(x=mce,y=EF1,color=fuel, size=fuel)) + xlab('MCE')+ ylab(xvar) + theme_classic()+
    theme(text = element_text(size=20))+ scale_color_brewer(palette="Dark2")
}

plotEFvsS = function(fire, xvar, yvar){
  ind = which(fire$variable == xvar )
  xvarfire = fire[ind,]
  ind2 = which(fire$variable == yvar)
  yvarfire = fire[ind2,]
  cc=colnames(fire)
  bothfire = merge(xvarfire, yvarfire, by='uniqueid', all.x=TRUE)
  
  ggplot(bothfire) +
    geom_point(aes(x=mce.5hz.x,y=ChosenEF.5hz.x,color=(round(ChosenEF.5hz.y, digits=1))), size=4) + 
    xlab('MCE')+ ylab(xvar) + theme_classic()+
    theme(text = element_text(size=20))+ scale_colour_gradientn(colours=rainbow(4))
  
#  ggplot(bothfire) +
##    geom_point(aes(x=mce.5Hz.x,y=ChosenEF.5hz.x,color=factor(round(ChosenEF.5hz.y, digits=1)), shape=fuel.x), size=4) + 
#    xlab('MCE')+ ylab(xvar) + theme_classic()+
#    theme(text = element_text(size=20))+ scale_color_brewer(palette="Dark2")
}

getfuel = function(fire, fuel, variable){
  ind = which(fire$variable == vv & fire$fuel == fuel)
  print(c("#Fires",length(unique(allBOTH.ag$fire[ind]))))
  print(c("#Plumes",length(unique(allBOTH.ag$uniqueid[ind]))))
  fueldata = fire$ChosenEF.5hz[ind]
  fueldata = cbind(fueldata, fire$mce.5hz[ind])
  fueldata = as.data.frame(fueldata)
  colnames(fueldata) = c("EF","mce")
  return(fueldata)
}
# 00000000000 OLD 0000000000000000000000
doBOTH5Hz = function(fire.ERstoCO2, fire, fuel){
  fire.ERstoCO2.ER = ERs5hz(fire.ERstoCO2,'CO2_ppb')
  fire.ERstoCO2.EF = getEFs5Hz(fire.ERstoCO2.ER)
  fire.ERstoCO2.EF$fire = fire
  fire.ERstoCO2.EF$fuel = fuel
  
  return(fire.ERstoCO2.EF)
}
getEFs5Hz = function(fire.ERstoCO2){
  # ------ Carbon content - for now just 50%
  cContentCorn = 500 # Assuming 50 % C, 500 g/kg
  
  # ------- # carbons -------
  
  nCCO2 = 1; nCCO = 1; nCCH4 = 1 ; nCOC=1
  fire.ERstoCO2$variable =c("mce" ,
                            # DISKIN
                            "ER_CO2_DISKIN", "ER_CO_DISKIN" ,"ER_CH4_DISKIN", 
                            # WARNEKE
                            'ER_HCN_WARNEKE','ER_CH2O_WARNEKE','ER_CH3OH_WARNEKE','ER_CH3CN_WARNEKE','ER_HNCO_WARNEKE','ER_CH3CHO_WARNEKE','ER_C2H5OH_WARNEKE','ER_HCOOH_WARNEKE','ER_C3H3N_WARNEKE',
                            'ER_C3H4O_WARNEKE','ER_C3H6O_WARNEKE','ER_C2H4O2_WARNEKE','ER_CH3NO2_WARNEKE','ER_C2H6S_WARNEKE','ER_C4H5N_WARNEKE','ER_C4H4O_WARNEKE','ER_C4H6O_WARNEKE','ER_C4H8O_WARNEKE',
                            'ER_C3H6O2_WARNEKE',	'ER_C6H6_WARNEKE','ER_C5H6O_WARNEKE',	'ER_C4H4O2_WARNEKE', 'ER_C4H6O2_WARNEKE',	'ER_C7H8_WARNEKE','ER_C6H6O_WARNEKE','ER_C5H4O2_WARNEKE',	'ER_C6H8O_WARNEKE',	
                            'ER_C4H2O3_WARNEKE',	'ER_C7H5N_WARNEKE',	'ER_C8H8_WARNEKE', 'ER_C7H6O_WARNEKE','ER_C8H10_WARNEKE',	'ER_C7H8O_WARNEKE',	'ER_C6H6O2_WARNEKE','ER_C8H6O_WARNEKE',	'ER_C9H12_WARNEKE','ER_C6H4O3_WARNEKE',
                            'ER_C7H8O2_WARNEKE',	'ER_C10H8_WARNEKE',	'ER_C10H16_WARNEKE','ER_C8H10O2_WARNEKE',	'ER_C8H10O3_WARNEKE',
                            # HANISCO
                            'ER_CH2O_HANISCO',
                            # ROLLINS
                            'ER_SO2_ROLLINS','ER_NO_ROLLINS',
                            # WENNBERG
                            'ER_BUTENE.HN_WENNBERG',  'ER_BUTENE.HP_WENNBERG', 'ER_ETHENE.HN_WENNBERG','ER_ETHENE.HP_WENNBERG',
                            'ER_PHENOL_WENNBERG', 'ER_H2O2_WENNBERG',# 'ER_MHP_WENNBERG', 'ER_GLYC_WENNBERG', 'ER_IEPOX_WENNBERG',
                            'ER_ISOPN_WENNBERG',# 'ER_ISOPOOH_WENNBERG','ER_PAA_WENNBERG',
                            'ER_PROPENE.HN_WENNBERG','ER_PROPENE.HP_WENNBERG',
                            'ER_HCN_WENNBERG','ER_HNO3_WENNBERG',
                            # HUEY
                            'ER_PAN_HUEY','ER_APAN_HUEY',
                            # MaleicAnhydride to Furan ratio
                            'MAtoF',
                            'MAtoF_rsq',
                            # size of plume used to calculate EFs
                            'npts')
  # For now, just include CO2, CO, CH4
  fire.ERstoCO2$TC = fire.ERstoCO2$ERtoCO2[2]* nCCO2 + fire.ERstoCO2$ERtoCO2[3]* nCCO + fire.ERstoCO2$ERtoCO2[4] *nCCH4 
  fire.ERstoCO2$EF = NaN
  fire.ERstoCO2$EF[2]   = cContentCorn * mWCO2/12 * fire.ERstoCO2$ERtoCO2[2] /fire.ERstoCO2$TC[2]
  fire.ERstoCO2$EF[3]   = cContentCorn * mWCO/12 * fire.ERstoCO2$ERtoCO2[3] /fire.ERstoCO2$TC[3]
  fire.ERstoCO2$EF[4]  = cContentCorn * mWCH4/12 * fire.ERstoCO2$ERtoCO2[4] /fire.ERstoCO2$TC[4]
  
  # WARNEKE VOCs
  wMWs = c(27.0,30.0,32.0,41.1,43.0,44.0,46.1,46.0,55.1,56.1,58.1,60.0,61.0,62.1,67.1,68.1,70.1,74.1,74.1,78.1,82.1,84.1,86.1,92.1,94.1,96.1,96.1,98.1,103.1,104.1,106.1,106.2,108.1,110.1,118.1,120.2,124.1,124.1,128.2,136.2,138.1,154.1 )
  for (jj in 5:46){
    fire.ERstoCO2$EF[jj]  = cContentCorn * wMWs[jj-4]/12 * fire.ERstoCO2$ERtoCO2[jj]/fire.ERstoCO2$TC[jj]
  }
  # HANISCO
  fire.ERstoCO2$EF[47] = cContentCorn * mWCH2O/12 * fire.ERstoCO2$ERtoCO2[47]/fire.ERstoCO2$TC[47]
  # ROLLINS
  fire.ERstoCO2$EF[48] = cContentCorn * mWSO2/12 * fire.ERstoCO2$ERtoCO2[48]/fire.ERstoCO2$TC[48]
  fire.ERstoCO2$EF[49]  = cContentCorn * mWNO/12  * fire.ERstoCO2$ERtoCO2[49]/fire.ERstoCO2$TC[49]
  # WENNBERG
  mWsW = c(  mWBUTENEHN, mWBUTENEHP, mWETHENEHN, mWETHENEHP,mWPHENOL,mWH2O2,# mWMHP, mWGLYC, mWIEPOX,
             mWISOPN,# mWISOPOOH,mWPAA, 
             mWPROPENEHN,mWPROPENEHP, mWHCN,mWHNO3)
  for (jj in 50:60){
    fire.ERstoCO2$EF[jj] = cContentCorn * mWsW[jj-49] * fire.ERstoCO2$ERtoCO2[jj]/fire.ERstoCO2$TC[jj]
  }
  fire.ERstoCO2$EF[61]  = cContentCorn * mWPAN/12  * fire.ERstoCO2$ERtoCO2[61]/fire.ERstoCO2$TC[61]
  fire.ERstoCO2$EF[62]  = cContentCorn * mWAPAN/12  * fire.ERstoCO2$ERtoCO2[62]/fire.ERstoCO2$TC[62]
  
  
  # ------datareturn
  return(fire.ERstoCO2)
}
getEFs = function(fire.ERstoCO2){
  mWCO2 = 44.01 ; mWCO = 28.01 ; mWCH4 = 16.04; mWNO = 30.01 ; mWNO2 = 46.0055 ; mWHONO = 47 
  mWOC = 12 ; mWCH2O = 30.031 ; mWC2H6 = 30.07 ; mWSO2 = 64.066 ; mWCHOCHO =58.04 ; mWCH3COCHO = 72.06
  mWN2O5 =108.01; mWBrCN = 105.921; mWBrCl =115.357
  mWBrO = 95.904; mWCH3COOCl = 78.5;  mWCl2 = 70.9 ;mWClNO2= 81.46
  mWHCN = 27.0253;  mWHCOOH= 46.03;  mWHNCO= 43.03 ;  mWHPMTF= 108.12
  mWNH3 = MolecularWeight(formula=list(N=1, H=3))
  mWHNO3 = MolecularWeight(formula=list( O=3,N=1, H=1))
  mWHCN = MolecularWeight(formula=list( C=1,N=1, H=1))
  mWPROPENEHP = MolecularWeight(formula=list( C=3,O=3, H=8))
  mWPROPENEHN = MolecularWeight(formula=list( C=3,O=4, H=7, N=1))
  mWPAA = MolecularWeight(formula=list( C=2,O=3, H=4))
  mWIEPOX = MolecularWeight(formula=list( C=5,O=3, H=10))
  mWGLYC = MolecularWeight(formula=list( C=2,O=2, H=4))
  mWETHENEHN = MolecularWeight(formula=list(C=2, O=4, H=5, N=1))
  mWETHENEHP = MolecularWeight(formula=list(C=2, O=3, H=6))
  mWBUTENEHP = MolecularWeight(formula=list(C=4, O=3, H=10))
  mWBUTENEHN = MolecularWeight(formula=list(C=4, O=4, H=9, N=1))
  mWPHENOL = MolecularWeight(formula=list(C=6, O=1, H=6))
  mWH2O2 = MolecularWeight(formula=list( O=2, H=2))
  mWMHP = MolecularWeight(formula=list( C=1,O=2, H=4))
  mWISOPN = MolecularWeight(formula=list( C=5,O=4, H=9, N=1))
  mWISOPOOH = MolecularWeight(formula=list( C=5,O=3, H=10))
  mWC6H10O5 = MolecularWeight(formula=list(C=6,H=10, O=5))
  mWNH4 = MolecularWeight(formula=list(N=1,H=4))
  mWSO4 = MolecularWeight(formula=list(S=1,O=4))
  mWNO3 = MolecularWeight(formula=list(N=1,O=3))
  
  nCCO2 = 1; nCCO = 1; nCCH4 = 1 ; nCOC=1
  
  colnames(fire.ERstoCO2) =c("mce" ,"mce_SCHWARZ","r_sq",
                             # DISKIN
                             "ER_CO2_DISKIN", "ER_CO_DISKIN" ,"ER_CH4_DISKIN", 
                             # RYERSON
                             "ER_NO2_RYERSON","ER_NO_RYERSON","ER_O3_RYERSON","ER_NOy_RYERSON",
                             # JIMENEZ
                             "ER_OC_JIMENEZ","ER_Sulfate_JIMENEZ",
                             "ER_Nitrate_JIMENEZ", "ER_Ammonium_JIMENEZ","ER_NR_Chloride_JIMENEZ",
                             "ER_Seasalt_JIMENEZ", "ER_MSA_JIMENEZ","ER_ClO4_JIMENEZ","ER_Bromine_JIMENEZ_ppb",
                             "ER_Iodine_JIMENEZ","ER_C6H10O5_JIMENEZ",
                             # SCHWARZ
                             "ER_BC_SCHWARZ", 
                             # WARNEKE
                             'ER_HCN_WARNEKE','ER_CH2O_WARNEKE','ER_CH3OH_WARNEKE','ER_CH3CN_WARNEKE','ER_HNCO_WARNEKE','ER_CH3CHO_WARNEKE','ER_C2H5OH_WARNEKE','ER_HCOOH_WARNEKE','ER_C3H3N_WARNEKE',
                             'ER_C3H4O_WARNEKE','ER_C3H6O_WARNEKE','ER_C2H4O2_WARNEKE','ER_CH3NO2_WARNEKE','ER_C2H6S_WARNEKE','ER_C4H5N_WARNEKE','ER_C4H4O_WARNEKE','ER_C4H6O_WARNEKE','ER_C4H8O_WARNEKE',
                             'ER_C3H6O2_WARNEKE',	'ER_C6H6_WARNEKE','ER_C5H6O_WARNEKE',	'ER_C4H4O2_WARNEKE', 'ER_C4H6O2_WARNEKE',	'ER_C7H8_WARNEKE','ER_C6H6O_WARNEKE','ER_C5H4O2_WARNEKE',	'ER_C6H8O_WARNEKE',	
                             'ER_C4H2O3_WARNEKE',	'ER_C7H5N_WARNEKE',	'ER_C8H8_WARNEKE', 'ER_C7H6O_WARNEKE','ER_C8H10_WARNEKE',	'ER_C7H8O_WARNEKE',	'ER_C6H6O2_WARNEKE','ER_C8H6O_WARNEKE',	'ER_C9H12_WARNEKE','ER_C6H4O3_WARNEKE',
                             'ER_C7H8O2_WARNEKE',	'ER_C10H8_WARNEKE',	'ER_C10H16_WARNEKE','ER_C8H10O2_WARNEKE',	'ER_C8H10O3_WARNEKE',
                             # FRIED
                             'ER_CH2O_CAMS','ER_C2H6_CAMS',
                             # HANISCO
                             'ER_CH2O_HANISCO',
                             # ROLLINS
                             'ER_SO2_ROLLINS','ER_NO_ROLLINS',
                             # WOMACK
                             'ER_NO2_WOMACK', 'ER_HNO2_WOMACK', 'ER_CH3COCHO_WOMACK', 'ER_CHOCHO_WOMACK',
                             # St Clair
                             'ER_NO2_STCLAIR',
                             # VERES
                             'ER_HNO2_NOAACIMS_VERES','ER_N2O5_NOAACIMS_VERES','ER_BrCN_NOAACIMS_VERES',
                             'ER_BrCl_NOAACIMS_VERES', 'ER_BrO_NOAACIMS_VERES', 'ER_CH3COOCl_NOAACIMS_VERES',
                             'ER_Cl2_NOAACIMS_VERES',  'ER_ClNO2_NOAACIMS_VERES', 'ER_HCN_NOAACIMS_VERES',
                             'ER_HCOOH_NOAACIMS_VERES', 'ER_HNCO_NOAACIMS_VERES','ER_HPMTF_NOAACIMS_VERES',
                             # WISTHALER
                             'ER_NH3_WISTHALER', 'ER_CH2O_WISTHALER',
                             # WENNBERG
                             'ER_BUTENE.HN_WENNBERG',  'ER_BUTENE.HP_WENNBERG', 'ER_ETHENE.HN_WENNBERG','ER_ETHENE.HP_WENNBERG',
                             'ER_PHENOL_WENNBERG', 'ER_H2O2_WENNBERG', 'ER_MHP_WENNBERG', 'ER_GLYC_WENNBERG', 'ER_IEPOX_WENNBERG',
                             'ER_ISOPN_WENNBERG', 'ER_ISOPOOH_WENNBERG','ER_PAA_WENNBERG', 'ER_PROPENE.HN_WENNBERG','ER_PROPENE.HP_WENNBERG',
                             'ER_HCN_WENNBERG','ER_HNO3_WENNBERG',
                             # HUEY
                             'ER_PAN_HUEY','ER_PPN_HUEY','ER_APAN_HUEY','ER_PBN_HUEY',
                             # BLAKE
                             'ER_OCS_WAS_BLAKE','ER_DMS_WAS_BLAKE','ER_CFC12_WAS_BLAKE','ER_CFC11_WAS_BLAKE','ER_CFC113_WAS_BLAKE','ER_CFC114_WAS_BLAKE', 
                             'ER_HFC152a_WAS_BLAKE','ER_HFC134a_WAS_BLAKE','ER_HFC365mfc_WAS_BLAKE','ER_HCFC22_WAS_BLAKE','ER_HCFC142b_WAS_BLAKE','ER_HCFC141b_WAS_BLAKE', 
                             'ER_H1301_WAS_BLAKE','ER_H2402_WAS_BLAKE','ER_H1211_WAS_BLAKE','ER_CH3CCl3_WAS_BLAKE','ER_CCl4_WAS_BLAKE','ER_CHCl3_WAS_BLAKE','ER_CH2Cl2_WAS_BLAKE', 
                             'ER_C2HCl3_WAS_BLAKE','ER_C2Cl4_WAS_BLAKE','ER_CH3Cl_WAS_BLAKE','ER_CH3Br_WAS_BLAKE','ER_CH3I_WAS_BLAKE','ER_CH2Br2_WAS_BLAKE','ER_CHBrCl2_WAS_BLAKE', 
                             'ER_CHBr2Cl_WAS_BLAKE','ER_CHBr3_WAS_BLAKE','ER_CH2ClCH2Cl_WAS_BLAKE','ER_C2H5Cl_WAS_BLAKE','ER_MeONO2_WAS_BLAKE','ER_EthONO2_WAS_BLAKE', 
                             'ER_iPropONO2_WAS_BLAKE','ER_nPropONO2_WAS_BLAKE','ER_x2ButONO2_WAS_BLAKE','ER_x3PentONO2_WAS_BLAKE','ER_x2PentONO2_WAS_BLAKE', 
                             'ER_x3Me2ButONO2_WAS_BLAKE','ER_Ethane_WAS_BLAKE','ER_Ethene_WAS_BLAKE','ER_Ethyne_WAS_BLAKE','ER_Propene_WAS_BLAKE','ER_Propane_WAS_BLAKE', 
                             'ER_Propadiene_WAS_BLAKE','ER_Propyne_WAS_BLAKE','ER_iButane_WAS_BLAKE','ER_nButane_WAS_BLAKE','ER_x1Butene_WAS_BLAKE','ER_iButene_WAS_BLAKE', 
                             'ER_t2Butene_WAS_BLAKE','ER_c2Butene_WAS_BLAKE','ER_x13Butadiene_WAS_BLAKE','ER_x12Butadiene_WAS_BLAKE','ER_x1Buten3yne_WAS_BLAKE', 
                             'ER_x13Butadyine_WAS_BLAKE','ER_x1Butyne_WAS_BLAKE','ER_x2Butyne_WAS_BLAKE','ER_iPentane_WAS_BLAKE','ER_nPentane_WAS_BLAKE', 
                             'ER_Isoprene_WAS_BLAKE','ER_x1Pentene_WAS_BLAKE','ER_t2Pentene_WAS_BLAKE','ER_c2Pentene_WAS_BLAKE','ER_3Me1Butene_WAS_BLAKE', 
                             'ER_2Me1Butene_WAS_BLAKE','ER_2Me2Butene_WAS_BLAKE','ER_x13Pentadienes_WAS_BLAKE','ER_x3Me1PenteneAnd4Me1Pentene_WAS_BLAKE', 
                             'ER_x1Hexene_WAS_BLAKE','ER_x1Heptene_WAS_BLAKE','ER_x1Octene_WAS_BLAKE','ER_x1Nonene_WAS_BLAKE','ER_x1Decene_WAS_BLAKE','ER_nHexane_WAS_BLAKE', 
                             'ER_nHeptane_WAS_BLAKE','ER_nOctane_WAS_BLAKE','ER_nNonane_WAS_BLAKE','ER_nDecane_WAS_BLAKE','ER_nUndecane_WAS_BLAKE','ER_x22Dimebutane_WAS_BLAKE', 
                             'ER_x23Dimebutane_WAS_BLAKE','ER_x2MePentane_WAS_BLAKE','ER_x3MePentane_WAS_BLAKE','ER_x2MeHexane_WAS_BLAKE','ER_x3MeHexane_WAS_BLAKE', 
                             'ER_x23DimePentane_BLAKE','ER_x224TrimePentane_WAS_BLAKE','ER_x234TrimePentane_WAS_BLAKE','ER_CycPentane_WAS_BLAKE', 
                             'ER_MeCycPentane_WAS_BLAKE','ER_CycHexane_WAS_BLAKE','ER_MeCycHexane_WAS_BLAKE','ER_CycPentene_WAS_BLAKE','ER_Benzene_WAS_BLAKE', 
                             'ER_Toluene_WAS_BLAKE','ER_EthBenzene_WAS_BLAKE','ER_mpXylene_WAS_BLAKE','ER_oXylene_WAS_BLAKE','ER_Styrene_WAS_BLAKE', 
                             'ER_EthynylBenzene_WAS_BLAKE','ER_iPropBenzene_WAS_BLAKE','ER_nPropBenzene_WAS_BLAKE','ER_x3EthToluene_WAS_BLAKE','ER_x4EthToluene_WAS_BLAKE', 
                             'ER_x2EthToluene_WAS_BLAKE','ER_x135rimeBenzene_WAS_BLAKE','ER_x124rimeBenzene_WAS_BLAKE','ER_ClBenzene_WAS_BLAKE','ER_aPinene_WAS_BLAKE', 
                             'ER_bPinene_WAS_BLAKE','ER_Tricyclene_WAS_BLAKE','ER_Camphene_WAS_BLAKE','ER_Myrcene_WAS_BLAKE','ER_Limonene_WAS_BLAKE','ER_Furan_WAS_BLAKE', 
                             'ER_x2MeFuran_WAS_BLAKE','ER_x3MeFuran_WAS_BLAKE','ER_BenzFuran_WAS_BLAKE','ER_iButanal_WAS_BLAKE','ER_Butanal_WAS_BLAKE', 
                             'ER_AcetonePropanal_WAS_BLAKE','ER_MEK_WAS_BLAKE','ER_MAC_WAS_BLAKE','ER_MVK_WAS_BLAKE','ER_Acrolein_WAS_BLAKE','ER_iPropanol_WAS_BLAKE', 
                             'ER_Nitromethane_WAS_BLAKE','ER_Acrylonitrile_WAS_BLAKE','ER_PropNitrile_WAS_BLAKE', 'ER_MeAcetate_WAS_BLAKE',
                             # APEL
                             'ER_HFC134a_TOGA_APEL', 'ER_HCFC141b_TOGA_APEL','ER_HCFC142b_TOGA_APEL','ER_HCFC22_TOGA_APEL','ER_CH2Cl2_TOGA_APEL', 
                             'ER_CHCl3_TOGA_APEL','ER_CH2ClCH2Cl_TOGA_APEL','ER_CH3CCl3_TOGA_APEL','ER_C2Cl4_TOGA_APEL','ER_ClBenzene_TOGA_APEL',
                             'ER_CHBrCl2_TOGA_APEL','ER_CHBr2Cl_TOGA_APEL','ER_CH3Br_TOGA_APEL','ER_CH2Br2_TOGA_APEL','ER_CHBr3_TOGA_APEL',
                             'ER_CH2ClI_TOGA_APEL','ER_CH3I_TOGA_APEL','ER_CS2_TOGA_APEL','ER_CH3SH_TOGA_APEL','ER_DMS_TOGA_APEL','ER_Propane_TOGA_APEL',
                             'ER_iButane_TOGA_APEL','ER_nButane_TOGA_APEL','ER_iPentane_TOGA_APEL','ER_nPentane_TOGA_APEL','ER_x2MePentane_TOGA_APEL',
                             'ER_x3MePentane_TOGA_APEL','ER_nHexane_TOGA_APEL','ER_x224TrimePentane_TOGA_APEL','ER_nHeptane_TOGA_APEL',
                             'ER_nOctane_TOGA_APEL','ER_Propene_TOGA_APEL','ER_iButene1Butene_TOGA_APEL','ER_Isoprene_TOGA_APEL', 
                             'ER_Tricyclene_TOGA_APEL','ER_aPinene_TOGA_APEL','ER_Camphene_TOGA_APEL','ER_bPineneMyrcene_TOGA_APEL', 
                             'ER_LimoneneD3Carene_TOGA_APEL','ER_Benzene_TOGA_APEL','ER_Toluene_TOGA_APEL','ER_EthBenzene_TOGA_APEL', 
                             'ER_mpXylene_TOGA_APEL','ER_oXylene_TOGA_APEL','ER_Styrene_TOGA_APEL','ER_EthynylBenzene_TOGA_APEL','ER_CH2O_TOGA_APEL', 
                             'ER_CH3CHO_TOGA_APEL','ER_Propanal_TOGA_APEL','ER_Butanal_TOGA_APEL','ER_iButanal_TOGA_APEL','ER_Acrolein_TOGA_APEL', 
                             'ER_x2Butenals_TOGA_APEL','ER_Acetone_TOGA_APEL','ER_MEK_TOGA_APEL','ER_CH3OH_TOGA_APEL','ER_C2H5OH_TOGA_APEL', 
                             'ER_iPropanol_TOGA_APEL','ER_MBO_TOGA_APEL','ER_MAC_TOGA_APEL','ER_MVK_TOGA_APEL','ER_MeFormate_TOGA_APEL', 
                             'ER_HCN_TOGA_APEL','ER_CH3CN_TOGA_APEL','ER_PropNitrile_TOGA_APEL','ER_Acrylonitrile_TOGA_APEL',
                             'ER_MeAcrylonitrile_TOGA_APEL','ER_Pyrrole_TOGA_APEL','ER_Nitromethane_TOGA_APEL','ER_MeONO2_TOGA_APEL', 
                             'ER_EthONO2_TOGA_APEL','ER_iPropONO2_TOGA_APEL','ER_x2ButONO2iButONO2_TOGA_APEL',
                             # Gilman
                             'ER_C2Cl4_NOAAiWAS_GILMAN','ER_CHCl3_NOAAiWAS_GILMAN','ER_Ethane_NOAAiWAS_GILMAN','ER_Propane_NOAAiWAS_GILMAN', 
                             'ER_nButane_NOAAiWAS_GILMAN','ER_iButane_NOAAiWAS_GILMAN','ER_nPentane_NOAAiWAS_GILMAN','ER_iPentane_NOAAiWAS_GILMAN',
                             'ER_nHexane_NOAAiWAS_GILMAN','ER_x2MePentane_NOAAiWAS_GILMAN','ER_x3MePentane_NOAAiWAS_GILMAN',
                             'ER_x22DiMeButane_NOAAiWAS_GILMAN','ER_x24DiMePentane_NOAAiWAS_GILMAN','ER_nOctane_NOAAiWAS_GILMAN',
                             'ER_x224TriMePentane_NOAAiWAS_GILMAN','ER_nNonane_NOAAiWAS_GILMAN','ER_nDecane_NOAAiWAS_GILMAN', 
                             'ER_MeCycPentane_NOAAiWAS_GILMAN','ER_CycHexane_NOAAiWAS_GILMAN','ER_MeCycHexane_NOAAiWAS_GILMAN', 
                             'ER_Ethyne_NOAAiWAS_GILMAN','ER_Ethene_NOAAiWAS_GILMAN','ER_Propene_NOAAiWAS_GILMAN','ER_x1Butene_NOAAiWAS_GILMAN',
                             'ER_c2Butene_NOAAiWAS_GILMAN','ER_t2butene_NOAAiWAS_GILMAN','ER_iButene_NOAAiWAS_GILMAN','ER_x1Pentene_NOAAiWAS_GILMAN', 
                             'ER_c2Pentene_NOAAiWAS_GILMAN','ER_t2Pentene_NOAAiWAS_GILMAN','ER_x2Me1Butene_NOAAiWAS_GILMAN','ER_x3Me1Butene_NOAAiWAS_GILMAN', 
                             'ER_t13Pentadiene_NOAAiWAS_GILMAN','ER_Isoprene_NOAAiWAS_GILMAN','ER_aPinene_NOAAiWAS_GILMAN','ER_Benzene_NOAAiWAS_GILMAN', 
                             'ER_Toluene_NOAAiWAS_GILMAN','ER_EthBenzene_NOAAiWAS_GILMAN','ER_oXylene_NOAAiWAS_GILMAN','ER_mpXylene_NOAAiWAS_GILMAN', 
                             'ER_Acetone_NOAAiWAS_GILMAN','ER_MEK_NOAAiWAS_GILMAN','ER_MeFormate_NOAAiWAS_GILMAN','ER_Furan_NOAAiWAS_GILMAN', 
                             'ER_CH3CN_NOAAiWAS_GILMAN','ER_Acrylonitrile_NOAAiWAS_GILMAN',
                             # MaleicAnhydride to Furan ratio
                             'MAtoF',
                             'MAtoF_rsq',
                             # NOxtoNOy
                             'NOxtoNOy',
                             'NOxtoNOy_rsq',
                             # size of plume used to calculate EFs
                             'npts',
                             # smoke age
                             'AGE_HOLMES',
                             'AGE_SCHWARZ')
  # For now, just include CO2, CO, CH4, and OC
  fire.ERstoCO2$TC = fire.ERstoCO2$ER_CO2_DISKIN* nCCO2 + fire.ERstoCO2$ER_CO_DISKIN* nCCO + fire.ERstoCO2$ER_CH4_DISKIN *nCCH4 
  
  #fire.ERstoCO2$TC = fire.ERstoCO2$TC  +
  #  fire.ERstoCO2$ER_OC_JIMENEZ * nCOC# + ... other VOCs eventually
  
  fire.ERstoCO2$EF_CO2_DISKIN   = cContentCorn * mWCO2/12 * fire.ERstoCO2$ER_CO2_DISKIN /fire.ERstoCO2$TC
  fire.ERstoCO2$EF_CO_DISKIN   = cContentCorn * mWCO/12 * fire.ERstoCO2$ER_CO_DISKIN /fire.ERstoCO2$TC
  fire.ERstoCO2$EF_CH4_DISKIN   = cContentCorn * mWCH4/12 * fire.ERstoCO2$ER_CH4_DISKIN /fire.ERstoCO2$TC
  
  fire.ERstoCO2$EF_BC_SCHWARZ   = cContentCorn * mWOC/12 * fire.ERstoCO2$ER_BC_SCHWARZ/fire.ERstoCO2$TC
  
  fire.ERstoCO2$EF_OC_JIMENEZ   = cContentCorn * mWOC/12 * fire.ERstoCO2$ER_OC_JIMENEZ/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_NH4_JIMENEZ   = cContentCorn * mWNH4/12 * fire.ERstoCO2$ER_Ammonium_JIMENEZ/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_SO4_JIMENEZ   = cContentCorn * mWSO4/12 * fire.ERstoCO2$ER_Sulfate_JIMENEZ/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_NO3_JIMENEZ   = cContentCorn * mWNO3/12 * fire.ERstoCO2$ER_Nitrate_JIMENEZ/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C6H10O5_JIMENEZ   = cContentCorn * mWC6H10O5/12 * fire.ERstoCO2$ER_C6H10O5_JIMENEZ/fire.ERstoCO2$TC
  
  fire.ERstoCO2$EF_NO_RYERSON   = cContentCorn * mWNO/12 * fire.ERstoCO2$ER_NO_RYERSON/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_NO2_RYERSON  = cContentCorn * mWNO2/12 * fire.ERstoCO2$ER_NO2_RYERSON/fire.ERstoCO2$TC
  #fire.ERstoCO2$EF_NOy  = cContentCorn * mWNOy/12 * fire.ERstoCO2$ER_NOy/fire.ERstoCO2$TC
  
  # WARNEKE VOCs
  wMWs = c(27.0,30.0,32.0,41.1,43.0,44.0,46.1,46.0,55.1,56.1,58.1,60.0,61.0,62.1,67.1,68.1,70.1,74.1,74.1,78.1,82.1,84.1,86.1,92.1,94.1,96.1,96.1,98.1,103.1,104.1,106.1,106.2,108.1,110.1,118.1,120.2,124.1,124.1,128.2,136.2,138.1,154.1 )
  
  fire.ERstoCO2$EF_HCN_WARNEKE    = cContentCorn * wMWs[1]/12 * fire.ERstoCO2$ER_HCN_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_CH2O_WARNEKE   = cContentCorn * wMWs[2]/12 * fire.ERstoCO2$ER_CH2O_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_CH3OH_WARNEKE  = cContentCorn * wMWs[3]/12 * fire.ERstoCO2$ER_CH3OH_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_CH3CN_WARNEKE  = cContentCorn * wMWs[4]/12 * fire.ERstoCO2$ER_CH3CN_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_HNCO_WARNEKE   = cContentCorn * wMWs[5]/12 * fire.ERstoCO2$ER_HNCO_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_CH3CHO_WARNEKE = cContentCorn * wMWs[6]/12 * fire.ERstoCO2$ER_CH3CHO_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C2H5OH_WARNEKE = cContentCorn * wMWs[7]/12 * fire.ERstoCO2$ER_C2H5OH_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_HCOOH_WARNEKE  = cContentCorn * wMWs[8]/12 * fire.ERstoCO2$ER_HCOOH_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C3H3N_WARNEKE  = cContentCorn * wMWs[9]/12 * fire.ERstoCO2$ER_C3H3N_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C3H4O_WARNEKE  = cContentCorn * wMWs[10]/12 * fire.ERstoCO2$ER_C3H4O_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C3H6O_WARNEKE  = cContentCorn * wMWs[11]/12 * fire.ERstoCO2$ER_C3H6O_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C2H4O2_WARNEKE = cContentCorn * wMWs[12]/12 * fire.ERstoCO2$ER_C2H4O2_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_CH3NO2_WARNEKE = cContentCorn * wMWs[13]/12 * fire.ERstoCO2$ER_CH3NO2_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C2H6S_WARNEKE  = cContentCorn * wMWs[14]/12 * fire.ERstoCO2$ER_C2H6S_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C4H5N_WARNEKE  = cContentCorn * wMWs[15]/12 * fire.ERstoCO2$ER_C4H5N_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C4H4O_WARNEKE  = cContentCorn * wMWs[16]/12 * fire.ERstoCO2$ER_C4H4O_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C4H6O_WARNEKE  = cContentCorn * wMWs[17]/12 * fire.ERstoCO2$ER_C4H6O_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C4H8O_WARNEKE  = cContentCorn * wMWs[18]/12 * fire.ERstoCO2$ER_C4H8O_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C3H6O2_WARNEKE = cContentCorn * wMWs[19]/12 * fire.ERstoCO2$ER_C3H6O2_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C6H6_WARNEKE   = cContentCorn * wMWs[20]/12 * fire.ERstoCO2$ER_C6H6_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C5H6O_WARNEKE  = cContentCorn * wMWs[21]/12 * fire.ERstoCO2$ER_C5H6O_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C4H4O2_WARNEKE = cContentCorn * wMWs[22]/12 * fire.ERstoCO2$ER_C4H4O2_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C4H6O2_WARNEKE = cContentCorn * wMWs[23]/12 * fire.ERstoCO2$ER_C4H6O2_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C7H8_WARNEKE  = cContentCorn * wMWs[24]/12 * fire.ERstoCO2$ER_C7H8_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C6H6O_WARNEKE  = cContentCorn * wMWs[25]/12 * fire.ERstoCO2$ER_C6H6O_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C5H4O2_WARNEKE = cContentCorn * wMWs[26]/12 * fire.ERstoCO2$ER_C5H4O2_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C6H8O_WARNEKE  = cContentCorn * wMWs[27]/12 * fire.ERstoCO2$ER_C6H8O_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C4H2O3_WARNEKE = cContentCorn * wMWs[28]/12 * fire.ERstoCO2$ER_C4H2O3_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C7H5N_WARNEKE  = cContentCorn * wMWs[29]/12 * fire.ERstoCO2$ER_C7H5N_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C8H8_WARNEKE   = cContentCorn * wMWs[30]/12 * fire.ERstoCO2$ER_C8H8_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C7H6O_WARNEKE  = cContentCorn * wMWs[31]/12 * fire.ERstoCO2$ER_C7H6O_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C8H10_WARNEKE  = cContentCorn * wMWs[32]/12 * fire.ERstoCO2$ER_C8H10_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C7H8O_WARNEKE  = cContentCorn * wMWs[33]/12 * fire.ERstoCO2$ER_C7H8O_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C6H6O2_WARNEKE = cContentCorn * wMWs[34]/12 * fire.ERstoCO2$ER_C6H6O2_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C8H6O_WARNEKE  = cContentCorn * wMWs[35]/12 * fire.ERstoCO2$ER_C8H6O_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C9H12_WARNEKE  = cContentCorn * wMWs[36]/12 * fire.ERstoCO2$ER_C9H12_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C6H4O3_WARNEKE = cContentCorn * wMWs[37]/12 * fire.ERstoCO2$ER_C6H4O3_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C7H8O2_WARNEKE = cContentCorn * wMWs[38]/12 * fire.ERstoCO2$ER_C7H8O2_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C10H8_WARNEKE  = cContentCorn * wMWs[39]/12 * fire.ERstoCO2$ER_C10H8_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C10H16_WARNEKE = cContentCorn * wMWs[40]/12 * fire.ERstoCO2$ER_C10H16_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C8H10O2_WARNEKE= cContentCorn * wMWs[41]/12 * fire.ERstoCO2$ER_C8H10O2_WARNEKE/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C8H10O3_WARNEKE= cContentCorn * wMWs[42]/12 * fire.ERstoCO2$ER_C8H10O3_WARNEKE/fire.ERstoCO2$TC
  # CAMS
  fire.ERstoCO2$EF_CH2O_CAMS = cContentCorn * mWCH2O/12 * fire.ERstoCO2$ER_CH2O_CAMS/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_C2H6_CAMS = cContentCorn * mWC2H6/12 * fire.ERstoCO2$ER_C2H6_CAMS/fire.ERstoCO2$TC
  # HANISCO
  fire.ERstoCO2$EF_CH2O_HANISCO = cContentCorn * mWCH2O/12 * fire.ERstoCO2$ER_CH2O_HANISCO/fire.ERstoCO2$TC
  # ROLLINS
  fire.ERstoCO2$EF_SO2_ROLLINS = cContentCorn * mWSO2/12 * fire.ERstoCO2$ER_SO2_ROLLINS/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_NO_ROLLINS  = cContentCorn * mWNO/12  * fire.ERstoCO2$ER_NO_ROLLINS/fire.ERstoCO2$TC
  #WOMACK
  fire.ERstoCO2$EF_NO2_WOMACK      = cContentCorn*mWNO2/12 * fire.ERstoCO2$ER_NO2_WOMACK/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_HNO2_WOMACK     = cContentCorn*mWHONO/12 * fire.ERstoCO2$ER_HNO2_WOMACK/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_CHOCHO_WOMACK   = cContentCorn*mWCHOCHO/12 * fire.ERstoCO2$ER_CHOCHO_WOMACK/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_CH3COCHO_WOMACK = cContentCorn*mWCH3COCHO/12 * fire.ERstoCO2$ER_CH3COCHO_WOMACK/fire.ERstoCO2$TC
  # St Clair
  fire.ERstoCO2$EF_NO2_STCLAIR     = cContentCorn*mWNO2/12 * fire.ERstoCO2$ER_NO2_STCLAIR/fire.ERstoCO2$TC
  # VERES
  fire.ERstoCO2$EF_HNO2_NOAACIMS_VERES    = cContentCorn*mWHONO/12 * fire.ERstoCO2$ER_HNO2_NOAACIMS_VERES/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_N2O5_NOAACIMS_VERES    = cContentCorn*mWN2O5/12  * fire.ERstoCO2$ER_N2O5_NOAACIMS_VERES/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_BrCN_NOAACIMS_VERES    = cContentCorn*mWBrCN/12 * fire.ERstoCO2$ER_BrCN_NOAACIMS_VERES/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_BrCl_NOAACIMS_VERES    = cContentCorn*mWBrCl/12 * fire.ERstoCO2$ER_BrCl_NOAACIMS_VERES/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_BrO_NOAACIMS_VERES     = cContentCorn*mWBrO/12  * fire.ERstoCO2$ER_BrO_NOAACIMS_VERES/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_CH3COOCl_NOAACIMS_VERES= cContentCorn*mWCH3COOCl/12 * fire.ERstoCO2$ER_CH3COOCl_NOAACIMS_VERES/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_Cl2_NOAACIMS_VERES    = cContentCorn*mWCl2/12   * fire.ERstoCO2$ER_Cl2_NOAACIMS_VERES/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_ClNO2_NOAACIMS_VERES  = cContentCorn*mWClNO2/12 * fire.ERstoCO2$ER_ClNO2_NOAACIMS_VERES/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_HCN_NOAACIMS_VERES    = cContentCorn*mWHCN/12   * fire.ERstoCO2$ER_HCN_NOAACIMS_VERES/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_HCOOH_NOAACIMS_VERES  = cContentCorn*mWHCOOH/12 * fire.ERstoCO2$ER_HCOOH_NOAACIMS_VERES/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_HNCO_NOAACIMS_VERES   = cContentCorn*mWHNCO/12  * fire.ERstoCO2$ER_HNCO_NOAACIMS_VERES/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_HPMTF_NOAACIMS_VERES  = cContentCorn*mWHPMTF/12 * fire.ERstoCO2$ER_HPMTF_NOAACIMS_VERES/fire.ERstoCO2$TC
  
  # WISTHALER
  fire.ERstoCO2$EF_NH3_WISTHALER = cContentCorn * mWNH3/12 * fire.ERstoCO2$ER_NH3_WISTHALER/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_CH2O_WISTHALER = cContentCorn * mWCH2O/12 * fire.ERstoCO2$ER_CH2O_WISTHALER/fire.ERstoCO2$TC
  
  # WENNBERG
  fire.ERstoCO2$EF_BUTENE.HN_WENNBERG = cContentCorn * mWBUTENEHN/12 * fire.ERstoCO2$ER_BUTENE.HN_WENNBERG/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_BUTENE.HP_WENNBERG = cContentCorn * mWBUTENEHP/12 * fire.ERstoCO2$ER_BUTENE.HP_WENNBERG/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_ETHENE.HN_WENNBERG = cContentCorn * mWETHENEHN/12 * fire.ERstoCO2$ER_ETHENE.HN_WENNBERG/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_ETHENE.HP_WENNBERG = cContentCorn * mWETHENEHP/12 * fire.ERstoCO2$ER_ETHENE.HP_WENNBERG/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_PHENOL_WENNBERG = cContentCorn * mWPHENOL/12 * fire.ERstoCO2$ER_PHENOL_WENNBERG/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_H2O2_WENNBERG = cContentCorn * mWH2O2/12 *  fire.ERstoCO2$ER_H2O2_WENNBERG/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_MHP_WENNBERG = cContentCorn * mWMHP/12 * fire.ERstoCO2$ER_MHP_WENNBERG/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_GLYC_WENNBERG = cContentCorn * mWGLYC/12 * fire.ERstoCO2$ER_GLYC_WENNBERG/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_IEPOX_WENNBERG = cContentCorn * mWIEPOX/12 * fire.ERstoCO2$ER_IEPOX_WENNBERG/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_ISOPN_WENNBERG = cContentCorn * mWISOPN/12 * fire.ERstoCO2$ER_ISOPN_WENNBERG/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_ISOPOOH_WENNBERG = cContentCorn * mWISOPOOH/12 * fire.ERstoCO2$ER_ISOPOOH_WENNBERG/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_PAA_WENNBERG = cContentCorn * mWPAA/12 * fire.ERstoCO2$ER_PAA_WENNBERG/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_PROPENE.HN_WENNBERG = cContentCorn * mWPROPENEHN/12 * fire.ERstoCO2$ER_PROPENE.HN_WENNBERG/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_PROPENE.HP_WENNBERG = cContentCorn * mWPROPENEHP/12 * fire.ERstoCO2$ER_PROPENE.HP_WENNBERG/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_HCN_WENNBERG = cContentCorn * mWHCN/12 * fire.ERstoCO2$ER_HCN_WENNBERG/fire.ERstoCO2$TC
  fire.ERstoCO2$EF_HNO3_WENNBERG = cContentCorn * mWHNO3/12 * fire.ERstoCO2$ER_HNO3_WENNBERG/fire.ERstoCO2$TC
  #fire.ERstoCO2[2,365:453] = NaN
  return(fire.ERstoCO2)
}

fuelvsfuel = function(allBOTH,fuel,variable){
  ind = which(allBOTH$fuel == fuel & allBOTH$variable == variable)
  ind2 = which(allBOTH$fuel == 'rice'& allBOTH$variable == variable)
  ttcornrice=t.test(allBOTH$mce.5hz[ind], allBOTH$mce.5hz[ind2])
  ind2 = which(allBOTH$fuel == 'soybean'& allBOTH$variable == variable)
  ttcornsoy=t.test(allBOTH$mce.5hz[ind], allBOTH$mce.5hz[ind2])
  ind2 = which(allBOTH$fuel == 'grass'& allBOTH$variable == variable)
  ttcorngrass=t.test(allBOTH$mce.5hz[ind], allBOTH$mce.5hz[ind2])
  ind2 = which(allBOTH$fuel == 'slash' & allBOTH$variable == variable)
  ttcornslash=t.test(allBOTH$mce.5hz[ind], allBOTH$mce.5hz[ind2])
  ind2 = which(allBOTH$fuel == 'pile'& allBOTH$variable == variable)
  ttcornpile=t.test(allBOTH$mce.5hz[ind], allBOTH$mce.5hz[ind2])
  ind2 = which(allBOTH$fuel == 'corn'& allBOTH$variable == variable)
  ttcorncorn=t.test(allBOTH$mce.5hz[ind], allBOTH$mce.5hz[ind2])
  cornvsfuel = as.data.frame(c(ttcornrice$p.value, ttcornsoy$p.value, ttcorngrass$p.value, 
                               ttcornslash$p.value, ttcornpile$p.value, ttcorncorn$p.value))
  colnames(cornvsfuel) = fuel
  rownames(cornvsfuel) = c("rice","soy","grass","slash","pile","corn")
  return(cornvsfuel)
}

makeEREFplot = function(alldata,var,specA,specAD){
  require(ggplot2)
  require(ggpubr)
  indA = which(akagi$Species == specA)
  tEF = ggplot(alldata) + 
    theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
    geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, size=fuel))+
    # geom_point(data=tmp.avg,aes(x=mce.5hz,y=FinalEF,  shape=Group.1),col='grey', size=8, stroke=1)+
    ylab("NO EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
    scale_color_manual(values = c(cbp1 ), limits=c("corn","soybean","rice","winter wheat","grass","pile","slash"))+
    scale_shape_manual(values = c(19,15,17,18,7,8,4), limits=c("corn","soybean","rice","winter wheat","grass","pile","slash"))+
    theme(legend.background=element_blank())+theme(legend.position = c(0.19, 0.7))+
    theme(text = element_text(size = sz)) +
    geom_point(data=xiaoxi, aes(x=MCE, y=NO), pch=0,col='green',stroke=2)+
    labs(col="",shape="")+
    # annotate(geom="text", x=0.83, y=0.8, label="Liu et al., 2016",size=6, col='green')+
    geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)+
    geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
    geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)
  
  
  tER = ggplot(alldata) + 
    theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
    geom_point(aes(x=mce.5hz, y=FinalERtoCO*1E3, col=fuel, size=fuel))+
    # geom_point(data=tmp.avg,aes(x=mce.5hz,y=FinalEF,  shape=Group.1),col='grey', size=8, stroke=1)+
    ylab("NO ER, ppt/ppb CO") + xlab("MCE")  + xlim(c(0.8, 1))+
    scale_color_manual(values = c(cbp1 ), limits=c("corn","soybean","rice","winter wheat","grass","pile","slash"))+
    scale_shape_manual(values = c(19,15,17,18,7,8,4), limits=c("corn","soybean","rice","winter wheat","grass","pile","slash"))+
    theme(legend.background=element_blank())+theme(legend.position = "none")+
    theme(text = element_text(size = sz)) +
    labs(col="",shape="")
  tmp = ggarrange(tEF,tER)
  return(tmp)
}
