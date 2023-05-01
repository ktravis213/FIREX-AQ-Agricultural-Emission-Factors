getICARTTdata = function(ff, year, month, day){
  tt = readLines(ff,n = 1)
  tt = strsplit(tt,",")
  tt = tt[[1]][1]
  tdata=read.csv(ff, header=TRUE, skip=as.numeric(tt)-1,sep=',')

  tdata$Year = year
  tdata$Month = month
  tdata$Day = day
  tmp = tdata$Time_Start/60/60
  tdata$Hour = floor(tmp)
  tmpM = (tmp - tdata$Hour)*60
  tdata$Min = floor(tmpM)
  tmpS = (tmpM - tdata$Min)*60
  tdata$Sec = round(tmpS)
  ind = which(tdata$Hour >= 24)
  if (length(ind) > 0){
    tdata$Hour[ind] = tdata$Hour[ind] - 24
    tdata$Day[ind] = tdata$Day[ind] + 1
  }
  #tdata$Date = getDates(tdata)
  tdata[tdata == -99999] = NaN
  tdata[tdata == -9999] = NaN
  tdata[tdata == -999] = NaN
  tdata[tdata ==-999.9] = NaN
  tdata[tdata == -999999] = NaN
  tdata[tdata ==-9999.99] = NaN
  tdata[tdata == -888888] = NaN
  namesF = colnames(tdata)
  # remove NaNs
  ind = which(namesF == "Smoke_flag_SCHWARZ")
  if (length(ind) > 0){
    ind = which(tdata$Smoke_flag_SCHWARZ == 1)
    ind2 = which(is.na(tdata$Smoke_flag_SCHWARZ))
    if (length(ind2) > 0){tdata$Smoke_flag_SCHWARZ[ind2] = 0}
  #  ind = which(mergedata.803$CO_DACOM_DISKIN > 500 & mergedata.803$Smoke_flag_SCHWARZ != 1) 
  #  ind = ind[19:133]
  #  mergedata.803$Smoke_flag_SCHWARZ[ind] = 1
    # count plume transects
    tdata$transect_plume_number = NaN
    start = 0
    out = 0
    for (i in 1:length(tdata$Smoke_flag_SCHWARZ)){
      if(tdata$Smoke_flag_SCHWARZ[i] == 1){
        tdata$transect_plume_number[i] = start
        if (tdata$Smoke_flag_SCHWARZ[i+1] != 1 ){ start = start + 1} # end
        i=i+1
        out = 0 # in a plume, set out to zero
      }
    }
  }
  doohr = 0
  if (doohr == 1){
    # ------------------- OH reactivity --------------------
    temp = tdata$Static_Air_Temp + 273.15 # convert to K
    nair =  (6.022E23 * tdata$Static_Pressure * 100)/(8.314 * temp)/1E6 # molecules air/cm3
    press = tdata$Static_Pressure
    # rates
    h2rate = 2.8E-12 * exp(-1800/temp)
  
    ch4rate = 2.45E-12 * exp(-1775.0/temp)
    tdata$ch4_hz = tdata$CH4_DACOM_DISKIN * nair*1E-9*ch4rate
    
    A0 = 1.5E-13
    B0 = 0.0
    C0 = 0.0
    R0 =  A0 * exp(C0/temp) * (300/temp)^(B0)
    R0 = R0 * (1 + 0.6*9.871E7*press)
    # new OH+CO rate from JPL2006.
    KLO1 = 5.9E-33*(300.0/temp)^(1.0)
    KHI1 = 1.1E-12*(300.0/temp)^(-1.3)
    XYRAT1 = KLO1*nair/KHI1
    BLOG1 = log10(XYRAT1)
    Fexp1 = 1.0/(1.0+BLOG1*BLOG1)
    KCO1=KLO1*nair * 0.6^Fexp1/(1.0+XYRAT1)
    KLO2=1.5E-13*(300.0/temp)^(0.0)
    KHI2=2.1E09 *(300.0/temp)^(-6.1)
    XYRAT2=KLO2*nair/KHI2
    BLOG2=log10(XYRAT2)
    Fexp2=1.0/(1.0+BLOG2*BLOG2)
    KCO2=KLO2*0.6^Fexp2/(1.0+XYRAT2)
    corate = KCO1+KCO2
    tdata$co_hz = tdata$CO_DACOM_DISKIN* nair*1E-9*corate
    
    o3rate =  1.7E-12 * exp(-940.0/temp)
    tdata$o3_hz = tdata$O3_CL_RYERSON* nair*1E-9*o3rate
    
    mohrate      =  2.9E-12 * exp(-345.0/temp)
    ind = which(namesF == "CH3OH_NOAAPTR_WARNEKE")
    if (length(ind) > 0){tdata$moh_hz = tdata$CH3OH_NOAAPTR_WARNEKE*nair*1e-9*mohrate}else{tdata$moh_hz = NaN}
    ald2rate     =  4.63E-12 * exp(350.0/temp)
    ind = which(namesF == "CH3CHO_NOAAPTR_WARNEKE")
    if (length(ind) > 0){tdata$ald2_hz = tdata$CH3CHO_NOAAPTR_WARNEKE*nair*1e-9*ald2rate}else{tdata$ald2_hz = NaN}
    hchorate     =  5.5E-12 * exp(125.0/temp)
    tdata$hcho_hz = tdata$CH2O_CAMS_pptv_FRIED*nair*1e-12*hchorate
    
    mprate =  3.8E-12 * exp(200.0/temp)  
    tdata$mp_hz = tdata$MHP_CIT_WENNBERG*nair*1e-12*mprate
    
    h2o2rate =  1.8E-12
    tdata$h2o2_hz = tdata$H2O2_CIT_firexaq.CIT.H2O2*nair*1E-12*h2o2rate
    
    #ho2rate =  4.8E-11 * exp(250./temp)
    #ohrate =  1.8E-12
    
    RLOW = ( 7.00E-31 * exp(0.0/temp) * (300/temp)^2.6) * nair
    RHIGH = (3.6E-11 * exp(0.0/temp) * (300/temp)^0.1)
    FV =  0.6
    XYRAT          = RLOW/RHIGH
    BLOG           = log10(XYRAT)
    Fexp           = 1.0/ (1.0 + BLOG * BLOG)
    norate         = RLOW*FV^Fexp/(1.0+XYRAT)
    tdata$no_hz = tdata$NO_CL_RYERSON*nair*1E-9*norate
    
    RLOW = ( 1.80E-30 * exp(0/temp) * (300/temp)^3.0) * nair
    RHIGH = (2.8E-11 * exp(0/temp) * (300/temp)^0.0)
    FV =  0.6
    XYRAT          = RLOW/RHIGH
    BLOG           = log10(XYRAT)
    Fexp           = 1.0/ (1.0 + BLOG * BLOG)
    no2rate        = RLOW*FV^Fexp/(1.0+XYRAT)
    tdata$no_hz = tdata$NO2_CL_RYERSON*nair*1E-9*no2rate
    
    A0 =  2.41E-14
    B0 =  0.0
    C0 =  460.0
    A1 =  2.69E-17
    B1 =  0.0
    C1 =  2199.0
    A2 =  6.51E-34
    B2 =  0.0
    C2 =  1335.
    R0 =  A0 * exp(C0/temp) * (300.0/temp)^(B0)
    R1 =  A1 * exp(C1/temp) * (300.0/temp)^(B1)
    R2 =  nair*(A2 * exp(C2/temp) * (300.0/temp)^(B2))
    hno3rate = R0 + R2/(1.0 + R2/R1)
    tdata$hno3_hz = tdata$HNO3_CIT_WENNBERG*nair*1E-12*hno3rate
    
    RLOW = ( 4.60E-27 * exp(0.0/temp) * (300/temp)^4.0) * nair
    RHIGH = (2.6E-11 * exp(0.0/temp) * (300/temp)^1.3)
    FV =  0.5
    XYRAT          = RLOW/RHIGH
    BLOG           = log10(XYRAT)
    Fexp           = 1.0/ (1.0 + BLOG * BLOG)
    prperate       = RLOW*FV^Fexp/(1.0+XYRAT)
    ind = which(namesF == "Propene_WAS_BLAKE")
    if (length(ind) > 0){tdata$prpe_hz= tdata$Propene_WAS_BLAKE*nair*1E-12*prperate}else{ tdata$prpe_hz=NaN}
    
    acetrate = 1.33E-13+3.82E-11*exp(-2000.0/temp)
    ind = which(namesF == "Acetone_TOGA_APEL")
    if (length(ind) > 0){tdata$acet_hz= tdata$Acetone_TOGA_APEL*nair*1E-12*acetrate} else{ tdata$acet_hz=NaN}
    
    c3h8rate =  7.6E-12 * exp(-585.0/temp)
    ind = which(namesF == "Propane_TOGA_APEL")
    if (length(ind) > 0){tdata$c3h8_hz= tdata$Propane_TOGA_APEL*nair*1E-12*c3h8rate} else{ tdata$c3h8_hz=NaN}
    
    c2h6rate =  7.66E-12 * exp(-1020.0/temp)
    tdata$c2h6_hz= tdata$C2H6_CAMS_pptv_FRIED*nair*1E-12*c2h6rate
    
    eohrate =  3.35E-12
    ind = which(namesF == "Ethanol_TOGA_APEL")
    if (length(ind) > 0){tdata$eoh_hz= tdata$Ethanol_TOGA_APEL*nair*1E-12*eohrate} else{ tdata$c3h8_hz=NaN}
    
    isoprate = 3.0E-11*exp(360.0/temp)
    ind = which(namesF == "Isoprene_NOAAPTR_WARNEKE")
    if (length(ind) > 0){tdata$isop_hz = tdata$Isoprene_NOAAPTR_WARNEKE*nair*1E-9*isoprate} else{tdata$isop_hz=NaN}
    
    rchorate =  6.0E-12 * exp(410.0/temp)
    
    alk4rate =  9.1E-12 * exp(-405.0/temp)
    
    actarate =  3.15E-14 * exp(920/temp)
    
    # dms
    dms1 =    1.2E-11 * exp(-280.0/temp)
    A0 = 8.20E-39
    B0 = 0.0E+00
    C0 = 5376.0
    A1 = 1.05E-5
    B1 = 0.0
    C1 = 3644.0
    R0 = A0 * exp((C0)/temp) * (300.0/temp)^(B0)
    R1 = A1 * exp((C1)/temp) * (300.0/temp)^(B1)
    dms2 = (R0*nair*0.2095)/(1+R1*0.2095)
    dmsrate =  c()
    for (i in 1:length(dms1) ){
      dmsrate =  c(dmsrate,max(c(dms1[i], dms2[i]), na.rm=TRUE))
    }
    ind = which(namesF == "DMS_WAS_BLAKE")
    if (length(ind)> 0){tdata$dms_hz =  tdata$DMS_WAS_BLAKE*nair*1E-12*dmsrate}else{tdata$dms_hz = NaN}
    
    hcoohrate = 4.00E-13
    tdata$hcooh_hz= tdata$HCOOH_NOAACIMS_VERES*nair*1E-12*hcoohrate 
    
    hno2rate = 1.8E-11 * exp(-390./temp)
    tdata$hno3_hz= tdata$HNO2_NOAACIMS_VERES*nair*1E-12*hno2rate
    
    hno4rate = 1.30E-12 * exp(380./temp)
    no3rate  = 2.2E-11
    isn1rate = 7.48E-12 * exp(410./temp)
    r4n2rate = 1.60E-12
    
    mekrate  = 1.3E-12 * exp(-25./temp)
    ind = which(namesF == "MEK_TOGA_APEL")
    if (length(ind)> 0){tdata$mek_hz= tdata$MEK_TOGA_APEL*nair*1E-12*mekrate} else{ tdata$c3h8_hz=NaN}
    
    npmnrate = 2.9E-11 
    ipmnrate = 3.0E-11
    imaerate = 3.48E-12
    glycrate = 8.00E-12
    tdata$glyc_hz = tdata$GLYC_CIT_WENNBERG*nair*1E-12*glycrate
    
    glyxrate = 3.10E-12 * exp(340./temp)
    
    mglyrate = 1.50E-11
    tdata$mgly_hz = tdata$CH3COCHO_ACES_ppbv_WOMACK*nair*1E-9*mglyrate
    
    hacrate  = 2.15E-12 * exp(305./temp) 
    tdata$hac_hz = tdata$HAC_CIT_WENNBERG * nair*1E-12*hacrate
    
    benzrate = 2.33E-12 * exp(-193.0/temp)
    ind = which(namesF == "Benzene_NOAAPTR_WARNEKE")
    if (length(ind) > 0){tdata$benz_hz= tdata$Benzene_NOAAPTR_WARNEKE*nair*1E-9*benzrate} else{tdata$benz_hz=NaN}
    
    tolurate = 1.81E-12 * exp(338.0/temp)
    ind = which(namesF == "Toluene_NOAAPTR_WARNEKE")
    if (length(ind) > 0){tdata$tolu_hz= tdata$Toluene_NOAAPTR_WARNEKE*nair*1E-9*tolurate} else{tdata$tolu_hz=NaN}
    
    xylerate =  2.31E-11
    ind = which(namesF == "oXylene_TOGA_APEL")
    if (length(ind) > 0){tdata$xyle_hz = (tdata$oXylene_TOGA_APEL + tdata$mpXylene_TOGA_APEL)*nair*1E-12*xylerate} else{tdata$xyle_hz = NaN}
    
    A0 =  3.30E-31
    B0 =  4.3
    C0 =  0.0
    A1 =  1.6E-012
    B1 = 0.0
    C1 = 0.0
    FV = 0.6
    XYRAT   = RLOW/RHIGH
    BLOG    = log10(XYRAT)
    Fexp    = 1.0 / (1.0 + BLOG * BLOG)
    so2rate = RLOW*FV^Fexp/(1.0+XYRAT)
    tdata$so2_hz = tdata$SO2_LIF_ROLLINS*nair*1E-12*so2rate
    
    riparate =   6.13E-12 * exp(200/temp)
    ripbrate =   4.14E-12 * exp(200/temp)
    ripdrate =   5.11E-12 * exp(200/temp)
    iepoxarate = 3.73E-11 * exp(-400/temp)
    iepoxbrate = 5.79E-11 * exp(-400/temp)
    iepoxdrate = 3.2E-11 * exp(-400/temp)
    tdata$rip_hz = tdata$C5O3H10_CIT_WENNBERG*nair*1E-12*iepoxarate # isoppooh + iepox
    
    lvocrate =   4.82E-011 * exp(-400/temp)
  
    hbrrate =  5.5E-12 * exp(200/temp)
    br2rate =  2.1E-11 * exp(240/temp)
    brorate =  1.70E-11 * exp(250/temp)
    chbr3rate  = 9.0E-13 * exp(-360/temp)
    ch2br2rate = 2.0E-12 * exp(-840./temp)
    ch3brrate  = 1.42E-12 * exp(-1150./temp)
    isopndrate = 1.20E-11 * exp(652./temp)
    isopnbrate = 2.40E-12 * exp(745.0/temp)
    
    tdata$isopn_hz = tdata$ISOPN_CIT_WENNBERG * nair*1E-12*isopndrate
    mvknrate   = 4.40E-13 * exp(380.0/temp)
    macrnrate  = 8.79E-13 * exp(380.0/temp)
    tdata$mvkn_hz = tdata$C4O5H7N_CIT_WENNBERG*nair*1E-12*macrnrate
    mobarate =  2.79E-11 * exp(380.0/temp)
    ethlnrate = 2.40E-12
    propnnrate = 6.7E-13
    mtparate = 1.21E-11 * exp(440.0/temp)
    limorate = 4.20E-11 * exp(401.0/temp)
    monitsrate =  4.80E-12
    moniturate =  7.29E-11
    mvkrate = 2.33E-12 * exp(-193.0/temp)
    macrrate = 1.81E-12 * exp(338.0/temp)
    ind = which(namesF == "MAC_TOGA_APEL")
    if (length(ind) > 0){tdata$macr_hz= tdata$MAC_TOGA_APEL*nair*1E-12*macrrate} else{tdata$macr_hz=NaN}
    
    A0 = 2.41E-14
    B0 = 0.0E+00
    C0 =  460.0
    A1 =  2.69E-17
    B1 =  0.0
    C1 =  2199.
    A2 =  6.51E-34
    B2 =  0.E0
    C2 =  1335.0
    R0 =  (A0) * exp(C0/temp) * (300.0/temp)^(B0)
    R1 =  (A1) * exp((C1)/temp) * (300.0/temp)^(B1)
    R2 =  nair*((A2) * exp((C2)/temp) * (300.0/temp)^(B2))
    honitrate = R0 + R2/(1.0 + R2/R1)
    
    hc187rate  =  1.40E-11
    hpaldrate  =  5.10E-11
    tdata$hpald_hz = tdata$C5O3H8_CIT_WENNBERG*nair*1e-12*hpaldrate
    
    phenolrate = 4.7E-13*exp(1220/temp)
    tdata$phenol_hz = tdata$PHENOL_CIT_WENNBERG*nair*1E-12*phenolrate
    
    cresolrate =  4.65E-11
    tdata$cresol_hz = tdata$CRESOL_CIT_WENNBERG*nair*1E-12*cresolrate
    
    #C2H2 GCKMT17      
    K170 = 5.0E-30*nair*(temp/300.)^(-1.5)
    K17I = 1.0E-12	
    KR17 = K170/K17I
    FC17 = 0.17*exp(-51.0/temp)+exp(-temp/204.)
    NC17 = 0.75-1.27*(log10(FC17))
    F17 = 10^(log10(FC17)/(1.0+(log10(KR17)/NC17)^2))
    c2h2rate = (K170*K17I*F17)/(K170+K17I)
    ind = which(namesF == "Ethyne_WAS_BLAKE")
    if (length(ind) > 0){tdata$c2h2_hz =  tdata$Ethyne_WAS_BLAKE*nair*1e-12*c2h2rate}else{tdata$c2h2_hz=NaN}
    
    #C2H4 GCKMT15()
    K150 = 8.6E-29*nair*(temp/300.0)^(-3.1)
    K15I = 9.0E-12*(temp/300)^(-0.85)
    KR15 = K150/K15I
    FC15 = 0.48
    NC15 = 0.75-1.27*(log10(FC15))
    F15 = 10^(log10(FC15)/(1+(log10(KR15)/NC15)^2.))
    c2h4rate = (K150*K15I)*F15/(K150+K15I)
    ind = which(namesF == "Ethene_WAS_BLAKE")
    if (length(ind) > 0){ tdata$c2h4_hz =  tdata$Ethene_WAS_BLAKE*nair*1e-12*c2h4rate}else{tdata$c2h2_hz=NaN}
    
    maprate    = 6.13E-13 * exp(200./temp)
    tdata$map_hz =  tdata$PAA_CIT_WENNBERG*nair*1e-12*maprate
      
    cl2rate    = 2.60E-12 * exp(-1100.0/temp)
    clorate    = 7.40E-12 * exp(270.0/temp)
    oclorate   = 1.40E-12 * exp(600.0/temp)
    hclrate    = 1.80E-12 * exp(-250.0/temp)
    hoclrate   = 3.00E-12 * exp(-500.0/temp)
    clno2rate  = 2.40E-12 * exp(-1250.0/temp)
    clno3rate  = 1.20E-12 * exp(-330.0/temp)
    ch3clrate  = 1.96E-12 * exp(-1200.0/temp)
    ch2cl2rate = 2.61E-12 * exp(-944.0/temp)
    chcl3rate  = 4.69E-12 * exp(-1134.0/temp)
    i2rate  = 1.80E-10
    hirate  = 3.00E-11
    hoirate = 5.00E-12
    ch3irate = 2.90E-12 * exp(-1100.0/temp)
    
    tmp = as.data.frame(cbind(tdata$co_hz, tdata$hac_hz,tdata$mp_hz,tdata$no_hz, tdata$ch4_hz, tdata$o3_hz,tdata$h2o2_hz, tdata$hcho_hz,
      tdata$hno3_hz,tdata$map_hz,tdata$mek_hz,tdata$rip_hz,tdata$so2_hz,tdata$acet_hz, tdata$ald2_hz, tdata$benz_hz, 
      tdata$eoh_hz,tdata$glyc_hz,tdata$hcooh_hz,tdata$hpald_hz,tdata$isop_hz, tdata$macr_hz, tdata$mgly_hz, 
      tdata$moh_hz, tdata$mvkn_hz,tdata$prpe_hz, tdata$tolu_hz, tdata$xyle_hz,tdata$c2h2_hz,tdata$c2h4_hz,
      tdata$c2h6_hz,tdata$c3h8_hz,tdata$isopn_hz,tdata$cresol_hz, tdata$phenol_hz, tdata$dms_hz))
    tdata$ohr = rowSums(tmp, na.rm=TRUE)
    tmp2 = cbind(tdata$hac_hz,tdata$hcho_hz,
                tdata$map_hz,tdata$mek_hz,tdata$rip_hz,tdata$acet_hz, tdata$ald2_hz, tdata$benz_hz, 
                tdata$eoh_hz,tdata$glyc_hz,tdata$hcooh_hz,tdata$hpald_hz,tdata$isop_hz, tdata$macr_hz, tdata$mgly_hz, 
                tdata$moh_hz, tdata$mvkn_hz,tdata$prpe_hz, tdata$tolu_hz, tdata$xyle_hz,tdata$c2h2_hz,tdata$c2h4_hz,
                tdata$c2h6_hz,tdata$c3h8_hz,tdata$isopn_hz,tdata$cresol_hz, tdata$phenol_hz)
    names = c("hac","hcho","map","mek","isopooh","acet","ald2","benz","eoh","glyc","hcooh",
              "hpald","isop","macr","mgly","moh","mvkn","prpe","tolu","xyle","c2h2","c2h4","c2h6","c3h8","isopn","cresol","phenol","dms")
    tmp3 = colMeans(tmp2,na.rm=TRUE)
    means = as.data.frame(cbind(names=names, means=round(tmp3,digits = 2)))
    tdata$nmvoc_ohr = rowSums(tmp2, na.rm=TRUE)
  }
  return(tdata)
}
getICARTTdataALL = function(ff){
  require(chron)
  tt = readLines(ff,n = 1)
  tt = strsplit(tt,",")
  tt = tt[[1]][1]
  tdata=read.csv(ff, header=TRUE, skip=as.numeric(tt)-1,sep=',')
  dDate= chron('2019-01-01', '00:00:00', format=c('y-m-d','h:m:s'))
  cc = colnames(tdata)
  ind = which( cc == 'Fractional_Day')
  if (length(ind) > 0){
    tdata$Date = dDate + tdata$Fractional_Day 
  } else{
    tdata$Date = dDate + tdata$Day_Of_Year + tdata$Time_Start/60/60/24
  }
  #tdata$Date = getDates(tdata)
  tdata[tdata == -99999] = NaN
  tdata[tdata == -9999] = NaN
  tdata[tdata == -999] = NaN
  tdata[tdata ==-999.9] = NaN
  tdata[tdata == -999999] = NaN
  tdata[tdata ==-9999.99] = NaN
  tdata[tdata == -888888] = NaN
  namesF = colnames(tdata)
  return(tdata)
}
getICARTTdataSIMPLE = function(ff){
  require(chron)
  tt = readLines(ff,n = 1)
  tt = strsplit(tt,",")
  tt = tt[[1]][1]
  tdata=read.csv(ff, header=TRUE, skip=as.numeric(tt)-1,sep=',')
  tdata[tdata == -99999] = NaN
  tdata[tdata == -9999] = NaN
  tdata[tdata == -999] = NaN
  tdata[tdata ==-999.9] = NaN
  tdata[tdata == -999999] = NaN
  tdata[tdata ==-9999.99] = NaN
  tdata[tdata == -888888] = NaN
  tdata[tdata == -8888.0] = NaN
  return(tdata)
}