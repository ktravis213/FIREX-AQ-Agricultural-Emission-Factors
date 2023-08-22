## ---------------------------
##
## Script name: CalculateEmissionFactorsv2.R 
##
## Purpose of script: Process FIREX-AQ DC8 planeflight files and calculate emission factors
##
## Author: Dr. Katherine Travis
##
## Date Created: 2022-09-29
##
## Copyright (c) Katherine Travis, 2022
## Email: katherine.travis@nasa.gov
##
## ---------------------------
##
## Notes:
##   Emission factors from biomass burning activities are calculated according to Equation 1 (Yokelson et al., 1999),
##   EF_i (g/kg)=F_c×1000(g/kg)×(〖MW〗_i (g))/(12(g))×C_i/C_T  ,								                       (1)
##       where EF_i is the mass (g) of a species (i) emitted per kg of fuel burned, F_c is the fuel carbon fraction, 
##      〖MW〗_i is the molecular weight of the species (i), C_i is the number of moles of the species (i), and
##       C_T is the total number of moles of emitted carbon. 
##   The value of C_T here is calculated according to Equation 2, where NC_i is the number of carbons in species i and (∆C_i)/(∆C_x ) is the emission ratio (ER) of species i to species X (here we use CO).  
##   C_i/C_T =((∆C_i)/(∆C_x ))/(∑(NC_i×(∆C_i)/(∆C_x )) ), 			                								        (2)
##  We assume that F_c is 41% for agricultural fuels, 46% for grass, and 51% for pile and slash fuels (Stockwell et al., 2014). 
##  The emitted carbon is assumed to be encompassed by CO2, CO, and CH4. 

##
## ---------------------------
# Location of FIREX Working directory
setwd('/Users/ktravis1/OneDrive - NASA/FIREX/FinalAnalysisForGithub/')
# ------------------------------------------------
# Required Packages ---------------------------------
require(ggmap) ; require(OrgMassSpecR); library(readxl) ;require(plyr) ; require(dplyr); require(ggpubr); require(ncdf4)
doclear=1
if (doclear == 1){ rm(list=ls()) } # clear all analysis
doFM =  0 # attempt to look at fuel moisture

source('getICARTTdata.R')
source('york_regression_function.R')
source('getGEO.R')
source('SpeciesProperties.R')
source('getERsv2.R')

doread =0; doplot = 0 # read in data or load Rdata?  Plot? 
# Species to calculate emission ratios to in addition to CO
xspecies='CO2_ppb'

# -------------------Read in flags ------
SLOW = 0 # don't relax R2 criteria, using fast data

# ------ Start stop times for the  Adapted from the original file, and modified to break apart
# ------- plumes as described in Travis et al.
file = 'InputFiles/FIREXAQ-FIREFLAG-TABULARDATA_Analysis_20190724_R9_thru20190905.xlsx'
flags  <- read_excel(file,  skip = 153)
#
if (doread == 1){
  # ----- [[[[[[[[[[[[[[[[[[[[ Aug 21st ]]]]]]]]]]]]]]]]]]]]] -------
  # --------#######--------- Get 1 Hz Data individual 8/21------#######---------
  # plume tags
  tags = getICARTTdataSIMPLE('InputFiles/firexaq-fire-Flags-1HZ_DC8_20190821_R9.ict') ;tags$Time_Start = tags$TIME_START
  # MET DATA
  met.821.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-MetNav_DC8_20190821_R1.ict')
  met.821.1hz = merge(met.821.1hz, tags, by='Time_Start')
  # CO2
  co2.821.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CO2-7000_DC8_20190821_R2.ict')
  # -------- DISKIN -----CO, CH4
  co.ch4.821.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-DACOM_DC8_20190821_R1.ict')
  # --------- WARNEKE ----  VOCs
  warneke.821.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-NOAAPTR-VOCs-1Hz_DC8_20190821_R3.ict')
  # ------ HANISCO - ISAF HCHO - merged to 5Hz from the online merge
  isaf.821.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-ISAF-CH2O-1Hz_DC8_20190821_R0.ict')
  #  ------- ROLLINS - SO2 and NO
  rollinsno.821.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-NO_DC8_20190821_R1.ict')
  rollinsno.821.1hz$Time_Start = rollinsno.821.1hz$time_mid
  rollinsso2.821.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-SO2_DC8_20190821_R1.ict')
  rollinsso2.821.1hz$Time_Start = rollinsso2.821.1hz$time_mid
  #  ----- WENNBERG - CIT VOCs - 
  cit.821.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190821_R0_CIT.ict')
   # ------ HUEY - GTCIMS PANs - not sure how to match up peaks here
  gtcims.821.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190821_R0_Huey.ict')
  # ------ RYERSON
  ryerson.A = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NO_DC8_20190821_R1.ict')
  ryerson.B = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NO2_DC8_20190821_R1.ict')
  ryerson.C = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NOy_DC8_20190821_R1.ict')
  ryerson.D = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-O3_DC8_20190821_R1.ict')
  ryerson.821.1hz = cbind(ryerson.A,ryerson.B,ryerson.C,ryerson.D) ; ryerson.821.1hz$Time_Start = ryerson.821.1hz$Time_start
  # ----- JIMENEZ ---
  jimenez.821.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-EESI_DC8_20190821_R1.ict')
  # ----- SCHWARZ ---
  schwarz.821.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-SP2-BC-1HZ_DC8_20190821_R2.ict')
  # ----- FREID ---
  freid.c2h6.821.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-C2H6_DC8_20190821_R3.ict')
  freid.ch2o.821.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CH2O_DC8_20190821_R3.ict')
  # ------ WOMACK ---
  womackA = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-CH3COCHO_DC8_20190821_R1.ict')
  womackB = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-CHOCHO_DC8_20190821_R1.ict')
  womackC = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-HNO2_DC8_20190821_R1.ict')
  womackD = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-NO2_DC8_20190821_R1.ict')
  womack.821.1hz = cbind(womackA, womackB, womackC, womackD)
  # -------St Clair
  stclair.821.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-CANOE-NO2_DC8_20190821_R0.ict')
  stclair.821.1hz$Time_Start = stclair.821.1hz$Time_start
  # ------- VERES
  veres.A = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-ClNO2_DC8_20190821_R0.ict')
  veres.B = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HCOOH_DC8_20190821_R1.ict')
  veres.C = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HNO2_DC8_20190821_R1.ict')
  veres.D = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-N2O5_DC8_20190821_R0.ict')
  veres.E = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HPMTF_DC8_20190821_R0.ict')
  veres.F = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-CH3COOCl_DC8_20190821_R0.ict')
  veres.G = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-Cl2_DC8_20190821_R0.ict')
  veres.H = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrCl_DC8_20190821_R0.ict')
  veres.I = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrCN_DC8_20190821_R0.ict')
  veres.J = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrO_DC8_20190821_R0.ict')
  veres.K = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HCN_DC8_20190821_R0.ict')
  veres.L = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HNCO_DC8_20190821_R0.ict')
  veres.821.1hz = cbind(veres.A,veres.B,veres.C,veres.D,veres.E,veres.F,veres.G,veres.H,veres.I,veres.J,veres.K,veres.L)
  
  # --- WISTHALER
  #wisthaler.821.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-PTRMS-NH3-1Hz_DC8_20190821_R1.ict')

  # ---- BLAKE
  blake.821.1hz = getICARTTdataSIMPLE('InputFiles/WAS-MERGE/firexaq-mrgWAS-dc8_merge_20190821_R1.ict')
  cc = colnames(blake.821.1hz)
  blake.821.merge = blake.821.1hz[,c(1,2,96:225)]
  blake.821.merge$CO_DACOM_DISKIN_BLAKE = blake.821.1hz$CO_DACOM_DISKIN
  blake.821.merge$CO2_7000_ppm_DISKIN_BLAKE = blake.821.1hz$CO2_7000_ppm_DISKIN
  
  # ------ APEL
  apel.821.1hz = getICARTTdataSIMPLE('InputFiles/TOGA-MERGE/firexaq-mrgTOGA-dc8_merge_20190821_R1.ict')
  cc = colnames(apel.821.1hz)
  apel.821.merge = apel.821.1hz[,c(1,2,226:315)]
  apel.821.merge$CO_DACOM_DISKIN_APEL = apel.821.1hz$CO_DACOM_DISKIN
  apel.821.merge$CO2_7000_ppm_DISKIN_APEL =apel.821.1hz$CO2_7000_ppm_DISKIN
  
  # Becky's better merge
  file = 'InputFiles/Hornbrook/FIREX-AQ weighted TOGA merge 2022-01-24_0821.xlsx'
  newTOGA.821 = readxl::read_xlsx(file); newTOGA.821[newTOGA.821==-999] = NaN; newTOGA.821[newTOGA.821==-888] = NaN
  newTOGA.821$CO_DACOM_DISKIN_BECKY = newTOGA.821$CO_DACOM_DISKIN
  newTOGA.821$CO2_7000_ppm_DISKIN_BECKY = NaN
  newTOGA.821$Time_Start=newTOGA.821$Time_Start...4
  # ----GILMAN
  gilman.821.1hz = getICARTTdataSIMPLE('InputFiles/iWAS-MERGE/firexaq-mrgiWAS-dc8_merge_20190821_R1.ict')
  cc = colnames(gilman.821.1hz)
  gilman.821.merge = gilman.821.1hz[,c(1,2,316:361)]
  gilman.821.merge$CO_DACOM_DISKIN_GILMAN = gilman.821.1hz$CO_DACOM_DISKIN
  gilman.821.merge$CO2_7000_ppm_DISKIN_GILMAN = gilman.821.1hz$CO2_7000_ppm_DISKIN
  #gilman.821.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAA-iWAS-VOCs_DC8_20190821_R0.ict')
  
  # ------ Moore 
  moore.821fast = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190821_R0_MOORE.ict')
  
  moore.821p1 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-AerosolCloudConc_DC8_20190821_R0.ict')
  moore.821p2 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-LAScold_DC8_20190821_R0.ict')
  moore.821p3 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-LAShot_DC8_20190821_R0.ict')
  moore.821p4 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-CPSPD_DC8_20190821_R0.ict')
  moore.821p5 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-CDP_DC8_20190821_R0.ict')
  moore.821p6 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-SMPS_DC8_20190821_R0.ict')
  moore.821 =merge(moore.821p1, moore.821p2, by='Time_mid', all = TRUE, incomparables = NA)
  moore.821 =merge(moore.821, moore.821p3, by='Time_mid', all = TRUE, incomparables = NA)
  moore.821 =merge(moore.821, moore.821p4, by='Time_mid', all = TRUE, incomparables = NA)
  moore.821 =merge(moore.821, moore.821p5, by='Time_mid', all = TRUE, incomparables = NA)
  moore.821 =merge(moore.821, moore.821p6, by='Time_mid', all = TRUE, incomparables = NA)
  # ------- append PI to colnames 1hz ----------
  cc = colnames(co2.821.1hz)
  cc[2] = paste(cc[2],'_DISKIN',sep='')
  colnames(co2.821.1hz) = cc 
  cc = colnames(co.ch4.821.1hz)
  cc[2:4] =  paste(cc[2:4],'_DISKIN',sep='')
  colnames(co.ch4.821.1hz) = cc
  cc=colnames(met.821.1hz)
  
  colnames(isaf.821.1hz) = c("Time_Start"  , "CH2O_ISAF_HANISCO" ,"CH2O_ISAF_precision_HANISCO")
  cc[2:36] = paste(cc[2:36],'_YANG',sep='')
  colnames(met.821.1hz) = cc
  cc = colnames(warneke.821.1hz)
  cc[2:43] = paste(cc[2:43],'_WARNEKE', sep='')
  colnames(warneke.821.1hz) = cc
  cc = colnames(rollinsno.821.1hz)
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsno.821.1hz) = cc
  cc = colnames(rollinsso2.821.1hz)
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsso2.821.1hz) = cc
  # make GTCIMS consistent with 1s merge
  cc=colnames(gtcims.821.1hz)
  cc[4:7] = c("PAN_GTCIMS_HUEY" , "PPN_GTCIMS_HUEY"  ,"APAN_GTCIMS_HUEY" ,"PBN_GTCIMS_HUEY")
  colnames(gtcims.821.1hz)=cc
  # since ISAF, CIT, and GTCIMS came from merge tool, alread has PI's appended.
  ryerson.821.1hz <- ryerson.821.1hz[, !duplicated(colnames(ryerson.821.1hz))]
  
  cc = colnames(ryerson.821.1hz)
  cc[2:9] = paste(cc[2:9],'_RYERSON',sep='')
  colnames(ryerson.821.1hz)=cc
  
  cc = colnames(schwarz.821.1hz)
  cc[2:3] = paste(cc[2:3],'_SCHWARZ', sep='')
  colnames(schwarz.821.1hz) =cc
  
  cc = colnames(freid.c2h6.821.1hz)
  cc[4:7] = paste(cc[4:7], '_FRIED',sep='')
  colnames(freid.c2h6.821.1hz) = cc
  
  cc = colnames(freid.ch2o.821.1hz)
  cc[4:7] = paste(cc[4:7], '_FRIED',sep='')
  colnames(freid.ch2o.821.1hz) = cc
  
  colnames(freid.ch2o.821.1hz) = cc
  womack.821.1hz <- womack.821.1hz[, !duplicated(colnames(womack.821.1hz))]
  cc = colnames(womack.821.1hz)
  cc[2:5] = paste(cc[2:5], '_WOMACK',sep='')
  colnames(womack.821.1hz) = cc
  
  veres.821.1hz <- veres.821.1hz[, !duplicated(colnames(veres.821.1hz))]
  cc = colnames(veres.821.1hz)
  cc[2:13] = paste(cc[2:13], '_VERES',sep='')
  colnames(veres.821.1hz) = cc
  
  #cc = colnames(wisthaler.821.1hz)
  #cc[4:5] = paste(cc[4:5], '_WISTHALER',sep='')
  #colnames(wisthaler.821.1hz) = cc
  
  cc = colnames(jimenez.821.1hz)
  cc[2:8] = paste(cc[2:8], '_JIMENEZ',sep='')
  colnames(jimenez.821.1hz) = cc
  
  # --------#######--------- Get 5 or 10 Hz Data 8/21------#######---------
       # MET DATA - BUI + YANG + DLH
  met.821.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190821_R0_met.ict')
  
       # CO2
  co2.821.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CO2-7000-5Hz_DC8_20190821_R1.ict')
       #CO, CH4
  co.ch4.821.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-DACOM-5Hz_DC8_20190821_R1.ict')
       # WARNEKE VOCs
  warneke.821.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-NOAAPTR-VOCs-5Hz_DC8_20190821_R3.ict')
      # ISAF HCHO - merged to 5Hz from the online merge
  isaf.821.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190821_R0_ISAF.ict')
     # ROLLINS SO2 and NO
  rollinsno.821.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-NO-5Hz_DC8_20190821_R0.ict')
  rollinsso2.821.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-SO2-5Hz_DC8_20190821_R1.ict')
    # CIT VOCs - 
  cit.821.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190821_R0_CIT.ict')
    # GTCIMS PANs - not sure how to match up peaks here
  gtcims.821.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190821_R0_huey.ict')
  
  # ----- Jimenez ---
  jimenez.821.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_AMS_20190821_R0_20230314T134024.ict')
  jimenez.821.5hz$OC_PM1_AMS_JIMENEZ = jimenez.821.5hz$OA_PM1_AMS_JIMENEZ/jimenez.821.5hz$OAtoOC_PM1_AMS
  
  # ------- append PI to colnames ----------
  cc = colnames(co2.821.5hz)
  cc[2] = paste(cc[2],'_DISKIN',sep='')
  colnames(co2.821.5hz) = cc 
  cc = colnames(co.ch4.821.5hz)
  cc[2:4] =  paste(cc[2:4],'_DISKIN',sep='')
  colnames(co.ch4.821.5hz) = cc
  cc = colnames(warneke.821.5hz)
  cc[2:43] = paste(cc[2:43],'_WARNEKE', sep='')
  colnames(warneke.821.5hz) = cc
  cc = colnames(rollinsno.821.5hz)
  cc[1] = "Time_Start" # really Time_Mid but need for merge
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsno.821.5hz) = cc
  cc = colnames(rollinsso2.821.5hz)
  cc[1] = "Time_Start" # really Time_Mid but need for merge
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsso2.821.5hz) = cc
  # make GTCIMS consistent with 1s merge
  cc=colnames(gtcims.821.5hz)
  cc[38:39] = c("APAN_GTCIMS_HUEY" , "PAN_GTCIMS_HUEY")
  colnames(gtcims.821.5hz)=cc
  # since ISAF, CIT, and GTCIMS came from merge tool, alread has PI's appended.
  
  # ------- get fuel moisture data --------
  if (doFM == 1){
    f1 = '/Users/ktravis1/Library/CloudStorage/Box-Box/FuelMoisture/fuel_moisture_content-20210715T1049Z/fmc_20190821_20Z.nc'
    fid = nc_open(f1)
    fuelMDead = ncvar_get(fid, varid = 'FMCG2D')
    fuelMLive = ncvar_get(fid, varid = 'FMCGLH2D')
    xlon = ncvar_get(fid, varid="XLONG_M")
    xlat = ncvar_get(fid, varid="XLAT_M")
    nc_close(fid)
  }  
  # ========================== Copper Breaks ======================
  fire="Copper Breaks"; fuel='forest'
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    # get plume start and stop times
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO   
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
       # Time align 1hz data

    # cut data to plume tags
    # shrink plume tags to relevant peaks
    if (i == 1){start=startO+143; stop=stopO-5; startB = startO + 0; stopB = startO + 2} # Smoke in picture, but poor correlation
    tmp = time_align(start,stop,co2.821.5hz,  co.ch4.821.5hz, 
                     warneke.821.5hz,  isaf.821.5hz,
                     rollinsno.821.5hz, rollinsso2.821.5hz, 
                     cit.821.5hz, gtcims.821.5hz,moore.821fast,
                     jimenez.821.5hz,met.821.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    aug21st.fire = time_alignSLOWNOCANS(start,stop, co2.821.1hz, co.ch4.821.1hz, 
                            warneke.821.1hz, isaf.821.1hz,rollinsno.821.1hz, rollinsso2.821.1hz, 
                            cit.821.1hz,gtcims.821.1hz,ryerson.821.1hz , jimenez.821.1hz, 
                            schwarz.821.1hz, freid.c2h6.821.1hz, freid.ch2o.821.1hz, 
                            womack.821.1hz, stclair.821.1hz,veres.821.1hz, #wisthaler.821.1hz,
                            moore.821, met.821.1hz)
    
    ind2 = which(aug21st.fire$Time_Start >= start & aug21st.fire$Time_Start <= stop	) 
    ind2B = which(aug21st.fire$Time_Start >= startB & aug21st.fire$Time_Start <= stopB	) 
    tmp$fire = ''; tmp$fuel = ''
    aug21st.fire$fire = ''; aug21st.fire$fuel = ''
    
    # -------- Blake cans 
    indBLAKE = which(blake.821.merge$Time_Start >= (start-5) & blake.821.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.821.merge[indBLAKE,]
      blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN; 
      blakeBG = blake.821.merge[1,]; blakeBG$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.821.merge$Time_Start >= (start-0) & gilman.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.821.merge[indGILMAN,]
      GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.821.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.821.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.821.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.821.merge$Time_Start >= (start-0) & apel.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.821.merge[indAPEL,]
      APELBG = apel.821.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.821.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.821.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.821.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.821.merge[1,]
    }
    
    # Becky
    indBECKY = which(newTOGA.821$Time_Start >= (start-0) & newTOGA.821$Time_Start <= stop & newTOGA.821$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.821[indBECKY,]
      BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.821$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.821[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.821[1,]
    }
    
    # plot plume?
    if (doplot == 1){
      plotpass5hzJUSTCO(aug21st.fire[ind2,])
      plotpass5hz(aug21st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      plotpass5hzJUSTOA(tmp[ind,])
      
      plotpass1hz(aug21st.fire[ind2,])
      
      plotpass1hzBLAKE(aug21st.fire[ind2,])
      
      # plot background?
      plotpass5hzJUSTCO(aug21st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    
    tmp$fire[ind] = fire
    aug21st.fire$fire[ind2] = fire
    tmp$fuel[ind] = fuel
    aug21st.fire$fuel[ind2] = fuel
    BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i
    if (i == 1){
      allfires.5hz =  tmp[ind,]
      allfires.1hz =  aug21st.fire[ind2,]
      toga.all = BECKY
      was.all = blake
      iwas.all = GILMAN
    } else{
      allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
      allfires.1hz = rbind.fill(allfires.1hz, aug21st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    }
    
    # ------------------ 1Hz ERs/EFs
    tmpEF1 = ERsEFsALLhzv2(aug21st.fire[ind2,],aug21st.fire[ind2B,],xspecies,fire,fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    # ------------------ 5Hz ERs/EFs
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire,fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # append start and stop times to EFs
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop
    # Get fuel moisture
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    # append to fire dataframe
    if (i == 1){CopperBreaks.1hz.EF = tmpEF1 ;CopperBreaks.5hz.EF = tmpEF5 }
    if (i > 1){CopperBreaks.1hz.EF = rbind(CopperBreaks.1hz.EF, tmpEF1) ; CopperBreaks.5hz.EF = rbind(CopperBreaks.5hz.EF, tmpEF5) }
    
  }
  # order by species and append fire ID
  CopperBreaks.5hz.EF = CopperBreaks.5hz.EF[order(CopperBreaks.5hz.EF$variable),]
  CopperBreaks.5hz.EF$transect_source_fire_ID = indA[1]
  CopperBreaks.1hz.EF = CopperBreaks.1hz.EF[order(CopperBreaks.1hz.EF$variable),]
  CopperBreaks.1hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== VIVIAN ======================
  fire = "Vivian"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    # get plume start and stop times
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO   
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    # Time align 1hz data
    
    # shrink plume tags to relevant peaks
    if (i == 1){start=startO+10 ; stop=stopO-7; startB = startO + 0; stopB = startO + 2} # Also poor correlation, CO < 200 ppb
    # cut data to plume tags
    tmp = time_align(start,stop,co2.821.5hz,  co.ch4.821.5hz, 
                     warneke.821.5hz,  isaf.821.5hz,
                     rollinsno.821.5hz, rollinsso2.821.5hz, 
                     cit.821.5hz, gtcims.821.5hz,moore.821fast, jimenez.821.5hz,met.821.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    aug21st.fire = time_alignSLOWNOCANS(start,stop, co2.821.1hz, co.ch4.821.1hz, 
                                        warneke.821.1hz, isaf.821.1hz,rollinsno.821.1hz, rollinsso2.821.1hz, 
                                        cit.821.1hz,gtcims.821.1hz,ryerson.821.1hz , jimenez.821.1hz, 
                                        schwarz.821.1hz, freid.c2h6.821.1hz, freid.ch2o.821.1hz, 
                                        womack.821.1hz, stclair.821.1hz,veres.821.1hz, #wisthaler.821.1hz,
                                        moore.821, met.821.1hz)
    ind2 = which(aug21st.fire$Time_Start >= start & aug21st.fire$Time_Start <= stop	) 
    ind2B = which(aug21st.fire$Time_Start >= startB & aug21st.fire$Time_Start <= stopB	) 
    
    tmp$fire = ''; tmp$fuel = ''
    aug21st.fire$fire = ''; aug21st.fire$fuel = ''
    tmp$fuel[ind] = 'savannah';tmp$fire[ind] = fire
    aug21st.fire$fuel[ind2]='savannah';aug21st.fire$fire[ind2] = fire
    
    # -------- Blake cans 
    indBLAKE = which(blake.821.merge$Time_Start >= (start-5) & blake.821.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.821.merge[indBLAKE,]
      blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.821.merge$Time_Start >= (start-0) & gilman.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.821.merge[indGILMAN,]
      GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.821.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.821.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.821.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.821.merge$Time_Start >= (start-0) & apel.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.821.merge[indAPEL,]
      APELBG = apel.821.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.821.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.821.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.821.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.821.merge[1,]
    }
    
    # Becky
    indBECKY = which(newTOGA.821$Time_Start >= (start-0) & newTOGA.821$Time_Start <= stop & newTOGA.821$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.821[indBECKY,]
      BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.821$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.821[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.821[1,]
    }
    
    # plot plume?
    if (doplot == 1){
      plotpass5hzJUSTCO(aug21st.fire[ind2,])
      plotpass5hz(aug21st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      plotpass1hzBLAKE(aug21st.fire[ind2,])
      
      # plot background?
      plotpass5hzJUSTCO(aug21st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    
   
   # Save out all data for analysis
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     ; BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug21st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}; 
    # ------------------ 1Hz ERs/EFs
    tmpEF1 = ERsEFsALLhzv2(aug21st.fire[ind2,],aug21st.fire[ind2B,],xspecies,fire, fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz ERs/EFs
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire, fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # append start and stop times to EFs
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop
    # Get fuel moisture
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    # append to fire dataframe
    if (i == 1){Vivian.1hz.EF = tmpEF1 ;Vivian.5hz.EF = tmpEF5 }
    if (i > 1){Vivian.1hz.EF = rbind(Vivian.1hz.EF, tmpEF1) ; Vivian.5hz.EF = rbind(Vivian.5hz.EF, tmpEF5) }
    
  }
  # order by species and append fire ID
  Vivian.5hz.EF = Vivian.5hz.EF[order(Vivian.5hz.EF$variable),]
  Vivian.5hz.EF$transect_source_fire_ID = indA[1]
  Vivian.1hz.EF = Vivian.1hz.EF[order(Vivian.1hz.EF$variable),]
  Vivian.1hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== HALF PINT ======================
  fire="Half pint"; fuel = "slash"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    # get plume start and stop times
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    
    # shrink plume tags to relevant peaks
    if (i == 1){start=startO+3.4;  stop=stopO-2.8; startB = startO + 0; stopB = startO + 2}
    if (i == 2){start=startO+9;  stop=stopO-.4;  startB = startO + 0; stopB = startO + 2}
    if (i == 3){start=startO+0;    stop=stopO-16;  startB = 66003; stopB = 66005}# split from above
    if (i == 4){start=startO+55;   stop= stopO-70;  startB = startO + 0; stopB = startO + 2}
    tmp = time_align(start,stop,co2.821.5hz,  co.ch4.821.5hz, 
                     warneke.821.5hz,  isaf.821.5hz,
                     rollinsno.821.5hz, rollinsso2.821.5hz, 
                     cit.821.5hz, gtcims.821.5hz,moore.821fast, jimenez.821.5hz,met.821.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i==2){start=start-1}
    aug21st.fire = time_alignSLOWNOCANS(start,stop, co2.821.1hz, co.ch4.821.1hz, 
                                        warneke.821.1hz, isaf.821.1hz,rollinsno.821.1hz, rollinsso2.821.1hz, 
                                        cit.821.1hz,gtcims.821.1hz,ryerson.821.1hz , jimenez.821.1hz, 
                                        schwarz.821.1hz, freid.c2h6.821.1hz, freid.ch2o.821.1hz, 
                                        womack.821.1hz, stclair.821.1hz,veres.821.1hz, #wisthaler.821.1hz,
                                        moore.821, met.821.1hz)
    
    ind2 = which(aug21st.fire$Time_Start >= start & aug21st.fire$Time_Start  <= stop	) 
    ind2B = which(aug21st.fire$Time_Start >= startB & aug21st.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug21st.fire$fire = ''; aug21st.fire$fuel = ''
    tmp$fuel[ind] = 'slash' ; tmp$fire = fire
    aug21st.fire$fuel[ind2] = 'slash' ; aug21st.fire$fire[ind2] = fire
    
    # -------- Blake cans 
    indBLAKE = which(blake.821.merge$Time_Start >= (start-5) & blake.821.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.821.merge[indBLAKE,]
      blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.821.merge$Time_Start >= (start-0) & gilman.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.821.merge[indGILMAN,]
      GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.821.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.821.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.821.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.821.merge$Time_Start >= (start-0) & apel.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.821.merge[indAPEL,]
      APELBG = apel.821.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.821.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.821.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.821.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.821.merge[1,]
    }
    
    # Becky
    indBECKY = which(newTOGA.821$Time_Start >= (start-0) & newTOGA.821$Time_Start <= stop & newTOGA.821$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.821[indBECKY,]
      BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.821$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.821[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.821[1,]
    }
    
    if (doplot == 1){
      plotpass5hzJUSTCO(aug21st.fire[ind2,])
      plotpass5hz(aug21st.fire[ind2,])
      plotpass1hz(aug21st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug21st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     ; BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug21st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz 
    tmpEF1 = ERsEFsALLhzv2(aug21st.fire[ind2,],aug21st.fire[ind2B,],xspecies,fire,fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire,fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Halfpint.1hz.EF = tmpEF1 ;Halfpint.5hz.EF = tmpEF5 }
    if (i > 1){Halfpint.1hz.EF = rbind(Halfpint.1hz.EF, tmpEF1) ; Halfpint.5hz.EF = rbind(Halfpint.5hz.EF, tmpEF5) }
  }
  Halfpint.1hz.EF = Halfpint.1hz.EF[order(Halfpint.1hz.EF$variable),]
  Halfpint.1hz.EF$transect_source_fire_ID = indA[1]
  Halfpint.5hz.EF = Halfpint.5hz.EF[order(Halfpint.5hz.EF$variable),]
  Halfpint.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== LORETTA ======================
  fire="Loretta"; fuel = "pile"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    # get plume start and stop times
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
  
    # shrink plume tags to relevant peaks
    if (i == 1){start=startO+5;    stop=stopO-17;   startB = startO + 0; stopB = startO + 2}
    if (i == 2){start=startO+6;    stop=stopO-6;    startB = startO + 0; stopB = startO + 2}
    if (i == 3){start=startO+61;   stop=stopO-16;  startB = 66832; stopB = 66834} # split from above
    if (i == 4){start=startO+49; stop=stopO-170;  startB = startO + 0; stopB = startO + 2} # bad pass, get rid of it?
    tmp = time_align(start,stop,co2.821.5hz,  co.ch4.821.5hz, 
                     warneke.821.5hz,  isaf.821.5hz,
                     rollinsno.821.5hz, rollinsso2.821.5hz, 
                     cit.821.5hz, gtcims.821.5hz,moore.821fast, jimenez.821.5hz,met.821.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    aug21st.fire = time_alignSLOWNOCANS(start,stop, co2.821.1hz, co.ch4.821.1hz, 
                                        warneke.821.1hz, isaf.821.1hz,rollinsno.821.1hz, rollinsso2.821.1hz, 
                                        cit.821.1hz,gtcims.821.1hz,ryerson.821.1hz , jimenez.821.1hz, 
                                        schwarz.821.1hz, freid.c2h6.821.1hz, freid.ch2o.821.1hz, 
                                        womack.821.1hz, stclair.821.1hz,veres.821.1hz, #wisthaler.821.1hz,
                                        moore.821, met.821.1hz)
    ind2 = which(aug21st.fire$Time_Start >= start & aug21st.fire$Time_Start  <= stop	) 
    ind2B = which(aug21st.fire$Time_Start >= startB & aug21st.fire$Time_Start <= stopB	) 

    tmp$fire = ''; tmp$fuel = ''
    aug21st.fire$fire = ''; aug21st.fire$fuel = ''
    tmp$fuel[ind] = fuel;tmp$fire[ind] =  'Loretta'
    aug21st.fire$fuel[ind2] = fuel;aug21st.fire$fire[ind2] =  'Loretta'
    
    # -------- Blake cans 
    indBLAKE = which(blake.821.merge$Time_Start >= (start-5) & blake.821.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.821.merge[indBLAKE,]
      blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.821.merge$Time_Start >= (start-0) & gilman.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.821.merge[indGILMAN,]
      GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.821.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.821.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.821.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.821.merge$Time_Start >= (start-0) & apel.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.821.merge[indAPEL,]
      APELBG = apel.821.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.821.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.821.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.821.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.821.merge[1,]
    }
    
    # Becky
    indBECKY = which(newTOGA.821$Time_Start >= (start-0) & newTOGA.821$Time_Start <= stop & newTOGA.821$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.821[indBECKY,]
      BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.821$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.821[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.821[1,]
    }
    
    if (doplot == 1){
      plotpass5hzJUSTCO(aug21st.fire[ind2,])
      plotpass5hz(aug21st.fire[ind2,])
      plotpass1hz(aug21st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug21st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     ; BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug21st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz is slash correct?
    tmpEF1 = ERsEFsALLhzv2(aug21st.fire[ind2,],aug21st.fire[ind2B,],xspecies,fire, fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire, fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Loretta.1hz.EF = tmpEF1 ;Loretta.5hz.EF = tmpEF5 }
    if (i > 1){Loretta.1hz.EF = rbind(Loretta.1hz.EF, tmpEF1) ; Loretta.5hz.EF = rbind(Loretta.5hz.EF, tmpEF5) }
  }
  Loretta.1hz.EF = Loretta.1hz.EF[order(Loretta.1hz.EF$variable),]
  Loretta.1hz.EF$transect_source_fire_ID = indA[1]
  Loretta.5hz.EF = Loretta.5hz.EF[order(Loretta.5hz.EF$variable),]
  Loretta.5hz.EF$transect_source_fire_ID = indA[1]
  # ========================== LIL DEBBIE =========================================
  fire="Lil Debby"; fuel = "corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  #plotfire(aug21st.fire, flags, indA,"LilDebbie", c(66E3, 729E2))  # i = 2 start =68994 stop= 69006
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    
    if (i == 1){start = startO +  6; stop = stopO -12.8; startB = startO -1; stopB = startO+1}
    if (i == 2){start = startO +  1; stop = stopO -0.2; startB = startO -2; stopB = startO}
    if (i == 3){start = startO +  0; stop = stopO -66; startB = 68992; stopB = 68994} # split from above
    if (i == 4){start = startO + 4; stop = stopO -8; startB = startO + 0; stopB = startO+2 }
    if (i == 5){start = startO + 24; stop = stopO - 0; startB = startO -1; stopB = startO+1} 
    if (i == 6){start = startO + 0; stop = stopO - 0; startB = 71330; stopB = 71332} #split from above
    if (i == 7){start = startO + 0; stop = stopO - 8.2; startB = 71330; stopB = 71332} # split from above
    
    tmp = time_align(start,stop,co2.821.5hz,  co.ch4.821.5hz, 
                       warneke.821.5hz,  isaf.821.5hz,
                       rollinsno.821.5hz, rollinsso2.821.5hz, 
                       cit.821.5hz, gtcims.821.5hz,moore.821fast, jimenez.821.5hz,met.821.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    # just for i=1 or7, need one more point for 1hz data
    if (i == 1){start = start - 1; startB=startB-1;stopB=stopB-1}
    if (i == 2){start = start - 1}
    if (i == 3){stop = stop + 1}
    if (i == 4){start = start - 1}
    if (i == 7){start = start - 1; stop=stop+1}
    
    aug21st.fire = time_alignSLOWNOCANS(start,stop, co2.821.1hz, co.ch4.821.1hz, 
                              warneke.821.1hz, isaf.821.1hz,rollinsno.821.1hz, rollinsso2.821.1hz, 
                              cit.821.1hz,gtcims.821.1hz,ryerson.821.1hz , jimenez.821.1hz, 
                              schwarz.821.1hz, freid.c2h6.821.1hz, freid.ch2o.821.1hz, 
                              womack.821.1hz, stclair.821.1hz,veres.821.1hz, #wisthaler.821.1hz,
                              moore.821, met.821.1hz)
    
    ind2 = which(aug21st.fire$Time_Start >= start & aug21st.fire$Time_Start  <= stop	) 
    ind2B = which(aug21st.fire$Time_Start >= startB & aug21st.fire$Time_Start <= stopB	) 

    if (doplot == 1){
      plotpass5hzJUSTCO(aug21st.fire[ind2,])
      plotpass5hz(aug21st.fire[ind2,])
      plotpass1hz(aug21st.fire[ind2,])
      plotpass1hzHNO2(aug21st.fire[ind2,])
      plotpass5hzJUSTCOHNO2(aug21st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass5hz(tmp[ind,])
      plotpass1hzBLAKE(aug21st.fire[ind2,])
      
      # plot background?
      plotpass5hzJUSTCO(aug21st.fire[ind2B,])
      plotpass5hzJUSTCOHNO2(aug21st.fire[ind2B,])
    }
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug21st.fire$fire = ''; aug21st.fire$fuel = ''
    tmp$fuel[ind] = fuel;tmp$fire[ind] =  'Lil Debbie'
    aug21st.fire$fuel[ind2] = fuel;aug21st.fire$fire[ind2] =  'Lil Debbie'
   
    # -------- Blake cans 
    indBLAKE = which(blake.821.merge$Time_Start >= (start-5) & blake.821.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.821.merge[indBLAKE,]
      blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.821.merge$Time_Start >= (start-0) & gilman.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.821.merge[indGILMAN,]
      GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.821.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.821.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.821.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.821.merge$Time_Start >= (start-0) & apel.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.821.merge[indAPEL,]
      APELBG = apel.821.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.821.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.821.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.821.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.821.merge[1,]
    }
    
    # Becky's new merge 
    indBECKY = which(newTOGA.821$Time_Start >= (start-0) & newTOGA.821$Time_Start <= stop & newTOGA.821$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.821[indBECKY,]
      BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.821$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.821[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.821[1,]
    }
    
    
    if (doplot == 1){
      plotpass5hzJUSTCO(aug21st.fire[ind2,])
      plotpass5hz(aug21st.fire[ind2,])
      plotpass1hz(aug21st.fire[ind2,])
      plotpass1hzHNO2(aug21st.fire[ind2,])
      plotpass5hzJUSTCOHNO2(aug21st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass5hz(tmp[ind,])
      plotpass1hzBLAKE(aug31st.fire[ind2,])
      
      # plot background?
      plotpass5hzJUSTCO(aug21st.fire[ind2B,])
      plotpass5hzJUSTCOHNO2(aug21st.fire[ind2B,])
      
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     ; BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug21st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug21st.fire[ind2,],aug21st.fire[ind2B,],xspecies,fire,fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire,fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){
      LilDebbie.1hz.EF = tmpEF1 
      LilDebbie.5hz.EF = tmpEF5 
    }
    if (i > 1){
      LilDebbie.1hz.EF = rbind(LilDebbie.1hz.EF, tmpEF1) 
      LilDebbie.5hz.EF = rbind(LilDebbie.5hz.EF, tmpEF5) 
      }
  }
  LilDebbie.5hz.EF = LilDebbie.5hz.EF[order(LilDebbie.5hz.EF$variable),]
  LilDebbie.5hz.EF$transect_source_fire_ID = indA[1]
  LilDebbie.1hz.EF = LilDebbie.1hz.EF[order(LilDebbie.1hz.EF$variable),]
  LilDebbie.1hz.EF$transect_source_fire_ID = indA[1]

  # ========================== RICEARONI =====================
  fire="Rice-a-Roni";fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
 
    if (i == 1){start = startO + 7.8; stop = stopO -24; startB = startO -1; stopB =startO + 1}
    if (i == 2){start = startO + 5; stop = stopO -40.4; startB = startO + 0; stopB = startO + 2} # bad pass, CO < 200 ppb
    if (i == 3){start = startO + 3.; stop = stopO -0; startB = startO -1; stopB = startO + 1} 
    if (i == 4){start = startO + 0; stop = stopO -9; startB = 69635; stopB = 69637} # split from above
    if (i == 5){start = startO + 3; stop = stopO -3; startB = startO + 0; stopB = startO + 2} 
    tmp = time_align(start,stop,co2.821.5hz,  co.ch4.821.5hz, 
                     warneke.821.5hz,  isaf.821.5hz,
                     rollinsno.821.5hz, rollinsso2.821.5hz, 
                     cit.821.5hz, gtcims.821.5hz,moore.821fast, jimenez.821.5hz,met.821.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	)
    
    if (i ==2){start=start-1;stop=stop+1}
    aug21st.fire = time_alignSLOWNOCANS(start,stop, co2.821.1hz, co.ch4.821.1hz, 
                                        warneke.821.1hz, isaf.821.1hz,rollinsno.821.1hz, rollinsso2.821.1hz, 
                                        cit.821.1hz,gtcims.821.1hz,ryerson.821.1hz , jimenez.821.1hz, 
                                        schwarz.821.1hz, freid.c2h6.821.1hz, freid.ch2o.821.1hz, 
                                        womack.821.1hz, stclair.821.1hz,veres.821.1hz, #wisthaler.821.1hz,
                                        moore.821, met.821.1hz)
    ind2 = which(aug21st.fire$Time_Start >= start & aug21st.fire$Time_Start  <= stop	) 
    ind2B = which(aug21st.fire$Time_Start >= startB & aug21st.fire$Time_Start <= stopB	) 

    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug21st.fire$fire = ''; aug21st.fire$fuel = ''
    tmp$fuel[ind] = fuel;tmp$fire[ind] =  'Ricearoni'
    aug21st.fire$fuel[ind2]= fuel;aug21st.fire$fire[ind2] =  'Ricearoni'
    
    # -------- Blake cans 
    indBLAKE = which(blake.821.merge$Time_Start >= (start-5) & blake.821.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.821.merge[indBLAKE,]
      blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.821.merge$Time_Start >= (start-0) & gilman.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.821.merge[indGILMAN,]
      GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.821.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.821.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.821.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.821.merge$Time_Start >= (start-0) & apel.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.821.merge[indAPEL,]
      APELBG = apel.821.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.821.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.821.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.821.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.821.merge[1,]
    }
    
    # Becky
    indBECKY = which(newTOGA.821$Time_Start >= (start-0) & newTOGA.821$Time_Start <= stop & newTOGA.821$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.821[indBECKY,]
      BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (as.numeric(BECKYBG$CO_DACOM_DISKIN_BECKY[1]) > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.821$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.821[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.821[1,]
    }
    
    if (doplot == 1){
      toplot = aug21st.fire[ind2,]
      par(mfrow=c(2,2), mar=c(5,5,5,5))
      plot(toplot$CO_DACOM_DISKIN, toplot$CNgt3nm_stdPT, pch=19, xlab='CO, ppb', ylab='CNgt3nm_stdPT')
      plot(toplot$CO_DACOM_DISKIN, toplot$CNgt6nm_stdPT, pch=19, xlab='CO, ppb', ylab='CNgt6nm_stdPT')
      plot(toplot$CO_DACOM_DISKIN, toplot$CNgt20nm_stdPT, pch=19, xlab='CO, ppb', ylab='CNgt20nm_stdPT')
      plot(toplot$CO_DACOM_DISKIN, toplot$CNgt20nm_nonvol_stdPT, pch=19, xlab='CO, ppb', ylab='CNgt20nm_non vol')
      
      plotpass5hzJUSTCO(aug21st.fire[ind2,])
      plotpass5hz(aug21st.fire[ind2,])
      plotpass1hz(aug21st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug21st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
      
      pass = tmp[ind,]
      ind = which(is.finite(pass$CH2O_ISAF_HANISCO) & is.finite(pass$CO_DACOM_DISKIN))
      pass = pass[ind,]
      pass$xKTE = 1#pass$CO_DACOM_DISKIN*0.02 # 2 pct
      pass$sKTE = 1#(pass$CH2O_ISAF_HANISCO*0.1+10)/1E3*
      pKT = york_regression(pass, x='CO_DACOM_DISKIN',y='CH2O_ISAF_HANISCO_ppb', variance_x = xKTE, variance_y = sKTE )
      #plot(pass$CO_DACOM_DISKIN, pass$CH2O_ISAF_HANISCO, xlab='CO, ppb', ylab='CH2O, ppt')
      #abline(a= pKT$intercept,b=pKT$slope*1E3)
      #points(2267.3,92672.9, col='red', pch=19)
      CH2O = c(92673,	75888	,96353)
      CO = c(1779.82,	1046.93,	2405.39)
      CH2O=c(92672.9,	75887.95,	96353.33)
      CO = c(2267.3	,1327.27,	2676.97525)
      #abline(lm(CH2O~CO), col='blue')
      tmp = lmodel2(pass$CH2O_ISAF_HANISCO ~ pass$CO_DACOM_DISKIN, range.x = 'relative',range.y='relative')
      #abline(tmp$regression.results$Intercept[1], tmp$regression.results$Slope[1], col='blue')
      #abline(tmp$regression.results$Intercept[2], tmp$regression.results$Slope[2], col='green')
      #abline(tmp$regression.results$Intercept[3], tmp$regression.results$Slope[3], col='orange')
      #abline(tmp$regression.results$Intercept[4], tmp$regression.results$Slope[4], col='pink')
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug21st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug21st.fire[ind2,],aug21st.fire[ind2B,],xspecies,fire, fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire, fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Ricearoni.1hz.EF = tmpEF1 ;Ricearoni.5hz.EF = tmpEF5 }
    if (i > 1){Ricearoni.1hz.EF = rbind(Ricearoni.1hz.EF, tmpEF1) ; Ricearoni.5hz.EF = rbind(Ricearoni.5hz.EF, tmpEF5) }
  }
  Ricearoni.1hz.EF = Ricearoni.1hz.EF[order(Ricearoni.1hz.EF$variable),]
  Ricearoni.1hz.EF$transect_source_fire_ID = indA[1]
  Ricearoni.5hz.EF = Ricearoni.5hz.EF[order(Ricearoni.5hz.EF$variable),]
  Ricearoni.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== CRAWBABY-HOUSE  =====================
  fire="Crawbaby house";fuel="house"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    
    if (i == 1){start = startO + 1; stop = stopO -4; startB = startO -1; stopB = startO + 1}
    if (i == 2){start = startO + 4; stop = stopO -6; startB = startO + 0; stopB = startO + 2}
    tmp = time_align(start,stop,co2.821.5hz,  co.ch4.821.5hz, 
                     warneke.821.5hz,  isaf.821.5hz,
                     rollinsno.821.5hz, rollinsso2.821.5hz, 
                     cit.821.5hz, gtcims.821.5hz,moore.821fast, jimenez.821.5hz,met.821.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    # just for 1hz data
    if (i ==1){start = start - 1}
    if (i ==2){start = start - 1}
    aug21st.fire = time_alignSLOWNOCANS(start,stop, co2.821.1hz, co.ch4.821.1hz, 
                                        warneke.821.1hz, isaf.821.1hz,rollinsno.821.1hz, rollinsso2.821.1hz, 
                                        cit.821.1hz,gtcims.821.1hz,ryerson.821.1hz , jimenez.821.1hz, 
                                        schwarz.821.1hz, freid.c2h6.821.1hz, freid.ch2o.821.1hz, 
                                        womack.821.1hz, stclair.821.1hz,veres.821.1hz, #wisthaler.821.1hz,
                                        moore.821, met.821.1hz)
    ind2 = which(aug21st.fire$Time_Start >= start & aug21st.fire$Time_Start  <= stop	) 
    ind2B = which(aug21st.fire$Time_Start >= startB & aug21st.fire$Time_Start <= stopB	)
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug21st.fire$fire = ''; aug21st.fire$fuel = ''
    
    tmp$fire[ind] =  'Crawbaby House'
    aug21st.fire$fire[ind2] =  'Crawbaby House'
    
    # -------- Blake cans 
    indBLAKE = which(blake.821.merge$Time_Start >= (start-5) & blake.821.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.821.merge[indBLAKE,]
      blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.821.merge$Time_Start >= (start-0) & gilman.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.821.merge[indGILMAN,]
      GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.821.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.821.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.821.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.821.merge$Time_Start >= (start-0) & apel.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.821.merge[indAPEL,]
      APELBG = apel.821.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.821.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.821.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.821.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.821.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.821$Time_Start >= (start-0) & newTOGA.821$Time_Start <= stop & newTOGA.821$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.821[indBECKY,]
      BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.821$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.821[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.821[1,]
    }
    
    if (doplot == 1){
      plotpass5hzJUSTCO(aug21st.fire[ind2,])
      plotpass5hz(aug21st.fire[ind2,])
      plotpass1hz(aug21st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug21st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug21st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug21st.fire[ind2,],aug21st.fire[ind2B,],xspecies,fire,fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire,fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){CrawbabyHouse.1hz.EF = tmpEF1 ;CrawbabyHouse.5hz.EF = tmpEF5 }
    if (i > 1){CrawbabyHouse.1hz.EF = rbind(CrawbabyHouse.1hz.EF, tmpEF1) ; CrawbabyHouse.5hz.EF = rbind(CrawbabyHouse.5hz.EF, tmpEF5) }
  }
  CrawbabyHouse.1hz.EF = CrawbabyHouse.1hz.EF[order(CrawbabyHouse.1hz.EF$variable),]
  CrawbabyHouse.1hz.EF$transect_source_fire_ID = indA[1]
  CrawbabyHouse.5hz.EF = CrawbabyHouse.5hz.EF[order(CrawbabyHouse.5hz.EF$variable),]
  CrawbabyHouse.5hz.EF$transect_source_fire_ID = indA[1]
  # ========================== CRAWDADDY  =====================
  fire="Crawdaddy";fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    
    if (i == 1){start = startO + 4.4; stop = stopO -9.8; startB = startO-1; stopB = startO + 1}
    if (i == 2){start = startO + 1; stop = stopO -7; startB = startO -1; stopB = startO + 0}
    if (i == 3){start = startO + 2.4; stop = stopO -0; startB = startO -1; stopB = startO +0.8}
    if (i == 4){start = startO + 0; stop = stopO -13; startB = 71177; stopB = 71178.8} # split from above
    tmp = time_align(start,stop,co2.821.5hz,  co.ch4.821.5hz, 
                     warneke.821.5hz,  isaf.821.5hz,
                     rollinsno.821.5hz, rollinsso2.821.5hz, 
                     cit.821.5hz, gtcims.821.5hz,moore.821fast, jimenez.821.5hz,met.821.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    # just for 1hz
    if (i ==1){start = start - 1}
    if (i ==2){start = start - 1}
    if (i ==3){start = start - 2; startB = startB - 1; stopB=stopB - 1}
    if (i ==4){start = start - 1}
    aug21st.fire = time_alignSLOWNOCANS(start,stop, co2.821.1hz, co.ch4.821.1hz, 
                                        warneke.821.1hz, isaf.821.1hz,rollinsno.821.1hz, rollinsso2.821.1hz, 
                                        cit.821.1hz,gtcims.821.1hz,ryerson.821.1hz , jimenez.821.1hz, 
                                        schwarz.821.1hz, freid.c2h6.821.1hz, freid.ch2o.821.1hz, 
                                        womack.821.1hz, stclair.821.1hz,veres.821.1hz, #wisthaler.821.1hz,
                                        moore.821, met.821.1hz)
    ind2 = which(aug21st.fire$Time_Start >= start & aug21st.fire$Time_Start  <= stop	) 
    ind2B = which(aug21st.fire$Time_Start >= startB & aug21st.fire$Time_Start <= stopB	) 

    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug21st.fire$fire = ''; aug21st.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = fire
    aug21st.fire$fuel[ind2] = fuel ; aug21st.fire$fire[ind2] = fire
    
    # -------- Blake cans 
    indBLAKE = which(blake.821.merge$Time_Start >= (start-5) & blake.821.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.821.merge[indBLAKE,]
      blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.821.merge$Time_Start >= (start-0) & gilman.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.821.merge[indGILMAN,]
      GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.821.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.821.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.821.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.821.merge$Time_Start >= (start-0) & apel.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.821.merge[indAPEL,]
      APELBG = apel.821.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.821.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.821.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.821.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.821.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.821$Time_Start >= (start-0) & newTOGA.821$Time_Start <= stop & newTOGA.821$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.821[indBECKY,]
      BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.821$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.821[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.821[1,]
    }
    
    if (doplot == 1){
      plotpass5hzJUSTCO(aug21st.fire[ind2,])
      plotpass5hz(aug21st.fire[ind2,])
      plotpass1hz(aug21st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug21st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug21st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug21st.fire[ind2,],aug21st.fire[ind2B,],xspecies,fire, fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire, fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop  
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Crawdaddy.1hz.EF = tmpEF1 ;Crawdaddy.5hz.EF = tmpEF5 }
    if (i > 1){Crawdaddy.1hz.EF = rbind(Crawdaddy.1hz.EF, tmpEF1) ; Crawdaddy.5hz.EF = rbind(Crawdaddy.5hz.EF, tmpEF5) }
  }
  
  Crawdaddy.1hz.EF = Crawdaddy.1hz.EF[order(Crawdaddy.1hz.EF$variable),]
  Crawdaddy.1hz.EF$transect_source_fire_ID = indA[1]
  Crawdaddy.5hz.EF = Crawdaddy.5hz.EF[order(Crawdaddy.5hz.EF$variable),]
  Crawdaddy.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== GUMBO =====================
  fire="Gumbo"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 3; stop = stopO -7.6; startB = startO + 0; stopB = startO + 2}
    if (i == 2){start = startO + 3.8; stop = stopO -1; startB = startO + 0; stopB = startO + 2}
    if (i == 3){start = startO + 2; stop = stopO -33; startB = 70900; stopB = 70902} # split from above
    tmp = time_align(start,stop,co2.821.5hz,  co.ch4.821.5hz, 
                     warneke.821.5hz,  isaf.821.5hz,
                     rollinsno.821.5hz, rollinsso2.821.5hz, 
                     cit.821.5hz, gtcims.821.5hz,moore.821fast, jimenez.821.5hz,met.821.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i == 1){ stop = stop -2}    # adjust for 1hz data
    if (i == 2){ start = start - 1}    # adjust for 1hz data
    if (i == 3){ start = start - 1}    # adjust for 1hz data
    aug21st.fire = time_alignSLOWNOCANS(start,stop, co2.821.1hz, co.ch4.821.1hz, 
                                        warneke.821.1hz, isaf.821.1hz,rollinsno.821.1hz, rollinsso2.821.1hz, 
                                        cit.821.1hz,gtcims.821.1hz,ryerson.821.1hz , jimenez.821.1hz, 
                                        schwarz.821.1hz, freid.c2h6.821.1hz, freid.ch2o.821.1hz, 
                                        womack.821.1hz, stclair.821.1hz,veres.821.1hz, #wisthaler.821.1hz,
                                        moore.821, met.821.1hz)
    ind2 = which(aug21st.fire$Time_Start >= start & aug21st.fire$Time_Start  <= stop	) 
    ind2B = which(aug21st.fire$Time_Start >= startB & aug21st.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug21st.fire$fire = ''; aug21st.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Gumbo'
    aug21st.fire$fuel[ind2] = fuel ; aug21st.fire$fire[ind2] = 'Gumbo'
    
    # -------- Blake cans 
    indBLAKE = which(blake.821.merge$Time_Start >= (start-5) & blake.821.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.821.merge[indBLAKE,]
      blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.821.merge$Time_Start >= (start-0) & gilman.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      # kludge
      if (length(indGILMAN) > 1){indGILMAN = indGILMAN[2]}
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.821.merge[indGILMAN,]
      GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.821.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.821.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.821.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.821.merge$Time_Start >= (start-0) & apel.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.821.merge[indAPEL,]
      APELBG = apel.821.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.821.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.821.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.821.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.821.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.821$Time_Start >= (start-0) & newTOGA.821$Time_Start <= stop & newTOGA.821$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.821[indBECKY,]
      BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.821$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.821[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.821[1,]
    }
    
    if (doplot == 1){
      plotpass5hzJUSTCO(aug21st.fire[ind2,])
      plotpass5hz(aug21st.fire[ind2,])
      plotpass1hz(aug21st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug21st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug21st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug21st.fire[ind2,],aug21st.fire[ind2B,],xspecies,'Gumbo', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Gumbo', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop  
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Gumbo.1hz.EF = tmpEF1 ;Gumbo.5hz.EF = tmpEF5 }
    if (i > 1){Gumbo.1hz.EF = rbind(Gumbo.1hz.EF, tmpEF1) ; Gumbo.5hz.EF = rbind(Gumbo.5hz.EF, tmpEF5) }
    
  }
  Gumbo.1hz.EF = Gumbo.1hz.EF[order(Gumbo.1hz.EF$variable),]
  Gumbo.1hz.EF$transect_source_fire_ID = indA[1]
  Gumbo.5hz.EF = Gumbo.5hz.EF[order(Gumbo.5hz.EF$variable),]
  Gumbo.5hz.EF$transect_source_fire_ID = indA[1]
  
  #========================== JUMBALAYA =====================
 
  fire="Jumbalaya"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 1; stop = stopO -7; startB = startO -2; stopB =  startO + 0}
    tmp = time_align(start,stop,co2.821.5hz,  co.ch4.821.5hz, 
                     warneke.821.5hz,  isaf.821.5hz,
                     rollinsno.821.5hz, rollinsso2.821.5hz, 
                     cit.821.5hz, gtcims.821.5hz,moore.821fast, jimenez.821.5hz,met.821.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i ==1){start=start-1; stop=stop+1}
    aug21st.fire = time_alignSLOWNOCANS(start,stop, co2.821.1hz, co.ch4.821.1hz, 
                                        warneke.821.1hz, isaf.821.1hz,rollinsno.821.1hz, rollinsso2.821.1hz, 
                                        cit.821.1hz,gtcims.821.1hz,ryerson.821.1hz , jimenez.821.1hz, 
                                        schwarz.821.1hz, freid.c2h6.821.1hz, freid.ch2o.821.1hz, 
                                        womack.821.1hz, stclair.821.1hz,veres.821.1hz, #wisthaler.821.1hz,
                                        moore.821, met.821.1hz)
    ind2 = which(aug21st.fire$Time_Start >= start & aug21st.fire$Time_Start  <= stop	) 
    ind2B = which(aug21st.fire$Time_Start >= startB & aug21st.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug21st.fire$fire = ''; aug21st.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Jumbalaya'
    aug21st.fire$fuel[ind2] = fuel ; aug21st.fire$fire[ind2] = 'Jumbalaya'
    
    # -------- Blake cans 
    indBLAKE = which(blake.821.merge$Time_Start >= (start-5) & blake.821.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.821.merge[indBLAKE,]
      blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.821.merge$Time_Start >= (start-0) & gilman.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.821.merge[indGILMAN,]
      GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.821.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.821.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.821.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.821.merge$Time_Start >= (start-0) & apel.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.821.merge[indAPEL,]
      APELBG = apel.821.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.821.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.821.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.821.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.821.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.821$Time_Start >= (start-0) & newTOGA.821$Time_Start <= stop & newTOGA.821$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.821[indBECKY,]
      BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.821$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.821[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.821[1,]
    }
    
    if (doplot == 1){
      plotpass5hzJUSTCO(aug21st.fire[ind2,])
      plotpass5hz(aug21st.fire[ind2,])
      plotpass1hz(aug21st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug21st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug21st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug21st.fire[ind2,],aug21st.fire[ind2B,],xspecies,'Jumbalaya', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Jumbalaya', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop     
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Jumbalaya.1hz.EF = tmpEF1 ;Jumbalaya.5hz.EF = tmpEF5 }
    if (i > 1){Jumbalaya.1hz.EF = rbind(Jumbalaya.1hz.EF, tmpEF1) ; Jumbalaya.5hz.EF = rbind(Jumbalaya.5hz.EF, tmpEF5) }
  }
  Jumbalaya.1hz.EF = Jumbalaya.1hz.EF[order(Jumbalaya.1hz.EF$variable),]
  Jumbalaya.1hz.EF$transect_source_fire_ID = indA[1]
  Jumbalaya.5hz.EF = Jumbalaya.5hz.EF[order(Jumbalaya.5hz.EF$variable),]
  Jumbalaya.5hz.EF$transect_source_fire_ID = indA[1]
  
  #========================== JUMBALAYA JR =====================
  fire="Jumbalaya Jr"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    
    if (i == 1){start = startO + 8; stop = stopO -27.4; startB = startO -2; stopB = startO +0} # bad pass
    tmp = time_align(start,stop,co2.821.5hz,  co.ch4.821.5hz, 
                     warneke.821.5hz,  isaf.821.5hz,
                     rollinsno.821.5hz, rollinsso2.821.5hz, 
                     cit.821.5hz, gtcims.821.5hz,moore.821fast, jimenez.821.5hz,met.821.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i ==1){start=start-1; stop=stop+1}
    aug21st.fire = time_alignSLOWNOCANS(start,stop, co2.821.1hz, co.ch4.821.1hz, 
                                        warneke.821.1hz, isaf.821.1hz,rollinsno.821.1hz, rollinsso2.821.1hz, 
                                        cit.821.1hz,gtcims.821.1hz,ryerson.821.1hz , jimenez.821.1hz, 
                                        schwarz.821.1hz, freid.c2h6.821.1hz, freid.ch2o.821.1hz, 
                                        womack.821.1hz, stclair.821.1hz,veres.821.1hz, #wisthaler.821.1hz,
                                        moore.821, met.821.1hz)
    
    ind2 = which(aug21st.fire$Time_Start >= start & aug21st.fire$Time_Start  <= stop	) 
    ind2B = which(aug21st.fire$Time_Start >= startB & aug21st.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug21st.fire$fire = ''; aug21st.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Jumbalaya Jr'
    aug21st.fire$fuel[ind2] = fuel ; aug21st.fire$fire[ind2] = 'Jumbalaya Jr'
    
    # -------- Blake cans 
    indBLAKE = which(blake.821.merge$Time_Start >= (start-5) & blake.821.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.821.merge[indBLAKE,]
      blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.821.merge$Time_Start >= (start-0) & gilman.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.821.merge[indGILMAN,]
      GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.821.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.821.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.821.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.821.merge$Time_Start >= (start-0) & apel.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.821.merge[indAPEL,]
      APELBG = apel.821.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.821.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.821.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.821.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.821.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.821$Time_Start >= (start-0) & newTOGA.821$Time_Start <= stop & newTOGA.821$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.821[indBECKY,]
      BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.821$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.821[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.821[1,]
    }
    
    if (doplot == 1){
      plotpass5hzJUSTCO(aug21st.fire[ind2,])
      plotpass5hz(aug21st.fire[ind2,])
      plotpass1hz(aug21st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug21st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug21st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug21st.fire[ind2,],aug21st.fire[ind2B,],xspecies,'Jumbalaya Jr', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Jumbalaya Jr', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop     
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){JumbalayaJr.1hz.EF = tmpEF1 ;JumbalayaJr.5hz.EF = tmpEF5 }
    if (i > 1){JumbalayaJr.1hz.EF = rbind(JumbalayaJr.1hz.EF, tmpEF1) ; JumbalayaJr.5hz.EF = rbind(JumbalayaJr.5hz.EF, tmpEF5) }
  }
  JumbalayaJr.1hz.EF = JumbalayaJr.1hz.EF[order(JumbalayaJr.1hz.EF$variable),]
  JumbalayaJr.1hz.EF$transect_source_fire_ID = indA[1]
  JumbalayaJr.5hz.EF = JumbalayaJr.5hz.EF[order(JumbalayaJr.5hz.EF$variable),]
  JumbalayaJr.5hz.EF$transect_source_fire_ID = indA[1]
  # ========================== CROUTON =====================

  fire="Crouton"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 121; stop = stopO -15; startB = startO +0; stopB = startO + 2}
    if (i == 2){start = startO + 18; stop = stopO -65; startB = startO + 0; stopB = startO + 2}
    tmp = time_align(start,stop,co2.821.5hz,  co.ch4.821.5hz, 
                     warneke.821.5hz,  isaf.821.5hz,
                     rollinsno.821.5hz, rollinsso2.821.5hz, 
                     cit.821.5hz, gtcims.821.5hz,moore.821fast, jimenez.821.5hz,met.821.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1}
    if (i ==2){start=start-1}
    aug21st.fire = time_alignSLOWNOCANS(start,stop, co2.821.1hz, co.ch4.821.1hz, 
                                        warneke.821.1hz, isaf.821.1hz,rollinsno.821.1hz, rollinsso2.821.1hz, 
                                        cit.821.1hz,gtcims.821.1hz,ryerson.821.1hz , jimenez.821.1hz, 
                                        schwarz.821.1hz, freid.c2h6.821.1hz, freid.ch2o.821.1hz, 
                                        womack.821.1hz, stclair.821.1hz,veres.821.1hz, #wisthaler.821.1hz,
                                        moore.821, met.821.1hz)
    ind2 = which(aug21st.fire$Time_Start >= start & aug21st.fire$Time_Start  <= stop	) 
    ind2B = which(aug21st.fire$Time_Start >= startB & aug21st.fire$Time_Start <= stopB	) 
    
        # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug21st.fire$fire = ''; aug21st.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = fire
    aug21st.fire$fuel[ind2] = fuel ; aug21st.fire$fire[ind2] = fire
    
    # -------- Blake cans 
    tmp$fire = ''; tmp$fuel = ''
    aug21st.fire$fire = ''; aug21st.fire$fuel = ''
    indBLAKE = which(blake.821.merge$Time_Start >= (start-5) & blake.821.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.821.merge[indBLAKE,]
      blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.821.merge$Time_Start >= (start-0) & gilman.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.821.merge[indGILMAN,]
      GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.821.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.821.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.821.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.821.merge$Time_Start >= (start-0) & apel.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.821.merge[indAPEL,]
      APELBG = apel.821.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.821.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.821.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.821.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.821.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.821$Time_Start >= (start-0) & newTOGA.821$Time_Start <= stop & newTOGA.821$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.821[indBECKY,]
      BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.821$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.821[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.821[1,]
    }
    
    if (doplot == 1){
      plotpass5hzJUSTCO(aug21st.fire[ind2,])
      plotpass5hz(aug21st.fire[ind2,])
      plotpass1hz(aug21st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug21st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug21st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug21st.fire[ind2,],aug21st.fire[ind2B,],xspecies,fire, fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire, fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Crouton.1hz.EF = tmpEF1 ;Crouton.5hz.EF = tmpEF5 }
    if (i > 1){Crouton.1hz.EF = rbind(Crouton.1hz.EF, tmpEF1) ; Crouton.5hz.EF = rbind(Crouton.5hz.EF, tmpEF5) }
    
  }
  Crouton.1hz.EF = Crouton.1hz.EF[order(Crouton.1hz.EF$variable),]
  Crouton.1hz.EF$transect_source_fire_ID = indA[1]
  Crouton.5hz.EF = Crouton.5hz.EF[order(Crouton.5hz.EF$variable),]
  Crouton.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== PO BOY =====================
  
  fire="Po_Boy"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    # Adjust endpoints
    if (i == 1){start = startO + 2; stop = stopO -10; startB = startO - 1; stopB = startO + 1}
    if (i == 2){start = startO + 14.6; stop = stopO -9; startB = startO + 0; stopB = startO + 2}
    # Time align  data
    tmp = time_align(start,stop,co2.821.5hz,  co.ch4.821.5hz, 
                       warneke.821.5hz,  isaf.821.5hz,
                       rollinsno.821.5hz, rollinsso2.821.5hz, 
                       cit.821.5hz, gtcims.821.5hz,moore.821fast, jimenez.821.5hz,met.821.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 

    if (i ==2){start=start-1}
    aug21st.fire = time_alignSLOWNOCANS(start,stop, co2.821.1hz, co.ch4.821.1hz, 
                            warneke.821.1hz, isaf.821.1hz,rollinsno.821.1hz, rollinsso2.821.1hz, 
                            cit.821.1hz,gtcims.821.1hz,ryerson.821.1hz , jimenez.821.1hz, 
                            schwarz.821.1hz, freid.c2h6.821.1hz, freid.ch2o.821.1hz, 
                            womack.821.1hz, stclair.821.1hz,veres.821.1hz, #wisthaler.821.1hz,
                            moore.821, met.821.1hz)
    ind2 = which(aug21st.fire$Time_Start >= start & aug21st.fire$Time_Start  <= stop	) 
    ind2B = which(aug21st.fire$Time_Start >= startB & aug21st.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug21st.fire$fire = ''; aug21st.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = fire
    aug21st.fire$fuel[ind2] = fuel ; aug21st.fire$fire[ind2] = fire
    # -------- Blake cans 
    indBLAKE = which(blake.821.merge$Time_Start >= (start-5) & blake.821.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.821.merge[indBLAKE,]
      blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.821.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.821.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.821.merge$Time_Start >= (start-0) & gilman.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.821.merge[indGILMAN,]
      GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.821.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.821.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.821.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.821.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.821.merge$Time_Start >= (start-0) & apel.821.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.821.merge[indAPEL,]
      APELBG = apel.821.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.821.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.821.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.821.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.821.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.821$Time_Start >= (start-0) & newTOGA.821$Time_Start <= stop & newTOGA.821$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.821[indBECKY,]
      BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.821$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.821[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.821[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.821[1,]
    }
    
    if (doplot == 1){
      plotpass5hzJUSTCO(aug21st.fire[ind2,])
      plotpass5hz(aug21st.fire[ind2,])
      plotpass1hz(aug21st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug21st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug21st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug21st.fire[ind2,],aug21st.fire[ind2B,],xspecies,fire, fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire, fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug21st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug21st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug21st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug21st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug21st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop     
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){PoBoy.1hz.EF = tmpEF1 ;PoBoy.5hz.EF = tmpEF5 }
    if (i > 1){PoBoy.1hz.EF = rbind(PoBoy.1hz.EF, tmpEF1) ; PoBoy.5hz.EF = rbind(PoBoy.5hz.EF, tmpEF5) }
    
  }
  PoBoy.1hz.EF = PoBoy.1hz.EF[order(PoBoy.1hz.EF$variable),]
  PoBoy.1hz.EF$transect_source_fire_ID = indA[1]
  PoBoy.5hz.EF = PoBoy.5hz.EF[order(PoBoy.5hz.EF$variable),]
  PoBoy.5hz.EF$transect_source_fire_ID = indA[1]
  
  # -------- all aug 8/21 together ------
   # -----****----- PLOT ALL FIRES: EFs for CO -----****-------
  par(mfrow=c(1,2))
  ind2 = which(LilDebbie.1hz.EF$variable =='CO_DACOM_DISKIN' & LilDebbie.1hz.EF$R2toX >= 0.9 )
  ind4 = which(Crouton.1hz.EF$variable =='CO_DACOM_DISKIN' & Crouton.1hz.EF$R2toX >= 0.9)
  ind6 = which(Ricearoni.1hz.EF$variable =='CO_DACOM_DISKIN' & Ricearoni.1hz.EF$R2toX >= 0.9)
  ind8 = which(Crawdaddy.1hz.EF$variable =='CO_DACOM_DISKIN' & Crawdaddy.1hz.EF$R2toX >= 0.9)
  ind10 = which(Gumbo.1hz.EF$variable =='CO_DACOM_DISKIN' & Gumbo.1hz.EF$R2toX >= 0.9)
  ind12 = which(Jumbalaya.1hz.EF$variable =='CO_DACOM_DISKIN' & Jumbalaya.1hz.EF$R2toX >= 0.9)
  ind14 = which(PoBoy.1hz.EF$variable =='CO_DACOM_DISKIN' & PoBoy.1hz.EF$R2toX >= 0.9)
  
  ind2B = which(LilDebbie.5hz.EF$variable =='CO_DACOM_DISKIN' & LilDebbie.5hz.EF$R2toX >= 0.9)
  ind4B = which(Crouton.5hz.EF$variable =='CO_DACOM_DISKIN' & Crouton.5hz.EF$R2toX >= 0.9)
  ind6B = which(Ricearoni.5hz.EF$variable =='CO_DACOM_DISKIN' & Ricearoni.5hz.EF$R2toX >= 0.9)
  ind8B = which(Crawdaddy.5hz.EF$variable =='CO_DACOM_DISKIN' & Crawdaddy.5hz.EF$R2toX >= 0.9)
  ind10B = which(Gumbo.5hz.EF$variable =='CO_DACOM_DISKIN' & Gumbo.5hz.EF$R2toX >= 0.9)
  ind12B = which(Jumbalaya.5hz.EF$variable =='CO_DACOM_DISKIN' & Jumbalaya.5hz.EF$R2toX >= 0.9)
  ind14B = which(PoBoy.5hz.EF$variable =='CO_DACOM_DISKIN' & PoBoy.5hz.EF$R2toX >= 0.9)
  # slope method #1
  if (doplot == 1){
    plot(LilDebbie.1hz.EF$EF1[ind2],LilDebbie.1hz.EF$mce[ind2], xlim=c(10,110), ylim=c(0.89, 1),
         main='8/03, corn: slope method, CO2 + CO + CH4', pch=19, xlab='CO EF, g/kg', ylab='mce')
    if (length(ind4)> 0) {points(Crouton.1hz.EF$EF1[ind4],Crouton.1hz.EF$mce[ind4], pch=2)}
    if (length(ind6)> 0) {points(Ricearoni.1hz.EF$EF1[ind6],Ricearoni.1hz.EF$mce[ind6],pch=4)}
    if (length(ind8)> 0) {points(Crawdaddy.1hz.EF$EF1[ind8],Crawdaddy.1hz.EF$mce[ind8],  pch=6)}
    if (length(ind10)> 0) {points(Gumbo.1hz.EF$EF1[ind10],Gumbo.1hz.EF$mce[ind10], pch=8)}
    if (length(ind12)> 0) {points(Jumbalaya.1hz.EF$EF1[ind12],Jumbalaya.1hz.EF$mce[ind12],  pch=10)}
    if (length(ind14)> 0) {points(PoBoy.1hz.EF$EF1[ind14],PoBoy.1hz.EF$mce[ind14], pch=12)}
  
    if (length(ind2B)> 0){ points(LilDebbie.5hz.EF$EF1[ind2B],LilDebbie.5hz.EF$mce[ind2B], col='red')}
    if (length(ind4B)> 0){ points(Crouton.5hz.EF$EF1[ind4B],Crouton.5hz.EF$mce[ind4B], col='red', pch=2)}
    if (length(ind6B)> 0){ points(Ricearoni.5hz.EF$EF1[ind6B],Ricearoni.5hz.EF$mce[ind6B], col='red', pch=4)}
    if (length(ind8B)> 0){ points(Crawdaddy.5hz.EF$EF1[ind8B],Crawdaddy.5hz.EF$mce[ind8B], col='red', pch=6)}
    if (length(ind10B)> 0){ points(Gumbo.5hz.EF$EF1[ind10B],Gumbo.5hz.EF$mce[ind10B], col='red', pch=8)}
    if (length(ind12B)> 0){ points(Jumbalaya.5hz.EF$EF1[ind12B],Jumbalaya.5hz.EF$mce[ind12B], col='red', pch=10)}
    if (length(ind14B)> 0){ points(PoBoy.5hz.EF$EF1[ind14B],PoBoy.5hz.EF$mce[ind14B], col='red', pch=12)}
    legend("bottomleft", legend=c("Lil Debbie","Crouton","Ricearoni",
                               "Crawdaddy","Gumbo","Jumbalaya","PoBoy"),
           pch=c(19,2,4,6,8,10,12), bty='n')
    legend("bottomright", legend=c("1Hz","5Hz"),pch=c(19,19), bty='n', col=c("black","red"))
  
      # integration method #1
    plot(LilDebbie.1hz.EF$EF1int[ind2],LilDebbie.1hz.EF$mce_int[ind2], xlim=c(10,110), ylim=c(0.89, 1),
         main='8/03, corn: integration method, CO2 + CO + CH4', pch=19, xlab='CO EF, g/kg', ylab='mce')
    if (length(ind4)> 0){ points(Crouton.1hz.EF$EF1int[ind4],Crouton.1hz.EF$mce_int[ind4], pch=2)}
    if (length(ind6)> 0){ points(Ricearoni.1hz.EF$EF1int[ind6],Ricearoni.1hz.EF$mce_int[ind6],pch=4)}
    if (length(ind8)> 0){ points(Crawdaddy.1hz.EF$EF1int[ind8],Crawdaddy.1hz.EF$mce_int[ind8],  pch=6)}
    if (length(ind10)> 0){ points(Gumbo.1hz.EF$EF1int[ind10],Gumbo.1hz.EF$mce_int[ind10], pch=8)}
    if (length(ind12)> 0){ points(Jumbalaya.1hz.EF$EF1int[ind12],Jumbalaya.1hz.EF$mce_int[ind12],  pch=10)}
    if (length(ind14)> 0){ points(PoBoy.1hz.EF$EF1int[ind14],PoBoy.1hz.EF$mce_int[ind14], pch=12)}
  
    if (length(ind2B)> 0){ points(LilDebbie.5hz.EF$EF1int[ind2B],LilDebbie.5hz.EF$mce_int[ind2B], col='red')}
    if (length(ind4B)> 0){ points(Crouton.5hz.EF$EF1int[ind4B],Crouton.5hz.EF$mce_int[ind4B], col='red', pch=2)}
    if (length(ind6B)> 0){ points(Ricearoni.5hz.EF$EF1int[ind6B],Ricearoni.5hz.EF$mce_int[ind6B], col='red', pch=4)}
    if (length(ind8B)> 0){ points(Crawdaddy.5hz.EF$EF1int[ind8B],Crawdaddy.5hz.EF$mce_int[ind8B], col='red', pch=6)}
    if (length(ind10B)> 0){ points(Gumbo.5hz.EF$EF1int[ind10B],Gumbo.5hz.EF$mce_int[ind10B], col='red', pch=8)}
    if (length(ind12B)> 0){ points(Jumbalaya.5hz.EF$EF1int[ind12B],Jumbalaya.5hz.EF$mce_int[ind12B], col='red', pch=10)}
    if (length(ind14B)> 0){ points(PoBoy.5hz.EF$EF1int[ind14B],PoBoy.5hz.EF$mce_int[ind14B], col='red', pch=12)}
  }    
  # ----- [[[[[[[[[[[[[[[[[[[[ Aug 23rd ]]]]]]]]]]]]]]]]]]]]] -------
  # --------#######--------- Get 1 Hz Data individual 8/23------#######---------
  # plume tags
  tags = getICARTTdataSIMPLE('InputFiles/firexaq-fire-Flags-1HZ_DC8_20190823_R9.ict') ;tags$Time_Start = tags$TIME_START
  # MET DATA
  met.823.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-MetNav_DC8_20190823_R1.ict')
  met.823.1hz = merge(met.823.1hz, tags, by='Time_Start')
  # CO2
  co2.823.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CO2-7000_DC8_20190823_R2.ict')
  # -------- DISKIN -----CO, CH4
  co.ch4.823.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-DACOM_DC8_20190823_R1.ict')
  # --------- WARNEKE ----  VOCs
  warneke.823.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-NOAAPTR-VOCs-1Hz_DC8_20190823_R3.ict')
  # ------ HANISCO - ISAF HCHO - merged to 5Hz from the online merge
  isaf.823.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-ISAF-CH2O-1Hz_DC8_20190823_R0.ict')
  
  #  ------- ROLLINS - SO2 and NO
  rollinsno.823.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-NO_DC8_20190823_R1.ict')
  rollinsno.823.1hz$Time_Start = rollinsno.823.1hz$time_mid
  rollinsso2.823.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-SO2_DC8_20190823_R1.ict')
  rollinsso2.823.1hz$Time_Start = rollinsso2.823.1hz$time_mid
  
  #  ----- WENNBERG - CIT VOCs - 
  cit.823.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190823_R0_CIT.ict')
  # ------ HUEY - GTCIMS PANs - not sure how to match up peaks here
  gtcims.823.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190823_R0_Huey.ict')
  # ------ RYERSON
  ryerson.A = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NO_DC8_20190823_R1.ict')
  ryerson.B = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NO2_DC8_20190823_R1.ict')
  ryerson.C = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NOy_DC8_20190823_R1.ict')
  ryerson.D = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-O3_DC8_20190823_R1.ict')
  ryerson.823.1hz = cbind(ryerson.A,ryerson.B,ryerson.C,ryerson.D) 
  ryerson.823.1hz$Time_Start = ryerson.823.1hz$Time_start
  # ----- JIMENEZ ---  # ----- JIMENEZ ---
  jimenez.823.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-EESI_DC8_20190823_R1.ict')
  
  # ----- SCHWARZ ---
  schwarz.823.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-SP2-BC-1HZ_DC8_20190823_R2.ict')
#  schwarz.823.5hz = read.table('InputFiles/SP2_5Hz_FIREXAQ_20190823_QL.txt', header=TRUE, sep=',')
#  schwarz.823.5hz$Time_Start = schwarz.823.5hz$bin_Htimes
  # ----- FREID ---
  freid.c2h6.823.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-C2H6_DC8_20190823_R3.ict')
  freid.ch2o.823.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CH2O_DC8_20190823_R3.ict')
  # ------ WOMACK ---
  womackA = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-CH3COCHO_DC8_20190823_R1.ict')
  womackB = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-CHOCHO_DC8_20190823_R1.ict')
  womackC = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-HNO2_DC8_20190823_R1.ict')
  womackD = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-NO2_DC8_20190823_R1.ict')
  womack.823.1hz = cbind(womackA, womackB, womackC, womackD)
  # -------St Clair
  stclair.823.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-CANOE-NO2_DC8_20190823_R0.ict')
  stclair.823.1hz$Time_Start = stclair.823.1hz$Time_start
  # data is all missing, kludge for now
  stclair.823.1hz = womack.823.1hz
  stclair.823.1hz$NO2_CANOE = womack.823.1hz$CH3COCHO_ACES*NaN
  stclair.823.1hz=stclair.823.1hz[,c(1,9)]
  # ------- VERES
  veres.A = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-ClNO2_DC8_20190823_R0.ict')
  veres.B = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HCOOH_DC8_20190823_R1.ict')
  veres.C = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HNO2_DC8_20190823_R1.ict')
  veres.D = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-N2O5_DC8_20190823_R0.ict')
  veres.E = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HPMTF_DC8_20190823_R0.ict')
  veres.F = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-CH3COOCl_DC8_20190823_R0.ict')
  veres.G = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-Cl2_DC8_20190823_R0.ict')
  veres.H = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrCl_DC8_20190823_R0.ict')
  veres.I = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrCN_DC8_20190823_R0.ict')
  veres.J = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrO_DC8_20190823_R0.ict')
  veres.K = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HCN_DC8_20190823_R0.ict')
  veres.L = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HNCO_DC8_20190823_R0.ict')
  veres.823.1hz = cbind(veres.A,veres.B,veres.C,veres.D,veres.E,veres.F,veres.G,veres.H,veres.I,veres.J,veres.K,veres.L)
  
  # --- WISTHALER
  #wisthaler.823.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190823_R0_Wisthaler.ict')
  
  # ---- BLAKE
  blake.823.1hz = getICARTTdataSIMPLE('InputFiles/WAS-MERGE/firexaq-mrgWAS-dc8_merge_20190823_R1.ict')
  cc = colnames(blake.823.1hz)
  blake.823.merge = blake.823.1hz[,c(1,2,96:225)]
  blake.823.merge$CO_DACOM_DISKIN_BLAKE = blake.823.1hz$CO_DACOM_DISKIN
  blake.823.merge$CO2_7000_ppm_DISKIN_BLAKE = blake.823.1hz$CO2_7000_ppm_DISKIN
  
  # ------ APEL
  apel.823.1hz = getICARTTdataSIMPLE('InputFiles/TOGA-MERGE/firexaq-mrgTOGA-dc8_merge_20190823_R1.ict')
  cc = colnames(apel.823.1hz)
  apel.823.merge = apel.823.1hz[,c(1,2,226:315)]
  apel.823.merge$CO_DACOM_DISKIN_APEL = apel.823.1hz$CO_DACOM_DISKIN
  apel.823.merge$CO2_7000_ppm_DISKIN_APEL =apel.823.1hz$CO2_7000_ppm_DISKIN
  
  # Becky's better merge
  file = 'InputFiles/Hornbrook/FIREX-AQ weighted TOGA merge 2022-01-24_0823.xlsx'
  newTOGA.823 = readxl::read_xlsx(file); newTOGA.823[newTOGA.823==-999] = NaN; newTOGA.823[newTOGA.823==-888] = NaN
  newTOGA.823$CO_DACOM_DISKIN_BECKY = newTOGA.823$CO_DACOM_DISKIN
  newTOGA.823$CO2_7000_ppm_DISKIN_BECKY = NaN
  newTOGA.823$Time_Start=newTOGA.823$Time_Start...4

  # ----GILMAN
  gilman.823.1hz = getICARTTdataSIMPLE('InputFiles/iWAS-MERGE/firexaq-mrgiWAS-dc8_merge_20190823_R1.ict')
  cc = colnames(gilman.823.1hz)
  gilman.823.merge = gilman.823.1hz[,c(1,2,316:361)]
  gilman.823.merge$CO_DACOM_DISKIN_GILMAN = gilman.823.1hz$CO_DACOM_DISKIN
  gilman.823.merge$CO2_7000_ppm_DISKIN_GILMAN = gilman.823.1hz$CO2_7000_ppm_DISKIN

  # ------ Moore 
  moore.823fast = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190823_R0_MOORE.ict')
  
  moore.823p1 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-AerosolCloudConc_DC8_20190823_R0.ict')
  moore.823p2 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-LAScold_DC8_20190823_R0.ict')
  moore.823p3 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-LAShot_DC8_20190823_R0.ict')
  moore.823p4 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-CPSPD_DC8_20190823_R0.ict')
  moore.823p5 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-CDP_DC8_20190823_R0.ict')
  moore.823p6 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-SMPS_DC8_20190823_R0.ict')
  moore.823 =merge(moore.823p1, moore.823p2, by='Time_mid', all = TRUE, incomparables = NA)
  moore.823 =merge(moore.823, moore.823p3, by='Time_mid', all = TRUE, incomparables = NA)
  moore.823 =merge(moore.823, moore.823p4, by='Time_mid', all = TRUE, incomparables = NA)
  moore.823 =merge(moore.823, moore.823p5, by='Time_mid', all = TRUE, incomparables = NA)
  moore.823 =merge(moore.823, moore.823p6, by='Time_mid', all = TRUE, incomparables = NA)
  
  # ------- append PI to colnames 1hz ----------
  cc = colnames(co2.823.1hz)
  cc[2] = paste(cc[2],'_DISKIN',sep='')
  colnames(co2.823.1hz) = cc 
  cc = colnames(co.ch4.823.1hz)
  cc[2:4] =  paste(cc[2:4],'_DISKIN',sep='')
  colnames(co.ch4.823.1hz) = cc
  cc=colnames(met.823.1hz)
  
  colnames(isaf.823.1hz) = c("Time_Start"  , "CH2O_ISAF_HANISCO" ,"CH2O_ISAF_precision_HANISCO")
  cc[2:36] = paste(cc[2:36],'_YANG',sep='')
  colnames(met.823.1hz) = cc
  cc = colnames(warneke.823.1hz)
  cc[2:43] = paste(cc[2:43],'_WARNEKE', sep='')
  colnames(warneke.823.1hz) = cc
  cc = colnames(rollinsno.823.1hz)
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsno.823.1hz) = cc
  cc = colnames(rollinsso2.823.1hz)
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsso2.823.1hz) = cc
  # make GTCIMS consistent with 1s merge
  cc=colnames(gtcims.823.1hz)
  cc[4:7] = c("PAN_GTCIMS_HUEY" , "PPN_GTCIMS_HUEY"  ,"APAN_GTCIMS_HUEY" ,"PBN_GTCIMS_HUEY")
  colnames(gtcims.823.1hz)=cc
  # since ISAF, CIT, and GTCIMS came from merge tool, alread has PI's appended.
  ryerson.823.1hz <- ryerson.823.1hz[, !duplicated(colnames(ryerson.823.1hz))]
  
  cc = colnames(ryerson.823.1hz)
  cc[2:9] = paste(cc[2:9],'_RYERSON',sep='')
  colnames(ryerson.823.1hz)=cc
  
  cc = colnames(schwarz.823.1hz)
  cc[2:3] = paste(cc[2:3],'_SCHWARZ', sep='')
  colnames(schwarz.823.1hz) =cc
  
  cc = colnames(freid.c2h6.823.1hz)
  cc[4:7] = paste(cc[4:7], '_FRIED',sep='')
  colnames(freid.c2h6.823.1hz) = cc
  
  cc = colnames(freid.ch2o.823.1hz)
  cc[4:7] = paste(cc[4:7], '_FRIED',sep='')
  colnames(freid.ch2o.823.1hz) = cc
  
  colnames(freid.ch2o.823.1hz) = cc
  womack.823.1hz <- womack.823.1hz[, !duplicated(colnames(womack.823.1hz))]
  cc = colnames(womack.823.1hz)
  cc[2:5] = paste(cc[2:5], '_WOMACK',sep='')
  colnames(womack.823.1hz) = cc

  veres.823.1hz <- veres.823.1hz[, !duplicated(colnames(veres.823.1hz))]
  cc = colnames(veres.823.1hz)
  cc[2:13] = paste(cc[2:13], '_VERES',sep='')
  colnames(veres.823.1hz) = cc
  
  #cc = colnames(wisthaler.823.1hz)
  #cc[4:5] = paste(cc[4:5], '_WISTHALER',sep='')
  #colnames(wisthaler.823.1hz) = cc
  cc = colnames(jimenez.823.1hz)
  cc[2:8] = paste(cc[2:8], '_JIMENEZ',sep='')
  colnames(jimenez.823.1hz) = cc
  # --------#######--------- Get 5 or 10 Hz Data 8/23 ------#######---------
  met.823.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190823_R0_met.ict')
  
  # CO2
  co2.823.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CO2-7000-5Hz_DC8_20190823_R1.ict')
  #CO, CH4
  co.ch4.823.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-DACOM-5Hz_DC8_20190823_R1.ict')
  # WARNEKE VOCs
  warneke.823.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-NOAAPTR-VOCs-5Hz_DC8_20190823_R3.ict')
  # ISAF HCHO - merged to 5Hz from the online merge
  isaf.823.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190823_R0_ISAF.ict')
  # ROLLINS SO2 and NO
  rollinsno.823.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-NO-5Hz_DC8_20190823_R0.ict')
  rollinsso2.823.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-SO2-5Hz_DC8_20190823_R1.ict')
  # CIT VOCs - 
  cit.823.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190823_R0_CIT.ict')
  # GTCIMS PANs - not sure how to match up peaks here
  gtcims.823.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190823_R0_huey.ict')
  
  # ----- Jimenez ---
  jimenez.823.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_AMS_20190823_R0_20230314T134036.ict')
  
  jimenez.823.5hz$OC_PM1_AMS_JIMENEZ = jimenez.823.5hz$OA_PM1_AMS_JIMENEZ/jimenez.823.5hz$OAtoOC_PM1_AMS
  
  # ------- append PI to colnames ----------
  cc = colnames(co2.823.5hz)
  cc[2] = paste(cc[2],'_DISKIN',sep='')
  colnames(co2.823.5hz) = cc 
  cc = colnames(co.ch4.823.5hz)
  cc[2:4] =  paste(cc[2:4],'_DISKIN',sep='')
  colnames(co.ch4.823.5hz) = cc
  cc=colnames(met.823.5hz)
  cc[2:36] = paste(cc[2:36],'_YANG',sep='')
  colnames(met.823.5hz) = cc
  cc = colnames(warneke.823.5hz)
  cc[2:43] = paste(cc[2:43],'_WARNEKE', sep='')
  colnames(warneke.823.5hz) = cc
  cc = colnames(rollinsno.823.5hz)
  cc[1] = 'Time_Start'
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsno.823.5hz) = cc
  cc = colnames(rollinsso2.823.5hz)
  cc[1] = 'Time_Start'
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsso2.823.5hz) = cc
  # make GTCIMS consistent with 1s merge
  cc=colnames(gtcims.823.5hz)
  cc[38:39] = c("APAN_GTCIMS_HUEY" , "PAN_GTCIMS_HUEY")
  colnames(gtcims.823.5hz)=cc
  # since ISAF, CIT, and GTCIMS came from merge tool, alread has PI's appended.
  # ------- get fuel moisture data --------
  if (doFM == 1){
    f1 = '/Users/ktravis1/Library/CloudStorage/Box-Box/FuelMoisture/fuel_moisture_content-20210715T1049Z/fmc_20190823_20Z.nc'
    fid = nc_open(f1)
    fuelMDead = ncvar_get(fid, varid = 'FMCG2D')
    fuelMLive = ncvar_get(fid, varid = 'FMCGLH2D')
    xlon = ncvar_get(fid, varid="XLONG_M")
    xlat = ncvar_get(fid, varid="XLAT_M")
    nc_close(fid)
  }
  # ========================== ANT  =====================
 
  fire="Ant"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 2;stop = stopO -4.6; startB = startO -1; stopB = startO + 1}
    if (i == 2){start = startO + 6.2;stop = stopO -12; startB = startO -1; stopB = startO +1}
    # Time align 1hz data
    tmp = time_align(start,stop,co2.823.5hz,  co.ch4.823.5hz, 
                     warneke.823.5hz,  isaf.823.5hz,
                     rollinsno.823.5hz, rollinsso2.823.5hz, 
                     cit.823.5hz, gtcims.823.5hz,moore.823fast, jimenez.823.5hz,met.823.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i ==1){ start=start-1;stop = stop+3; startB=startB-1; stopB=stopB-1}
    if (i ==2){ start=start-2;stop = stop+0}
    aug23rd.fire = time_alignSLOWNOCANS(start,stop, co2.823.1hz, co.ch4.823.1hz, 
                            warneke.823.1hz, isaf.823.1hz,rollinsno.823.1hz, rollinsso2.823.1hz, 
                            cit.823.1hz,gtcims.823.1hz,ryerson.823.1hz , jimenez.823.1hz, 
                            schwarz.823.1hz, freid.c2h6.823.1hz, freid.ch2o.823.1hz, 
                            womack.823.1hz, stclair.823.1hz,veres.823.1hz, #wisthaler.823.1hz,
                            moore.823, met.823.1hz)
    
    ind2 = which(aug23rd.fire$Time_Start >= start & aug23rd.fire$Time_Start  <= stop	) 
    ind2B = which(aug23rd.fire$Time_Start >= startB & aug23rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug23rd.fire$fire = ''; aug23rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Ant'
    aug23rd.fire$fuel[ind2] = fuel ; aug23rd.fire$fire[ind2] = 'Ant'
    
    # -------- Blake cans 
    indBLAKE = which(blake.823.merge$Time_Start >= (start-5) & blake.823.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.823.merge[indBLAKE,]
      blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.823.merge$Time_Start >= (start-0) & gilman.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.823.merge[indGILMAN,]
      GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.823.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.823.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.823.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.823.merge$Time_Start >= (start-0) & apel.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.823.merge[indAPEL,]
      APELBG = apel.823.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.823.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG =apel.823.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.823.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.823.merge[1,]
    }
    
    # Becky
    indBECKY = which(newTOGA.823$Time_Start >= (start-0) & newTOGA.823$Time_Start <= stop & newTOGA.823$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.823[indBECKY,]
      BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.823[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.823[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug23rd.fire[ind2,])
      plotpass5hz(aug23rd.fire[ind2,])
      plotpass1hz(aug23rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass1hzBLAKE(aug31st.fire[ind2,])
      
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug23rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug23rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug23rd.fire[ind2,],aug23rd.fire[ind2B,],xspecies,'Ant', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Ant', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Ant.1hz.EF = tmpEF1 ;Ant.5hz.EF = tmpEF5 }
    if (i > 1){Ant.1hz.EF = rbind(Ant.1hz.EF, tmpEF1) ; Ant.5hz.EF = rbind(Ant.5hz.EF, tmpEF5) }
  }
  Ant.1hz.EF = Ant.1hz.EF[order(Ant.1hz.EF$variable),]
  Ant.1hz.EF$transect_source_fire_ID = indA[1]
  Ant.5hz.EF = Ant.5hz.EF[order(Ant.5hz.EF$variable),]
  Ant.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== BLANKET  =====================

  fire="Blanket"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 2.4;stop = stopO -4; startB = startO -2; stopB = startO +0} # low CO
    if (i == 2){start = startO + 3;stop = stopO -10; startB = startO + 0; stopB = startO+2}
    if (i == 3){start = startO + 7;stop = stopO -7; startB = startO -1; stopB = startO+1}
    # Time align 1hz data
    tmp = time_align(start,stop,co2.823.5hz,  co.ch4.823.5hz, 
                     warneke.823.5hz,  isaf.823.5hz,
                     rollinsno.823.5hz, rollinsso2.823.5hz, 
                     cit.823.5hz, gtcims.823.5hz,moore.823fast, jimenez.823.5hz,met.823.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i ==1){ start = start - 1; startB=startB-1;stopB=stopB-1}
    if (i ==3){ start = start - 1; startB=startB-1;stopB=stopB-1}
    aug23rd.fire = time_alignSLOWNOCANS(start,stop, co2.823.1hz, co.ch4.823.1hz, 
                                        warneke.823.1hz, isaf.823.1hz,rollinsno.823.1hz, rollinsso2.823.1hz, 
                                        cit.823.1hz,gtcims.823.1hz,ryerson.823.1hz , jimenez.823.1hz, 
                                        schwarz.823.1hz, freid.c2h6.823.1hz, freid.ch2o.823.1hz, 
                                        womack.823.1hz, stclair.823.1hz,veres.823.1hz, #wisthaler.823.1hz,
                                        moore.823, met.823.1hz)
    ind2 = which(aug23rd.fire$Time_Start >= start & aug23rd.fire$Time_Start  <= stop	) 
    ind2B = which(aug23rd.fire$Time_Start >= startB & aug23rd.fire$Time_Start <= stopB	)
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug23rd.fire$fire = ''; aug23rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Blanket'
    aug23rd.fire$fuel[ind2] = fuel ; aug23rd.fire$fire[ind2] = 'Blanket'
    # -------- Blake cans 
    indBLAKE = which(blake.823.merge$Time_Start >= (start-5) & blake.823.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.823.merge[indBLAKE,]
      blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.823.merge$Time_Start >= (start-0) & gilman.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.823.merge[indGILMAN,]
      GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.823.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.823.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.823.merge[1,]
    }
    
    # Apel  
    indAPEL = which( apel.823.merge$Time_Start >= (start-0) &  apel.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.823.merge[indAPEL,]
      APELBG = apel.823.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.823.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.823.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.823.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.823$Time_Start >= (start-0) & newTOGA.823$Time_Start <= stop & newTOGA.823$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.823[indBECKY,]
      BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.823[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.823[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug23rd.fire[ind2,])
      plotpass5hz(aug23rd.fire[ind2,])
      plotpass1hz(aug23rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug23rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug23rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug23rd.fire[ind2,],aug23rd.fire[ind2B,],xspecies,'Blanket', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Blanket', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Blanket.1hz.EF = tmpEF1 ;Blanket.5hz.EF = tmpEF5 }
    if (i > 1){Blanket.1hz.EF = rbind(Blanket.1hz.EF, tmpEF1) ; Blanket.5hz.EF = rbind(Blanket.5hz.EF, tmpEF5) }
  }
  Blanket.1hz.EF = Blanket.1hz.EF[order(Blanket.1hz.EF$variable),]
  Blanket.1hz.EF$transect_source_fire_ID = indA[1]
  Blanket.5hz.EF = Blanket.5hz.EF[order(Blanket.5hz.EF$variable),]
  Blanket.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== CHIPS  =====================
  # example for paper
  start = 68225
  stop= 68232
  ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
  plotpass5hzJUSTCO(tmp[ind,])
  fire="Chips"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    
    if (i == 1){start = startO + 8;stop = stopO -8; startB = startO + 0; stopB = startO + 2}
    if (i == 2){start = startO + 2;stop = stopO -8.4; startB = startO -1; stopB = startO + 1}
    if (i == 3){start = startO + .8;stop = stopO -5; startB = 68225; stopB = 68226} # split from above
    if (i == 4){start = startO + 0;stop = stopO -0; startB = startO + 0; stopB = stopO - 0} # bad pass, CO < 100
    if (i == 5){start = startO + 19;stop = stopO -0; startB = startO + 0; stopB =  startO+2} #
    if (i == 6){start = startO + 0;stop = stopO -5; startB = 68689; stopB =  68691} # split from above
    if (i == 7){start = startO + 1;stop = stopO -34; startB = startO -2; stopB = startO +0}
    if (i == 8){start = startO + 34.4;stop = stopO -6; startB = startO + 0; stopB = startO + 2}
    if (i == 9){start = startO +9;stop = stopO -0; startB = startO + 0; stopB = startO + 2}
    if (i == 10){start = startO +.2;stop = stopO -30; startB = 69632; stopB = 69634} # split from above
    #if (i == 11){start = startO + 1;stop = stopO -5; startB = 69127; stopB = 69129}# split from above
    #if (i == 11){start = startO + 8;stop = stopO -66; startB = startO + 0; stopB = startO + 2}
    tmp = time_align(start,stop,co2.823.5hz,  co.ch4.823.5hz, 
                     warneke.823.5hz,  isaf.823.5hz,
                     rollinsno.823.5hz, rollinsso2.823.5hz, 
                     cit.823.5hz, gtcims.823.5hz,moore.823fast, jimenez.823.5hz,met.823.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){ start = start -1; stop=stop+0}
    if (i ==2){ start = start -2; stop=stop+1}
    if (i ==6){ start = start -0; stop=stop-1}
    if (i ==7){ start = start -2; stop=stop-0}
    if (i ==8){ start=start-1;startB = startB -1; stopB=stopB-1}
    if (i ==9){ start = start -0; stop=stop-0}
    #if (i ==11){ start = start -1; stop=stop+0}
    aug23rd.fire = time_alignSLOWNOCANS(start,stop, co2.823.1hz, co.ch4.823.1hz, 
                                        warneke.823.1hz, isaf.823.1hz,rollinsno.823.1hz, rollinsso2.823.1hz, 
                                        cit.823.1hz,gtcims.823.1hz,ryerson.823.1hz , jimenez.823.1hz, 
                                        schwarz.823.1hz, freid.c2h6.823.1hz, freid.ch2o.823.1hz, 
                                        womack.823.1hz, stclair.823.1hz,veres.823.1hz, #wisthaler.823.1hz,
                                        moore.823, met.823.1hz)
    ind2 = which(aug23rd.fire$Time_Start >= start & aug23rd.fire$Time_Start  <= stop	) 
    ind2B = which(aug23rd.fire$Time_Start >= startB & aug23rd.fire$Time_Start <= stopB	) 
    # invetsigate for i=8
    #plot(aug23rd.fire$Time_Start[c(ind2[1]-7,ind2[1]-6,ind2[1]-5,ind2[1]-4,ind2[1]-3,ind2[1]-2,ind2[1]-1,ind2)], aug23rd.fire$CO_DACOM_DISKIN[c(ind2[1]-7,ind2[1]-6,ind2[1]-5,ind2[1]-4,ind2[1]-3,ind2[1]-2,ind2[1]-1,ind2)], type='o')
    #points(aug23rd.fire$Time_Start[ind2B], aug23rd.fire$CO_DACOM_DISKIN[ind2B], col='red', pch=19)
    #points(GILMAN$Time_Start, GILMAN$CO_DACOM_DISKIN_GILMAN, col='green')
    #points(aug23rd.fire$Time_Start[ind2B[1]], GILMANBG$CO_DACOM_DISKIN_GILMAN, col='green')
     
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug23rd.fire$fire = ''; aug23rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Chips'
    aug23rd.fire$fuel[ind2] = fuel ; aug23rd.fire$fire[ind2] = 'Chips'
    # try and get a background as close to the plume as possible?
    # -------- Blake cans 
    indBLAKE = which(blake.823.merge$Time_Start >= (start-5) & blake.823.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.823.merge[indBLAKE,]
      blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary, obs can't have BG greater than value
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 | blakeBG$CO_DACOM_DISKIN_BLAKE[1] > blake$CO_DACOM_DISKIN_BLAKE[1] ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.823.merge$Time_Start >= (start-0) & gilman.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.823.merge[indGILMAN,]
      GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.823.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.823.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.823.merge[1,]
    }
    
    # Apel  
    indAPEL = which( apel.823.merge$Time_Start >= (start-0) &  apel.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.823.merge[indAPEL,]
      APELBG = apel.823.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.823.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.823.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.823.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.823$Time_Start >= (start-0) & newTOGA.823$Time_Start <= stop & newTOGA.823$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.823[indBECKY,]
      BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.823[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.823[1,]
    }
   
    if (doplot == 1){
      plotpass5hzJUSTCO(aug23rd.fire[ind2,])
      plotpass5hz(aug23rd.fire[ind2,])
      plotpass1hz(aug23rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug23rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug23rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug23rd.fire[ind2,],aug23rd.fire[ind2B,],xspecies,'Chips', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG,doBECKY)
    tmpEF1$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Chips', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG,doBECKY)
    tmpEF5$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop 
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Chips.1hz.EF = tmpEF1 ;Chips.5hz.EF = tmpEF5 }
    if (i > 1){Chips.1hz.EF = rbind(Chips.1hz.EF, tmpEF1) ; Chips.5hz.EF = rbind(Chips.5hz.EF, tmpEF5) }
  }
  Chips.1hz.EF = Chips.1hz.EF[order(Chips.1hz.EF$variable),]
  Chips.1hz.EF$transect_source_fire_ID = indA[1]
  Chips.5hz.EF = Chips.5hz.EF[order(Chips.5hz.EF$variable),]
  Chips.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== DIP =====================
 
  fire="Dip";fuel="pile"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) ) 
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 1;stop = stop -5; startB = startO -1; stopB = startO + 1}
    if (i == 2){start = startO + 1;stop = stop -0; startB = startO -1; stopB = startO + 1}
    if (i == 3){start = startO + 1;stop = stop -0; startB = startO -1; stopB = startO + 1}
    if (i == 4){start = startO + 14.4;stop = stop -15.8; startB = startO -1; stopB = startO + 1}
    tmp = time_align(start,stop,co2.823.5hz,  co.ch4.823.5hz, 
                     warneke.823.5hz,  isaf.823.5hz,
                     rollinsno.823.5hz, rollinsso2.823.5hz, 
                     cit.823.5hz, gtcims.823.5hz,moore.823fast, jimenez.823.5hz,met.823.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1}
    if (i ==4){start=start-1}
    
    aug23rd.fire = time_alignSLOWNOCANS(start,stop, co2.823.1hz, co.ch4.823.1hz, 
                                        warneke.823.1hz, isaf.823.1hz,rollinsno.823.1hz, rollinsso2.823.1hz, 
                                        cit.823.1hz,gtcims.823.1hz,ryerson.823.1hz , jimenez.823.1hz, 
                                        schwarz.823.1hz, freid.c2h6.823.1hz, freid.ch2o.823.1hz, 
                                        womack.823.1hz, stclair.823.1hz,veres.823.1hz, #wisthaler.823.1hz,
                                        moore.823, met.823.1hz)
    ind2 = which(aug23rd.fire$Time_Start >= start & aug23rd.fire$Time_Start  <= stop	) 
    ind2B = which(aug23rd.fire$Time_Start >= startB & aug23rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug23rd.fire$fire = ''; aug23rd.fire$fuel = ''
    tmp$fuel[ind] = 'pile' ; tmp$fire[ind] = 'Dip'
    aug23rd.fire$fuel[ind2] = 'pile' ; aug23rd.fire$fire[ind2] = 'Dip'
    # -------- Blake cans 
    indBLAKE = which(blake.823.merge$Time_Start >= (start-5) & blake.823.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.823.merge[indBLAKE,]
      blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.823.merge$Time_Start >= (start-0) & gilman.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      if (length(indGILMAN) > 1 & i == 2){indGILMAN = indGILMAN[1]}
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.823.merge[indGILMAN,]
      GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.823.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.823.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.823.merge[1,]
    }
    
    # Apel  
    indAPEL = which( apel.823.merge$Time_Start >= (start-0) &  apel.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.823.merge[indAPEL,]
      APELBG = apel.823.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.823.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.823.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.823.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.823$Time_Start >= (start-0) & newTOGA.823$Time_Start <= stop & newTOGA.823$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.823[indBECKY,]
      BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.823[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.823[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug23rd.fire[ind2,])
      plotpass5hz(aug23rd.fire[ind2,])
      plotpass1hz(aug23rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug23rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug23rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug23rd.fire[ind2,],aug23rd.fire[ind2B,],xspecies,'Dip', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Dip', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Dip.1hz.EF = tmpEF1 ;Dip.5hz.EF = tmpEF5 }
    if (i > 1){Dip.1hz.EF = rbind(Dip.1hz.EF, tmpEF1) ; Dip.5hz.EF = rbind(Dip.5hz.EF, tmpEF5) }
  }
  Dip.1hz.EF = Dip.1hz.EF[order(Dip.1hz.EF$variable),]
  Dip.1hz.EF$transect_source_fire_ID = indA[1]
  Dip.5hz.EF = Dip.5hz.EF[order(Dip.5hz.EF$variable),]
  Dip.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== ESCARGOT =====================
  fire="Escargot"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) ) 
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 1;stop = stopO -.2; startB = startO -1; stopB =  startO + 1}
    if (i == 2){start = startO + 0;stop = stopO -7; startB = 70003; stopB = 70005} #split from above
    if (i == 3){start = startO + 1;stop = stopO -3; startB = startO + 0; stopB = startO + 2}
    if (i == 4){start = startO + 6.2;stop = stopO -10; startB = startO + 0; stopB = startO + 2}
    if (i == 5){start = startO + 24.4;stop = stopO -55.2; startB = startO + 0; stopB = startO + 2}
    if (i == 6){start = startO + 3.8;stop = stopO -42; startB = startO + 0; stopB = startO+2}
    tmp = time_align(start,stop,co2.823.5hz,  co.ch4.823.5hz, 
                     warneke.823.5hz,  isaf.823.5hz,
                     rollinsno.823.5hz, rollinsso2.823.5hz, 
                     cit.823.5hz, gtcims.823.5hz,moore.823fast, jimenez.823.5hz,met.823.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i ==1){start=start-1}
    if (i ==2){ start = start - 1}
    if (i == 3){stop=stop+2}
    if (i ==4){ start = start - 1}
    if (i ==6){ start = start - 1}
    aug23rd.fire = time_alignSLOWNOCANS(start,stop, co2.823.1hz, co.ch4.823.1hz, 
                                        warneke.823.1hz, isaf.823.1hz,rollinsno.823.1hz, rollinsso2.823.1hz, 
                                        cit.823.1hz,gtcims.823.1hz,ryerson.823.1hz , jimenez.823.1hz, 
                                        schwarz.823.1hz, freid.c2h6.823.1hz, freid.ch2o.823.1hz, 
                                        womack.823.1hz, stclair.823.1hz,veres.823.1hz, #wisthaler.823.1hz,
                                        moore.823, met.823.1hz)
    ind2 = which(aug23rd.fire$Time_Start >= start & aug23rd.fire$Time_Start  <= stop	) 
    ind2B = which(aug23rd.fire$Time_Start >= startB & aug23rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug23rd.fire$fire = ''; aug23rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Escargot'
    aug23rd.fire$fuel[ind2] = fuel ; aug23rd.fire$fire[ind2] = 'Escargot'
    # -------- Blake cans 
    indBLAKE = which(blake.823.merge$Time_Start >= (start-5) & blake.823.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.823.merge[indBLAKE,]
      blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.823.merge$Time_Start >= (start-0) & gilman.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.823.merge[indGILMAN,]
      GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.823.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.823.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.823.merge[1,]
    }
    
    # Apel  
    indAPEL = which( apel.823.merge$Time_Start >= (start-0) &  apel.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.823.merge[indAPEL,]
      APELBG = apel.823.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.823.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.823.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.823.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.823$Time_Start >= (start-0) & newTOGA.823$Time_Start <= stop & newTOGA.823$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.823[indBECKY,]
      BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.823[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.823[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug23rd.fire[ind2,])
      plotpass5hz(aug23rd.fire[ind2,])
      plotpass1hz(aug23rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug23rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug23rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug23rd.fire[ind2,],aug23rd.fire[ind2B,],xspecies,'Escargot', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Escargot', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop   
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Escargot.1hz.EF = tmpEF1 ;Escargot.5hz.EF = tmpEF5 }
    if (i > 1){Escargot.1hz.EF = rbind(Escargot.1hz.EF, tmpEF1) ; Escargot.5hz.EF = rbind(Escargot.5hz.EF, tmpEF5) }
    
  }
  Escargot.1hz.EF = Escargot.1hz.EF[order(Escargot.1hz.EF$variable),]
  Escargot.1hz.EF$transect_source_fire_ID = indA[1]
  Escargot.5hz.EF = Escargot.5hz.EF[order(Escargot.5hz.EF$variable),]
  Escargot.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== FRISBEE =====================
  fire="Frisbee"; fuel="soybean"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) ) 
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 4;stop = stopO -24; startB = startO + 0; stopB = startO + 2}
    if (i == 2){start = startO + 2;stop = stopO -2; startB = startO + 0; stopB =startO + 2}
    if (i == 3){start = startO + 0.2;stop = stopO -10; startB = 71508; stopB = 71510} # split from above
    tmp = time_align(start,stop,co2.823.5hz,  co.ch4.823.5hz, 
                     warneke.823.5hz,  isaf.823.5hz,
                     rollinsno.823.5hz, rollinsso2.823.5hz, 
                     cit.823.5hz, gtcims.823.5hz,moore.823fast, jimenez.823.5hz,met.823.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i == 1){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i == 2){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i == 3){start=start-1;startB=startB-1;stopB=stopB-1}
    aug23rd.fire = time_alignSLOWNOCANS(start,stop, co2.823.1hz, co.ch4.823.1hz, 
                                        warneke.823.1hz, isaf.823.1hz,rollinsno.823.1hz, rollinsso2.823.1hz, 
                                        cit.823.1hz,gtcims.823.1hz,ryerson.823.1hz , jimenez.823.1hz, 
                                        schwarz.823.1hz, freid.c2h6.823.1hz, freid.ch2o.823.1hz, 
                                        womack.823.1hz, stclair.823.1hz,veres.823.1hz, #wisthaler.823.1hz,
                                        moore.823, met.823.1hz)
    ind2 = which(aug23rd.fire$Time_Start >= start & aug23rd.fire$Time_Start  <= stop	) 
    ind2B = which(aug23rd.fire$Time_Start >= startB & aug23rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug23rd.fire$fire = ''; aug23rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Frisbee'
    aug23rd.fire$fuel[ind2] = fuel ; aug23rd.fire$fire[ind2] = 'Frisbee'
    # -------- Blake cans 
    indBLAKE = which(blake.823.merge$Time_Start >= (start-5) & blake.823.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.823.merge[indBLAKE,]
      blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.823.merge$Time_Start >= (start-0) & gilman.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.823.merge[indGILMAN,]
      GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.823.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.823.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.823.merge[1,]
    }
    
    # Apel  
    indAPEL = which( apel.823.merge$Time_Start >= (start-0) &  apel.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.823.merge[indAPEL,]
      APELBG = apel.823.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.823.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.823.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.823.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.823$Time_Start >= (start-0) & newTOGA.823$Time_Start <= stop & newTOGA.823$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.823[indBECKY,]
      BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.823[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.823[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug23rd.fire[ind2,])
      plotpass5hz(aug23rd.fire[ind2,])
      plotpass1hz(aug23rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug23rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug23rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug23rd.fire[ind2,],aug23rd.fire[ind2B,],xspecies,'Frisbee', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Frisbee', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop  
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Frisbee.1hz.EF = tmpEF1 ;Frisbee.5hz.EF = tmpEF5 }
    if (i > 1){Frisbee.1hz.EF = rbind(Frisbee.1hz.EF, tmpEF1) ; Frisbee.5hz.EF = rbind(Frisbee.5hz.EF, tmpEF5) }
  }
  Frisbee.1hz.EF = Frisbee.1hz.EF[order(Frisbee.1hz.EF$variable),]
  Frisbee.1hz.EF$transect_source_fire_ID = indA[1]
  Frisbee.5hz.EF = Frisbee.5hz.EF[order(Frisbee.5hz.EF$variable),]
  Frisbee.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== GUAC =====================
  fire="Guac";fuel="shrub"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){    
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 8.6;stop = stopO -115; startB = startO + 0; stopB = startO + 2} # bad pass, CO
    #if (i == 2){start = startO + 91;stop = stopO -68; startB = startO + 0; stopB = startO + 2}
    if (i == 2){start = startO + 94;stop = stopO-0.2; startB = startO + 0; stopB = startO + 2}
    if (i == 3){start = startO+0.4 ;stop = stopO-70.4; startB = startO + 0; stopB = startO + 2}
    tmp = time_align(start,stop,co2.823.5hz,  co.ch4.823.5hz, 
                     warneke.823.5hz,  isaf.823.5hz,
                     rollinsno.823.5hz, rollinsso2.823.5hz, 
                     cit.823.5hz, gtcims.823.5hz,moore.823fast, jimenez.823.5hz,met.823.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1}
    
    aug23rd.fire = time_alignSLOWNOCANS(start,stop, co2.823.1hz, co.ch4.823.1hz, 
                                        warneke.823.1hz, isaf.823.1hz,rollinsno.823.1hz, rollinsso2.823.1hz, 
                                        cit.823.1hz,gtcims.823.1hz,ryerson.823.1hz , jimenez.823.1hz, 
                                        schwarz.823.1hz, freid.c2h6.823.1hz, freid.ch2o.823.1hz, 
                                        womack.823.1hz, stclair.823.1hz,veres.823.1hz, #wisthaler.823.1hz,
                                        moore.823, met.823.1hz)
    ind2 = which(aug23rd.fire$Time_Start >= start & aug23rd.fire$Time_Start  <= stop	) 
    ind2B = which(aug23rd.fire$Time_Start >= startB & aug23rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug23rd.fire$fire = ''; aug23rd.fire$fuel = ''
    tmp$fuel[ind] = 'shrub' ; tmp$fire[ind] = 'Guac'
    aug23rd.fire$fuel[ind2] = 'shrub' ; aug23rd.fire$fire[ind2] = 'Guac'
    # -------- Blake cans 
    indBLAKE = which(blake.823.merge$Time_Start >= (start-5) & blake.823.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.823.merge[indBLAKE,]
      blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.823.merge$Time_Start >= (start-0) & gilman.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.823.merge[indGILMAN,]
      GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.823.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.823.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.823.merge[1,]
    }
    
    # Apel  
    indAPEL = which( apel.823.merge$Time_Start >= (start-0) &  apel.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.823.merge[indAPEL,]
      APELBG = apel.823.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.823.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.823.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.823.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.823$Time_Start >= (start-0) & newTOGA.823$Time_Start <= stop & newTOGA.823$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.823[indBECKY,]
      BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.823[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.823[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug23rd.fire[ind2,])
      plotpass5hz(aug23rd.fire[ind2,])
      plotpass1hz(aug23rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug23rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug23rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug23rd.fire[ind2,],aug23rd.fire[ind2B,],xspecies,'Guac', 'shrub',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Guac', 'shrub',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Guac.1hz.EF = tmpEF1 ;Guac.5hz.EF = tmpEF5 }
    if (i > 1){Guac.1hz.EF = rbind(Guac.1hz.EF, tmpEF1) ; Guac.5hz.EF = rbind(Guac.5hz.EF, tmpEF5) }
  }
  Guac.1hz.EF = Guac.1hz.EF[order(Guac.1hz.EF$variable),]
  Guac.1hz.EF$transect_source_fire_ID = indA[1]
  Guac.5hz.EF = Guac.5hz.EF[order(Guac.5hz.EF$variable),]
  Guac.5hz.EF$transect_source_fire_ID = indA[1]
  # ========================== HAMBURGER =====================
  fire="Hamburger"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) ) 
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 16;stop = stopO -6.2; startB = startO + 0; stopB = startO +2} # bad pass
    if (i == 2){start = startO + 30.8;  stop = stopO -0; startB = startB = startO + 0; stopB = startO +2}
    if (i == 3){start = startO + 1;  stop = stopO -1; startB = startB = 72613; stopB = 72615} # split from above
    if (i == 4){start = startO + 0;  stop = stopO -0; startB = startB = 72613; stopB =72615} # split from above
    if (i == 5){start = startO + 1;  stop = stopO -99.2; startB = startB = 72613; stopB =72615} # split from above
    if (i == 6){start = startO + 32.6;  stop = stopO -17; tartB = startO + 0; stopB = startO +2} # low CO
    if (i == 7){start = startO + 3;stop = stopO -1; startB = startO + 0; stopB = startO +2} #
    if (i == 8){start = startO + 1;stop = stopO -7; startB = 72870; stopB = 72872} # split from above
    if (i == 9){start = startO + 27;stop = stopO -1; startB = startO + 0; stopB = startO +2} #
    if (i == 10){start = startO + 2;stop = stopO -0; startB = 73304; stopB = 73306} # split from above
    if (i == 11){start = startO + 0;stop = stopO -47; startB = 73304; stopB = 73306} # split from above
    if (i == 12){start = startO + 3;stop = stopO -29; startB = startO + 0; stopB = startO +2} #
    tmp = time_align(start,stop,co2.823.5hz,  co.ch4.823.5hz, 
                     warneke.823.5hz,  isaf.823.5hz,
                     rollinsno.823.5hz, rollinsso2.823.5hz, 
                     cit.823.5hz, gtcims.823.5hz,moore.823fast, jimenez.823.5hz,met.823.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i == 1){start = start -1}
    if (i == 2){start=start-1;startB=startB-1; stopB=stopB-1}
    if (i == 3){start=start-1;startB=startB-1; stopB=stopB-1}
    if (i == 4){start = start -1; startB=startB-1; stopB=stopB-1}
    if (i == 5){start = start -1; startB=startB-1; stopB=stopB-1}
    if (i == 6){start = start -1; startB=startB-1; stopB=stopB-1}
    #if (i == 7){start = start -1; startB=startB-1; stopB=stopB-1}
    #if (i == 8){start = start -1; startB=startB-1; stopB=stopB-1}
    #if (i == 9){start = start -1; startB=startB-1; stopB=stopB-1}
    #if (i == 10){start = start -1; startB=startB-1; stopB=stopB-1}
    #if (i == 11){start = start -1; stop=stop+2;startB=startB-1; stopB=stopB-1}
    #if (i == 12){start = start -1; startB=startB-1; stopB=stopB-1}
    
    aug23rd.fire = time_alignSLOWNOCANS(start,stop, co2.823.1hz, co.ch4.823.1hz, 
                                        warneke.823.1hz, isaf.823.1hz,rollinsno.823.1hz, rollinsso2.823.1hz, 
                                        cit.823.1hz,gtcims.823.1hz,ryerson.823.1hz , jimenez.823.1hz, 
                                        schwarz.823.1hz, freid.c2h6.823.1hz, freid.ch2o.823.1hz, 
                                        womack.823.1hz, stclair.823.1hz,veres.823.1hz, #wisthaler.823.1hz,
                                        moore.823, met.823.1hz)
    ind2 = which(aug23rd.fire$Time_Start >= start & aug23rd.fire$Time_Start  <= stop	) 
    ind2B = which(aug23rd.fire$Time_Start >= startB & aug23rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug23rd.fire$fire = ''; aug23rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Hamburger'
    aug23rd.fire$fuel[ind2] = fuel ; aug23rd.fire$fire[ind2] = 'Hamburger'
     # -------- Blake cans 
    indBLAKE = which(blake.823.merge$Time_Start >= (start-5) & blake.823.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.823.merge[indBLAKE,]
      blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.823.merge$Time_Start >= (start-0) & gilman.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.823.merge[indGILMAN,]
      GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.823.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.823.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.823.merge[1,]
    }
    
    # Apel  
    indAPEL = which( apel.823.merge$Time_Start >= (start-0) &  apel.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.823.merge[indAPEL,]
      APELBG = apel.823.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.823.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.823.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.823.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.823$Time_Start >= (start-0) & newTOGA.823$Time_Start <= stop & newTOGA.823$Percent_plume > 0.5) # check at least 5 sec before?
    # kludge for i == 7 #72854
    if (i == 7){indBECKY = which(newTOGA.823$Time_Start >= (start-17) & newTOGA.823$Time_Start <= stop)}
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.823[indBECKY,]
      BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.823[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.823[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug23rd.fire[ind2,])
      plotpass5hz(aug23rd.fire[ind2,])
      plotpass1hz(aug23rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug23rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug23rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug23rd.fire[ind2,],aug23rd.fire[ind2B,],xspecies,'Hamburger', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Hamburger', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Hamburger.1hz.EF = tmpEF1 ;Hamburger.5hz.EF = tmpEF5 }
    if (i > 1){Hamburger.1hz.EF = rbind(Hamburger.1hz.EF, tmpEF1) ; Hamburger.5hz.EF = rbind(Hamburger.5hz.EF, tmpEF5) }
    
    #iiq = which(tmpEF5$variable == 'CO_DACOM_DISKIN')

    #if (i >1 & tmpEF5$maxval[iiq] > 400 & tmpEF5$R2toX[iiq] > 0.75){
    #  print(c('Here2',i))
    #  tmp$pass = i
    #  if (i == 2){HamburgerData.5hz =  tmp[ind,]}
    #  if (i > 2){HamburgerData.5hz = rbind.fill(HamburgerData.5hz, tmp[ind,])}
    #  
    #  if (doBECKY == 1 & i == 2){
    #    benz = BECKY$Benzene_ppt; benz_BG = BECKYBG$Benzene_ppt
    #    co = BECKY$CO_X; co_BG = BECKYBG$CO_X
    #    timeB = BECKY$Time_Start
    #  }
    #  if (doBECKY == 1 & i > 2){
    #    benz = c(benz,BECKY$Benzene_ppt); benz_BG = c(benz_BG,BECKYBG$Benzene_ppt)
    #    co = c(co,BECKY$CO_X); co_BG=c(co_BG, BECKYBG$CO_X)
    #    timeB = c(timeB,BECKY$Time_Start)
    #  }
    #}
  }
  Hamburger.1hz.EF = Hamburger.1hz.EF[order(Hamburger.1hz.EF$variable),]
  Hamburger.1hz.EF$transect_source_fire_ID = indA[1]
  Hamburger.5hz.EF = Hamburger.5hz.EF[order(Hamburger.5hz.EF$variable),]
  Hamburger.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== IPA =====================
  fire="IPA"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 1.6;stop = stopO -9.2; startB = startO -2; stopB = startO + 0}
    if (i == 2){start = startO + 5.2;stop = stopO -30; startB = startO -1; stopB = startO + 1}
    tmp = time_align(start,stop,co2.823.5hz,  co.ch4.823.5hz, 
                     warneke.823.5hz,  isaf.823.5hz,
                     rollinsno.823.5hz, rollinsso2.823.5hz, 
                     cit.823.5hz, gtcims.823.5hz,moore.823fast, jimenez.823.5hz,met.823.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i == 1){start = start-2;startB=startB-1;stopB=stopB-1}
    if (i == 2){start = start-2;startB=startB-1;stopB=stopB-1}
    aug23rd.fire = time_alignSLOWNOCANS(start,stop, co2.823.1hz, co.ch4.823.1hz, 
                                        warneke.823.1hz, isaf.823.1hz,rollinsno.823.1hz, rollinsso2.823.1hz, 
                                        cit.823.1hz,gtcims.823.1hz,ryerson.823.1hz , jimenez.823.1hz, 
                                        schwarz.823.1hz, freid.c2h6.823.1hz, freid.ch2o.823.1hz, 
                                        womack.823.1hz, stclair.823.1hz,veres.823.1hz, #wisthaler.823.1hz,
                                        moore.823, met.823.1hz)
    ind2 = which(aug23rd.fire$Time_Start >= start & aug23rd.fire$Time_Start  <= stop	) 
    ind2B = which(aug23rd.fire$Time_Start >= startB & aug23rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug23rd.fire$fire = ''; aug23rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'IPA'
    aug23rd.fire$fuel[ind2] = fuel ; aug23rd.fire$fire[ind2] = 'IPA'
    # -------- Blake cans 
    indBLAKE = which(blake.823.merge$Time_Start >= (start-5) & blake.823.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.823.merge[indBLAKE,]
      blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.823.merge$Time_Start >= (start-0) & gilman.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.823.merge[indGILMAN,]
      GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.823.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.823.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.823.merge[1,]
    }
    
    # Apel  
    indAPEL = which( apel.823.merge$Time_Start >= (start-0) &  apel.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.823.merge[indAPEL,]
      APELBG = apel.823.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.823.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.823.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.823.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.823$Time_Start >= (start-0) & newTOGA.823$Time_Start <= stop & newTOGA.823$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.823[indBECKY,]
      BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.823[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.823[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug23rd.fire[ind2,])
      plotpass5hz(aug23rd.fire[ind2,])
      plotpass1hz(aug23rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug23rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug23rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug23rd.fire[ind2,],aug23rd.fire[ind2B,],xspecies,'IPA', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'IPA', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){IPA.1hz.EF = tmpEF1 ;IPA.5hz.EF = tmpEF5 }
    if (i > 1){IPA.1hz.EF = rbind(IPA.1hz.EF, tmpEF1) ; IPA.5hz.EF = rbind(IPA.5hz.EF, tmpEF5) }
    
  }
  IPA.1hz.EF = IPA.1hz.EF[order(IPA.1hz.EF$variable),]
  IPA.1hz.EF$transect_source_fire_ID = indA[1]
  IPA.5hz.EF = IPA.5hz.EF[order(IPA.5hz.EF$variable),]
  IPA.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== JELLO =====================
  fire="Jello"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 2;stop = stop -4; startB = startO + 0; stopB = startO + 2}
    if (i == 2){start = startO + 1.6;stop = stop -21.6; startB = startO -2; stopB = startO -1}
    tmp = time_align(start,stop,co2.823.5hz,  co.ch4.823.5hz, 
                     warneke.823.5hz,  isaf.823.5hz,
                     rollinsno.823.5hz, rollinsso2.823.5hz, 
                     cit.823.5hz, gtcims.823.5hz,moore.823fast, jimenez.823.5hz,met.823.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i==1){start=start-1;stop=stop+1;startB=startB-1;stopB=stopB-1}
    if (i==2){start=start-1;startB=startB-1;stopB=stopB-1}
    aug23rd.fire = time_alignSLOWNOCANS(start,stop, co2.823.1hz, co.ch4.823.1hz, 
                                        warneke.823.1hz, isaf.823.1hz,rollinsno.823.1hz, rollinsso2.823.1hz, 
                                        cit.823.1hz,gtcims.823.1hz,ryerson.823.1hz , jimenez.823.1hz, 
                                        schwarz.823.1hz, freid.c2h6.823.1hz, freid.ch2o.823.1hz, 
                                        womack.823.1hz, stclair.823.1hz,veres.823.1hz, #wisthaler.823.1hz,
                                        moore.823, met.823.1hz)
    ind2 = which(aug23rd.fire$Time_Start >= start & aug23rd.fire$Time_Start  <= stop	) 
    ind2B = which(aug23rd.fire$Time_Start >= startB & aug23rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug23rd.fire$fire = ''; aug23rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Jello'
    aug23rd.fire$fuel[ind2] = fuel ; aug23rd.fire$fire[ind2] = 'Jello'
    # -------- Blake cans 
    indBLAKE = which(blake.823.merge$Time_Start >= (start-5) & blake.823.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.823.merge[indBLAKE,]
      blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.823.merge$Time_Start >= (start-0) & gilman.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.823.merge[indGILMAN,]
      GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.823.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.823.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.823.merge[1,]
    }
    
    # Apel  
    indAPEL = which( apel.823.merge$Time_Start >= (start-0) &  apel.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.823.merge[indAPEL,]
      APELBG = apel.823.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.823.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.823.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.823.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.823$Time_Start >= (start-0) & newTOGA.823$Time_Start <= stop & newTOGA.823$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.823[indBECKY,]
      BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.823[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.823[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug23rd.fire[ind2,])
      plotpass5hz(aug23rd.fire[ind2,])
      plotpass1hz(aug23rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug23rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug23rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug23rd.fire[ind2,],aug23rd.fire[ind2B,],xspecies,'Jello', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Jello', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop  
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Jello.1hz.EF = tmpEF1 ;Jello.5hz.EF = tmpEF5 }
    if (i > 1){Jello.1hz.EF = rbind(Jello.1hz.EF, tmpEF1) ; Jello.5hz.EF = rbind(Jello.5hz.EF, tmpEF5) }
  }
  Jello.1hz.EF = Jello.1hz.EF[order(Jello.1hz.EF$variable),]
  Jello.1hz.EF$transect_source_fire_ID = indA[1]
  Jello.5hz.EF = Jello.5hz.EF[order(Jello.5hz.EF$variable),]
  Jello.5hz.EF$transect_source_fire_ID = indA[1]
  # ========================== KEBAB =====================
  fire='Kebab'; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == 'Kebab')
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 31;stop = stop -3; startB = startO + 0; stopB = startO + 2}
    if (i == 2){start = startO + 24;stop = stop -55; startB = 78490; stopB = 78492} # split from above
    tmp = time_align(start,stop,co2.823.5hz,  co.ch4.823.5hz, 
                     warneke.823.5hz,  isaf.823.5hz,
                     rollinsno.823.5hz, rollinsso2.823.5hz, 
                     cit.823.5hz, gtcims.823.5hz,moore.823fast, jimenez.823.5hz,met.823.5hz)
   
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i==1){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i==2){start=start-1;startB=startB-1;stopB=stopB-1}
    
    aug23rd.fire = time_alignSLOWNOCANS(start,stop, co2.823.1hz, co.ch4.823.1hz, 
                                        warneke.823.1hz, isaf.823.1hz,rollinsno.823.1hz, rollinsso2.823.1hz, 
                                        cit.823.1hz,gtcims.823.1hz,ryerson.823.1hz , jimenez.823.1hz, 
                                        schwarz.823.1hz, freid.c2h6.823.1hz, freid.ch2o.823.1hz, 
                                        womack.823.1hz, stclair.823.1hz,veres.823.1hz, #wisthaler.823.1hz,
                                        moore.823, met.823.1hz)
    ind2 = which(aug23rd.fire$Time_Start >= start & aug23rd.fire$Time_Start  <= stop	) 
    ind2B = which(aug23rd.fire$Time_Start >= startB & aug23rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug23rd.fire$fire = ''; aug23rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Kebab'
    aug23rd.fire$fuel[ind2] = fuel ; aug23rd.fire$fire[ind2] = 'Kebab'
    # -------- Blake cans 
    indBLAKE = which(blake.823.merge$Time_Start >= (start-5) & blake.823.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.823.merge[indBLAKE,]
      blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.823.merge$Time_Start >= (start-0) & gilman.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.823.merge[indGILMAN,]
      GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.823.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.823.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.823.merge[1,]
    }
    
    # Apel  
    indAPEL = which( apel.823.merge$Time_Start >= (start-0) &  apel.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.823.merge[indAPEL,]
      APELBG = apel.823.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.823.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.823.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.823.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.823$Time_Start >= (start-0) & newTOGA.823$Time_Start <= stop & newTOGA.823$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.823[indBECKY,]
      BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.823[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.823[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug23rd.fire[ind2,])
      plotpass5hz(aug23rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass5hzJUSTCH4(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug23rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug23rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug23rd.fire[ind2,],aug23rd.fire[ind2B,],xspecies,'Kebab', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Kebab', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Kebab.1hz.EF = tmpEF1 ;Kebab.5hz.EF = tmpEF5 }
    if (i > 1){Kebab.1hz.EF = rbind(Kebab.1hz.EF, tmpEF1) ; Kebab.5hz.EF = rbind(Kebab.5hz.EF, tmpEF5) }
  }
  Kebab.1hz.EF = Kebab.1hz.EF[order(Kebab.1hz.EF$variable),]
  Kebab.1hz.EF$transect_source_fire_ID = indA[1]
  Kebab.5hz.EF = Kebab.5hz.EF[order(Kebab.5hz.EF$variable),]
  Kebab.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== LIMONCELLO =====================
  fire="Limoncello"; fuel="rice"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 4.6;stop = stop -.6; startB = startO + 0; stopB = startO + 2}
    if (i == 2){start = startO + 1;stop = stop -25; startB = 79301; stopB = 79303} # split from above
    if (i == 3){start = startO + 16.4;stop = stop -10.6; startB = startO -1; stopB = startO +1}
    if (i == 4){start = startO + 18.6;stop = stop -8.4; startB = startO -2; stopB =  startO +0}
    if (i == 5){start = startO + 15;stop = stop -13.2; startB = startO + 0; stopB = startO +1} # bad pass?
    tmp = time_align(start,stop,co2.823.5hz,  co.ch4.823.5hz, 
                     warneke.823.5hz,  isaf.823.5hz,
                     rollinsno.823.5hz, rollinsso2.823.5hz, 
                     cit.823.5hz, gtcims.823.5hz,moore.823fast, jimenez.823.5hz,met.823.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i == 1){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i == 2){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i == 3){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i == 4){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i == 5){start=start-1;startB=startB-1;stopB=stopB-1}
    aug23rd.fire = time_alignSLOWNOCANS(start,stop, co2.823.1hz, co.ch4.823.1hz, 
                                        warneke.823.1hz, isaf.823.1hz,rollinsno.823.1hz, rollinsso2.823.1hz, 
                                        cit.823.1hz,gtcims.823.1hz,ryerson.823.1hz , jimenez.823.1hz, 
                                        schwarz.823.1hz, freid.c2h6.823.1hz, freid.ch2o.823.1hz, 
                                        womack.823.1hz, stclair.823.1hz,veres.823.1hz, #wisthaler.823.1hz,
                                        moore.823, met.823.1hz)
    ind2 = which(aug23rd.fire$Time_Start >= start & aug23rd.fire$Time_Start  <= stop	) 
    ind2B = which(aug23rd.fire$Time_Start >= startB & aug23rd.fire$Time_Start <= stopB	) 

    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug23rd.fire$fire = ''; aug23rd.fire$fuel = ''
    tmp$fuel[ind] = 'rice' ; tmp$fire[ind] = 'Limoncello'
    aug23rd.fire$fuel[ind2] = 'rice' ; aug23rd.fire$fire[ind2] = 'Limoncello'
    # -------- Blake cans 
    indBLAKE = which(blake.823.merge$Time_Start >= (start-5) & blake.823.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.823.merge[indBLAKE,]
      blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.823.merge$Time_Start >= (start-0) & gilman.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.823.merge[indGILMAN,]
      GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.823.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.823.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.823.merge[1,]
    }
    
    # Apel  
    indAPEL = which( apel.823.merge$Time_Start >= (start-0) &  apel.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.823.merge[indAPEL,]
      APELBG = apel.823.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.823.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.823.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.823.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.823$Time_Start >= (start-0) & newTOGA.823$Time_Start <= stop & newTOGA.823$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.823[indBECKY,]
      BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.823[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.823[1,]
    }

    if (doplot == 1){
      plotpass5hzJUSTCO(aug23rd.fire[ind2,])
      plotpass5hz(aug23rd.fire[ind2,])
      plotpass1hz(aug23rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug23rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug23rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug23rd.fire[ind2,],aug23rd.fire[ind2B,],xspecies,'Limoncello', 'rice',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Limoncello', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop  
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Limoncello.1hz.EF = tmpEF1 ;Limoncello.5hz.EF = tmpEF5 }
    if (i > 1){Limoncello.1hz.EF = rbind(Limoncello.1hz.EF, tmpEF1) ; Limoncello.5hz.EF = rbind(Limoncello.5hz.EF, tmpEF5) }
  }
  Limoncello.1hz.EF = Limoncello.1hz.EF[order(Limoncello.1hz.EF$variable),]
  Limoncello.1hz.EF$transect_source_fire_ID = indA[1]
  Limoncello.5hz.EF = Limoncello.5hz.EF[order(Limoncello.5hz.EF$variable),]
  Limoncello.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== MUSTARD both very low CO =====================
  fire="Mustard"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 55.8;stop = stopO -7.4; startB = startO + 0; stopB = startO + 2}
    if (i == 2){start = startO + 14;stop = stopO -26; startB = startO + 0; stopB = startO + 2} # bad pass
    tmp = time_align(start,stop,co2.823.5hz,  co.ch4.823.5hz, 
                     warneke.823.5hz,  isaf.823.5hz,
                     rollinsno.823.5hz, rollinsso2.823.5hz, 
                     cit.823.5hz, gtcims.823.5hz,moore.823fast, jimenez.823.5hz,met.823.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i == 1){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i == 2){start=start-1;startB=startB-1;stopB=stopB-1}
    aug23rd.fire = time_alignSLOWNOCANS(start,stop, co2.823.1hz, co.ch4.823.1hz, 
                                        warneke.823.1hz, isaf.823.1hz,rollinsno.823.1hz, rollinsso2.823.1hz, 
                                        cit.823.1hz,gtcims.823.1hz,ryerson.823.1hz , jimenez.823.1hz, 
                                        schwarz.823.1hz, freid.c2h6.823.1hz, freid.ch2o.823.1hz, 
                                        womack.823.1hz, stclair.823.1hz,veres.823.1hz, #wisthaler.823.1hz,
                                        moore.823, met.823.1hz)
    ind2 = which(aug23rd.fire$Time_Start >= start & aug23rd.fire$Time_Start  <= stop	) 
    ind2B = which(aug23rd.fire$Time_Start >= startB & aug23rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug23rd.fire$fire = ''; aug23rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Mustard'
    aug23rd.fire$fuel[ind2] = fuel ; aug23rd.fire$fire[ind2] = 'Mustard'
    # -------- Blake cans 
    indBLAKE = which(blake.823.merge$Time_Start >= (start-5) & blake.823.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.823.merge[indBLAKE,]
      blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.823.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.823.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.823.merge$Time_Start >= (start-0) & gilman.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.823.merge[indGILMAN,]
      GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.823.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.823.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.823.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.823.merge[1,]
    }
    
    # Apel  
    indAPEL = which( apel.823.merge$Time_Start >= (start-0) &  apel.823.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.823.merge[indAPEL,]
      APELBG = apel.823.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.823.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.823.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.823.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.823$Time_Start >= (start-0) & newTOGA.823$Time_Start <= stop & newTOGA.823$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.823[indBECKY,]
      BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.823$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.823[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.823[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.823[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug23rd.fire[ind2,])
      plotpass5hz(aug23rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug23rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug23rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug23rd.fire[ind2,],aug23rd.fire[ind2B,],xspecies,'Mustard', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Mustard', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug23rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug23rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug23rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug23rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug23rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop   
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Mustard.1hz.EF = tmpEF1 ;Mustard.5hz.EF = tmpEF5 }
    if (i > 1){Mustard.1hz.EF = rbind(Mustard.1hz.EF, tmpEF1) ; Mustard.5hz.EF = rbind(Mustard.5hz.EF, tmpEF5) }
  }
  Mustard.1hz.EF = Mustard.1hz.EF[order(Mustard.1hz.EF$variable),]
  Mustard.1hz.EF$transect_source_fire_ID = indA[1]
  Mustard.5hz.EF = Mustard.5hz.EF[order(Mustard.5hz.EF$variable),]
  Mustard.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ----- [[[[[[[[[[[[[[[[[[[[ Aug 26th ]]]]]]]]]]]]]]]]]]]]] -------
  # --------#######--------- Get 1 Hz Data individual 8/26------#######---------
  # plume tags
  tags = getICARTTdataSIMPLE('InputFiles/firexaq-fire-Flags-1HZ_DC8_20190826_R9.ict') ;tags$Time_Start = tags$TIME_START
  # MET DATA
  met.826.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-MetNav_DC8_20190826_R1.ict')
  met.826.1hz = merge(met.826.1hz, tags, by='Time_Start')
  # CO2
  co2.826.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CO2-7000_DC8_20190826_R2.ict')
  # -------- DISKIN -----CO, CH4
  co.ch4.826.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-DACOM_DC8_20190826_R1.ict')
  # --------- WARNEKE ----  VOCs
  warneke.826.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-NOAAPTR-VOCs-1Hz_DC8_20190826_R3.ict')
  # ------ HANISCO - ISAF HCHO - merged to 5Hz from the online merge
  isaf.826.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-ISAF-CH2O-1Hz_DC8_20190826_R0.ict')
  #  ------- ROLLINS - SO2 and NO
  rollinsno.826.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-NO_DC8_20190826_R1.ict')
  rollinsno.826.1hz$Time_Start = rollinsno.826.1hz$time_mid
  rollinsso2.826.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-SO2_DC8_20190826_R1.ict')
  rollinsso2.826.1hz$Time_Start = rollinsso2.826.1hz$time_mid
  
  #  ----- WENNBERG - CIT VOCs - 
  cit.826.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190826_R0_CIT.ict')
  # ------ HUEY - GTCIMS PANs - not sure how to match up peaks here
  gtcims.826.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190826_R0_Huey.ict')
  
  # ------ RYERSON
  ryerson.A = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NO_DC8_20190826_R1.ict')
  ryerson.B = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NO2_DC8_20190826_R1.ict')
  ryerson.C = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NOy_DC8_20190826_R1.ict')
  ryerson.D = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-O3_DC8_20190826_R1.ict')
  ryerson.826.1hz = cbind(ryerson.A,ryerson.B,ryerson.C,ryerson.D) 
  ryerson.826.1hz$Time_Start = ryerson.826.1hz$Time_start
  # ----- JIMENEZ ---
  jimenez.826.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-EESI_DC8_20190826_R1.ict')
  
  # ----- SCHWARZ ---
  schwarz.826.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-SP2-BC-1HZ_DC8_20190826_R2.ict')
  # ----- FREID ---
  freid.c2h6.826.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-C2H6_DC8_20190826_R3.ict')
  freid.ch2o.826.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CH2O_DC8_20190826_R3.ict')
  # ------ WOMACK ---
  womackA = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-CH3COCHO_DC8_20190826_R1.ict')
  womackB = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-CHOCHO_DC8_20190826_R1.ict')
  womackC = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-HNO2_DC8_20190826_R1.ict')
  womackD = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-NO2_DC8_20190826_R1.ict')
  womack.826.1hz = cbind(womackA, womackB, womackC, womackD)
  # -------St Clair
  stclair.826.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-CANOE-NO2_DC8_20190826_R0.ict')
  stclair.826.1hz$Time_Start = stclair.826.1hz$Time_start
  # data is all missing, kludge for now
  stclair.826.1hz = womack.826.1hz
  stclair.826.1hz$NO2_CANOE = womack.826.1hz$CH3COCHO_ACES*NaN
  stclair.826.1hz=stclair.826.1hz[,c(1,9)]
  # ------- VERES
  veres.A = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-ClNO2_DC8_20190826_R0.ict')
  veres.B = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HCOOH_DC8_20190826_R1.ict')
  veres.C = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HNO2_DC8_20190826_R1.ict')
  veres.D = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-N2O5_DC8_20190826_R0.ict')
  veres.E = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HPMTF_DC8_20190826_R0.ict')
  veres.F = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-CH3COOCl_DC8_20190826_R0.ict')
  veres.G = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-Cl2_DC8_20190826_R0.ict')
  veres.H = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrCl_DC8_20190826_R0.ict')
  veres.I = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrCN_DC8_20190826_R0.ict')
  veres.J = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrO_DC8_20190826_R0.ict')
  veres.K = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HCN_DC8_20190826_R0.ict')
  veres.L = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HNCO_DC8_20190826_R0.ict')
  veres.826.1hz = cbind(veres.A,veres.B,veres.C,veres.D,veres.E,veres.F,veres.G,veres.H,veres.I,veres.J,veres.K,veres.L)
  
  # --- WISTHALER
  #wisthaler.826.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190826_R0_Wisthaler.ict')
  
  # ---- BLAKE
  blake.826.1hz = getICARTTdataSIMPLE('InputFiles/WAS-MERGE/firexaq-mrgWAS-dc8_merge_20190826_R1.ict')
  cc = colnames(blake.826.1hz)
  blake.826.merge = blake.826.1hz[,c(1,2,96:225)]
  blake.826.merge$CO_DACOM_DISKIN_BLAKE = blake.826.1hz$CO_DACOM_DISKIN
  blake.826.merge$CO2_7000_ppm_DISKIN_BLAKE = blake.826.1hz$CO2_7000_ppm_DISKIN
  
  # ------ APEL
  apel.826.1hz = getICARTTdataSIMPLE('InputFiles/TOGA-MERGE/firexaq-mrgTOGA-dc8_merge_20190826_R1.ict')
  cc = colnames(apel.826.1hz)
  apel.826.merge = apel.826.1hz[,c(1,2,226:315)]
  apel.826.merge$CO_DACOM_DISKIN_APEL = apel.826.1hz$CO_DACOM_DISKIN
  apel.826.merge$CO2_7000_ppm_DISKIN_APEL =apel.826.1hz$CO2_7000_ppm_DISKIN
  # Becky's better merge
  file = 'InputFiles/Hornbrook/FIREX-AQ weighted TOGA merge 2022-01-24_0826.xlsx'
  newTOGA.826 = readxl::read_xlsx(file); newTOGA.826[newTOGA.826==-999] = NaN; newTOGA.826[newTOGA.826==-888] = NaN
  newTOGA.826$CO_DACOM_DISKIN_BECKY = newTOGA.826$CO_DACOM_DISKIN
  newTOGA.826$CO2_7000_ppm_DISKIN_BECKY = NaN
  newTOGA.826$Time_Start=newTOGA.826$Time_Start...4
  
  # ----GILMAN
  gilman.826.1hz = getICARTTdataSIMPLE('InputFiles/iWAS-MERGE/firexaq-mrgiWAS-dc8_merge_20190826_R1.ict')
  cc = colnames(gilman.826.1hz)
  gilman.826.merge = gilman.826.1hz[,c(1,2,316:361)]
  gilman.826.merge$CO_DACOM_DISKIN_GILMAN = gilman.826.1hz$CO_DACOM_DISKIN
  gilman.826.merge$CO2_7000_ppm_DISKIN_GILMAN = gilman.826.1hz$CO2_7000_ppm_DISKIN
  #gilman.826.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAA-iWAS-VOCs_DC8_20190826_R0.ict')
  
  # ------ Moore 
  moore.826fast = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190826_R0_MOORE.ict')
  
  moore.826p1 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-AerosolCloudConc_DC8_20190826_R0.ict')
  moore.826p2 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-LAScold_DC8_20190826_R0.ict')
  moore.826p3 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-LAShot_DC8_20190826_R0.ict')
  moore.826p4 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-CPSPD_DC8_20190826_R0.ict')
  moore.826p5 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-CDP_DC8_20190826_R0.ict')
  moore.826p6 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-SMPS_DC8_20190826_R0.ict')
  moore.826 =merge(moore.826p1, moore.826p2, by='Time_mid', all = TRUE, incomparables = NA)
  moore.826 =merge(moore.826, moore.826p3, by='Time_mid', all = TRUE, incomparables = NA)
  moore.826 =merge(moore.826, moore.826p4, by='Time_mid', all = TRUE, incomparables = NA)
  moore.826 =merge(moore.826, moore.826p5, by='Time_mid', all = TRUE, incomparables = NA)
  moore.826 =merge(moore.826, moore.826p6, by='Time_mid', all = TRUE, incomparables = NA)
  # ------- append PI to colnames 1hz ----------
  cc = colnames(co2.826.1hz)
  cc[2] = paste(cc[2],'_DISKIN',sep='')
  colnames(co2.826.1hz) = cc 
  cc = colnames(co.ch4.826.1hz)
  cc[2:4] =  paste(cc[2:4],'_DISKIN',sep='')
  colnames(co.ch4.826.1hz) = cc
  cc=colnames(met.826.1hz)
  
  colnames(isaf.826.1hz) = c("Time_Start"  , "CH2O_ISAF_HANISCO" ,"CH2O_ISAF_precision_HANISCO")
  cc[2:36] = paste(cc[2:36],'_YANG',sep='')
  colnames(met.826.1hz) = cc
  cc = colnames(warneke.826.1hz)
  cc[2:43] = paste(cc[2:43],'_WARNEKE', sep='')
  colnames(warneke.826.1hz) = cc
  cc = colnames(rollinsno.826.1hz)
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsno.826.1hz) = cc
  cc = colnames(rollinsso2.826.1hz)
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsso2.826.1hz) = cc
  # make GTCIMS consistent with 1s merge
  cc=colnames(gtcims.826.1hz)
  cc[4:7] = c("PAN_GTCIMS_HUEY" , "PPN_GTCIMS_HUEY"  ,"APAN_GTCIMS_HUEY" ,"PBN_GTCIMS_HUEY")
  colnames(gtcims.826.1hz)=cc
  # since ISAF, CIT, and GTCIMS came from merge tool, alread has PI's appended.
  ryerson.826.1hz <- ryerson.826.1hz[, !duplicated(colnames(ryerson.826.1hz))]
  
  cc = colnames(ryerson.826.1hz)
  cc[2:9] = paste(cc[2:9],'_RYERSON',sep='')
  colnames(ryerson.826.1hz)=cc
  
  cc = colnames(schwarz.826.1hz)
  cc[2:3] = paste(cc[2:3],'_SCHWARZ', sep='')
  colnames(schwarz.826.1hz) =cc
  
  cc = colnames(freid.c2h6.826.1hz)
  cc[4:7] = paste(cc[4:7], '_FRIED',sep='')
  colnames(freid.c2h6.826.1hz) = cc
  
  cc = colnames(freid.ch2o.826.1hz)
  cc[4:7] = paste(cc[4:7], '_FRIED',sep='')
  colnames(freid.ch2o.826.1hz) = cc
  
  womack.826.1hz <- womack.826.1hz[, !duplicated(colnames(womack.826.1hz))]
  cc = colnames(womack.826.1hz)
  cc[2:5] = paste(cc[2:5], '_WOMACK',sep='')
  colnames(womack.826.1hz) = cc

  veres.826.1hz <- veres.826.1hz[, !duplicated(colnames(veres.826.1hz))]
  cc = colnames(veres.826.1hz)
  cc[2:13] = paste(cc[2:13], '_VERES',sep='')
  colnames(veres.826.1hz) = cc
  
  #cc = colnames(wisthaler.826.1hz)
  #cc[4:5] = paste(cc[4:5], '_WISTHALER',sep='')
  #colnames(wisthaler.826.1hz) = cc
  cc = colnames(jimenez.826.1hz)
  cc[2:8] = paste(cc[2:8], '_JIMENEZ',sep='')
  colnames(jimenez.826.1hz) = cc
  # --------#######--------- Get 5 or 10 Hz Data 8/26 ------#######---------
  # MET DATA - BUI + YANG + DLH
  met.826.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190826_R0_met.ict')
  #CO, CH4
  co.ch4.826.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-DACOM-5Hz_DC8_20190826_R1.ict')
  # CO2
  co2.826.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CO2-7000-5Hz_DC8_20190826_R1.ict')
   # WARNEKE VOCs
  warneke.826.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-NOAAPTR-VOCs-5Hz_DC8_20190826_R3.ict')
  # ISAF HCHO - merged to 5Hz from the online merge
  isaf.826.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190826_R0_ISAF.ict')
  # ROLLINS SO2 and NO
  rollinsno.826.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-NO-5Hz_DC8_20190826_R0.ict')
  rollinsso2.826.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-SO2-5Hz_DC8_20190826_R1.ict')
  # CIT VOCs - 
  cit.826.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190826_R0_CIT.ict')
  # GTCIMS PANs - not sure how to match up peaks here
  gtcims.826.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190826_R0_huey.ict')
  # ----- Jimenez ---
  jimenez.826.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_AMS_20190826_R0_20230314T134049.ict')
  jimenez.826.5hz$OC_PM1_AMS_JIMENEZ = jimenez.826.5hz$OA_PM1_AMS_JIMENEZ/jimenez.826.5hz$OAtoOC_PM1_AMS
  
  # ------- append PI to colnames ----------
  cc = colnames(co2.826.5hz)
  cc[2] = paste(cc[2],'_DISKIN',sep='')
  colnames(co2.826.5hz) = cc 
  cc = colnames(co.ch4.826.5hz)
  cc[2:4] =  paste(cc[2:4],'_DISKIN',sep='')
  colnames(co.ch4.826.5hz) = cc
  cc=colnames(met.826.5hz)
  cc[2:36] = paste(cc[2:36],'_YANG',sep='')
  colnames(met.826.5hz) = cc
  cc = colnames(warneke.826.5hz)
  cc[2:43] = paste(cc[2:43],'_WARNEKE', sep='')
  colnames(warneke.826.5hz) = cc
  cc = colnames(rollinsno.826.5hz)
  cc[1] = "Time_Start" # really Time_Mid but need for merge
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsno.826.5hz) = cc
  cc = colnames(rollinsso2.826.5hz)
  cc[1] = "Time_Start" # really Time_Mid but need for merge
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsso2.826.5hz) = cc
  # make GTCIMS consistent with 1s merge
  cc=colnames(gtcims.826.5hz)
  cc[38:39] = c("APAN_GTCIMS_HUEY" , "PAN_GTCIMS_HUEY")
  colnames(gtcims.826.5hz)=cc
  # since ISAF, CIT, and GTCIMS came from merge tool, alread has PI's appended.
  
  # ------- get fuel moisture data --------
  if (doFM == 1){
    f1 = '/Users/ktravis1/Library/CloudStorage/Box-Box/FuelMoisture/fuel_moisture_content-20210715T1049Z/fmc_20190826_20Z.nc'
    fid = nc_open(f1)
    fuelMDead = ncvar_get(fid, varid = 'FMCG2D')
    fuelMLive = ncvar_get(fid, varid = 'FMCGLH2D')
    xlon = ncvar_get(fid, varid="XLONG_M")
    xlat = ncvar_get(fid, varid="XLAT_M")
    nc_close(fid)
  }  
  # ========================== DaysLater=====================
  fire="28DaysLater"; fuel="slash"
  indA = which(flags$transec_source_fire_namestr ==fire)
  for (i in 1:length(indA)){
    print(c(fire,i) ) # - no pictures
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 115;stop = stopO - 8; startB = startO -1; stopB = startO + 1}
    if (i == 2){start = startO + 27;stop = stopO - 3.8; startB = startO; stopB =  startO+2} 
    if (i == 3){start = startO + 0;stop = stopO - 18; startB = 64092; stopB =  64094} #  split from above
    if (i == 4){start = startO + 27;stop = stopO - 114; startB = 64092; stopB =  64094} #  split from above
    if (i == 5){start = startO + 38;stop = stopO - 12; startB = startO -1; stopB = startO + 1}
    if (i == 6){start = startO + 21;stop = stopO - 14; startB = startO + 0; stopB = startO + 2}
    if (i == 7){start = startO + 8.4;stop = stopO - 242; startB = startO + 0; stopB = startO + 2}
    if (i == 8){start = startO + 27.2;stop = stopO - 57; startB = startO + 0; stopB = startO + 2}
    if (i == 9){start = startO + 15;stop = stopO - 65; startB = startO + 0; stopB = startO + 2} # bad pass
    if (i == 10){start = startO + 23;stop = stopO - 44; startB = startO + 0; stopB = startO + 2}
    if (i == 11){start = startO + 46;stop = stopO - 13; startB = startO +1; stopB = startO +2} # bad background
    if (i == 12){start = startO + 15;stop = stopO - 28; startB = startO + 0; stopB = startO + 2} # bad pass, CO < 200 ppb
    tmp = time_align(start,stop,co2.826.5hz,  co.ch4.826.5hz, 
                     warneke.826.5hz,  isaf.826.5hz,
                     rollinsno.826.5hz, rollinsso2.826.5hz, 
                     cit.826.5hz, gtcims.826.5hz,moore.826fast, jimenez.826.5hz,met.826.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i>=1){start=start-1; startB=startB-1;stopB=stopB-1} # apply to all?
    aug26th.fire= time_alignSLOWNOCANS(start,stop, co2.826.1hz, co.ch4.826.1hz, 
                                        warneke.826.1hz, isaf.826.1hz,rollinsno.826.1hz, rollinsso2.826.1hz, 
                                        cit.826.1hz,gtcims.826.1hz,ryerson.826.1hz , jimenez.826.1hz, 
                                        schwarz.826.1hz, freid.c2h6.826.1hz, freid.ch2o.826.1hz, 
                                        womack.826.1hz, stclair.826.1hz,veres.826.1hz, #wisthaler.826.1hz,
                                        moore.826, met.826.1hz)
    ind2 = which(aug26th.fire$Time_Start >= start & aug26th.fire$Time_Start  <= stop	) 
    ind2B = which(aug26th.fire$Time_Start >= startB & aug26th.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug26th.fire$fire = ''; aug26th.fire$fuel = ''
    tmp$fuel[ind] = 'slash' ; tmp$fire[ind] = 'DaysLater'
    aug26th.fire$fuel[ind2] = 'slash' ; aug26th.fire$fire[ind2] = 'DaysLater'
    # -------- Blake cans 
    indBLAKE = which(blake.826.merge$Time_Start >= (start-5) & blake.826.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.826.merge[indBLAKE,]
      blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.826.merge[1,]
      blakeBG = blake.826.merge[1,]
    }
    
    # Gilman cans
    indGILMAN = which(gilman.826.merge$Time_Start >= (start-0) & gilman.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.826.merge[indGILMAN,]
      GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.826.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.826.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.826.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.826.merge$Time_Start >= (start-0) & apel.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.826.merge[indAPEL,]
      APELBG = apel.826.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.826.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.826.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.826.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.826.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.826$Time_Start >= (start-0) & newTOGA.826$Time_Start <= stop & newTOGA.826$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.826[indBECKY,]
      BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.826$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.826[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.826[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug26th.fire[ind2,])
      plotpass5hz(aug26th.fire[ind2,])
      plotpass1hz(aug26th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug26th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug26th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug26th.fire[ind2,],aug26th.fire[ind2B,],xspecies,'DaysLater', 'slash',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'DaysLater', 'slash',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop  
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){DaysLater.1hz.EF = tmpEF1 ;DaysLater.5hz.EF = tmpEF5 }
    if (i > 1){DaysLater.1hz.EF = rbind(DaysLater.1hz.EF, tmpEF1) ; DaysLater.5hz.EF = rbind(DaysLater.5hz.EF, tmpEF5) }
  }
  DaysLater.1hz.EF = DaysLater.1hz.EF[order(DaysLater.1hz.EF$variable),]
  DaysLater.1hz.EF$transect_source_fire_ID = indA[1]
  DaysLater.5hz.EF = DaysLater.5hz.EF[order(DaysLater.5hz.EF$variable),]
  DaysLater.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== ALIEN  =====================
  fire="Alien"; fuel="slash"
  indA = which(flags$transec_source_fire_namestr == fire)
  indA = indA[1:2]# the third one seems to be a duplicate
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO +12.2;stop = stopO - 73.6; startB = startO + 0; stopB = startO +2} # lets just do second peak
    if (i == 2){start = startO + 10.6;stop = stopO - 177.6; startB = startO + 0; stopB = startO + 2}
    if (i == 3){start = startO + 1;stop = stopO - 184; startB = 67711; stopB = 67713} # use background for #2
    tmp = time_align(start,stop,co2.826.5hz,  co.ch4.826.5hz, 
                     warneke.826.5hz,  isaf.826.5hz,
                     rollinsno.826.5hz, rollinsso2.826.5hz, 
                     cit.826.5hz, gtcims.826.5hz,moore.826fast, jimenez.826.5hz,met.826.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-2; startB=startB-1; stopB=stopB-1}
    if (i ==2){start=start-2; startB=startB-1; stopB=stopB-1}
    if (i ==3){start=start-2; startB=startB-1; stopB=stopB-1}
    aug26th.fire= time_alignSLOWNOCANS(start,stop, co2.826.1hz, co.ch4.826.1hz, 
                                       warneke.826.1hz, isaf.826.1hz,rollinsno.826.1hz, rollinsso2.826.1hz, 
                                       cit.826.1hz,gtcims.826.1hz,ryerson.826.1hz , jimenez.826.1hz, 
                                       schwarz.826.1hz, freid.c2h6.826.1hz, freid.ch2o.826.1hz, 
                                       womack.826.1hz, stclair.826.1hz,veres.826.1hz, #wisthaler.826.1hz,
                                       moore.826, met.826.1hz)
    ind2 = which(aug26th.fire$Time_Start >= start & aug26th.fire$Time_Start  <= stop	) 
    ind2B = which(aug26th.fire$Time_Start >= startB & aug26th.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug26th.fire$fire = ''; aug26th.fire$fuel = ''
    tmp$fuel[ind] = 'slash' ; tmp$fire[ind] = 'Alien'
    aug26th.fire$fuel[ind2] = 'slash' ; aug26th.fire$fire[ind2] = 'Alien'
    # -------- Blake cans 
    indBLAKE = which(blake.826.merge$Time_Start >= (start-5) & blake.826.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.826.merge[indBLAKE,]
      blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.826.merge[1,]
      blakeBG = blake.826.merge[1,]
    }
    
    # Gilman cans
    indGILMAN = which(gilman.826.merge$Time_Start >= (start-0) & gilman.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.826.merge[indGILMAN,]
      GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.826.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.826.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.826.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.826.merge$Time_Start >= (start-0) & apel.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.826.merge[indAPEL,]
      APELBG = apel.826.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.826.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.826.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.826.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.826.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.826$Time_Start >= (start-0) & newTOGA.826$Time_Start <= stop & newTOGA.826$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.826[indBECKY,]
      BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.826$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.826[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.826[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug26th.fire[ind2,])
      plotpass5hz(aug26th.fire[ind2,])
      plotpass1hz(aug26th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug26th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug26th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug26th.fire[ind2,],aug26th.fire[ind2B,],xspecies,'Alien', 'slash',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Alien', 'slash',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop   
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Alien.1hz.EF = tmpEF1 ;Alien.5hz.EF = tmpEF5 }
    if (i > 1){Alien.1hz.EF = rbind(Alien.1hz.EF, tmpEF1) ; Alien.5hz.EF = rbind(Alien.5hz.EF, tmpEF5) }
  }
  Alien.1hz.EF = Alien.1hz.EF[order(Alien.1hz.EF$variable),]
  Alien.1hz.EF$transect_source_fire_ID = indA[1]
  Alien.5hz.EF = Alien.5hz.EF[order(Alien.5hz.EF$variable),]
  Alien.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== BAMBI  =====================
  fire="Bambi"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 0.4;stop = stopO - 0; startB = startO -1; stopB = startO +.6}
    if (i == 2){start = startO + 0;stop = stopO - 0; startB =  69384; stopB = 69385.6}# split from above
    if (i == 3){start = startO + 0;stop = stopO - 3; startB =  69384; stopB = 69385.6} # split from above
    if (i == 4){start = startO + 1.6;stop = stopO - 1.2; startB = 69384; stopB = 69385.6} # split from above
    if (i == 5){start = startO + 1;stop = stopO -54.6; startB = 69384; stopB = 69385.6} # split from above
    #@Jim 70024.2 to 70120.6
    tmp = time_align(start,stop,co2.826.5hz,  co.ch4.826.5hz, 
                     warneke.826.5hz,  isaf.826.5hz,
                     rollinsno.826.5hz, rollinsso2.826.5hz, 
                     cit.826.5hz, gtcims.826.5hz,moore.826fast, jimenez.826.5hz,met.826.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i == 1){start = start-2; startB=startB-1;stopB=stopB-1}
    if (i == 2){start = start-1; startB=startB-1;stopB=stopB-1}
    if (i == 3){start = start-0;stop=stop+1; startB=startB-1;stopB=stopB-1}
    if (i == 4){start = start-1;startB=startB-1;stopB=stopB-1}
    if (i == 5){start = start-1;startB=startB-1;stopB=stopB-1}
    aug26th.fire= time_alignSLOWNOCANS(start,stop, co2.826.1hz, co.ch4.826.1hz, 
                                       warneke.826.1hz, isaf.826.1hz,rollinsno.826.1hz, rollinsso2.826.1hz, 
                                       cit.826.1hz,gtcims.826.1hz,ryerson.826.1hz , jimenez.826.1hz, 
                                       schwarz.826.1hz, freid.c2h6.826.1hz, freid.ch2o.826.1hz, 
                                       womack.826.1hz, stclair.826.1hz,veres.826.1hz, #wisthaler.826.1hz,
                                       moore.826, met.826.1hz)
    ind2 = which(aug26th.fire$Time_Start >= start & aug26th.fire$Time_Start  <= stop	) 
    ind2B = which(aug26th.fire$Time_Start >= startB & aug26th.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug26th.fire$fire = ''; aug26th.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Bambi'
    aug26th.fire$fuel[ind2] = fuel ; aug26th.fire$fire[ind2] = 'Bambi'
    # -------- Blake cans 
    indBLAKE = which(blake.826.merge$Time_Start >= (start-5) & blake.826.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.826.merge[indBLAKE,]
      blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.826.merge[1,]
      blakeBG = blake.826.merge[1,]
    }
    
    # Gilman cans
    indGILMAN = which(gilman.826.merge$Time_Start >= (start-0) & gilman.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.826.merge[indGILMAN,]
      GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.826.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.826.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.826.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.826.merge$Time_Start >= (start-0) & apel.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.826.merge[indAPEL,]
      APELBG = apel.826.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.826.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.826.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.826.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.826.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.826$Time_Start >= (start-0) & newTOGA.826$Time_Start <= stop & newTOGA.826$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.826[indBECKY,]
      BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.826$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.826[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.826[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug26th.fire[ind2,])
      plotpass5hz(aug26th.fire[ind2,])
      plotpass1hz(aug26th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug26th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug26th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug26th.fire[ind2,],aug26th.fire[ind2B,],xspecies,'Bambi', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Bambi', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop    
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Bambi.1hz.EF = tmpEF1 ;Bambi.5hz.EF = tmpEF5 }
    if (i > 1){Bambi.1hz.EF = rbind(Bambi.1hz.EF, tmpEF1) ; Bambi.5hz.EF = rbind(Bambi.5hz.EF, tmpEF5) }
  }
  Bambi.1hz.EF = Bambi.1hz.EF[order(Bambi.1hz.EF$variable),]
  Bambi.1hz.EF$transect_source_fire_ID = indA[1]
  Bambi.5hz.EF = Bambi.5hz.EF[order(Bambi.5hz.EF$variable),]
  Bambi.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== BAMBI JR =====================
  indA = which(flags$transec_source_fire_namestr == 'Bambi Jr')
  fire="Bambi Jr"; fuel="pile"
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 1;stop = stopO -7; startB = startO-1; stopB =startO+0} # 
    if (i == 2){start = startO + 0.4;stop = stopO -4.6; startB = startO-2; stopB =startO+0} # 
    if (i == 3){start = startO + 2;stop = stopO -116; startB = startO+0; stopB =startO+2} # 
    #@Jim 70024.2 to 70120.6
    tmp = time_align(start,stop,co2.826.5hz,  co.ch4.826.5hz, 
                     warneke.826.5hz,  isaf.826.5hz,
                     rollinsno.826.5hz, rollinsso2.826.5hz, 
                     cit.826.5hz, gtcims.826.5hz,moore.826fast, jimenez.826.5hz,met.826.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i == 1){start = start-1;startB=startB-1;stopB=stopB-1}
    if (i == 2){start = start-1;startB=startB-1;stopB=stopB-1}
    if (i == 3){start = start-1;startB=startB-1;stopB=stopB-1}
    aug26th.fire= time_alignSLOWNOCANS(start,stop, co2.826.1hz, co.ch4.826.1hz, 
                                       warneke.826.1hz, isaf.826.1hz,rollinsno.826.1hz, rollinsso2.826.1hz, 
                                       cit.826.1hz,gtcims.826.1hz,ryerson.826.1hz , jimenez.826.1hz, 
                                       schwarz.826.1hz, freid.c2h6.826.1hz, freid.ch2o.826.1hz, 
                                       womack.826.1hz, stclair.826.1hz,veres.826.1hz, #wisthaler.826.1hz,
                                       moore.826, met.826.1hz)
    ind2 = which(aug26th.fire$Time_Start >= start & aug26th.fire$Time_Start  <= stop	) 
    ind2B = which(aug26th.fire$Time_Start >= startB & aug26th.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug26th.fire$fire = ''; aug26th.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Bambi Jr'
    aug26th.fire$fuel[ind2] = fuel ; aug26th.fire$fire[ind2] = 'Bambi Jr'
    # -------- Blake cans 
    indBLAKE = which(blake.826.merge$Time_Start >= (start-5) & blake.826.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.826.merge[indBLAKE,]
      blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.826.merge[1,]
      blakeBG = blake.826.merge[1,]
    }
    
    # Gilman cans
    indGILMAN = which(gilman.826.merge$Time_Start >= (start-0) & gilman.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.826.merge[indGILMAN,]
      GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.826.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.826.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.826.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.826.merge$Time_Start >= (start-0) & apel.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.826.merge[indAPEL,]
      APELBG = apel.826.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.826.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.826.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.826.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.826.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.826$Time_Start >= (start-0) & newTOGA.826$Time_Start <= stop & newTOGA.826$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.826[indBECKY,]
      BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.826$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.826[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.826[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug26th.fire[ind2,])
      plotpass5hz(aug26th.fire[ind2,])
      plotpass1hz(aug26th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug26th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug26th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug26th.fire[ind2,],aug26th.fire[ind2B,],xspecies,'Bambi Jr', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Bambi Jr', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop    
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){BambiJr.1hz.EF = tmpEF1 ;BambiJr.5hz.EF = tmpEF5 }
    if (i > 1){BambiJr.1hz.EF = rbind(BambiJr.1hz.EF, tmpEF1) ; BambiJr.5hz.EF = rbind(BambiJr.5hz.EF, tmpEF5) }
  }
  BambiJr.1hz.EF = BambiJr.1hz.EF[order(BambiJr.1hz.EF$variable),]
  BambiJr.1hz.EF$transect_source_fire_ID = indA[1]
  BambiJr.5hz.EF = BambiJr.5hz.EF[order(BambiJr.5hz.EF$variable),]
  BambiJr.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== DEADPOOL  =====================
  
  fire="Deadpool"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 1;stop = stopO -13; startB = startO -1; stopB = startO + 1}
    if (i == 2){start = startO +12.6;stop = stopO -13.2; startB = startO + 0; stopB = startO + 2}
    if (i == 3){start = startO + 11.8;stop = stopO - 18.0; startB = startO + 0; stopB = startO + 2} # need to look at pictures on this one
    if (i == 4){start = startO + 11.2;stop = stopO - 13; startB = startO + 0; stopB =startO + 2} # need to look at pictures on this one
    tmp = time_align(start,stop,co2.826.5hz,  co.ch4.826.5hz, 
                     warneke.826.5hz,  isaf.826.5hz,
                     rollinsno.826.5hz, rollinsso2.826.5hz, 
                     cit.826.5hz, gtcims.826.5hz,moore.826fast, jimenez.826.5hz,met.826.5hz)
    #@Jim#2 74012.4 to 74029.6
    #@Jim#4 74573.6 to 74591.2
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1; startB=startB-1;stopB=stopB-1}
    if (i ==2){start=start-1; startB=startB-1;stopB=stopB-1}
    if (i ==3){start=start-1; startB=startB-1;stopB=stopB-1}
    if (i ==4){start=start-1; startB=startB-1;stopB=stopB-1}
    aug26th.fire= time_alignSLOWNOCANS(start,stop, co2.826.1hz, co.ch4.826.1hz, 
                                       warneke.826.1hz, isaf.826.1hz,rollinsno.826.1hz, rollinsso2.826.1hz, 
                                       cit.826.1hz,gtcims.826.1hz,ryerson.826.1hz , jimenez.826.1hz, 
                                       schwarz.826.1hz, freid.c2h6.826.1hz, freid.ch2o.826.1hz, 
                                       womack.826.1hz, stclair.826.1hz,veres.826.1hz, #wisthaler.826.1hz,
                                       moore.826, met.826.1hz)
    ind2 = which(aug26th.fire$Time_Start >= start & aug26th.fire$Time_Start  <= stop	) 
    ind2B = which(aug26th.fire$Time_Start >= startB & aug26th.fire$Time_Start <= stopB	) 
    
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug26th.fire$fire = ''; aug26th.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Deadpool'
    aug26th.fire$fuel[ind2] = fuel ; aug26th.fire$fire[ind2] = 'Deadpool'
    # -------- Blake cans 
    indBLAKE = which(blake.826.merge$Time_Start >= (start-5) & blake.826.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.826.merge[indBLAKE,]
      blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.826.merge[1,]
      blakeBG = blake.826.merge[1,]
    }
    
    # Gilman cans
    indGILMAN = which(gilman.826.merge$Time_Start >= (start-0) & gilman.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.826.merge[indGILMAN,]
      GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.826.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.826.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.826.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.826.merge$Time_Start >= (start-0) & apel.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.826.merge[indAPEL,]
      APELBG = apel.826.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.826.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.826.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.826.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.826.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.826$Time_Start >= (start-0) & newTOGA.826$Time_Start <= stop & newTOGA.826$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.826[indBECKY,]
      BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.826$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.826[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.826[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug26th.fire[ind2,])
      plotpass5hz(aug26th.fire[ind2,])
      plotpass1hz(aug26th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug26th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug26th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug26th.fire[ind2,],aug26th.fire[ind2B,],xspecies,'Deadpool', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Deadpool', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop 
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Deadpool.1hz.EF = tmpEF1 ;Deadpool.5hz.EF = tmpEF5 }
    if (i > 1){Deadpool.1hz.EF = rbind(Deadpool.1hz.EF, tmpEF1) ; Deadpool.5hz.EF = rbind(Deadpool.5hz.EF, tmpEF5) }
  }
  Deadpool.1hz.EF = Deadpool.1hz.EF[order(Deadpool.1hz.EF$variable),]
  Deadpool.1hz.EF$transect_source_fire_ID = indA[1]
  Deadpool.5hz.EF = Deadpool.5hz.EF[order(Deadpool.5hz.EF$variable),]
  Deadpool.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== ELF  =====================
  
  fire="Elf"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 2.2;stop = stop -4.2; startB = startO + 0; stopB = startO + 2}
    if (i == 2){start = startO + 0;stop = stop -4; startB = 74280; stopB = 74282} # split from above
    if (i == 3){start = startO + 3;stop = stop -11; startB = startO + 0; stopB = startO + 2} # get rid of second peak
    if (i == 4){start = startO + 5;stop = stop -14.4; startB = startO -1; stopB = startO +.8}
    tmp = time_align(start,stop,co2.826.5hz,  co.ch4.826.5hz, 
                     warneke.826.5hz,  isaf.826.5hz,
                     rollinsno.826.5hz, rollinsso2.826.5hz, 
                     cit.826.5hz, gtcims.826.5hz,moore.826fast, jimenez.826.5hz,met.826.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i ==2){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i ==3){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i ==4){start=start-1; startB=startB-1;stopB=stopB-1}
    aug26th.fire= time_alignSLOWNOCANS(start,stop, co2.826.1hz, co.ch4.826.1hz, 
                                       warneke.826.1hz, isaf.826.1hz,rollinsno.826.1hz, rollinsso2.826.1hz, 
                                       cit.826.1hz,gtcims.826.1hz,ryerson.826.1hz , jimenez.826.1hz, 
                                       schwarz.826.1hz, freid.c2h6.826.1hz, freid.ch2o.826.1hz, 
                                       womack.826.1hz, stclair.826.1hz,veres.826.1hz, #wisthaler.826.1hz,
                                       moore.826, met.826.1hz)
    ind2 = which(aug26th.fire$Time_Start >= start & aug26th.fire$Time_Start  <= stop	) 
    ind2B = which(aug26th.fire$Time_Start >= startB & aug26th.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug26th.fire$fire = ''; aug26th.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Elf'
    aug26th.fire$fuel[ind2] = fuel ; aug26th.fire$fire[ind2] = 'Elf'
    # -------- Blake cans 
    indBLAKE = which(blake.826.merge$Time_Start >= (start-5) & blake.826.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.826.merge[indBLAKE,]
      blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.826.merge[1,]
      blakeBG = blake.826.merge[1,]
    }
    
    # Gilman cans
    indGILMAN = which(gilman.826.merge$Time_Start >= (start-0) & gilman.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.826.merge[indGILMAN,]
      GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.826.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.826.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.826.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.826.merge$Time_Start >= (start-0) & apel.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.826.merge[indAPEL,]
      APELBG = apel.826.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.826.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.826.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.826.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.826.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.826$Time_Start >= (start-0) & newTOGA.826$Time_Start <= stop & newTOGA.826$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.826[indBECKY,]
      BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.826$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.826[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.826[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug26th.fire[ind2,])
      plotpass5hz(aug26th.fire[ind2,])
      plotpass1hz(aug26th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug26th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug26th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug26th.fire[ind2,],aug26th.fire[ind2B,],xspecies,'Elf', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Elf', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop  
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Elf.1hz.EF = tmpEF1 ;Elf.5hz.EF = tmpEF5 }
    if (i > 1){Elf.1hz.EF = rbind(Elf.1hz.EF, tmpEF1) ; Elf.5hz.EF = rbind(Elf.5hz.EF, tmpEF5) }
  }
  Elf.1hz.EF = Elf.1hz.EF[order(Elf.1hz.EF$variable),]
  Elf.1hz.EF$transect_source_fire_ID = indA[1]
  Elf.5hz.EF = Elf.5hz.EF[order(Elf.5hz.EF$variable),]
  Elf.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== FARGO  =====================
  fire="Fargo"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 5.2;   stop = stopO - 9; startB = startO + 0; stopB = startO + 2}
    if (i == 2){start = startO + 96;  stop = stopO - 14; startB = startO -1; stopB = startO + .6}
    if (i == 3){start = startO + 60.2 ;stop = stopO - 5; startB = startO -1; stopB = startO + 1}
    if (i == 4){start = startO + 1;   stop = stopO - 55.2;startB = startO -2; stopB = startO + 0}
    if (i == 6){start = startO + 141;stop = stopO - 73.2; startB = startO -1; stopB = startO + 1}
    tmp = time_align(start,stop,co2.826.5hz,  co.ch4.826.5hz, 
                     warneke.826.5hz,  isaf.826.5hz,
                     rollinsno.826.5hz, rollinsso2.826.5hz, 
                     cit.826.5hz, gtcims.826.5hz,moore.826fast, jimenez.826.5hz,met.826.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-2; stopB=stopB-1; startB=startB-1}
    if (i ==2){start=start-1; stopB=stopB-1; startB=startB-1}
    if (i ==3){start=start-1; stopB=stopB-1; startB=startB-1}
    if (i ==4){start=start-1;stop=stop+1;stopB=stopB-1; startB=startB-1}
    if (i ==5){start=start-1;stop=stop+1;stopB=stopB-1; startB=startB-1}
    if (i ==6){start=start-1;stop=stop+1;stopB=stopB-1; startB=startB-1}
    aug26th.fire= time_alignSLOWNOCANS(start,stop, co2.826.1hz, co.ch4.826.1hz, 
                                       warneke.826.1hz, isaf.826.1hz,rollinsno.826.1hz, rollinsso2.826.1hz, 
                                       cit.826.1hz,gtcims.826.1hz,ryerson.826.1hz , jimenez.826.1hz, 
                                       schwarz.826.1hz, freid.c2h6.826.1hz, freid.ch2o.826.1hz, 
                                       womack.826.1hz, stclair.826.1hz,veres.826.1hz, #wisthaler.826.1hz,
                                       moore.826, met.826.1hz)
    ind2 = which(aug26th.fire$Time_Start >= start & aug26th.fire$Time_Start  <= stop	)
    ind2B = which(aug26th.fire$Time_Start >= startB & aug26th.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug26th.fire$fire = ''; aug26th.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Fargo'
    aug26th.fire$fuel[ind2] = fuel ; aug26th.fire$fire[ind2] = 'Fargo'
    # -------- Blake cans 
    indBLAKE = which(blake.826.merge$Time_Start >= (start-5) & blake.826.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.826.merge[indBLAKE,]
      blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.826.merge[1,]
      blakeBG = blake.826.merge[1,]
    }
    
    # Gilman cans
    indGILMAN = which(gilman.826.merge$Time_Start >= (start-0) & gilman.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      if (length(indGILMAN) == 2 & i == 2){indGILMAN = indGILMAN[1]}
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.826.merge[indGILMAN,]
      GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.826.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.826.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.826.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.826.merge$Time_Start >= (start-0) & apel.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.826.merge[indAPEL,]
      APELBG = apel.826.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.826.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.826.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.826.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.826.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.826$Time_Start >= (start-0) & newTOGA.826$Time_Start <= stop & newTOGA.826$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.826[indBECKY,]
      BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.826$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.826[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.826[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug26th.fire[ind2,])
      plotpass5hz(aug26th.fire[ind2,])
      plotpass1hz(aug26th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug26th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug26th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug26th.fire[ind2,],aug26th.fire[ind2B,],xspecies,'Fargo', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Fargo', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Fargo.1hz.EF = tmpEF1 ;Fargo.5hz.EF = tmpEF5 }
    if (i > 1){Fargo.1hz.EF = rbind(Fargo.1hz.EF, tmpEF1) ; Fargo.5hz.EF = rbind(Fargo.5hz.EF, tmpEF5) }
    
  }
  Fargo.1hz.EF = Fargo.1hz.EF[order(Fargo.1hz.EF$variable),]
  Fargo.1hz.EF$transect_source_fire_ID = indA[1]
  Fargo.5hz.EF = Fargo.5hz.EF[order(Fargo.5hz.EF$variable),]
  Fargo.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== HELLBOY  =====================
  fire="hellboy"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c( fire,i) ) 
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 7;stop = stop - 66.6; startB = startO + 0; stopB = startO + 2} # 
    if (i == 2){start = startO + 2;stop = stop - 28; startB = startO -1; stopB = startO + 1}
    if (i == 3){start = startO + 1;stop = stop - 248; startB = startO -2; stopB = startO +0}
    if (i == 4){start = startO + 17.6;stop = stop - 24; startB = startO -1; stopB = startO +1}
    if (i == 5){start = startO + 6.4;stop = stop - 315; startB = startO + 0; stopB = startO +2}
    tmp = time_align(start,stop,co2.826.5hz,  co.ch4.826.5hz, 
                     warneke.826.5hz,  isaf.826.5hz,
                     rollinsno.826.5hz, rollinsso2.826.5hz, 
                     cit.826.5hz, gtcims.826.5hz,moore.826fast, jimenez.826.5hz,met.826.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1; startB=startB-1;stopB=stopB-1}
    if (i ==2){start=start-1; startB=startB-1;stopB=stopB-1}
    if (i ==3){start=start-1; startB=startB-1;stopB=stopB-1}
    if (i ==4){start=start-1; startB=startB-1;stopB=stopB-1}
    if (i ==5){start=start-1; startB=startB-1;stopB=stopB-1}
    aug26th.fire= time_alignSLOWNOCANS(start,stop, co2.826.1hz, co.ch4.826.1hz, 
                                       warneke.826.1hz, isaf.826.1hz,rollinsno.826.1hz, rollinsso2.826.1hz, 
                                       cit.826.1hz,gtcims.826.1hz,ryerson.826.1hz , jimenez.826.1hz, 
                                       schwarz.826.1hz, freid.c2h6.826.1hz, freid.ch2o.826.1hz, 
                                       womack.826.1hz, stclair.826.1hz,veres.826.1hz, #wisthaler.826.1hz,
                                       moore.826, met.826.1hz)
    ind2 = which(aug26th.fire$Time_Start >= start & aug26th.fire$Time_Start  <= stop	) 
    ind2B = which(aug26th.fire$Time_Start >= startB & aug26th.fire$Time_Start <= stopB	) 

    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug26th.fire$fire = ''; aug26th.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Hellboy'
    aug26th.fire$fuel[ind2] = fuel ; aug26th.fire$fire[ind2] = 'Hellboy'
    # -------- Blake cans 
    indBLAKE = which(blake.826.merge$Time_Start >= (start-5) & blake.826.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.826.merge[indBLAKE,]
      blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.826.merge[1,]
      blakeBG = blake.826.merge[1,]
    }
    
    # Gilman cans
    indGILMAN = which(gilman.826.merge$Time_Start >= (start-0) & gilman.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.826.merge[indGILMAN,]
      GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.826.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.826.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.826.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.826.merge$Time_Start >= (start-0) & apel.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.826.merge[indAPEL,]
      APELBG = apel.826.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.826.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.826.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.826.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.826.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.826$Time_Start >= (start-0) & newTOGA.826$Time_Start <= stop & newTOGA.826$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.826[indBECKY,]
      BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.826$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.826[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.826[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug26th.fire[ind2,])
      plotpass5hz(aug26th.fire[ind2,])
      plotpass1hz(aug26th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCH4(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug26th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug26th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug26th.fire[ind2,],aug26th.fire[ind2B,],xspecies,'Hellboy', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Hellboy', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop  
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Hellboy.1hz.EF = tmpEF1 ;Hellboy.5hz.EF = tmpEF5 }
    if (i > 1){Hellboy.1hz.EF = rbind(Hellboy.1hz.EF, tmpEF1) ; Hellboy.5hz.EF = rbind(Hellboy.5hz.EF, tmpEF5) }
    
  }
  Hellboy.1hz.EF = Hellboy.1hz.EF[order(Hellboy.1hz.EF$variable),]
  Hellboy.1hz.EF$transect_source_fire_ID = indA[1]
  Hellboy.5hz.EF = Hellboy.5hz.EF[order(Hellboy.5hz.EF$variable),]
  Hellboy.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== INVICTUS  =====================
  fire="Invictus"; fuel="forest"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) ) # check plumes again unless it is really coniferous/deciduous
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO  
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 3;stop = stop - 0; startB = startO + 0; stopB = startO +2}
    if (i == 2){start = startO + 1;stop = stop -3; startB = 85110; stopB =  85112} # split from above
    if (i == 3){start = startO + 6;stop = stop - 31.4; startB = 85110; stopB =  85112} # split from above
    if (i == 4){start = startO + 18;stop = stop - 14.6; startB = startO + 0; stopB = startO+1}
    if (i == 5){start = startO + 4;stop = stop - 5.6; startB = startO + 0; stopB = startO +2} # 
    if (i == 6){start = startO + 2;stop = stop - 2; startB = startO + 0; stopB = startO +2} # 
    if (i == 7){start = startO + .6;stop = stop - 40; startB = 86001; stopB = 86003} # split from above
    if (i == 8){start = startO + 8;stop = stop - 13; startB = startO-1; stopB = startO+1} # bad pass
    tmp = time_align(start,stop,co2.826.5hz,  co.ch4.826.5hz, 
                     warneke.826.5hz,  isaf.826.5hz,
                     rollinsno.826.5hz, rollinsso2.826.5hz, 
                     cit.826.5hz, gtcims.826.5hz,moore.826fast, jimenez.826.5hz,met.826.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1; startB=startB-1;stopB=stopB-1}
    if (i ==2){start=start-1; startB=startB-1;stopB=stopB-1}
    if (i ==3){start=start-1; startB=startB-1;stopB=stopB-1}
    if (i ==4){start=start-1; startB=startB-1;stopB=stopB-1}
    if (i ==5){start=start-1; startB=startB-1;stopB=stopB-1}
    if (i ==6){start=start-1; startB=startB-1;stopB=stopB-1}
    if (i ==7){start=start-2; startB=startB-1;stopB=stopB-1}
    aug26th.fire= time_alignSLOWNOCANS(start,stop, co2.826.1hz, co.ch4.826.1hz, 
                                       warneke.826.1hz, isaf.826.1hz,rollinsno.826.1hz, rollinsso2.826.1hz, 
                                       cit.826.1hz,gtcims.826.1hz,ryerson.826.1hz , jimenez.826.1hz, 
                                       schwarz.826.1hz, freid.c2h6.826.1hz, freid.ch2o.826.1hz, 
                                       womack.826.1hz, stclair.826.1hz,veres.826.1hz, #wisthaler.826.1hz,
                                       moore.826, met.826.1hz)
    ind2 = which(aug26th.fire$Time_Start >= start & aug26th.fire$Time_Start  <= stop	) 
    ind2B = which(aug26th.fire$Time_Start >= startB & aug26th.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug26th.fire$fire = ''; aug26th.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] =fire
    aug26th.fire$fuel[ind2] =fuel; aug26th.fire$fire[ind2] = fire
    # -------- Blake cans 
    indBLAKE = which(blake.826.merge$Time_Start >= (start-5) & blake.826.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.826.merge[indBLAKE,]
      blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.826.merge[1,]
      blakeBG = blake.826.merge[1,]
    }
    
    # Gilman cans
    indGILMAN = which(gilman.826.merge$Time_Start >= (start-0) & gilman.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.826.merge[indGILMAN,]
      GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.826.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.826.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.826.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.826.merge$Time_Start >= (start-0) & apel.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.826.merge[indAPEL,]
      APELBG = apel.826.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.826.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.826.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.826.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.826.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.826$Time_Start >= (start-0) & newTOGA.826$Time_Start <= stop & newTOGA.826$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.826[indBECKY,]
      BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.826$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.826[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.826[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug26th.fire[ind2,])
      plotpass5hz(aug26th.fire[ind2,])
      plotpass1hz(aug26th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug26th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug26th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug26th.fire[ind2,],aug26th.fire[ind2B,],xspecies,fire,fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Invictus', 'coniferous/decidous',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF5$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF5$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop     
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Invictus.1hz.EF = tmpEF1 ;Invictus.5hz.EF = tmpEF5 }
    if (i > 1){Invictus.1hz.EF = rbind(Invictus.1hz.EF, tmpEF1) ; Invictus.5hz.EF = rbind(Invictus.5hz.EF, tmpEF5) }
    
  }
  Invictus.1hz.EF = Invictus.1hz.EF[order(Invictus.1hz.EF$variable),]
  Invictus.1hz.EF$transect_source_fire_ID = indA[1]
  Invictus.5hz.EF = Invictus.5hz.EF[order(Invictus.5hz.EF$variable),]
  Invictus.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== INVICTUS-Unknown  =====================
  indA = which(flags$transec_source_fire_namestr == 'Invictus_unknown') # - bad fire
  indA = which(flags$transec_source_fire_namestr ==  'Invictus_unknown2')# check plumes again unless it is really coniferous/deciduous
  fire="Invictus_unknown2"; fuel="coniferous/decidous"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 7;stop = stop - 17; startB = startO + 0; stopB = startO + 2}
    tmp = time_align(start,stop,co2.826.5hz,  co.ch4.826.5hz, 
                     warneke.826.5hz,  isaf.826.5hz,
                     rollinsno.826.5hz, rollinsso2.826.5hz, 
                     cit.826.5hz, gtcims.826.5hz,moore.826fast, jimenez.826.5hz,met.826.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    aug26th.fire= time_alignSLOWNOCANS(start,stop, co2.826.1hz, co.ch4.826.1hz, 
                                       warneke.826.1hz, isaf.826.1hz,rollinsno.826.1hz, rollinsso2.826.1hz, 
                                       cit.826.1hz,gtcims.826.1hz,ryerson.826.1hz , jimenez.826.1hz, 
                                       schwarz.826.1hz, freid.c2h6.826.1hz, freid.ch2o.826.1hz, 
                                       womack.826.1hz, stclair.826.1hz,veres.826.1hz, #wisthaler.826.1hz,
                                       moore.826, met.826.1hz)
    ind2 = which(aug26th.fire$Time_Start >= start & aug26th.fire$Time_Start  <= stop	) 
    ind2B = which(aug26th.fire$Time_Start >= startB & aug26th.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug26th.fire$fire = ''; aug26th.fire$fuel = ''
    tmp$fuel[ind] = 'coniferous/decidous' ; tmp$fire[ind] = 'InvictusU'
    aug26th.fire$fuel[ind2] = 'coniferous/decidous' ; aug26th.fire$fire[ind2] = 'InvictusU'
    # -------- Blake cans 
    indBLAKE = which(blake.826.merge$Time_Start >= (start-5) & blake.826.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.826.merge[indBLAKE,]
      blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.826.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.826.merge[1,]
      blakeBG = blake.826.merge[1,]
    }
    
    # Gilman cans
    indGILMAN = which(gilman.826.merge$Time_Start >= (start-0) & gilman.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.826.merge[indGILMAN,]
      GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.826.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.826.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.826.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.826.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.826.merge$Time_Start >= (start-0) & apel.826.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.826.merge[indAPEL,]
      APELBG = apel.826.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.826.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.826.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.826.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.826.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.826$Time_Start >= (start-0) & newTOGA.826$Time_Start <= stop & newTOGA.826$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.826[indBECKY,]
      BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.826$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.826[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.826[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.826[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug26th.fire[ind2,])
      plotpass5hz(aug26th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      plotslopes5hz(aug26th.fire[ind2,])
      plotslopes5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug26th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug26th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug26th.fire[ind2,],aug26th.fire[ind2B,],xspecies,'InvictusU', 'coniferous/decidous',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug26th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug26th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug26th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'InvictusU', 'coniferous/decidous',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug26th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug26th.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop   
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){InvictusU.1hz.EF = tmpEF1 ;InvictusU.5hz.EF = tmpEF5 }
    if (i > 1){InvictusU.1hz.EF = rbind(InvictusU.1hz.EF, tmpEF1) ; InvictusU.5hz.EF = rbind(InvictusU.5hz.EF, tmpEF5) }
  }
  InvictusU.1hz.EF = InvictusU.1hz.EF[order(InvictusU.1hz.EF$variable),]
  InvictusU.1hz.EF$transect_source_fire_ID = indA[1]
  InvictusU.5hz.EF = InvictusU.5hz.EF[order(InvictusU.5hz.EF$variable),]
  InvictusU.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ----- [[[[[[[[[[[[[[[[[[[[ Aug 29th ]]]]]]]]]]]]]]]]]]]]] -------
  # ------------- Get 1 Hz Data -------------
  # --------#######--------- Get 1 Hz Data individual 8/29------#######---------
  # plume tags
  tags = getICARTTdataSIMPLE('InputFiles/firexaq-fire-Flags-1HZ_DC8_20190829_R9.ict') ;tags$Time_Start = tags$TIME_START
  # MET DATA
  met.829.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-MetNav_DC8_20190829_R1.ict')
  met.829.1hz = merge(met.829.1hz, tags, by='Time_Start')
  # CO2
  co2.829.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CO2-7000_DC8_20190829_R2.ict')
  # -------- DISKIN -----CO, CH4
  co.ch4.829.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-DACOM_DC8_20190829_R1.ict')
  # --------- WARNEKE ----  VOCs
  warneke.829.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-NOAAPTR-VOCs-1Hz_DC8_20190829_R3.ict')
  # ------ HANISCO - ISAF HCHO - merged to 5Hz from the online merge
  isaf.829.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-ISAF-CH2O-1Hz_DC8_20190829_R0.ict')
  #  ------- ROLLINS - SO2 and NO
  rollinsno.829.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-NO_DC8_20190829_R1.ict')
  rollinsno.829.1hz$Time_Start = rollinsno.829.1hz$time_mid
  rollinsso2.829.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-SO2_DC8_20190829_R1.ict')
  rollinsso2.829.1hz$Time_Start = rollinsso2.829.1hz$time_mid
  
  #  ----- WENNBERG - CIT VOCs - 
  cit.829.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190829_R0_CIT.ict')
  # ------ HUEY - GTCIMS PANs - not sure how to match up peaks here
  gtcims.829.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190829_R0_Huey.ict')
  
  # ------ RYERSON
  ryerson.A = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NO_DC8_20190829_R1.ict')
  ryerson.B = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NO2_DC8_20190829_R1.ict')
  ryerson.C = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NOy_DC8_20190829_R1.ict')
  ryerson.D = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-O3_DC8_20190829_R1.ict')
  ryerson.829.1hz = cbind(ryerson.A,ryerson.B,ryerson.C,ryerson.D) 
  ryerson.829.1hz$Time_Start = ryerson.829.1hz$Time_start
  # ----- JIMENEZ ---
  jimenez.829.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-EESI_DC8_20190829_R1.ict')
  
  # ----- SCHWARZ ---
  schwarz.829.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-SP2-BC-1HZ_DC8_20190829_R2.ict')
  # ----- FREID ---
  freid.c2h6.829.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-C2H6_DC8_20190829_R3.ict')
  freid.ch2o.829.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CH2O_DC8_20190829_R3.ict')
  # ------ WOMACK ---
  womackA = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-CH3COCHO_DC8_20190829_R1.ict')
  womackB = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-CHOCHO_DC8_20190829_R1.ict')
  womackC = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-HNO2_DC8_20190829_R1.ict')
  womackD = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-NO2_DC8_20190829_R1.ict')
  womack.829.1hz = cbind(womackA, womackB, womackC, womackD)
  # -------St Clair
  stclair.829.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-CANOE-NO2_DC8_20190829_R0.ict')
  stclair.829.1hz$Time_Start = stclair.829.1hz$Time_start

  # ------- VERES
  veres.A = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-ClNO2_DC8_20190829_R0.ict')
  veres.B = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HCOOH_DC8_20190829_R2.ict')
  veres.C = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HNO2_DC8_20190829_R1.ict')
  veres.D = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-N2O5_DC8_20190829_R0.ict')
  veres.E = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HPMTF_DC8_20190829_R0.ict')
  veres.F = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-CH3COOCl_DC8_20190829_R1.ict')
  veres.G = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-Cl2_DC8_20190829_R1.ict')
  veres.H = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrCl_DC8_20190829_R1.ict')
  veres.I = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrCN_DC8_20190829_R1.ict')
  veres.J = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrO_DC8_20190829_R1.ict')
  veres.K = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HCN_DC8_20190829_R1.ict')
  veres.L = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HNCO_DC8_20190829_R0.ict')
  veres.829.1hz = cbind(veres.A,veres.B,veres.C,veres.D,veres.E,veres.F,veres.G,veres.H,veres.I,veres.J,veres.K,veres.L)
  
  # --- WISTHALER
  #wisthaler.829.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190829_R0_Wisthaler.ict')

  # ---- BLAKE
  blake.829.1hz = getICARTTdataSIMPLE('InputFiles/WAS-MERGE/firexaq-mrgWAS-dc8_merge_20190829_R1.ict')
  cc = colnames(blake.829.1hz)
  blake.829.merge = blake.829.1hz[,c(1,2,96:225)]
  blake.829.merge$CO_DACOM_DISKIN_BLAKE = blake.829.1hz$CO_DACOM_DISKIN
  blake.829.merge$CO2_7000_ppm_DISKIN_BLAKE = blake.829.1hz$CO2_7000_ppm_DISKIN
  
  # ------ APEL
  apel.829.1hz = getICARTTdataSIMPLE('InputFiles/TOGA-MERGE/firexaq-mrgTOGA-dc8_merge_20190829_R1.ict')
  cc = colnames(apel.829.1hz)
  apel.829.merge = apel.829.1hz[,c(1,2,226:315)]
  apel.829.merge$CO_DACOM_DISKIN_APEL = apel.829.1hz$CO_DACOM_DISKIN
  apel.829.merge$CO2_7000_ppm_DISKIN_APEL =apel.829.1hz$CO2_7000_ppm_DISKIN
  # Becky's better merge
  file = 'InputFiles/Hornbrook/FIREX-AQ weighted TOGA merge 2022-01-24_0829.xlsx'
  newTOGA.829 = readxl::read_xlsx(file); newTOGA.829[newTOGA.829==-999] = NaN; newTOGA.829[newTOGA.829==-888] = NaN
  newTOGA.829$CO_DACOM_DISKIN_BECKY = newTOGA.829$CO_DACOM_DISKIN
  newTOGA.829$CO2_7000_ppm_DISKIN_BECKY = NaN
  newTOGA.829$Time_Start=newTOGA.829$Time_Start...4
  
  # ----GILMAN
  gilman.829.1hz = getICARTTdataSIMPLE('InputFiles/iWAS-MERGE/firexaq-mrgiWAS-dc8_merge_20190829_R1.ict')
  cc = colnames(gilman.829.1hz)
  gilman.829.merge = gilman.829.1hz[,c(1,2,316:361)]
  gilman.829.merge$CO_DACOM_DISKIN_GILMAN = gilman.829.1hz$CO_DACOM_DISKIN
  gilman.829.merge$CO2_7000_ppm_DISKIN_GILMAN = gilman.829.1hz$CO2_7000_ppm_DISKIN

  # ------ Moore 
  moore.829fast = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190829_R0_MOORE.ict')
  
  moore.829p1 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-AerosolCloudConc_DC8_20190829_R0.ict')
  moore.829p2 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-LAScold_DC8_20190829_R0.ict')
  moore.829p3 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-LAShot_DC8_20190829_R0.ict')
  moore.829p4 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-CPSPD_DC8_20190829_R0.ict')
  moore.829p5 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-CDP_DC8_20190829_R0.ict')
  moore.829p6 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-SMPS_DC8_20190829_R0.ict')
  moore.829 =merge(moore.829p1, moore.829p2, by='Time_mid', all = TRUE, incomparables = NA)
  moore.829 =merge(moore.829, moore.829p3, by='Time_mid', all = TRUE, incomparables = NA)
  moore.829 =merge(moore.829, moore.829p4, by='Time_mid', all = TRUE, incomparables = NA)
  moore.829 =merge(moore.829, moore.829p5, by='Time_mid', all = TRUE, incomparables = NA)
  moore.829 =merge(moore.829, moore.829p6, by='Time_mid', all = TRUE, incomparables = NA)
  
  # ------- append PI to colnames 1hz ----------
  cc = colnames(co2.829.1hz)
  cc[2] = paste(cc[2],'_DISKIN',sep='')
  colnames(co2.829.1hz) = cc 
  cc = colnames(co.ch4.829.1hz)
  cc[2:4] =  paste(cc[2:4],'_DISKIN',sep='')
  colnames(co.ch4.829.1hz) = cc
  cc=colnames(met.829.1hz)
  
  colnames(isaf.829.1hz) = c("Time_Start"  , "CH2O_ISAF_HANISCO" ,"CH2O_ISAF_precision_HANISCO")
  cc[2:36] = paste(cc[2:36],'_YANG',sep='')
  colnames(met.829.1hz) = cc
  cc = colnames(warneke.829.1hz)
  cc[2:43] = paste(cc[2:43],'_WARNEKE', sep='')
  colnames(warneke.829.1hz) = cc
  cc = colnames(rollinsno.829.1hz)
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsno.829.1hz) = cc
  cc = colnames(rollinsso2.829.1hz)
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsso2.829.1hz) = cc
  # make GTCIMS consistent with 1s merge
  cc=colnames(gtcims.829.1hz)
  cc[4:7] = c("PAN_GTCIMS_HUEY" , "PPN_GTCIMS_HUEY"  ,"APAN_GTCIMS_HUEY" ,"PBN_GTCIMS_HUEY")
  colnames(gtcims.829.1hz)=cc
  # since ISAF, CIT, and GTCIMS came from merge tool, alread has PI's appended.
  ryerson.829.1hz <- ryerson.829.1hz[, !duplicated(colnames(ryerson.829.1hz))]
  
  cc = colnames(ryerson.829.1hz)
  cc[2:9] = paste(cc[2:9],'_RYERSON',sep='')
  colnames(ryerson.829.1hz)=cc
  
  cc = colnames(schwarz.829.1hz)
  cc[2:3] = paste(cc[2:3],'_SCHWARZ', sep='')
  colnames(schwarz.829.1hz) =cc
  
  cc = colnames(freid.c2h6.829.1hz)
  cc[4:7] = paste(cc[4:7], '_FRIED',sep='')
  colnames(freid.c2h6.829.1hz) = cc
  
  cc = colnames(freid.ch2o.829.1hz)
  cc[4:7] = paste(cc[4:7], '_FRIED',sep='')
  colnames(freid.ch2o.829.1hz) = cc
  
  womack.829.1hz <- womack.829.1hz[, !duplicated(colnames(womack.829.1hz))]
  cc = colnames(womack.829.1hz)
  cc[2:5] = paste(cc[2:5], '_WOMACK',sep='')
  colnames(womack.829.1hz) = cc

  veres.829.1hz <- veres.829.1hz[, !duplicated(colnames(veres.829.1hz))]
  cc = colnames(veres.829.1hz)
  cc[2:14] = paste(cc[2:14], '_VERES',sep='')
  colnames(veres.829.1hz) = cc
  
  #cc = colnames(wisthaler.829.1hz)
  #cc[4:5] = paste(cc[4:5], '_WISTHALER',sep='')
  #colnames(wisthaler.829.1hz) = cc
  cc = colnames(jimenez.829.1hz)
  cc[2:8] = paste(cc[2:8], '_JIMENEZ',sep='')
  colnames(jimenez.829.1hz) = cc
  # --------#######--------- Get 5 or 10 Hz Data 8/29 ------#######---------
  # MET DATA - BUI + YANG + DLH
  met.829.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190829_R0_met.ict')
  #CO, CH4
  co.ch4.829.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-DACOM-5Hz_DC8_20190829_R1.ict')
  ind = which(co.ch4.829.5hz$CO_DACOM == max(co.ch4.829.5hz$CO_DACOM , na.rm=TRUE))
  maxco = (co.ch4.829.5hz$Time_Start[ind])
  # CO2
  co2.829.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CO2-7000-5Hz_DC8_20190829_R1.ict')
  ind = which(co2.829.5hz$CO2_7000_ppm == max(co2.829.5hz$CO2_7000_ppm , na.rm=TRUE))
  print(co2.829.5hz$Time_Start[ind] - maxco)
  # WARNEKE VOCs
  warneke.829.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-NOAAPTR-VOCs-5Hz_DC8_20190829_R3.ict')
  ind = which(warneke.829.5hz$Furan_NOAAPTR_ppbv== max(warneke.829.5hz$Furan_NOAAPTR_ppbv, na.rm=TRUE))
  print(warneke.829.5hz$Time_Start[ind] - maxco)
  # ISAF HCHO - merged to 5Hz from the online merge
  isaf.829.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190829_R0_ISAF.ict')
  # ROLLINS SO2 and NO
  rollinsno.829.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-NO-5Hz_DC8_20190829_R0.ict')
  ind = which(rollinsno.829.5hz$NO_LIF == max(rollinsno.829.5hz$NO_LIF, na.rm=TRUE))
  rollinsso2.829.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-SO2-5Hz_DC8_20190829_R1.ict')
  ind = which(rollinsso2.829.5hz$SO2_LIF == max(rollinsso2.829.5hz$SO2_LIF, na.rm=TRUE))
  # CIT VOCs - 
  cit.829.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190829_R0_CIT.ict')
  ind = which(cit.829.5hz$PHENOL.1Hz_CIT_WENNBERG == max(cit.829.5hz$PHENOL.1Hz_CIT_WENNBERG, na.rm=TRUE))
  # GTCIMS PANs - not sure how to match up peaks here
  gtcims.829.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190829_R0_huey.ict')
  # ----- Jimenez ---
  jimenez.829.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_AMS_20190829_R0_20230314T134101.ict')
  jimenez.829.5hz$OC_PM1_AMS_JIMENEZ = jimenez.829.5hz$OA_PM1_AMS_JIMENEZ/jimenez.829.5hz$OAtoOC_PM1_AMS
  
  # ------- append PI to colnames ----------
  cc = colnames(co2.829.5hz)
  cc[2] = paste(cc[2],'_DISKIN',sep='')
  colnames(co2.829.5hz) = cc 
  cc = colnames(co.ch4.829.5hz)
  cc[2:4] =  paste(cc[2:4],'_DISKIN',sep='')
  colnames(co.ch4.829.5hz) = cc
  cc=colnames(met.829.5hz)
  cc[2:36] = paste(cc[2:36],'_YANG',sep='')
  colnames(met.829.5hz) = cc
  cc = colnames(warneke.829.5hz)
  cc[2:43] = paste(cc[2:43],'_WARNEKE', sep='')
  colnames(warneke.829.5hz) = cc
  cc = colnames(rollinsno.829.5hz)
  cc[1] = "Time_Start" # really Time_Mid but need for merge
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsno.829.5hz) = cc
  cc = colnames(rollinsso2.829.5hz)
  cc[1] = "Time_Start" # really Time_Mid but need for merge
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsso2.829.5hz) = cc
  # make GTCIMS consistent with 1s merge
  cc=colnames(gtcims.829.5hz)
  cc[38:39] = c("APAN_GTCIMS_HUEY" , "PAN_GTCIMS_HUEY")
  colnames(gtcims.829.5hz)=cc
  # since ISAF, CIT, and GTCIMS came from merge tool, alread has PI's appended.
  
  # ------- get fuel moisture data --------
  if (doFM == 1){
    f1 = '/Users/ktravis1/Library/CloudStorage/Box-Box/FuelMoisture/fuel_moisture_content-20210715T1049Z/fmc_20190829_20Z.nc'
    fid = nc_open(f1)
    fuelMDead = ncvar_get(fid, varid = 'FMCG2D')
    fuelMLive = ncvar_get(fid, varid = 'FMCGLH2D')
    xlon = ncvar_get(fid, varid="XLONG_M")
    xlat = ncvar_get(fid, varid="XLAT_M")
    nc_close(fid)
  }  
  # ========================== HICKORY RIDGE  =====================

  fire = "Hickory Ridge"; fuel="grass"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 2;stop = stop - 2; startB = startO -1; stopB = startO+1}
    if (i == 2){start = startO + 3.2;stop = stop - 42.4; startB = 64415; stopB = 64417}# split from above
    if (i == 3){start = startO + 1.8;stop = stop - 1.8; startB = startO + 0; stopB = startO +2}
    if (i == 4){start = startO + 1;stop = stop - 2; startB = startO -2; stopB = startO +0}
    if (i == 5){start = startO + 13;stop = stop - 7.8; startB = startO + 0; stopB = startO +2}
    if (i == 6){start = startO + 25;stop = stop - 0; startB = startO + 0; stopB = startO +1} # bad pass
    if (i == 7){start = startO + 13;stop = stop - 1; startB = startO + 0; stopB = startO +1} # bad pass?
    if (i == 8){start = startO + 8.4;stop = stop - 170; startB = startO + 0; stopB = startO +1}
    tmp = time_align(start,stop,co2.829.5hz,  co.ch4.829.5hz, 
                     warneke.829.5hz,  isaf.829.5hz,
                     rollinsno.829.5hz, rollinsso2.829.5hz, 
                     cit.829.5hz, gtcims.829.5hz,moore.829fast, jimenez.829.5hz,met.829.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i==1){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i==2){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i==3){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i==4){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i==5){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i==6){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i==8){start=start-1;startB=startB-1;stopB=stopB-1}
    aug29th.fire= time_alignSLOWNOCANS(start,stop, co2.829.1hz, co.ch4.829.1hz, 
                                       warneke.829.1hz, isaf.829.1hz,rollinsno.829.1hz, rollinsso2.829.1hz, 
                                       cit.829.1hz,gtcims.829.1hz,ryerson.829.1hz , jimenez.829.1hz, 
                                       schwarz.829.1hz, freid.c2h6.829.1hz, freid.ch2o.829.1hz, 
                                       womack.829.1hz, stclair.829.1hz,veres.829.1hz, #wisthaler.829.1hz,
                                       moore.829, met.829.1hz)
    
    ind2 = which(aug29th.fire$Time_Start >= start & aug29th.fire$Time_Start  <= stop	) 
    ind2B = which(aug29th.fire$Time_Start >= startB & aug29th.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug29th.fire$fire = ''; aug29th.fire$fuel = ''
    tmp$fuel[ind] = 'grass' ; tmp$fire[ind] = 'HickoryRidge'
    aug29th.fire$fuel[ind2] = 'grass' ; aug29th.fire$fire[ind2] = 'HickoryRidge'
    # -------- Blake cans 
    indBLAKE = which(blake.829.merge$Time_Start >= (start-5) & blake.829.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.829.merge[indBLAKE,]
      blakeBG = blake.829.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.829.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.829.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.829.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.829.merge$Time_Start >= (start-0) & gilman.829.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.829.merge[indGILMAN,]
      GILMANBG = gilman.829.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.829.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.829.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.829.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.829.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.829.merge$Time_Start >= (start-0) & apel.829.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.829.merge[indAPEL,]
      APELBG = apel.829.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.829.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.829.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.829.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.829.merge[1,]
    }
    # Becky
    # Filter for >200 per Becky
    indBECKY = which(newTOGA.829$Time_Start >= (start-0) & newTOGA.829$Time_Start <= stop & newTOGA.829$Percent_plume > 0.5 &
                       newTOGA.829$CO_DACOM_DISKIN_BECKY > 200) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.829[indBECKY,]
      BECKYBG = newTOGA.829[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.829$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.829[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.829[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.829[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug29th.fire[ind2,])
      plotpass5hz(aug29th.fire[ind2,])
      plotpass1hz(aug29th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug29th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug29th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug29th.fire[ind2,],aug29th.fire[ind2B,],xspecies,'HickoryRidge', 'grass',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug29th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug29th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug29th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug29th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug29th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'HickoryRidge', 'grass',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug29th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug29th.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop  
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){HickoryRidge.1hz.EF = tmpEF1 ;HickoryRidge.5hz.EF = tmpEF5 }
    if (i > 1){HickoryRidge.1hz.EF = rbind(HickoryRidge.1hz.EF, tmpEF1) ; HickoryRidge.5hz.EF = rbind(HickoryRidge.5hz.EF, tmpEF5) }
    
  }
  HickoryRidge.1hz.EF = HickoryRidge.1hz.EF[order(HickoryRidge.1hz.EF$variable),]
  HickoryRidge.1hz.EF$transect_source_fire_ID = indA[1]
  HickoryRidge.5hz.EF = HickoryRidge.5hz.EF[order(HickoryRidge.5hz.EF$variable),]
  HickoryRidge.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== TALLGRASS  =====================
  fire="Tallgrass"; fuel="grass"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 6;stop = stop - 7; startB = startO - 1; stopB = startO +.6}
    if (i == 2){start = startO + 1.6;stop = stop - 8; startB = startO - 1; stopB = startO + 1}
    if (i == 3){start = startO + 1;stop = stop - 5; startB = startO - 1; stopB = startO + 1}
    if (i == 4){start = startO + 14;stop = stop - 4.4;  startB = startO + 0; stopB = startO +2}
    if (i == 5){start = startO + 1;stop = stop - 11.2; startB = startO -1; stopB = startO +1}
    if (i == 6){start = startO + 23;stop = stop - 28.2; startB = startO -1; stopB = startO +1}
    if (i == 7){start = startO + 204;stop = stop - 5; startB = startO + 0; stopB = startO +2}
    tmp = time_align(start,stop,co2.829.5hz,  co.ch4.829.5hz, 
                     warneke.829.5hz,  isaf.829.5hz,
                     rollinsno.829.5hz, rollinsso2.829.5hz, 
                     cit.829.5hz, gtcims.829.5hz,moore.829fast, jimenez.829.5hz,met.829.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i == 1){start=start-1; startB=startB-1; stopB=stopB-1}
    if (i == 2){start=start-1; startB=startB-1; stopB=stopB-1}
    if (i == 3){ start=start-1; startB=startB-1; stopB=stopB-1}
    if (i == 4){ start=start-1; startB=startB-1; stopB=stopB-1}
    if (i == 5){start=start-1; startB=startB-1;stopB=stopB-1}
    if (i == 6){start=start-1; startB=startB-1;stopB=stopB-1}
    if (i == 7){start=start-1; startB=startB-1;stopB=stopB-1}
    aug29th.fire= time_alignSLOWNOCANS(start,stop, co2.829.1hz, co.ch4.829.1hz, 
                                       warneke.829.1hz, isaf.829.1hz,rollinsno.829.1hz, rollinsso2.829.1hz, 
                                       cit.829.1hz,gtcims.829.1hz,ryerson.829.1hz , jimenez.829.1hz, 
                                       schwarz.829.1hz, freid.c2h6.829.1hz, freid.ch2o.829.1hz, 
                                       womack.829.1hz, stclair.829.1hz,veres.829.1hz, #wisthaler.829.1hz,
                                       moore.829, met.829.1hz)
    
    ind2 = which(aug29th.fire$Time_Start >= start & aug29th.fire$Time_Start  <= stop	) 
    ind2B = which(aug29th.fire$Time_Start >= startB & aug29th.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug29th.fire$fire = ''; aug29th.fire$fuel = ''
    tmp$fuel[ind] = 'grass' ; tmp$fire[ind] = 'Tallgrass'
    aug29th.fire$fuel[ind2] = 'grass' ; aug29th.fire$fire[ind2] = 'Tallgrass'
    # -------- Blake cans 
    indBLAKE = which(blake.829.merge$Time_Start >= (start-5) & blake.829.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.829.merge[indBLAKE,]
      blakeBG = blake.829.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.829.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.829.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.829.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.829.merge$Time_Start >= (start-0) & gilman.829.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.829.merge[indGILMAN,]
      GILMANBG = gilman.829.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.829.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.829.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.829.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.829.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.829.merge$Time_Start >= (start-0) & apel.829.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.829.merge[indAPEL,]
      APELBG = apel.829.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.829.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.829.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.829.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.829.merge[1,]
    }
    # Becky
    # Exclude <0.3ppm per Becky
    indBECKY = which(newTOGA.829$Time_Start >= (start-0) & newTOGA.829$Time_Start <= stop & 
                       newTOGA.829$Percent_plume > 0.5 & newTOGA.829$CO_DACOM_DISKIN > 300) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.829[indBECKY,]
      BECKYBG = newTOGA.829[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.829$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.829[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.829[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.829[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug29th.fire[ind2,])
      plotpass5hz(aug29th.fire[ind2,])
      plotpass1hz(aug29th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      plotslowstuff(tmp,startO,stopO,"Tallgrass",i,blake,blakeBG,GILMAN, GILMANBG, BECKY, BECKYBG, doBECKY)
      # plot background?
      plotpass5hzJUSTCO(aug29th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug29th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug29th.fire[ind2,],aug29th.fire[ind2B,],xspecies,'Tallgrass', 'grass',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug29th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug29th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug29th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug29th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug29th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Tallgrass', 'grass',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug29th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug29th.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop 
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Tallgrass.1hz.EF = tmpEF1 ;Tallgrass.5hz.EF = tmpEF5 }
    if (i > 1){Tallgrass.1hz.EF = rbind(Tallgrass.1hz.EF, tmpEF1) ; Tallgrass.5hz.EF = rbind(Tallgrass.5hz.EF, tmpEF5) }
  }
  Tallgrass.1hz.EF = Tallgrass.1hz.EF[order(Tallgrass.1hz.EF$variable),]
  Tallgrass.1hz.EF$transect_source_fire_ID = indA[1]
  Tallgrass.5hz.EF = Tallgrass.5hz.EF[order(Tallgrass.5hz.EF$variable),]
  Tallgrass.5hz.EF$transect_source_fire_ID = indA[1]
  # ========================== AKITO  - unidentified fuel/low CO =====================
  fire="Akito"; fuel="?"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 6;stop = stop - 7; startB = startO + 0; stopB = startO+1}
    if (i == 2){start = startO + 2;stop = stop - 1; startB = startO + 0; stopB = startO+2}
    if (i == 3){start = startO + 4;stop = stop - 7; startB = startO + 0; stopB = startO+2}
    tmp = time_align(start,stop,co2.829.5hz,  co.ch4.829.5hz, 
                     warneke.829.5hz,  isaf.829.5hz,
                     rollinsno.829.5hz, rollinsso2.829.5hz, 
                     cit.829.5hz, gtcims.829.5hz,moore.829fast, jimenez.829.5hz,met.829.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    aug29th.fire= time_alignSLOWNOCANS(start,stop, co2.829.1hz, co.ch4.829.1hz, 
                                       warneke.829.1hz, isaf.829.1hz,rollinsno.829.1hz, rollinsso2.829.1hz, 
                                       cit.829.1hz,gtcims.829.1hz,ryerson.829.1hz , jimenez.829.1hz, 
                                       schwarz.829.1hz, freid.c2h6.829.1hz, freid.ch2o.829.1hz, 
                                       womack.829.1hz, stclair.829.1hz,veres.829.1hz, #wisthaler.829.1hz,
                                       moore.829, met.829.1hz)
    ind2 = which(aug29th.fire$Time_Start >= start & aug29th.fire$Time_Start  <= stop	) 
    ind2B = which(aug29th.fire$Time_Start >= startB & aug29th.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug29th.fire$fire = ''; aug29th.fire$fuel = ''
    tmp$fuel[ind] = '?' ; tmp$fire[ind] = 'Akito'
    aug29th.fire$fuel[ind2] = '?' ; aug29th.fire$fire[ind2] = 'Akito'
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug29th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # -------- Blake cans 
    indBLAKE = which(blake.829.merge$Time_Start >= (start-5) & blake.829.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.829.merge[indBLAKE,]
      blakeBG = blake.829.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.829.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.829.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.829.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.829.merge$Time_Start >= (start-0) & gilman.829.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.829.merge[indGILMAN,]
      GILMANBG = gilman.829.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.829.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.829.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.829.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.829.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.829.merge$Time_Start >= (start-0) & apel.829.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.829.merge[indAPEL,]
      APELBG = apel.829.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.829.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.829.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.829.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.829.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.829$Time_Start >= (start-0) & newTOGA.829$Time_Start <= stop & newTOGA.829$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.829[indBECKY,]
      BECKYBG = newTOGA.829[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.829$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.829[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.829[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.829[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug29th.fire[ind2,])
      plotpass5hz(aug29th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug29th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug29th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug29th.fire[ind2,],aug29th.fire[ind2B,],xspecies,'Akito', '?',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug29th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug29th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug29th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug29th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug29th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Akito', '?',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug29th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug29th.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop    
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Akito.1hz.EF = tmpEF1 ;Akito.5hz.EF = tmpEF5 }
    if (i > 1){Akito.1hz.EF = rbind(Akito.1hz.EF, tmpEF1) ; Akito.5hz.EF = rbind(Akito.5hz.EF, tmpEF5) }
  }
  Akito.1hz.EF = Akito.1hz.EF[order(Akito.1hz.EF$variable),]
  Akito.1hz.EF$transect_source_fire_ID = indA[1]
  Akito.5hz.EF = Akito.5hz.EF[order(Akito.5hz.EF$variable),]
  Akito.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== BOXER  =====================
  fire="Boxer"; fuel="rice"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 3.4;stop = stop - 4.8; startB = startO + 0; stopB = startO+1}
    
    tmp = time_align(start,stop,co2.829.5hz,  co.ch4.829.5hz, 
                     warneke.829.5hz,  isaf.829.5hz,
                     rollinsno.829.5hz, rollinsso2.829.5hz, 
                     cit.829.5hz, gtcims.829.5hz,moore.829fast, jimenez.829.5hz,met.829.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-2;stop=stop+1;startB=startB-1;stopB=stopB-1}
    
    aug29th.fire= time_alignSLOWNOCANS(start,stop, co2.829.1hz, co.ch4.829.1hz, 
                                       warneke.829.1hz, isaf.829.1hz,rollinsno.829.1hz, rollinsso2.829.1hz, 
                                       cit.829.1hz,gtcims.829.1hz,ryerson.829.1hz , jimenez.829.1hz, 
                                       schwarz.829.1hz, freid.c2h6.829.1hz, freid.ch2o.829.1hz, 
                                       womack.829.1hz, stclair.829.1hz,veres.829.1hz, #wisthaler.829.1hz,
                                       moore.829, met.829.1hz)
    ind2 = which(aug29th.fire$Time_Start >= start & aug29th.fire$Time_Start  <= stop	) 
    ind2B = which(aug29th.fire$Time_Start >= startB & aug29th.fire$Time_Start <= stopB	) 
    
    if (doplot == 1){
      plotpass5hzJUSTCO(aug29th.fire[ind2,])
      plotpass5hz(aug29th.fire[ind2,])
      plotpass1hz(aug29th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug29th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug29th.fire$fire = ''; aug29th.fire$fuel = ''
    tmp$fuel[ind] = 'rice' ; tmp$fire[ind] = 'Boxer'
    aug29th.fire$fuel[ind2] = 'rice' ; aug29th.fire$fire[ind2] = 'Boxer'
    
    # -------- Blake cans 
    indBLAKE = which(blake.829.merge$Time_Start >= (start-5) & blake.829.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.829.merge[indBLAKE,]
      blakeBG = blake.829.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.829.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.829.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.829.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }

    # Gilman cans
    indGILMAN = which(gilman.829.merge$Time_Start >= (start-0) & gilman.829.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.829.merge[indGILMAN,]
      GILMANBG = gilman.829.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.829.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.829.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.829.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.829.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.829.merge$Time_Start >= (start-0) & apel.829.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.829.merge[indAPEL,]
      APELBG = apel.829.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.829.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.829.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.829.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.829.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.829$Time_Start >= (start-0) & newTOGA.829$Time_Start <= stop & newTOGA.829$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.829[indBECKY,]
      BECKYBG = newTOGA.829[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.829$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.829[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.829[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.829[1,]
    }
    
   
    if (doplot == 1){
      plotslowstuff(tmp[ind,], startO, stopO, 'Boxer',i, blake, blakeBG, GILMAN, GILMANBG, BECKY, BECKYBG)
      plotpass5hzJUSTCO(aug29th.fire[ind2,])
      plotpass5hz(aug29th.fire[ind2,])
      plotpass1hz(aug29th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug29th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug29th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug29th.fire[ind2,],aug29th.fire[ind2B,],xspecies,'Boxer', 'rice',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug29th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug29th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug29th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug29th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug29th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Boxer', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug29th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug29th.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop    
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Boxer.1hz.EF = tmpEF1 ;Boxer.5hz.EF = tmpEF5 }
    if (i > 1){Boxer.1hz.EF = rbind(Boxer.1hz.EF, tmpEF1) ; Boxer.5hz.EF = rbind(Boxer.5hz.EF, tmpEF5) }
    
  }
  Boxer.1hz.EF = Boxer.1hz.EF[order(Boxer.1hz.EF$variable),]
  Boxer.1hz.EF$transect_source_fire_ID = indA[1]
  Boxer.5hz.EF = Boxer.5hz.EF[order(Boxer.5hz.EF$variable),]
  Boxer.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== CHESSIE  =====================
  fire="Chessie"; fuel="winter wheat"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    tmp = time_align(start,stop,co2.829.5hz,  co.ch4.829.5hz, 
                     warneke.829.5hz,  isaf.829.5hz,
                     rollinsno.829.5hz, rollinsso2.829.5hz, 
                     cit.829.5hz, gtcims.829.5hz,moore.829fast, jimenez.829.5hz,met.829.5hz)
    if (i == 1){start = startO + 1.2;stop = stop - 2.5; startB = startO -1; stopB = startO+1}
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-2; startB=startB-2;stopB=stopB-2}
    aug29th.fire= time_alignSLOWNOCANS(start,stop, co2.829.1hz, co.ch4.829.1hz, 
                                       warneke.829.1hz, isaf.829.1hz,rollinsno.829.1hz, rollinsso2.829.1hz, 
                                       cit.829.1hz,gtcims.829.1hz,ryerson.829.1hz , jimenez.829.1hz, 
                                       schwarz.829.1hz, freid.c2h6.829.1hz, freid.ch2o.829.1hz, 
                                       womack.829.1hz, stclair.829.1hz,veres.829.1hz, #wisthaler.829.1hz,
                                       moore.829, met.829.1hz)
    ind2 = which(aug29th.fire$Time_Start >= start & aug29th.fire$Time_Start  <= stop	) 
    ind2B = which(aug29th.fire$Time_Start >= startB & aug29th.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug29th.fire$fire = ''; aug29th.fire$fuel = ''
    tmp$fuel[ind] = 'winter wheat' ; tmp$fire[ind] = 'Chessie'
    aug29th.fire$fuel[ind2] = 'winter wheat' ; aug29th.fire$fire[ind2] = 'Chessie'
    # -------- Blake cans 
    indBLAKE = which(blake.829.merge$Time_Start >= (start-5) & blake.829.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.829.merge[indBLAKE,]
      blakeBG = blake.829.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.829.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.829.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.829.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.829.merge$Time_Start >= (start-0) & gilman.829.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.829.merge[indGILMAN,]
      GILMANBG = gilman.829.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.829.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.829.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.829.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.829.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.829.merge$Time_Start >= (start-0) & apel.829.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.829.merge[indAPEL,]
      APELBG = apel.829.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.829.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.829.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.829.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.829.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.829$Time_Start >= (start-0) & newTOGA.829$Time_Start <= stop & newTOGA.829$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.829[indBECKY,]
      BECKYBG = newTOGA.829[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.829$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.829[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.829[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.829[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug29th.fire[ind2,])
      plotpass5hz(aug29th.fire[ind2,])
      plotpass1hz(aug29th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug29th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug29th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug29th.fire[ind2,],aug29th.fire[ind2B,],xspecies,'Chessie', 'winter wheat',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug29th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug29th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug29th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug29th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug29th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Chessie', 'winter wheat',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug29th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug29th.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop 
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Chessie.1hz.EF = tmpEF1 ;Chessie.5hz.EF = tmpEF5 }
    if (i > 1){Chessie.1hz.EF = rbind(Chessie.1hz.EF, tmpEF1) ; Chessie.5hz.EF = rbind(Chessie.5hz.EF, tmpEF5) }
    
  }
  Chessie.1hz.EF = Chessie.1hz.EF[order(Chessie.1hz.EF$variable),]
  Chessie.1hz.EF$transect_source_fire_ID = indA[1]
  Chessie.5hz.EF = Chessie.5hz.EF[order(Chessie.5hz.EF$variable),]
  Chessie.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== DINGO  =====================
  fire="Dingo";fuel="grass"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ 
    print(c(fire,i) )
    
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 0;stop = stopO - 4; startB = startO - 2.2; stopB = startO+.2}
    if (i == 2){start = startO + 1;stop = stopO - 16.4; startB =  77965.8; stopB = 77968.2} # split from above
    if (i == 3){start = startO + 3.2;stop = stopO - 195; startB =  77965.8; stopB =77968.2} # split from above
    tmp = time_align(start,stop,co2.829.5hz,  co.ch4.829.5hz, 
                     warneke.829.5hz,  isaf.829.5hz,
                     rollinsno.829.5hz, rollinsso2.829.5hz, 
                     cit.829.5hz, gtcims.829.5hz,moore.829fast, jimenez.829.5hz,met.829.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1; startB=startB-1;stopB=stopB-1}
    if (i ==2){start=start-1; stop=stop+1; startB=startB-1;stopB=stopB-1}
    if (i ==3){start=start-1; startB=startB-1;stopB=stopB-1}
    aug29th.fire= time_alignSLOWNOCANS(start,stop, co2.829.1hz, co.ch4.829.1hz, 
                                       warneke.829.1hz, isaf.829.1hz,rollinsno.829.1hz, rollinsso2.829.1hz, 
                                       cit.829.1hz,gtcims.829.1hz,ryerson.829.1hz , jimenez.829.1hz, 
                                       schwarz.829.1hz, freid.c2h6.829.1hz, freid.ch2o.829.1hz, 
                                       womack.829.1hz, stclair.829.1hz,veres.829.1hz, #wisthaler.829.1hz,
                                       moore.829, met.829.1hz)
    ind2 = which(aug29th.fire$Time_Start >= start & aug29th.fire$Time_Start  <= stop	) 
    ind2B = which(aug29th.fire$Time_Start >= startB & aug29th.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug29th.fire$fire = ''; aug29th.fire$fuel = ''
    tmp$fuel[ind] = 'grass' ; tmp$fire[ind] = 'Dingo'
    aug29th.fire$fuel[ind2] = 'grass' ; aug29th.fire$fire[ind2] = 'Dingo'
    # -------- Blake cans 
    indBLAKE = which(blake.829.merge$Time_Start >= (start-5) & blake.829.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.829.merge[indBLAKE,]
      blakeBG = blake.829.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.829.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.829.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.829.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.829.merge$Time_Start >= (start-0) & gilman.829.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.829.merge[indGILMAN,]
      GILMANBG = gilman.829.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.829.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.829.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.829.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.829.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.829.merge$Time_Start >= (start-0) & apel.829.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.829.merge[indAPEL,]
      APELBG = apel.829.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.829.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.829.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.829.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.829.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.829$Time_Start >= (start-0) & newTOGA.829$Time_Start <= stop & newTOGA.829$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.829[indBECKY,]
      BECKYBG = newTOGA.829[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.829$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.829[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.829[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.829[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug29th.fire[ind2,])
      plotpass5hz(aug29th.fire[ind2,])
      plotpass1hz(aug29th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug29th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug29th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug29th.fire[ind2,],aug29th.fire[ind2B,],xspecies,'Dingo', 'grass',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug29th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug29th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug29th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug29th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug29th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Dingo', 'grass',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug29th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug29th.fire$transect_type[ind2], na.rm=TRUE)
  
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop   
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop    
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Dingo.1hz.EF = tmpEF1 ;Dingo.5hz.EF = tmpEF5 }
    if (i > 1){Dingo.1hz.EF = rbind(Dingo.1hz.EF, tmpEF1) ; Dingo.5hz.EF = rbind(Dingo.5hz.EF, tmpEF5) }
  }
  Dingo.1hz.EF = Dingo.1hz.EF[order(Dingo.1hz.EF$variable),]
  Dingo.1hz.EF$transect_source_fire_ID = indA[1]
  Dingo.5hz.EF = Dingo.5hz.EF[order(Dingo.5hz.EF$variable),]
  Dingo.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== ELKHOUND  =====================
  fire="Elkhound";fuel="grass"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    # several consective plumes but they seem to behave about the same
    if (i == 1){start = startO + 2;stop = stopO - 338; startB = startO + 0; stopB = startO +2}
    
    tmp = time_align(start,stop,co2.829.5hz,  co.ch4.829.5hz, 
                     warneke.829.5hz,  isaf.829.5hz,
                     rollinsno.829.5hz, rollinsso2.829.5hz, 
                     cit.829.5hz, gtcims.829.5hz,moore.829fast, jimenez.829.5hz,met.829.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    aug29th.fire= time_alignSLOWNOCANS(start,stop, co2.829.1hz, co.ch4.829.1hz, 
                                       warneke.829.1hz, isaf.829.1hz,rollinsno.829.1hz, rollinsso2.829.1hz, 
                                       cit.829.1hz,gtcims.829.1hz,ryerson.829.1hz , jimenez.829.1hz, 
                                       schwarz.829.1hz, freid.c2h6.829.1hz, freid.ch2o.829.1hz, 
                                       womack.829.1hz, stclair.829.1hz,veres.829.1hz, #wisthaler.829.1hz,
                                       moore.829, met.829.1hz)
    ind2 = which(aug29th.fire$Time_Start >= start & aug29th.fire$Time_Start  <= stop	) 
    ind2B = which(aug29th.fire$Time_Start >= startB & aug29th.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug29th.fire$fire = ''; aug29th.fire$fuel = ''
    tmp$fuel[ind] = 'grass' ; tmp$fire[ind] = 'Elkhound'
    aug29th.fire$fuel[ind2] = 'grass' ; aug29th.fire$fire[ind2] = 'Elkhound'
    # -------- Blake cans 
    indBLAKE = which(blake.829.merge$Time_Start >= (start-5) & blake.829.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.829.merge[indBLAKE,]
      blakeBG = blake.829.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.829.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.829.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.829.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.829.merge$Time_Start >= (start-0) & gilman.829.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      #kludge
      if (length(indGILMAN) > 1 & i == 1){indGILMAN = indGILMAN[1]}
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.829.merge[indGILMAN,]
      GILMANBG = gilman.829.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.829.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.829.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.829.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.829.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.829.merge$Time_Start >= (start-0) & apel.829.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.829.merge[indAPEL,]
      APELBG = apel.829.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.829.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.829.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.829.merge[1,]; APEL$CO_DACOM_DISKIN_APEL = NaN
      APELBG = apel.829.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.829$Time_Start >= (start-0) & newTOGA.829$Time_Start <= stop & newTOGA.829$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.829[indBECKY,]
      BECKYBG = newTOGA.829[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.829$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.829[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.829[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.829[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug29th.fire[ind2,])
      plotpass5hz(aug29th.fire[ind2,])
      plotpass1hz(aug29th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug29th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug29th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug29th.fire[ind2,],aug29th.fire[ind2B,],xspecies,'Elkhound', 'grass',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug29th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug29th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug29th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug29th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug29th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Elkhound', 'grass',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug29th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug29th.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop 
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Elkhound.1hz.EF = tmpEF1 ;Elkhound.5hz.EF = tmpEF5 }
    if (i > 1){Elkhound.1hz.EF = rbind(Elkhound.1hz.EF, tmpEF1) ; Elkhound.5hz.EF = rbind(Elkhound.5hz.EF, tmpEF5) }
    
  }
  Elkhound.1hz.EF = Elkhound.1hz.EF[order(Elkhound.1hz.EF$variable),]
  Elkhound.1hz.EF$transect_source_fire_ID = indA[1]
  Elkhound.5hz.EF = Elkhound.5hz.EF[order(Elkhound.5hz.EF$variable),]
  Elkhound.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ----- [[[[[[[[[[[[[[[[[[[[ Aug 30rd ]]]]]]]]]]]]]]]]]]]]] -------
  # --------#######--------- Get 1 Hz Data individual 8/30------#######---------
  # plume tags
  tags = getICARTTdataSIMPLE('InputFiles/firexaq-fire-Flags-1HZ_DC8_20190830_R9.ict') ;tags$Time_Start = tags$TIME_START
  # MET DATA
  met.830.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-MetNav_DC8_20190830_R1.ict')
  met.830.1hz = merge(met.830.1hz, tags, by='Time_Start')
  # CO2
  co2.830.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CO2-7000_DC8_20190830_R2.ict')
  # -------- DISKIN -----CO, CH4
  co.ch4.830.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-DACOM_DC8_20190830_R1.ict')
  # --------- WARNEKE ----  VOCs
  warneke.830.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-NOAAPTR-VOCs-1Hz_DC8_20190830_R3.ict')
  # ------ HANISCO - ISAF HCHO - merged to 5Hz from the online merge
  isaf.830.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-ISAF-CH2O-1Hz_DC8_20190830_R0.ict')
  #  ------- ROLLINS - SO2 and NO
  rollinsno.830.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-NO_DC8_20190830_R1.ict')
  rollinsno.830.1hz$Time_Start = rollinsno.830.1hz$time_mid
  rollinsso2.830.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-SO2_DC8_20190830_R1.ict')
  rollinsso2.830.1hz$Time_Start = rollinsso2.830.1hz$time_mid
  
  #  ----- WENNBERG - CIT VOCs - 
  cit.830.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190830_R0_CIT.ict')
  # ------ HUEY - GTCIMS PANs - not sure how to match up peaks here
  gtcims.830.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190830_R0_Huey.ict')
  # ------ RYERSON
  ryerson.A = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NO_DC8_20190830_R1.ict')
  ryerson.B = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NO2_DC8_20190830_R1.ict')
  ryerson.C = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NOy_DC8_20190830_R1.ict')
  ryerson.D = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-O3_DC8_20190830_R1.ict')
  ryerson.830.1hz = cbind(ryerson.A,ryerson.B,ryerson.C,ryerson.D) 
  ryerson.830.1hz$Time_Start = ryerson.830.1hz$Time_start
  # ----- JIMENEZ ---
  jimenez.830.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-EESI_DC8_20190830_R1.ict')
  
  # ----- SCHWARZ ---
  schwarz.830.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-SP2-BC-1HZ_DC8_20190830_R2.ict')
  # ----- FREID ---
  freid.c2h6.830.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-C2H6_DC8_20190830_R3.ict')
  freid.ch2o.830.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CH2O_DC8_20190830_R3.ict')
  # ------ WOMACK ---
  womackA = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-CH3COCHO_DC8_20190830_R1.ict')
  womackB = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-CHOCHO_DC8_20190830_R1.ict')
  womackC = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-HNO2_DC8_20190830_R1.ict')
  womackD = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-NO2_DC8_20190830_R1.ict')
  womack.830.1hz = cbind(womackA, womackB, womackC, womackD)
  # -------St Clair
  stclair.830.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-CANOE-NO2_DC8_20190830_R0.ict')
  stclair.830.1hz$Time_Start = stclair.830.1hz$Time_start
  
  # ------- VERES
  veres.A = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-ClNO2_DC8_20190830_R0.ict')
  veres.B = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HCOOH_DC8_20190830_R1.ict')
  veres.C = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HNO2_DC8_20190830_R1.ict')
  veres.D = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-N2O5_DC8_20190830_R0.ict')
  veres.E = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HPMTF_DC8_20190830_R0.ict')
  veres.F = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-CH3COOCl_DC8_20190830_R0.ict')
  veres.G = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-Cl2_DC8_20190830_R0.ict')
  veres.H = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrCl_DC8_20190830_R0.ict')
  veres.I = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrCN_DC8_20190830_R0.ict')
  veres.J = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrO_DC8_20190830_R0.ict')
  veres.K = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HCN_DC8_20190830_R0.ict')
  veres.L = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HNCO_DC8_20190830_R0.ict')
  veres.830.1hz = cbind(veres.A,veres.B,veres.C,veres.D,veres.E,veres.F,veres.G,veres.H,veres.I,veres.J,veres.K,veres.L)
  
  # --- WISTHALER
  #wisthaler.830.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190830_R0_Wisthaler.ict')

  # ---- BLAKE
  blake.830.1hz = getICARTTdataSIMPLE('InputFiles/WAS-MERGE/firexaq-mrgWAS-dc8_merge_20190830_R1.ict')
  cc = colnames(blake.830.1hz)
  blake.830.merge = blake.830.1hz[,c(1,2,96:225)]
  blake.830.merge$CO_DACOM_DISKIN_BLAKE = blake.830.1hz$CO_DACOM_DISKIN
  blake.830.merge$CO2_7000_ppm_DISKIN_BLAKE = blake.830.1hz$CO2_7000_ppm_DISKIN
  
  # ------ APEL
  apel.830.1hz = getICARTTdataSIMPLE('InputFiles/TOGA-MERGE/firexaq-mrgTOGA-dc8_merge_20190830_R1.ict')
  cc = colnames(apel.830.1hz)
  apel.830.merge = apel.830.1hz[,c(1,2,226:315)]
  apel.830.merge$CO_DACOM_DISKIN_APEL = apel.830.1hz$CO_DACOM_DISKIN
  apel.830.merge$CO2_7000_ppm_DISKIN_APEL =apel.830.1hz$CO2_7000_ppm_DISKIN
  
  # Becky's better merge
  file = 'InputFiles/Hornbrook/FIREX-AQ weighted TOGA merge 2022-01-24_0830.xlsx'
  newTOGA.830 = readxl::read_xlsx(file); newTOGA.830[newTOGA.830==-999] = NaN; newTOGA.830[newTOGA.830==-888] = NaN
  newTOGA.830$CO_DACOM_DISKIN_BECKY = newTOGA.830$CO_DACOM_DISKIN
  newTOGA.830$CO2_7000_ppm_DISKIN_BECKY = NaN
  newTOGA.830$Time_Start=newTOGA.830$Time_Start...4
  
  # ----GILMAN
  gilman.830.1hz = getICARTTdataSIMPLE('InputFiles/iWAS-MERGE/firexaq-mrgiWAS-dc8_merge_20190830_R1.ict')
  cc = colnames(gilman.830.1hz)
  gilman.830.merge = gilman.830.1hz[,c(1,2,316:361)]
  gilman.830.merge$CO_DACOM_DISKIN_GILMAN = gilman.830.1hz$CO_DACOM_DISKIN
  gilman.830.merge$CO2_7000_ppm_DISKIN_GILMAN = gilman.830.1hz$CO2_7000_ppm_DISKIN
  
  # ------ Moore 
  moore.830fast = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190830_R0_MOORE.ict')
  
  moore.830p1 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-AerosolCloudConc_DC8_20190830_R0.ict')
  moore.830p2 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-LAScold_DC8_20190830_R0.ict')
  moore.830p3 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-LAShot_DC8_20190830_R0.ict')
  moore.830p4 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-CPSPD_DC8_20190830_R0.ict')
  moore.830p5 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-CDP_DC8_20190830_R0.ict')
  moore.830p6 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-SMPS_DC8_20190830_R0.ict')
  moore.830 =merge(moore.830p1, moore.830p2, by='Time_mid', all = TRUE, incomparables = NA)
  moore.830 =merge(moore.830, moore.830p3, by='Time_mid', all = TRUE, incomparables = NA)
  moore.830 =merge(moore.830, moore.830p4, by='Time_mid', all = TRUE, incomparables = NA)
  moore.830 =merge(moore.830, moore.830p5, by='Time_mid', all = TRUE, incomparables = NA)
  moore.830 =merge(moore.830, moore.830p6, by='Time_mid', all = TRUE, incomparables = NA)
  
  # ------- append PI to colnames 1hz ----------
  cc = colnames(co2.830.1hz)
  cc[2] = paste(cc[2],'_DISKIN',sep='')
  colnames(co2.830.1hz) = cc 
  cc = colnames(co.ch4.830.1hz)
  cc[2:4] =  paste(cc[2:4],'_DISKIN',sep='')
  colnames(co.ch4.830.1hz) = cc
  cc=colnames(met.830.1hz)
  
  colnames(isaf.830.1hz) = c("Time_Start"  , "CH2O_ISAF_HANISCO" ,"CH2O_ISAF_precision_HANISCO")
  cc[2:36] = paste(cc[2:36],'_YANG',sep='')
  colnames(met.830.1hz) = cc
  cc = colnames(warneke.830.1hz)
  cc[2:43] = paste(cc[2:43],'_WARNEKE', sep='')
  colnames(warneke.830.1hz) = cc
  cc = colnames(rollinsno.830.1hz)
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsno.830.1hz) = cc
  cc = colnames(rollinsso2.830.1hz)
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsso2.830.1hz) = cc
  # make GTCIMS consistent with 1s merge
  cc=colnames(gtcims.830.1hz)
  cc[4:7] = c("PAN_GTCIMS_HUEY" , "PPN_GTCIMS_HUEY"  ,"APAN_GTCIMS_HUEY" ,"PBN_GTCIMS_HUEY")
  colnames(gtcims.830.1hz)=cc
  # since ISAF, CIT, and GTCIMS came from merge tool, alread has PI's appended.
  ryerson.830.1hz <- ryerson.830.1hz[, !duplicated(colnames(ryerson.830.1hz))]
  
  cc = colnames(ryerson.830.1hz)
  cc[2:9] = paste(cc[2:9],'_RYERSON',sep='')
  colnames(ryerson.830.1hz)=cc
  
  cc = colnames(schwarz.830.1hz)
  cc[2:3] = paste(cc[2:3],'_SCHWARZ', sep='')
  colnames(schwarz.830.1hz) =cc
  
  cc = colnames(freid.c2h6.830.1hz)
  cc[4:7] = paste(cc[4:7], '_FRIED',sep='')
  colnames(freid.c2h6.830.1hz) = cc
  
  cc = colnames(freid.ch2o.830.1hz)
  cc[4:7] = paste(cc[4:7], '_FRIED',sep='')
  colnames(freid.ch2o.830.1hz) = cc
  
  colnames(freid.ch2o.830.1hz) = cc
  womack.830.1hz <- womack.830.1hz[, !duplicated(colnames(womack.830.1hz))]
  cc = colnames(womack.830.1hz)
  cc[2:5] = paste(cc[2:5], '_WOMACK',sep='')
  colnames(womack.830.1hz) = cc

  veres.830.1hz <- veres.830.1hz[, !duplicated(colnames(veres.830.1hz))]
  cc = colnames(veres.830.1hz)
  cc[2:14] = paste(cc[2:14], '_VERES',sep='')
  colnames(veres.830.1hz) = cc
  
  #cc = colnames(wisthaler.830.1hz)
  #cc[4:5] = paste(cc[4:5], '_WISTHALER',sep='')
  #colnames(wisthaler.830.1hz) = cc
  cc = colnames(jimenez.830.1hz)
  cc[2:8] = paste(cc[2:8], '_JIMENEZ',sep='')
  colnames(jimenez.830.1hz) = cc
  # --------#######--------- Get 5 or 10 Hz Data 8/30 ------#######---------
  # MET DATA - BUI + YANG + DLH
  met.830.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190830_R0_met.ict')
  
  #CO, CH4
  co.ch4.830.5hz = getICARTTdataSIMPLE('InputFiles//FIREXAQ-DACOM-5Hz_DC8_20190830_R1.ict')
  # CO2
  co2.830.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CO2-7000-5Hz_DC8_20190830_R1.ict')
  # WARNEKE VOCs
  warneke.830.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-NOAAPTR-VOCs-5Hz_DC8_20190830_R3.ict')
  # ISAF HCHO - merged to 5Hz from the online merge
  isaf.830.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190830_R0_ISAF.ict')
  # ROLLINS SO2 and NO
  rollinsno.830.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-NO-5Hz_DC8_20190830_R0.ict')
  rollinsso2.830.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-SO2-5Hz_DC8_20190830_R1.ict')
  # CIT VOCs - 
  cit.830.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190830_R0_CIT.ict')
  # GTCIMS PANs - not sure how to match up peaks here
  gtcims.830.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190830_R0_huey.ict')
  # ----- Jimenez ---
  jimenez.830.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_AMS_20190830_R0_20230314T134114.ict')
  jimenez.830.5hz$OC_PM1_AMS_JIMENEZ = jimenez.830.5hz$OA_PM1_AMS_JIMENEZ/jimenez.830.5hz$OAtoOC_PM1_AMS
  
  # ------- append PI to colnames ----------
  cc = colnames(co2.830.5hz)
  cc[2] = paste(cc[2],'_DISKIN',sep='')
  colnames(co2.830.5hz) = cc 
  cc = colnames(co.ch4.830.5hz)
  cc[2:4] =  paste(cc[2:4],'_DISKIN',sep='')
  colnames(co.ch4.830.5hz) = cc
  cc=colnames(met.830.5hz)
  cc[2:36] = paste(cc[2:36],'_YANG',sep='')
  colnames(met.830.5hz) = cc
  cc = colnames(warneke.830.5hz)
  cc[2:43] = paste(cc[2:43],'_WARNEKE', sep='')
  colnames(warneke.830.5hz) = cc
  cc = colnames(rollinsno.830.5hz)
  cc[1] = "Time_Start" # really Time_Mid but need for merge
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsno.830.5hz) = cc
  cc = colnames(rollinsso2.830.5hz)
  cc[1] = "Time_Start" # really Time_Mid but need for merge
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsso2.830.5hz) = cc
  # make GTCIMS consistent with 1s merge
  cc=colnames(gtcims.830.5hz)
  cc[38:39] = c("APAN_GTCIMS_HUEY" , "PAN_GTCIMS_HUEY")
  colnames(gtcims.830.5hz)=cc
  # since ISAF, CIT, and GTCIMS came from merge tool, alread has PI's appended.
  
  # ------- get fuel moisture data --------
  if (doFM == 1){
    f1 = '/Users/ktravis1/Library/CloudStorage/Box-Box/FuelMoisture/fuel_moisture_content-20210715T1049Z/fmc_20190830_20Z.nc'
    fid = nc_open(f1)
    fuelMDead = ncvar_get(fid, varid = 'FMCG2D')
    fuelMLive = ncvar_get(fid, varid = 'FMCGLH2D')
    xlon = ncvar_get(fid, varid="XLONG_M")
    xlat = ncvar_get(fid, varid="XLAT_M")
    nc_close(fid)
  }  
  # ========================== BLACKWATER RIVER  =====================
  
  types = c()
  
  # use ethane for Blake and Gilman
  #plot(aug30th.fire1$Ethane_WAS_BLAKE, aug30th.fire1$C2H6_CAMS_pptv_FRIED)
  #cor.test(aug30th.fire1$Ethane_WAS_BLAKE, aug30th.fire1$C2H6_CAMS_pptv_FRIED)
  #plot(aug30th.fire1$Ethane_NOAAiWAS_GILMAN, aug30th.fire1$C2H6_CAMS_pptv_FRIED)
  #cor.test(aug30th.fire1$Ethane_NOAAiWAS_GILMAN, aug30th.fire1$C2H6_CAMS_pptv_FRIED)
  # use benzene for Apel
  #plot(aug30th.fire1$Benzene_TOGA_APEL, aug30th.fire1$Benzene_NOAAPTR_ppbv_WARNEKE)
  #cor.test(aug30th.fire1$Benzene_TOGA_APEL, aug30th.fire1$Benzene_NOAAPTR_ppbv_WARNEKE)
  
  coBW = c();  lonBW=c();  latBW=c()
  fire="Blackwater River";fuel='forest'
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    types=c(types,flags$transect_type[indA[i]])
    if (i == 1){start = startO + 2;stop = stopO - 19; startB = startO -1; stopB = startO +1}
    if (i == 2){start = startO + 2;stop = stopO - 8.2; startB = 59821; stopB = 59823} # split from above
    if (i == 3){start = startO + 2;stop = stopO - 5; startB = startO + 0; stopB = startO +2}
    if (i == 4){start = startO + 1.6;stop = stopO - 5.8; startB = startO + 0; stopB = startO +1}
    #if (i == 5){start = startO + 0.2;stop = stopO - 5.8; startB = 60330; stopB = 60331} # split from above
    if (i == 5){start = startO + 6;stop = stopO - 0; startB = startO + 0; stopB = startO +2}
    if (i == 6){start = startO + 0;stop = stopO - 0; startB = 60598; stopB = 60600}# split from above
    if (i == 7){start = startO + 0;stop = stopO - 4.6; startB = 60598; stopB = 60600}# slit from #6
    if (i == 8){start = startO + 9;stop = stopO - 9; startB = startO -2; stopB = startO + 0} # tough background
    if (i == 9){start = startO + 7;stop = stopO - 15; startB = startO-1; stopB = startO +1}
    if (i == 10){start = startO + 2.2;stop = stopO - 3; startB = startO-1; stopB = startO +.6}
    if (i == 11){start = startO + 1;stop = stopO - 0; startB = 61456; stopB = 61457.6}# split from above
    if (i == 12){start = startO + 1.6;stop = stopO - 5.2; startB = 61456; stopB = 61457.6} # slit from #11
    if (i == 13){start = startO + 21;stop = stopO - 34; startB = startO -1; stopB = startO +1}
    if (i == 14){start = startO + 4;stop = stopO - 4; startB = startO + 0; stopB = startO +2}
    if (i == 15){start = startO + 24;stop = stopO - 50; startB = startO -1; stopB = startO +1}
    if (i == 16){start = startO + 32;stop = stopO - 4; startB = startO -1; stopB = startO +1}
    if (i == 17){start = startO + 8;stop = stopO - 3; startB = startO + 0; stopB = startO +2}
    if (i == 18){start = startO -2 ;stop = stopO - 86.2; startB = 62493; stopB = 62495} # split from above
    if (i == 19){start = startO + 30;stop = stopO - 23; startB = startO ; stopB = startO }
    if (i == 20){start = startO + 20;stop = stopO - 36; startB = startO + 0; stopB = startO +2} # very low CO
    if (i == 21){start = startO + 52.4;stop = stopO - 10.6; startB = startO ; stopB = startO +.8}
    if (i == 22){start = startO + 63;stop = stopO - 18; startB = startO + 0; stopB = startO +1} # low CO
    if (i == 23){start = startO + 2;stop = stopO - 2; startB = startO -1; stopB = startO +1}
    if (i == 24){start = startO + 0.8;stop = stopO - 0; startB = 63844; stopB = 63846}# split from above
    if (i == 25){start = startO + 0.0;stop = stopO - 219; startB = 63844; stopB = 63846}# split from above
    if (i == 26){start = startO + 8;stop = stopO - 3; startB = startO -1; stopB = startO + .4}
    if (i == 27){start = startO + 13;stop = stopO - 210; startB = 64586; stopB = 64587.4}# split from above
    if (i == 28){start = startO + 3;stop = stopO - 0; startB = startO + 0; stopB = startO +2}
    if (i == 29){start = startO + 0.6;stop = stopO - 0; startB = 65463; stopB = 65465}# split from above
    if (i == 30){start = startO + 17;stop = stopO - 0; startB = 65463; stopB = 65465}  # split from above
    if (i == 31){start = startO + 3.8;stop = stopO - 283; startB = 65463; stopB = 65465}  # split from above
    if (i == 32){start = startO + 1.4;stop = stopO - 4.4; startB = startO -1; stopB = startO +1}
    
    if (i == 33){start = startO + 32;stop = stopO - 34; startB = startO -1; stopB = startO +1}
    if (i == 34){start = startO + 0;stop = stopO - 1; startB = 66465; stopB =66467} # split from above
    if (i == 35){start = startO -1;stop = stopO - 93; startB = 66465; stopB =66467} # split from #34
    tmp = time_align(start,stop,co2.830.5hz,  co.ch4.830.5hz, 
                     warneke.830.5hz,  isaf.830.5hz,
                     rollinsno.830.5hz, rollinsso2.830.5hz, 
                     cit.830.5hz, gtcims.830.5hz,moore.830fast, jimenez.830.5hz,met.830.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
   
    if (i == 4){start= start-1; stop = stop-1; startB=startB-1;stopB=stopB-1}
    if (i == 5){start= start-2; stop = stop-1; startB=startB-1;stopB=stopB-1}
    if (i == 10){start= start-1; stop = stop-0; startB=startB-1;stopB=stopB-1}
    if (i == 21){start= start-1; stop = stop-0; startB=startB-1;stopB=stopB-1}
    if (i == 26){start= start-0; stop = stop-0; startB=startB-1;stopB=stopB-1}
    if (i == 29){start= start-1; stop = stop-0}
    if (i == 31){start= start-1; stop = stop-0}
    if (i == 32){start= start-1; stop = stop-0}
    if (i == 33){start= start-0; stop = stop-0; startB=startB-1;stopB=stopB-1}
    if (i == 34){start= start-0; stop = stop-0; startB=startB-1;stopB=stopB-1}
    if (i == 35){start= start-0; stop = stop-0; startB=startB-1;stopB=stopB-1}
    aug30th.fire= time_alignSLOWNOCANS(start,stop, co2.830.1hz, co.ch4.830.1hz, 
                                       warneke.830.1hz, isaf.830.1hz,rollinsno.830.1hz, rollinsso2.830.1hz, 
                                       cit.830.1hz,gtcims.830.1hz,ryerson.830.1hz , jimenez.830.1hz, 
                                       schwarz.830.1hz, freid.c2h6.830.1hz, freid.ch2o.830.1hz, 
                                       womack.830.1hz, stclair.830.1hz,veres.830.1hz, #wisthaler.830.1hz,
                                       moore.830, met.830.1hz)
    ind2 = which(aug30th.fire$Time_Start >= start & aug30th.fire$Time_Start  <= stop	) 
    ind2B = which(aug30th.fire$Time_Start >= startB & aug30th.fire$Time_Start <= stopB	) 
    print(c(i,' ','Start =',start))
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug30th.fire$fire = ''; aug30th.fire$fuel = ''
    tmp$fuel[ind] = 'forest' ; tmp$fire[ind] = 'BlackwaterRiver'
    aug30th.fire$fuel[ind2] = 'forest' ; aug30th.fire$fire[ind2] = 'BlackwaterRiver'
    # -------- Blake cans 
    indBLAKE = which(blake.830.merge$Time_Start >= (start-5) & blake.830.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.830.merge[indBLAKE,]
      blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.830.merge$Time_Start >= (start-0) & gilman.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.830.merge[indGILMAN,]
      GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      #ind = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])
      # don't use GILMAN if didn't capture peak - dont do this for Blackwater
      #if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} 
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.830.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.830.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.830.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.830.merge$Time_Start >= (start-0) & apel.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.830.merge[indAPEL,]
      APELBG = apel.830.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.830.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.830.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.830.merge[1,]
      APELBG = apel.830.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.830$Time_Start >= (start-0) & newTOGA.830$Time_Start <= stop & newTOGA.830$Percent_plume > 0.5) # check at least 5 sec before?
    
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.830[indBECKY,]
      BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.830$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.830[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.830[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug30th.fire[ind2,])
      plotpass5hz(aug30th.fire[ind2,])
      plotpass1hz(aug30th.fire[ind2,])
      plotpass1hzHNO2(aug30th.fire[ind2,])
      
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug30th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug30th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug30th.fire[ind2,],aug30th.fire[ind2B,],xspecies,'BlackwaterRiver', 'forest',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug30th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug30th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug30th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    tmpEF1$alt = median(aug30th.fire$Pressure_Altitude_YANG[ind2])/1E3
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'BlackwaterRiver', 'forest',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF5$alt = median(aug30th.fire$Pressure_Altitude_YANG[ind2])/1E3
    if ( median(aug30th.fire$transect_type[ind2], na.rm=TRUE) == 1){
      coBW = c(coBW, tmp$CO_DACOM_DISKIN[ind])
      lonBW = c(lonBW, tmp$Longitude_YANG[ind])
      latBW = c(latBW, tmp$Latitude_YANG[ind])
    }    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop   
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){Blackwaterriver.1hz.EF = tmpEF1 ;Blackwaterriver.5hz.EF = tmpEF5 }
    if (i > 1){Blackwaterriver.1hz.EF = rbind(Blackwaterriver.1hz.EF, tmpEF1) ; Blackwaterriver.5hz.EF = rbind(Blackwaterriver.5hz.EF, tmpEF5) }
    
  }
  Blackwaterriver.1hz.EF = Blackwaterriver.1hz.EF[order(Blackwaterriver.1hz.EF$variable),]
  Blackwaterriver.1hz.EF$transect_source_fire_ID = indA[1]
  Blackwaterriver.5hz.EF = Blackwaterriver.5hz.EF[order(Blackwaterriver.5hz.EF$variable),]
  Blackwaterriver.5hz.EF$transect_source_fire_ID = indA[1]
  
  #bwdata = as.data.frame(cbind(coBW, lonBW, latBW))
  #map.us <- get_map(c(-87.3,30.25,-86.65,31.2))
  #ind = which(aug30th.fire$Latitude_YANG <=  31.2 & aug30th.fire$Longitude_YANG <= -86.65 & aug30th.fire$Longitude_YANG >= -87.3 &
  #              aug30th.fire$Radar_Altitude_YANG < 5E3)
  #mapfire = ggmap(map.us) +
   #geom_point(data = bwdata, aes(x = lonBW,y = latBW,colour =coBW),  size = 4)+ 
   # ggtitle("Blackwater") +scale_color_viridis_c(trans = 'log10')
  
  # Analyze Blackwater
  # ------ Merge 1hz and 5hz together ------
  Blackwaterriver.1hz.EF$uniqueid = Blackwaterriver.1hz.EF$pass*Blackwaterriver.1hz.EF$transect_source_fire_ID
  Blackwaterriver.5hz.EF$uniqueid = Blackwaterriver.5hz.EF$pass*Blackwaterriver.5hz.EF$transect_source_fire_ID
  #Blackwaterriver.EF= merge(Blackwaterriver.1hz.EF, Blackwaterriver.5hz.EF, by=c('variable','fire','fuel;, 'transect_source_fire_ID','mWs',
  #                                                                               'nCs',"Start","Stop","StartO","StopO","pass","uniqueid"),
  #                          all = TRUE, suffixes = c(".1hz", ".5hz"))#, incomparables = NA) # x = 1Hz, y=5hz
  # ----- get rid of longitudonal passes?
  #ind = which(Blackwaterriver.EF$transect_type.1hz  == 1)
  #Blackwaterriver.EF = Blackwaterriver.EF[ind,]
  
  #ind = which(Blackwaterriver.EF$variable == 'CO_DACOM_DISKIN' & Blackwaterriver.EF$R2toX.5hz >= 0.8 & Blackwaterriver.EF$maxval.1hz > 400 ) 
  #goodpasses = Blackwaterriver.EF$uniqueid[ind]
  #ind = which(Blackwaterriver.EF$uniqueid %in% goodpasses)
  #Blackwaterriver.EF = Blackwaterriver.EF[ind,]
  
  #ind = which(Blackwaterriver.EF$variable == 'CO_DACOM_DISKIN' & Blackwaterriver.EF$R2toX.5hz >= 0.8 & Blackwaterriver.EF$maxval.1hz > 400 ) 
  #hist(Blackwaterriver.EF$mce.5hz[ind], xlab='MCE', main='Blackwater River (n=21)')
  #hist(Blackwaterriver.EF$age.5hz[ind]/60, xlab='Age, min',main='Blackwater River (n=21)')
 
  #plot(Blackwaterriver.EF$mce.5hz[ind], Blackwaterriver.EF$EF1.5hz[ind], xlab='MCE', ylab='CO, g/kg', pch=19)
  #plot(Blackwaterriver.EF$age.5hz[ind]/60, Blackwaterriver.EF$EF1.5hz[ind], xlab='Age, min', ylab='CO, g/kg', pch=19)
  
  #plot(Blackwaterriver.EF$mce.1hz[ind], Blackwaterriver.EF$EF1.1hz[ind], xlab='MCE', ylab='CO, g/kg', pch=19)
  #plot(Blackwaterriver.EF$age.5hz[ind]/60, Blackwaterriver.EF$EF1.5hz[ind], xlab='Age, min', ylab='CO, g/kg', pch=19)
  
  #par(mfrow=c(2,2), cex=1.2, oma=c(0,0,0,0))
  #ind = which(Blackwaterriver.EF$variable == 'CH4_DACOM_DISKIN' & Blackwaterriver.EF$R2toX.5hz >= 0.8 ) 
  ##plot(Blackwaterriver.EF$mce.5hz[ind], Blackwaterriver.EF$EF1.5hz[ind], xlab='MCE', ylab='CH4, g/kg', pch=19)
  #ind = which(Blackwaterriver.EF$variable == 'C2H6_CAMS_pptv_FRIED' & Blackwaterriver.EF$R2toX.1hz >= 0.8 ) 
  #plot(Blackwaterriver.EF$mce.1hz[ind], Blackwaterriver.EF$EF1.1hz[ind], xlab='MCE', ylab='C2H6, g/kg', pch=19)
  #ind = which(Blackwaterriver.EF$variable == 'Benzene_NOAAPTR_ppbv_WARNEKE' & Blackwaterriver.EF$R2toX.5hz >= 0.8 ) 
  #plot(Blackwaterriver.EF$mce.5hz[ind], Blackwaterriver.EF$EF1.5hz[ind], xlab='MCE', ylab='Benzene, g/kg', pch=19)
  #ind = which(Blackwaterriver.EF$variable == 'Toluene_NOAAPTR_ppbv_WARNEKE' & Blackwaterriver.EF$R2toX.5hz >= 0.8 ) 
  #plot(Blackwaterriver.EF$mce.5hz[ind], Blackwaterriver.EF$EF1.5hz[ind], xlab='MCE', ylab='Toluene, g/kg', pch=19)
#  
  #ind = which(Blackwaterriver.EF$variable == 'C2H5OH_NOAAPTR_ppbv_WARNEKE')# & Blackwaterriver.EF$R2toX.5hz >= 0.8 ) 
  #ind2= which(Blackwaterriver.EF$variable == 'CH4_DACOM_DISKIN' ) 
  #plot(Blackwaterriver.EF$EF1.5hz[ind2], Blackwaterriver.EF$EF1.5hz[ind], xlab='CH4, g/kg', ylab='HCN, g/kg', pch=19)
  #cor.test(Blackwaterriver.EF$EF1.5hz[ind2], Blackwaterriver.EF$EF1.5hz[ind])
  
  #ind = which(Blackwaterriver.EF$variable == 'OC_JIMENEZ' & Blackwaterriver.EF$R2toX.1hz >= 0.8 ) 
  #plot(Blackwaterriver.EF$mce.1hz[ind], Blackwaterriver.EF$EF1.1hz[ind], xlab='MCE', ylab='OC, g/kg', pch=19)
  
  #ind = which(Blackwaterriver.EF$variable == 'Nitrate_JIMENEZ' & Blackwaterriver.EF$R2toX.1hz >= 0.8 ) 
  #plot(Blackwaterriver.EF$mce.1hz[ind], Blackwaterriver.EF$EF1.1hz[ind], xlab='MCE', ylab='SO4, g/kg', pch=19)
  
  #ind = which(Blackwaterriver.EF$variable == 'BC_SCHWARZ' & Blackwaterriver.EF$R2toX.1hz >= 0.85 ) 
  #plot(Blackwaterriver.EF$mce.1hz[ind], Blackwaterriver.EF$EF1.1hz[ind], xlab='MCE', ylab='BC, g/kg', pch=19)
  
  #ind = which(Blackwaterriver.EF$variable == 'CH4_DACOM_DISKIN' & Blackwaterriver.EF$R2toX.5hz >= 0.8 ) 
  #plot(Blackwaterriver.EF$MAtoF.5hz[ind], Blackwaterriver.EF$EF1.5hz[ind], xlab='MAtoF', ylab='CH4, g/kg', pch=19)
  #plot(Blackwaterriver.EF$alt.5hz[ind], Blackwaterriver.EF$EF1.5hz[ind], xlab='Altitude, km', ylab='CH4, g/kg', pch=19)
  #plot(Blackwaterriver.EF$lat.5hz[ind], Blackwaterriver.EF$EF1.5hz[ind], xlab='Latitude', ylab='CH4, g/kg', pch=19)
  
  #ind = which(Blackwaterriver.EF$variable == 'CH4_DACOM_DISKIN' & Blackwaterriver.EF$R2toX.5hz >= 0.8 ) 
  #plot(Blackwaterriver.EF$MAtoF.5hz[ind], Blackwaterriver.EF$EF1.5hz[ind], xlab='MAtoF', ylab='CH4, g/kg', pch=19)
  #plot(Blackwaterriver.EF$alt.5hz[ind], Blackwaterriver.EF$EF1.5hz[ind], xlab='Altitude, km', ylab='CH4, g/kg', pch=19)
  #plot(Blackwaterriver.EF$lat.5hz[ind], Blackwaterriver.EF$EF1.5hz[ind], xlab='Latitude', ylab='CH4, g/kg', pch=19)
  
  #ind = which(Blackwaterriver.EF$variable == 'NO_RYERSON' )#& Blackwaterriver.EF$R2toX.1hz >= 0.85 ) 
  #ind2 = which(Blackwaterriver.EF$variable == 'NO2_RYERSON' )#& Blackwaterriver.EF$R2toX.1hz >= 0.85 ) 
  #ind3= which(Blackwaterriver.EF$variable == 'CH4_DACOM_DISKIN' ) 
  #plot(Blackwaterriver.EF$mce.1hz[ind2],  Blackwaterriver.EF$EF1.1hz[ind]+Blackwaterriver.EF$EF1.1hz[ind2], xlab='MCE', ylab='NOx, g/kg', pch=19)
  
  #get_VOC_EF =function(fireEF){
  #  ind = which(fireEF$variable == 'CH4_DACOM_DISKIN')
  #  VOC_EF = fireEF$EF1.5hz
   # 
  #}
  #plot(Blackwaterriver.EF$mce_int.5hz[ind], Blackwaterriver.EF$EF1int.5hz[ind], xlab='MCE', ylab='CH4, g/kg', pch=19)
  # 
  # ind = which(Blackwaterriver.EF$variable == 'CO_DACOM_DISKIN' & Blackwaterriver.EF$R2toX.5hz >= 0.8 & Blackwaterriver.EF$R2toX.1hz >= 0.8 & 
  #               Blackwaterriver.EF$maxval.1hz > 400 ) 
  # plot(Blackwaterriver.EF$EF1.1hz[ind], Blackwaterriver.EF$EF1.5hz[ind], xlab='1hz', ylab='5hz', pch=19, xlim=c(40,75), ylim=c(40,75))
  # 
  # 
  # ind = which(Blackwaterriver.EF$variable == 'CO_DACOM_DISKIN' & Blackwaterriver.EF$R2toX.5hz >= 0.8 & Blackwaterriver.EF$R2toX.1hz >= 0.8 & 
  #               Blackwaterriver.EF$maxval.1hz > 400 ) 
  # plot(Blackwaterriver.EF$EF1.1hz[ind], Blackwaterriver.EF$EF1.5hz[ind], xlab='1hz', ylab='5hz', pch=19, xlim=c(40,75), ylim=c(40,75))
  # 
  # 
  # 
  # # write out Blackwater for Holly
  # writeBW = 1
  # if (writeBW == 1){
  # 
  #   # ----- get rid of longitudonal passes?
  #   ind = which(Blackwaterriver.EF$pass != 1 & Blackwaterriver.EF$pass !=2 & Blackwaterriver.EF$pass <= 23)
  #   Blackwaterriver.EF = Blackwaterriver.EF[ind,]
  #   
  #    # ---- Cut data to > 400 ppb where CO to CO2 R2 > 0.8
  #   ind = which(Blackwaterriver.EF$variable == 'CO_DACOM_DISKIN' & Blackwaterriver.EF$R2toX.5hz >= 0.8 & Blackwaterriver.EF$maxval.1hz > 400 ) 
  #   goodpasses = Blackwaterriver.EF$uniqueid[ind]
  #   ind = which(Blackwaterriver.EF$uniqueid %in% goodpasses)
  #   Blackwaterriver.EF = Blackwaterriver.EF[ind,]
  #   # ---- Check aging?
  #    
  #   # ----- Cut out all other data where 1hz R2 to CO is < 0.25, assume bad data, don't get rid of NaNs in case they are 1hz but not 5hz
  #   ind = which(Blackwaterriver.EF$R2toX.1hz >= 0.50 | is.na(Blackwaterriver.EF$R2toX.1hz) |  is.nan(Blackwaterriver.EF$R2toX.1hz))
  #   Blackwaterriver.EF = Blackwaterriver.EF[ind,]
  #   ind = which(Blackwaterriver.EF$R2toX.5hz >= 0.50 | is.na(Blackwaterriver.EF$R2toX.5hz) |  is.nan(Blackwaterriver.EF$R2toX.1hz))
  #   Blackwaterriver.EF = Blackwaterriver.EF[ind,]
  #   # get rid of NaN ERs for species
  #   ind = which(is.finite(Blackwaterriver.EF$ERtoX.5hz) | is.finite(Blackwaterriver.EF$ERtoX.1hz))
  #   Blackwaterriver.EF = Blackwaterriver.EF[ind,]
  # 
  #   ind = which(Blackwaterriver.EF$variable == "CO_DACOM_DISKIN" )
  #   plot(Blackwaterriver.EF$EF1.1hz[ind], Blackwaterriver.EF$EF1int.1hz[ind], xlab='EF CO, slope', ylab='EF CO, integration method', main='1hz')
  #   plot(Blackwaterriver.EF$EF1.5hz[ind], Blackwaterriver.EF$EF1int.5hz[ind], xlab='EF CO, slope', ylab='EF CO, integration method', main='5hz')
  #   plot(Blackwaterriver.EF$EF1CO.5hz[ind], Blackwaterriver.EF$EF1COint.5hz[ind], xlab='EF CO, slope', ylab='EF CO, integration method', main='5hz')
  #   plot(Blackwaterriver.EF$EF1int.1hz[ind], Blackwaterriver.EF$EF1int.5hz[ind], xlab='EF CO, int 1hz', ylab='EF CO, int 5hz', main='1hz vs 5hz')
  #   write.csv(Blackwaterriver.EF, file='BlackwaterriverEF.csv')
  # }
  # 
  # ========================== WIGGINS  =====================
  # use ethane for Blake and Gilman
  #plot(aug30th.fire$Ethane_WAS_BLAKE, aug30th.fire$C2H6_CAMS_pptv_FRIED)
  #cor.test(aug30th.fire$Ethane_WAS_BLAKE, aug30th.fire$C2H6_CAMS_pptv_FRIED)
  #plot(aug30th.fire$Ethane_NOAAiWAS_GILMAN, aug30th.fire$C2H6_CAMS_pptv_FRIED)
  #cor.test(aug30th.fire$Ethane_NOAAiWAS_GILMAN, aug30th.fire$C2H6_CAMS_pptv_FRIED)
  # use benzene for Apel
  #plot(aug30th.fire$Benzene_TOGA_APEL, aug30th.fire$Benzene_NOAAPTR_ppbv_WARNEKE)
  #cor.test(aug30th.fire$Benzene_TOGA_APEL, aug30th.fire$Benzene_NOAAPTR_ppbv_WARNEKE)
  fire="Wiggins";fuel="pile"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = start+17 ;stop = stop - 12; startB = startO + 0; stopB = startO +2}
    tmp = time_align(start,stop,co2.830.5hz,  co.ch4.830.5hz, 
                     warneke.830.5hz,  isaf.830.5hz,
                     rollinsno.830.5hz, rollinsso2.830.5hz, 
                     cit.830.5hz, gtcims.830.5hz,moore.830fast, jimenez.830.5hz,met.830.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){ start=start-1; startB=startB-1; stopB=stopB-1}
    aug30th.fire= time_alignSLOWNOCANS(start,stop, co2.830.1hz, co.ch4.830.1hz, 
                                       warneke.830.1hz, isaf.830.1hz,rollinsno.830.1hz, rollinsso2.830.1hz, 
                                       cit.830.1hz,gtcims.830.1hz,ryerson.830.1hz , jimenez.830.1hz, 
                                       schwarz.830.1hz, freid.c2h6.830.1hz, freid.ch2o.830.1hz, 
                                       womack.830.1hz, stclair.830.1hz,veres.830.1hz, #wisthaler.830.1hz,
                                       moore.830, met.830.1hz)
    ind2 = which(aug30th.fire$Time_Start >= start & aug30th.fire$Time_Start  <= stop	) 
    ind2B = which(aug30th.fire$Time_Start >= startB & aug30th.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug30th.fire$fire = ''; aug30th.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Wiggins'
    aug30th.fire$fuel[ind2] = fuel ; aug30th.fire$fire[ind2] = 'Wiggins'
    # -------- Blake cans 
    indBLAKE = which(blake.830.merge$Time_Start >= (start-5) & blake.830.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.830.merge[indBLAKE,]
      blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.830.merge$Time_Start >= (start-0) & gilman.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.830.merge[indGILMAN,]
      GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.830.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.830.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.830.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.830.merge$Time_Start >= (start-0) & apel.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.830.merge[indAPEL,]
      APELBG = apel.830.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.830.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.830.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.830.merge[1,]
      APELBG = apel.830.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.830$Time_Start >= (start-0) & newTOGA.830$Time_Start <= stop & newTOGA.830$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.830[indBECKY,]
      BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.830$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.830[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.830[1,]
    }
    # Becky
    indBECKY = which(newTOGA.830$Time_Start >= (start-0) & newTOGA.830$Time_Start <= stop & newTOGA.830$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.830[indBECKY,]
      BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.830$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.830[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.830[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug30th.fire[ind2,])
      plotpass5hz(aug30th.fire[ind2,])
      plotpass1hz(aug30th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug30th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug30th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug30th.fire[ind2,],aug30th.fire[ind2B,],xspecies,'Wiggins', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug30th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug30th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug30th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Wiggins', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    if (i == 1){WIGGINS.1hz.EF = tmpEF1 ;WIGGINS.5hz.EF = tmpEF5 }
    if (i > 1){WIGGINS.1hz.EF = rbind(WIGGINS.1hz.EF, tmpEF1) ; WIGGINS.5hz.EF = rbind(WIGGINS.5hz.EF, tmpEF5) }
    
  }
  WIGGINS.1hz.EF = WIGGINS.1hz.EF[order(WIGGINS.1hz.EF$variable),]
  WIGGINS.1hz.EF$transect_source_fire_ID = indA[1]
  WIGGINS.5hz.EF = WIGGINS.5hz.EF[order(WIGGINS.5hz.EF$variable),]
  WIGGINS.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== WIGGINS-NEIGHBORS  =====================
  fire="WIggins_neighbors"; fuel="pile"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 2.8;stop = stop -62.4; startB = startO + 0; stopB = startO+1.6}
    if (i == 2){start = startO + 1.6;stop = stop -9.8; startB = startO + 0; stopB = startO+2}
    tmp = time_align(start,stop,co2.830.5hz,  co.ch4.830.5hz, 
                     warneke.830.5hz,  isaf.830.5hz,
                     rollinsno.830.5hz, rollinsso2.830.5hz, 
                     cit.830.5hz, gtcims.830.5hz,moore.830fast, jimenez.830.5hz,met.830.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1; startB=startB-1; stopB=stopB-1}
    if (i ==2){start=start-2; startB=startB-1; stopB=stopB-1}
    ind2 = which(aug30th.fire$Time_Start >= start & aug30th.fire$Time_Start  <= stop	) 
    ind2B = which(aug30th.fire$Time_Start >= startB & aug30th.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug30th.fire$fire = ''; aug30th.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = fire
    aug30th.fire$fuel[ind2] = fuel ; aug30th.fire$fire[ind2] = fire
    # -------- Blake cans 
    indBLAKE = which(blake.830.merge$Time_Start >= (start-5) & blake.830.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.830.merge[indBLAKE,]
      blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.830.merge$Time_Start >= (start-0) & gilman.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.830.merge[indGILMAN,]
      GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.830.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.830.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.830.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.830.merge$Time_Start >= (start-0) & apel.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.830.merge[indAPEL,]
      APELBG = apel.830.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.830.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.830.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.830.merge[1,]
      APELBG = apel.830.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.830$Time_Start >= (start-0) & newTOGA.830$Time_Start <= stop & newTOGA.830$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.830[indBECKY,]
      BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.830$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.830[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.830[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug30th.fire[ind2,])
      plotpass5hz(aug30th.fire[ind2,])
      plotpass1hz(aug30th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug30th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug30th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug30th.fire[ind2,],aug30th.fire[ind2B,],xspecies,'Wiggins-neighbor', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug30th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug30th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug30th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Wiggins-neighbor', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop   
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){WIGGINSNEIGHBORS.1hz.EF = tmpEF1 ;WIGGINSNEIGHBORS.5hz.EF = tmpEF5 }
    if (i > 1){WIGGINSNEIGHBORS.1hz.EF = rbind(WIGGINSNEIGHBORS.1hz.EF, tmpEF1) ; WIGGINSNEIGHBORS.5hz.EF = rbind(WIGGINSNEIGHBORS.5hz.EF, tmpEF5) }
    
  }
  WIGGINSNEIGHBORS.1hz.EF = WIGGINSNEIGHBORS.1hz.EF[order(WIGGINSNEIGHBORS.1hz.EF$variable),]
  WIGGINSNEIGHBORS.1hz.EF$transect_source_fire_ID = indA[1]
  WIGGINSNEIGHBORS.5hz.EF = WIGGINSNEIGHBORS.5hz.EF[order(WIGGINSNEIGHBORS.5hz.EF$variable),]
  WIGGINSNEIGHBORS.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== ZZTOP  =====================
  fire="ZZTop";fuel="pile"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 4;stop = stop -3.2; startB = startO -1; stopB = startO+1}
    if (i == 2){start = startO + 2;stop = stop -21.6; startB = startO -.6; stopB = startO+1}
    tmp = time_align(start,stop,co2.830.5hz,  co.ch4.830.5hz, 
                     warneke.830.5hz,  isaf.830.5hz,
                     rollinsno.830.5hz, rollinsso2.830.5hz, 
                     cit.830.5hz, gtcims.830.5hz,moore.830fast, jimenez.830.5hz,met.830.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i==1){ start = start-1; startB=startB-1;stopB=stopB-1}
    if (i==2){ start = start-2; stop=stop+1;startB=startB-1;stopB=stopB-1}
    aug30th.fire= time_alignSLOWNOCANS(start,stop, co2.830.1hz, co.ch4.830.1hz, 
                                       warneke.830.1hz, isaf.830.1hz,rollinsno.830.1hz, rollinsso2.830.1hz, 
                                       cit.830.1hz,gtcims.830.1hz,ryerson.830.1hz , jimenez.830.1hz, 
                                       schwarz.830.1hz, freid.c2h6.830.1hz, freid.ch2o.830.1hz, 
                                       womack.830.1hz, stclair.830.1hz,veres.830.1hz, #wisthaler.830.1hz,
                                       moore.830, met.830.1hz)
    ind2 = which(aug30th.fire$Time_Start >= start & aug30th.fire$Time_Start  <= stop	) 
    ind2B = which(aug30th.fire$Time_Start >= startB & aug30th.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug30th.fire$fire = ''; aug30th.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'ZZTop'
    aug30th.fire$fuel[ind2] = fuel ; aug30th.fire$fire[ind2] = 'ZZTop'
    # -------- Blake cans 
    indBLAKE = which(blake.830.merge$Time_Start >= (start-5) & blake.830.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.830.merge[indBLAKE,]
      blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.830.merge$Time_Start >= (start-0) & gilman.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.830.merge[indGILMAN,]
      GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.830.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.830.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.830.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.830.merge$Time_Start >= (start-0) & apel.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.830.merge[indAPEL,]
      APELBG = apel.830.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.830.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.830.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.830.merge[1,]
      APELBG = apel.830.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.830$Time_Start >= (start-0) & newTOGA.830$Time_Start <= stop & newTOGA.830$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.830[indBECKY,]
      BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.830$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.830[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.830[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug30th.fire[ind2,])
      plotpass5hz(aug30th.fire[ind2,])
      plotpass1hz(aug30th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug30th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug30th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug30th.fire[ind2,],aug30th.fire[ind2B,],xspecies,'ZZTop', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug30th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug30th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug30th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'ZZTop', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop  
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){ZZTop.1hz.EF = tmpEF1 ;ZZTop.5hz.EF = tmpEF5 }
    if (i > 1){ZZTop.1hz.EF = rbind(ZZTop.1hz.EF, tmpEF1) ; ZZTop.5hz.EF = rbind(ZZTop.5hz.EF, tmpEF5) }
    
  }
  ZZTop.1hz.EF = ZZTop.1hz.EF[order(ZZTop.1hz.EF$variable),]
  ZZTop.1hz.EF$transect_source_fire_ID = indA[1]
  ZZTop.5hz.EF = ZZTop.5hz.EF[order(ZZTop.5hz.EF$variable),]
  ZZTop.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== YOUNGMC  =====================
  fire="YoungMC"; fuel="pile"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 2.6;stop = stop -33.8; startB = startO -2; stopB = startO + 0}
    tmp = time_align(start,stop,co2.830.5hz,  co.ch4.830.5hz, 
                     warneke.830.5hz,  isaf.830.5hz,
                     rollinsno.830.5hz, rollinsso2.830.5hz, 
                     cit.830.5hz, gtcims.830.5hz,moore.830fast, jimenez.830.5hz,met.830.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1; startB=startB-1;stopB=stopB-1}
    aug30th.fire= time_alignSLOWNOCANS(start,stop, co2.830.1hz, co.ch4.830.1hz, 
                                       warneke.830.1hz, isaf.830.1hz,rollinsno.830.1hz, rollinsso2.830.1hz, 
                                       cit.830.1hz,gtcims.830.1hz,ryerson.830.1hz , jimenez.830.1hz, 
                                       schwarz.830.1hz, freid.c2h6.830.1hz, freid.ch2o.830.1hz, 
                                       womack.830.1hz, stclair.830.1hz,veres.830.1hz, #wisthaler.830.1hz,
                                       moore.830, met.830.1hz)
    ind2 = which(aug30th.fire$Time_Start >= start & aug30th.fire$Time_Start  <= stop	) 
    ind2B = which(aug30th.fire$Time_Start >= startB & aug30th.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug30th.fire$fire = ''; aug30th.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'YoungMC'
    aug30th.fire$fuel[ind2] = fuel ; aug30th.fire$fire[ind2] = 'YoungMC'
    # -------- Blake cans 
    indBLAKE = which(blake.830.merge$Time_Start >= (start-5) & blake.830.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.830.merge[indBLAKE,]
      blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.830.merge$Time_Start >= (start-0) & gilman.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.830.merge[indGILMAN,]
      GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      if (GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.830.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.830.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.830.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.830.merge$Time_Start >= (start-0) & apel.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.830.merge[indAPEL,]
      APELBG = apel.830.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.830.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.830.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.830.merge[1,]
      APELBG = apel.830.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.830$Time_Start >= (start-0) & newTOGA.830$Time_Start <= stop & newTOGA.830$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.830[indBECKY,]
      BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.830$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.830[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.830[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug30th.fire[ind2,])
      plotpass5hz(aug30th.fire[ind2,])
      plotpass1hz(aug30th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug30th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug30th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug30th.fire[ind2,],aug30th.fire[ind2B,],xspecies,'YoungMC', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug30th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug30th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug30th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'YoungMC', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop  
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){YoungMC.1hz.EF = tmpEF1 ;YoungMC.5hz.EF = tmpEF5 }
    if (i > 1){YoungMC.1hz.EF = rbind(YoungMC.1hz.EF, tmpEF1) ; YoungMC.5hz.EF = rbind(YoungMC.5hz.EF, tmpEF5) }
    
  }
  YoungMC.1hz.EF = YoungMC.1hz.EF[order(YoungMC.1hz.EF$variable),]
  YoungMC.1hz.EF$transect_source_fire_ID = indA[1]
  YoungMC.5hz.EF = YoungMC.5hz.EF[order(YoungMC.5hz.EF$variable),]
  YoungMC.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== XTC  =====================
  fire="XTC"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    
    if (i == 1){ start = startO + 3; stop = stopO -49; startB = startO + 0; stopB = startO +1.8}
    if (i == 2){ start = startO + 3; stop = stopO -3; startB = startO + 0; stopB = startO +2}
    if (i == 3){ start = startO + 2; stop = stopO -3; startB = startO + 0.2; stopB = startO +2}
    if (i == 4){ start = startO + 0; stop = stopO -3; startB = 72544.2; stopB = 72546} # split from above
    tmp = time_align(start,stop,co2.830.5hz,  co.ch4.830.5hz, 
                     warneke.830.5hz,  isaf.830.5hz,
                     rollinsno.830.5hz, rollinsso2.830.5hz, 
                     cit.830.5hz, gtcims.830.5hz,moore.830fast, jimenez.830.5hz,met.830.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-2;startB=startB-1;stopB=stopB-1}
    if (i ==2){start=start-2;startB=startB-1;stopB=stopB-1}
    if (i ==3){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i ==4){start=start-2;startB=startB-1;stopB=stopB-1}
    aug30th.fire= time_alignSLOWNOCANS(start,stop, co2.830.1hz, co.ch4.830.1hz, 
                                       warneke.830.1hz, isaf.830.1hz,rollinsno.830.1hz, rollinsso2.830.1hz, 
                                       cit.830.1hz,gtcims.830.1hz,ryerson.830.1hz , jimenez.830.1hz, 
                                       schwarz.830.1hz, freid.c2h6.830.1hz, freid.ch2o.830.1hz, 
                                       womack.830.1hz, stclair.830.1hz,veres.830.1hz, #wisthaler.830.1hz,
                                       moore.830, met.830.1hz)
    ind2 = which(aug30th.fire$Time_Start >= start & aug30th.fire$Time_Start  <= stop	) 
    ind2B = which(aug30th.fire$Time_Start >= startB & aug30th.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug30th.fire$fire = ''; aug30th.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'XTC'
    aug30th.fire$fuel[ind2] = fuel ; aug30th.fire$fire[ind2] = 'XTC'
    # -------- Blake cans 
    indBLAKE = which(blake.830.merge$Time_Start >= (start-5) & blake.830.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.830.merge[indBLAKE,]
      blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.830.merge$Time_Start >= (start-0) & gilman.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.830.merge[indGILMAN,]
      GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.830.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.830.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.830.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.830.merge$Time_Start >= (start-0) & apel.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.830.merge[indAPEL,]
      APELBG = apel.830.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.830.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.830.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.830.merge[1,]
      APELBG = apel.830.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.830$Time_Start >= (start-0) & newTOGA.830$Time_Start <= stop & newTOGA.830$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.830[indBECKY,]
      BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.830$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.830[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.830[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug30th.fire[ind2,])
      plotpass5hz(aug30th.fire[ind2,])
      plotpass1hz(aug30th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug30th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug30th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug30th.fire[ind2,],aug30th.fire[ind2B,],xspecies,'XTC', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug30th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug30th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug30th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'XTC', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop    
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){XTC.1hz.EF = tmpEF1 ;XTC.5hz.EF = tmpEF5 }
    if (i > 1){XTC.1hz.EF = rbind(XTC.1hz.EF, tmpEF1) ; XTC.5hz.EF = rbind(XTC.5hz.EF, tmpEF5) }
    
  }
  XTC.1hz.EF = XTC.1hz.EF[order(XTC.1hz.EF$variable),]
  XTC.1hz.EF$transect_source_fire_ID = indA[1]
  XTC.5hz.EF = XTC.5hz.EF[order(XTC.5hz.EF$variable),]
  XTC.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== WEEZER  =====================
  fire="Weezer";fuel="pile"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 0.6;stop = stop -9; startB = startO-1; stopB = startO + .8}
    if (i == 2){start = startO + 3.2;stop = stop -27.8; startB = startO + 0; stopB = startO + 2}
    tmp = time_align(start,stop,co2.830.5hz,  co.ch4.830.5hz, 
                     warneke.830.5hz,  isaf.830.5hz,
                     rollinsno.830.5hz, rollinsso2.830.5hz, 
                     cit.830.5hz, gtcims.830.5hz,moore.830fast, jimenez.830.5hz,met.830.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1; stop = stop + 1}
    if (i ==2){start=start-2; stop = stop+1}
    aug30th.fire= time_alignSLOWNOCANS(start,stop, co2.830.1hz, co.ch4.830.1hz, 
                                       warneke.830.1hz, isaf.830.1hz,rollinsno.830.1hz, rollinsso2.830.1hz, 
                                       cit.830.1hz,gtcims.830.1hz,ryerson.830.1hz , jimenez.830.1hz, 
                                       schwarz.830.1hz, freid.c2h6.830.1hz, freid.ch2o.830.1hz, 
                                       womack.830.1hz, stclair.830.1hz,veres.830.1hz, #wisthaler.830.1hz,
                                       moore.830, met.830.1hz)
    ind2 = which(aug30th.fire$Time_Start >= start & aug30th.fire$Time_Start  <= stop	) 
    ind2B = which(aug30th.fire$Time_Start >= startB & aug30th.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug30th.fire$fire = ''; aug30th.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Weezer'
    aug30th.fire$fuel[ind2] = fuel ; aug30th.fire$fire[ind2] = 'Weezer'
    # -------- Blake cans 
    indBLAKE = which(blake.830.merge$Time_Start >= (start-5) & blake.830.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.830.merge[indBLAKE,]
      blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.830.merge$Time_Start >= (start-0) & gilman.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.830.merge[indGILMAN,]
      GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.830.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.830.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.830.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.830.merge$Time_Start >= (start-0) & apel.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.830.merge[indAPEL,]
      APELBG = apel.830.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.830.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.830.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.830.merge[1,]
      APELBG = apel.830.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.830$Time_Start >= (start-0) & newTOGA.830$Time_Start <= stop & newTOGA.830$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.830[indBECKY,]
      BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.830$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.830[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.830[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug30th.fire[ind2,])
      plotpass5hz(aug30th.fire[ind2,])
      plotpass1hz(aug30th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug30th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug30th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug30th.fire[ind2,],aug30th.fire[ind2B,],xspecies,'Weezer', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug30th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug30th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug30th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Weezer', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop 
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Weezer.1hz.EF = tmpEF1 ;Weezer.5hz.EF = tmpEF5 }
    if (i > 1){Weezer.1hz.EF = rbind(Weezer.1hz.EF, tmpEF1) ; Weezer.5hz.EF = rbind(Weezer.5hz.EF, tmpEF5) }
  }
  Weezer.1hz.EF = Weezer.1hz.EF[order(Weezer.1hz.EF$variable),]
  Weezer.1hz.EF$transect_source_fire_ID = indA[1]
  Weezer.5hz.EF = Weezer.5hz.EF[order(Weezer.5hz.EF$variable),]
  Weezer.5hz.EF$transect_source_fire_ID = indA[1]
  # ========================== XTC-WEEZER-UNKNOWN  - bad fire =====================
  # ========================== W-V  - bad fire =====================
  # ========================== VIOLENT FEMMES  =====================
  fire="Violent Femmes"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 5.6;stop = stop - 0.8; startB = startO + 0; stopB = startO+2} # 
    if (i == 2){start = startO + 6;stop = stop -2; startB = 74376; stopB = 74378} # split from above
    if (i == 3){start = startO + 4;stop = stop - 80; startB = startO + 0.8; stopB = startO+2}
    tmp = time_align(start,stop,co2.830.5hz,  co.ch4.830.5hz, 
                     warneke.830.5hz,  isaf.830.5hz,
                     rollinsno.830.5hz, rollinsso2.830.5hz, 
                     cit.830.5hz, gtcims.830.5hz,moore.830fast, jimenez.830.5hz,met.830.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i==1){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i==2){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i==3){start=start-1;startB=startB-1;stopB=stopB-1}
    aug30th.fire= time_alignSLOWNOCANS(start,stop, co2.830.1hz, co.ch4.830.1hz, 
                                       warneke.830.1hz, isaf.830.1hz,rollinsno.830.1hz, rollinsso2.830.1hz, 
                                       cit.830.1hz,gtcims.830.1hz,ryerson.830.1hz , jimenez.830.1hz, 
                                       schwarz.830.1hz, freid.c2h6.830.1hz, freid.ch2o.830.1hz, 
                                       womack.830.1hz, stclair.830.1hz,veres.830.1hz, #wisthaler.830.1hz,
                                       moore.830, met.830.1hz)
    ind2 = which(aug30th.fire$Time_Start >= start & aug30th.fire$Time_Start  <= stop	) 
    ind2B = which(aug30th.fire$Time_Start >= startB & aug30th.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug30th.fire$fire = ''; aug30th.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'ViolentFemmes'
    aug30th.fire$fuel[ind2] = fuel ; aug30th.fire$fire[ind2] = 'ViolentFemmes'
    # -------- Blake cans 
    indBLAKE = which(blake.830.merge$Time_Start >= (start-5) & blake.830.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.830.merge[indBLAKE,]
      blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.830.merge$Time_Start >= (start-0) & gilman.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.830.merge[indGILMAN,]
      GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])   
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.830.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.830.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.830.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.830.merge$Time_Start >= (start-0) & apel.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.830.merge[indAPEL,]
      APELBG = apel.830.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.830.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.830.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.830.merge[1,]
      APELBG = apel.830.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.830$Time_Start >= (start-0) & newTOGA.830$Time_Start <= stop & newTOGA.830$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.830[indBECKY,]
      BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.830$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.830[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.830[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug30th.fire[ind2,])
      plotpass5hz(aug30th.fire[ind2,])
      plotpass1hz(aug30th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug30th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug30th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug30th.fire[ind2,],aug30th.fire[ind2B,],xspecies,'ViolentFemmes', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug30th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug30th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug30th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'ViolentFemmes', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop   
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){ViolentFemmes.1hz.EF = tmpEF1 ;ViolentFemmes.5hz.EF = tmpEF5 }
    if (i > 1){ViolentFemmes.1hz.EF = rbind(ViolentFemmes.1hz.EF, tmpEF1) ; ViolentFemmes.5hz.EF = rbind(ViolentFemmes.5hz.EF, tmpEF5) }
    
  }
  ViolentFemmes.1hz.EF = ViolentFemmes.1hz.EF[order(ViolentFemmes.1hz.EF$variable),]
  ViolentFemmes.1hz.EF$transect_source_fire_ID = indA[1]
  ViolentFemmes.5hz.EF = ViolentFemmes.5hz.EF[order(ViolentFemmes.5hz.EF$variable),]
  ViolentFemmes.5hz.EF$transect_source_fire_ID = indA[1]
  # ========================== U2  =====================
  fire='U2'; fuel="slash"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 8;  stop = stop - 43; startB = startO + 0; stopB = startO + 2}
    if (i == 2){start = startO + 1.6;stop = stop - 6;  startB = startO-1; stopB = startO + 1}
    if (i == 3){start = startO + 3;  stop = stop - 84; startB = startO + 0; stopB = startO + 2}
    tmp = time_align(start,stop,co2.830.5hz,  co.ch4.830.5hz, 
                     warneke.830.5hz,  isaf.830.5hz,
                     rollinsno.830.5hz, rollinsso2.830.5hz, 
                     cit.830.5hz, gtcims.830.5hz,moore.830fast, jimenez.830.5hz,met.830.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i ==2){start=start-2;stop=stop+1;startB=startB-1;stopB=stopB-1}
    if (i ==3){start=start-1;stop=stop+1;startB=startB-1;stopB=stopB-1}
    aug30th.fire= time_alignSLOWNOCANS(start,stop, co2.830.1hz, co.ch4.830.1hz, 
                                       warneke.830.1hz, isaf.830.1hz,rollinsno.830.1hz, rollinsso2.830.1hz, 
                                       cit.830.1hz,gtcims.830.1hz,ryerson.830.1hz , jimenez.830.1hz, 
                                       schwarz.830.1hz, freid.c2h6.830.1hz, freid.ch2o.830.1hz, 
                                       womack.830.1hz, stclair.830.1hz,veres.830.1hz, #wisthaler.830.1hz,
                                       moore.830, met.830.1hz)
    ind2 = which(aug30th.fire$Time_Start >= start & aug30th.fire$Time_Start  <= stop	) 
    ind2B = which(aug30th.fire$Time_Start >= startB & aug30th.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug30th.fire$fire = ''; aug30th.fire$fuel = ''
    tmp$fuel[ind] = 'slash' ; tmp$fire[ind] = 'U2'
    aug30th.fire$fuel[ind2] = 'slash' ; aug30th.fire$fire[ind2] = 'U2'
    # -------- Blake cans 
    indBLAKE = which(blake.830.merge$Time_Start >= (start-5) & blake.830.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.830.merge[indBLAKE,]
      blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.830.merge$Time_Start >= (start-0) & gilman.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.830.merge[indGILMAN,]
      GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])      
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.830.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.830.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.830.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.830.merge$Time_Start >= (start-0) & apel.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.830.merge[indAPEL,]
      APELBG = apel.830.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.830.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.830.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.830.merge[1,]
      APELBG = apel.830.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.830$Time_Start >= (start-0) & newTOGA.830$Time_Start <= stop & newTOGA.830$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.830[indBECKY,]
      BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.830$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.830[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.830[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug30th.fire[ind2,])
      plotpass5hz(aug30th.fire[ind2,])
      plotpass1hz(aug30th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug30th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug30th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug30th.fire[ind2,],aug30th.fire[ind2B,],xspecies,'U2', 'slash',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug30th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug30th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug30th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'U2', 'slash',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){U2.1hz.EF = tmpEF1 ;U2.5hz.EF = tmpEF5 }
    if (i > 1){U2.1hz.EF = rbind(U2.1hz.EF, tmpEF1) ; U2.5hz.EF = rbind(U2.5hz.EF, tmpEF5) }
    
  }
  U2.1hz.EF = U2.1hz.EF[order(U2.1hz.EF$variable),]
  U2.1hz.EF$transect_source_fire_ID = indA[1]
  U2.5hz.EF = U2.5hz.EF[order(U2.5hz.EF$variable),]
  U2.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== TOTO  =====================
  fire="toto"; fuel="slash"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 1;stop = stopO - 52; startB = startO + 0; stopB = startO +1.6}
    tmp = time_align(start,stop,co2.830.5hz,  co.ch4.830.5hz, 
                     warneke.830.5hz,  isaf.830.5hz,
                     rollinsno.830.5hz, rollinsso2.830.5hz, 
                     cit.830.5hz, gtcims.830.5hz,moore.830fast, jimenez.830.5hz,met.830.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i ==1){start=start-1; startB=startB-1;stopB=stopB-1}
    aug30th.fire= time_alignSLOWNOCANS(start,stop, co2.830.1hz, co.ch4.830.1hz, 
                                       warneke.830.1hz, isaf.830.1hz,rollinsno.830.1hz, rollinsso2.830.1hz, 
                                       cit.830.1hz,gtcims.830.1hz,ryerson.830.1hz , jimenez.830.1hz, 
                                       schwarz.830.1hz, freid.c2h6.830.1hz, freid.ch2o.830.1hz, 
                                       womack.830.1hz, stclair.830.1hz,veres.830.1hz, #wisthaler.830.1hz,
                                       moore.830, met.830.1hz)
    ind2 = which(aug30th.fire$Time_Start >= start & aug30th.fire$Time_Start  <= stop	) 
    ind2B = which(aug30th.fire$Time_Start >= startB & aug30th.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug30th.fire$fire = ''; aug30th.fire$fuel = ''
    tmp$fuel[ind] = 'slash' ; tmp$fire[ind] = fire
    aug30th.fire$fuel[ind2] = 'slash' ; aug30th.fire$fire[ind2] = fire
    # -------- Blake cans 
    indBLAKE = which(blake.830.merge$Time_Start >= (start-5) & blake.830.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.830.merge[indBLAKE,]
      blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.830.merge$Time_Start >= (start-0) & gilman.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.830.merge[indGILMAN,]
      GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])     
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.830.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.830.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.830.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.830.merge$Time_Start >= (start-0) & apel.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.830.merge[indAPEL,]
      APELBG = apel.830.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.830.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.830.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.830.merge[1,]
      APELBG = apel.830.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.830$Time_Start >= (start-0) & newTOGA.830$Time_Start <= stop & newTOGA.830$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.830[indBECKY,]
      BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.830$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.830[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.830[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug30th.fire[ind2,])
      plotpass5hz(aug30th.fire[ind2,])
      plotpass1hz(aug30th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug30th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug30th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug30th.fire[ind2,],aug30th.fire[ind2B,],xspecies,fire,fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug30th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug30th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug30th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire,fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop    
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Toto.1hz.EF = tmpEF1 ;Toto.5hz.EF = tmpEF5 }
    if (i > 1){Toto.1hz.EF = rbind(Toto.1hz.EF, tmpEF1) ; Toto.5hz.EF = rbind(Toto.5hz.EF, tmpEF5) }
    
  }
  Toto.1hz.EF = Toto.1hz.EF[order(Toto.1hz.EF$variable),]
  Toto.1hz.EF$transect_source_fire_ID = indA[1]
  Toto.5hz.EF = Toto.5hz.EF[order(Toto.5hz.EF$variable),]
  Toto.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== SUPERTRAMP  =====================
  fire="Supertramp"; fuel="slash"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 1;stop = stopO - 0.6; startB = startO + 0; stopB = startO +1} # 
    if (i == 2){start = startO + 0.4;stop = stopO - 3.6; startB = 77004; stopB = 77005} # split from above
    if (i == 3){start = startO + 0;stop = stopO - 82; startB = 77004; stopB = 77005} # split from above
    if (i == 4){start = startO + 4;stop = stopO - 12; startB = startO -1; stopB = startO +.8}
    if (i == 5){start = startO + 2;stop = stopO - 33; startB = startO -1; stopB = startO +1}
    if (i == 6){start = startO + 0.8;stop = stopO - 384; startB =77456; stopB = 77458} # take out trailing CO/CO2, use background from #5, background here is confusing
    tmp = time_align(start,stop,co2.830.5hz,  co.ch4.830.5hz, 
                     warneke.830.5hz,  isaf.830.5hz,
                     rollinsno.830.5hz, rollinsso2.830.5hz, 
                     cit.830.5hz, gtcims.830.5hz,moore.830fast, jimenez.830.5hz,met.830.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1;startB=startB-1; stopB=stopB-1}
    if (i ==2){start=start-1;stop=stop+2;startB=startB-1; stopB=stopB-1}
    if (i ==3){start=start-1;startB=startB-1; stopB=stopB-1}
    if (i ==4){start=start-1;startB=startB-1; stopB=stopB-1}
    if (i ==5){start=start-1;startB=startB-1; stopB=stopB-1}
    if (i ==6){start=start-1;startB=startB-1; stopB=stopB-1}
    aug30th.fire= time_alignSLOWNOCANS(start,stop, co2.830.1hz, co.ch4.830.1hz, 
                                       warneke.830.1hz, isaf.830.1hz,rollinsno.830.1hz, rollinsso2.830.1hz, 
                                       cit.830.1hz,gtcims.830.1hz,ryerson.830.1hz , jimenez.830.1hz, 
                                       schwarz.830.1hz, freid.c2h6.830.1hz, freid.ch2o.830.1hz, 
                                       womack.830.1hz, stclair.830.1hz,veres.830.1hz, #wisthaler.830.1hz,
                                       moore.830, met.830.1hz)
    ind2 = which(aug30th.fire$Time_Start >= start & aug30th.fire$Time_Start  <= stop	) 
    ind2B = which(aug30th.fire$Time_Start >= startB & aug30th.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug30th.fire$fire = ''; aug30th.fire$fuel = ''
    tmp$fuel[ind] = 'slash' ; tmp$fire[ind] = 'Supertramp'
    aug30th.fire$fuel[ind2] = 'slash' ; aug30th.fire$fire[ind2] = 'Supertramp'
    # -------- Blake cans 
    indBLAKE = which(blake.830.merge$Time_Start >= (start-5) & blake.830.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.830.merge[indBLAKE,]
      blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.830.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.830.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.830.merge$Time_Start >= (start-0) & gilman.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.830.merge[indGILMAN,]
      GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])     
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.830.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.830.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.830.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.830.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.830.merge$Time_Start >= (start-0) & apel.830.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.830.merge[indAPEL,]
      APELBG = apel.830.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.830.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.830.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.830.merge[1,]
      APELBG = apel.830.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.830$Time_Start >= (start-0) & newTOGA.830$Time_Start <= stop & newTOGA.830$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.830[indBECKY,]
      BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.830$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.830[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.830[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.830[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug30th.fire[ind2,])
      plotpass5hz(aug30th.fire[ind2,])
      plotpass1hz(aug30th.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug30th.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug30th.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug30th.fire[ind2,],aug30th.fire[ind2B,],xspecies,'Supertramp', 'slash',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug30th.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug30th.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug30th.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Supertramp', 'slash',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug30th.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug30th.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop 
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Supertramp.1hz.EF = tmpEF1 ;Supertramp.5hz.EF = tmpEF5 }
    if (i > 1){Supertramp.1hz.EF = rbind(Supertramp.1hz.EF, tmpEF1) ; Supertramp.5hz.EF = rbind(Supertramp.5hz.EF, tmpEF5) }
    
  }
  Supertramp.1hz.EF = Supertramp.1hz.EF[order(Supertramp.1hz.EF$variable),]
  Supertramp.1hz.EF$transect_source_fire_ID = indA[1]
  Supertramp.5hz.EF = Supertramp.5hz.EF[order(Supertramp.5hz.EF$variable),]
  Supertramp.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== NEAR SUPERTRAMP  - not a good fire=====================
  # ----- [[[[[[[[[[[[[[[[[[[[ Aug 31st ]]]]]]]]]]]]]]]]]]]]] -------
  # --------#######--------- Get 1 Hz Data individual 8/31------#######---------
  # plume tags
  tags = getICARTTdataSIMPLE('InputFiles/firexaq-fire-Flags-1HZ_DC8_20190831_R9.ict') ;tags$Time_Start = tags$TIME_START
  # MET DATA
  met.831.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-MetNav_DC8_20190831_R1.ict')
  met.831.1hz = merge(met.831.1hz, tags, by='Time_Start')
  # CO2
  co2.831.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CO2-7000_DC8_20190831_R2.ict')
  # -------- DISKIN -----CO, CH4
  co.ch4.831.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-DACOM_DC8_20190831_R1.ict')
  # --------- WARNEKE ----  VOCs
  warneke.831.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-NOAAPTR-VOCs-1Hz_DC8_20190831_R3.ict')
  # ------ HANISCO - ISAF HCHO - merged to 5Hz from the online merge
  isaf.831.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-ISAF-CH2O-1Hz_DC8_20190831_R0.ict')
  #  ------- ROLLINS - SO2 and NO
  rollinsno.831.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-NO_DC8_20190831_R1.ict')
  rollinsno.831.1hz$Time_Start = rollinsno.831.1hz$time_mid
  rollinsso2.831.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-SO2_DC8_20190831_R1.ict')
  rollinsso2.831.1hz$Time_Start = rollinsso2.831.1hz$time_mid
  
  #  ----- WENNBERG - CIT VOCs - 
  cit.831.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190831_R0_CIT.ict')
  # ------ HUEY - GTCIMS PANs - not sure how to match up peaks here
  gtcims.831.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190831_R0_Huey.ict')
  # ------ RYERSON
  ryerson.A = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NO_DC8_20190831_R1.ict')
  ryerson.B = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NO2_DC8_20190831_R1.ict')
  ryerson.C = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NOy_DC8_20190831_R1.ict')
  ryerson.D = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-O3_DC8_20190831_R1.ict')
  ryerson.831.1hz = cbind(ryerson.A,ryerson.B,ryerson.C,ryerson.D) 
  ryerson.831.1hz$Time_Start = ryerson.831.1hz$Time_start
  # ----- JIMENEZ ---
  jimenez.831.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-EESI_DC8_20190831_R1.ict')
  
  # ----- SCHWARZ ---
  schwarz.831.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-SP2-BC-1HZ_DC8_20190831_R2.ict')
  # ----- FREID ---
  freid.c2h6.831.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-C2H6_DC8_20190831_R3.ict')
  freid.ch2o.831.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CH2O_DC8_20190831_R3.ict')
  # ------ WOMACK ---
  womackA = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-CH3COCHO_DC8_20190831_R1.ict')
  womackB = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-CHOCHO_DC8_20190831_R1.ict')
  womackC = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-HNO2_DC8_20190831_R1.ict')
  womackD = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-NO2_DC8_20190831_R1.ict')
  womack.831.1hz = cbind(womackA, womackB, womackC, womackD)
  # -------St Clair
  stclair.831.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-CANOE-NO2_DC8_20190831_R0.ict')
  stclair.831.1hz$Time_Start = stclair.831.1hz$Time_start
  
  # ------- VERES
  veres.A = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-ClNO2_DC8_20190831_R0.ict')
  veres.B = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HCOOH_DC8_20190831_R1.ict')
  veres.C = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HNO2_DC8_20190831_R1.ict')
  veres.D = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-N2O5_DC8_20190831_R0.ict')
  veres.E = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HPMTF_DC8_20190831_R0.ict')
  veres.F = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-CH3COOCl_DC8_20190831_R0.ict')
  veres.G = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-Cl2_DC8_20190831_R0.ict')
  veres.H = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrCl_DC8_20190831_R0.ict')
  veres.I = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrCN_DC8_20190831_R0.ict')
  veres.J = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrO_DC8_20190831_R0.ict')
  veres.K = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HCN_DC8_20190831_R0.ict')
  veres.L = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HNCO_DC8_20190831_R0.ict')
  veres.831.1hz = cbind(veres.A,veres.B,veres.C,veres.D,veres.E,veres.F,veres.G,veres.H,veres.I,veres.J,veres.K,veres.L)
  
  # --- WISTHALER
  #wisthaler.831.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190831_R0_Wisthaler.ict')
  
  # ---- BLAKE
  blake.831.1hz = getICARTTdataSIMPLE('InputFiles/WAS-MERGE/firexaq-mrgWAS-dc8_merge_20190831_R1.ict')
  cc = colnames(blake.831.1hz)
  blake.831.merge = blake.831.1hz[,c(1,2,96:225)]
  blake.831.merge$CO_DACOM_DISKIN_BLAKE = blake.831.1hz$CO_DACOM_DISKIN
  blake.831.merge$CO2_7000_ppm_DISKIN_BLAKE = blake.831.1hz$CO2_7000_ppm_DISKIN
  
  # ------ APEL
  apel.831.1hz = getICARTTdataSIMPLE('InputFiles/TOGA-MERGE/firexaq-mrgTOGA-dc8_merge_20190831_R1.ict')
  cc = colnames(apel.831.1hz)
  apel.831.merge = apel.831.1hz[,c(1,2,226:315)]
  apel.831.merge$CO_DACOM_DISKIN_APEL = apel.831.1hz$CO_DACOM_DISKIN
  apel.831.merge$CO2_7000_ppm_DISKIN_APEL =apel.831.1hz$CO2_7000_ppm_DISKIN
  
  # Becky's better merge
  file = 'InputFiles/Hornbrook/FIREX-AQ weighted TOGA merge 2022-01-24_0831.xlsx'
  newTOGA.831 = readxl::read_xlsx(file); newTOGA.831[newTOGA.831==-999] = NaN; newTOGA.831[newTOGA.831==-888] = NaN
  newTOGA.831$CO_DACOM_DISKIN_BECKY = newTOGA.831$CO_DACOM_DISKIN
  newTOGA.831$CO2_7000_ppm_DISKIN_BECKY = NaN
  newTOGA.831$Time_Start=newTOGA.831$Time_Start...4
  
  # ----GILMAN
  gilman.831.1hz = getICARTTdataSIMPLE('InputFiles/iWAS-MERGE/firexaq-mrgiWAS-dc8_merge_20190831_R1.ict')
  cc = colnames(gilman.831.1hz)
  gilman.831.merge = gilman.831.1hz[,c(1,2,316:361)]
  gilman.831.merge$CO_DACOM_DISKIN_GILMAN = gilman.831.1hz$CO_DACOM_DISKIN
  gilman.831.merge$CO2_7000_ppm_DISKIN_GILMAN = gilman.831.1hz$CO2_7000_ppm_DISKIN

  # ------ Moore 
  moore.831fast = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190831_R0_MOORE.ict')
  
  moore.831p1 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-AerosolCloudConc_DC8_20190831_R0.ict')
  moore.831p2 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-LAScold_DC8_20190831_R0.ict')
  moore.831p3 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-LAShot_DC8_20190831_R0.ict')
  moore.831p4 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-CPSPD_DC8_20190831_R0.ict')
  moore.831p5 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-CDP_DC8_20190831_R0.ict')
  moore.831p6 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-SMPS_DC8_20190831_R0.ict')
  moore.831 =merge(moore.831p1, moore.831p2, by='Time_mid', all = TRUE, incomparables = NA)
  moore.831 =merge(moore.831, moore.831p3, by='Time_mid', all = TRUE, incomparables = NA)
  moore.831 =merge(moore.831, moore.831p4, by='Time_mid', all = TRUE, incomparables = NA)
  moore.831 =merge(moore.831, moore.831p5, by='Time_mid', all = TRUE, incomparables = NA)
  moore.831 =merge(moore.831, moore.831p6, by='Time_mid', all = TRUE, incomparables = NA)
  # ------- append PI to colnames 1hz ----------
  cc = colnames(co2.831.1hz)
  cc[2] = paste(cc[2],'_DISKIN',sep='')
  colnames(co2.831.1hz) = cc 
  cc = colnames(co.ch4.831.1hz)
  cc[2:4] =  paste(cc[2:4],'_DISKIN',sep='')
  colnames(co.ch4.831.1hz) = cc
  cc=colnames(met.831.1hz)
  
  colnames(isaf.831.1hz) = c("Time_Start"  , "CH2O_ISAF_HANISCO" ,"CH2O_ISAF_precision_HANISCO")
  cc[2:36] = paste(cc[2:36],'_YANG',sep='')
  colnames(met.831.1hz) = cc
  cc = colnames(warneke.831.1hz)
  cc[2:43] = paste(cc[2:43],'_WARNEKE', sep='')
  colnames(warneke.831.1hz) = cc
  cc = colnames(rollinsno.831.1hz)
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsno.831.1hz) = cc
  cc = colnames(rollinsso2.831.1hz)
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsso2.831.1hz) = cc
  # make GTCIMS consistent with 1s merge
  cc=colnames(gtcims.831.1hz)
  cc[4:7] = c("PAN_GTCIMS_HUEY" , "PPN_GTCIMS_HUEY"  ,"APAN_GTCIMS_HUEY" ,"PBN_GTCIMS_HUEY")
  colnames(gtcims.831.1hz)=cc
  # since ISAF, CIT, and GTCIMS came from merge tool, alread has PI's appended.
  ryerson.831.1hz <- ryerson.831.1hz[, !duplicated(colnames(ryerson.831.1hz))]
  
  cc = colnames(ryerson.831.1hz)
  cc[2:9] = paste(cc[2:9],'_RYERSON',sep='')
  colnames(ryerson.831.1hz)=cc
  
  cc = colnames(schwarz.831.1hz)
  cc[2:3] = paste(cc[2:3],'_SCHWARZ', sep='')
  colnames(schwarz.831.1hz) =cc
  
  cc = colnames(freid.c2h6.831.1hz)
  cc[4:7] = paste(cc[4:7], '_FRIED',sep='')
  colnames(freid.c2h6.831.1hz) = cc
  
  cc = colnames(freid.ch2o.831.1hz)
  cc[4:7] = paste(cc[4:7], '_FRIED',sep='')
  colnames(freid.ch2o.831.1hz) = cc
  
  colnames(freid.ch2o.831.1hz) = cc
  womack.831.1hz <- womack.831.1hz[, !duplicated(colnames(womack.831.1hz))]
  cc = colnames(womack.831.1hz)
  cc[2:5] = paste(cc[2:5], '_WOMACK',sep='')
  colnames(womack.831.1hz) = cc

  veres.831.1hz <- veres.831.1hz[, !duplicated(colnames(veres.831.1hz))]
  cc = colnames(veres.831.1hz)
  cc[2:13] = paste(cc[2:13], '_VERES',sep='')
  colnames(veres.831.1hz) = cc
  
  #c = colnames(wisthaler.831.1hz)
  #cc[4:5] = paste(cc[4:5], '_WISTHALER',sep='')
  #colnames(wisthaler.831.1hz) = cc
  cc = colnames(jimenez.831.1hz)
  cc[2:8] = paste(cc[2:8], '_JIMENEZ',sep='')
  colnames(jimenez.831.1hz) = cc
  # --------#######--------- Get 5 or 10 Hz Data 8/31 ------#######---------
  # MET DATA - BUI + YANG + DLH
  met.831.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190831_R0_met.ict')
  #CO, CH4
  co.ch4.831.5hz = getICARTTdataSIMPLE('InputFiles//FIREXAQ-DACOM-5Hz_DC8_20190831_R1.ict')
  # CO2
  co2.831.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CO2-7000-5Hz_DC8_20190831_R1.ict')
  # WARNEKE VOCs
  warneke.831.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-NOAAPTR-VOCs-5Hz_DC8_20190831_R3.ict')
  # ISAF HCHO - merged to 5Hz from the online merge
  isaf.831.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190831_R0_ISAF.ict')
  # ROLLINS SO2 and NO
  rollinsno.831.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-NO-5Hz_DC8_20190831_R0.ict')
  rollinsso2.831.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-SO2-5Hz_DC8_20190831_R1.ict')
  # CIT VOCs - 
  cit.831.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190831_R0_CIT.ict')
  # GTCIMS PANs - not sure how to match up peaks here
  gtcims.831.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190831_R0_huey.ict')
  # ----- Jimenez ---
  jimenez.831.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_AMS_20190831_R0_20230314T134128.ict')
  jimenez.831.5hz$OC_PM1_AMS_JIMENEZ = jimenez.831.5hz$OA_PM1_AMS_JIMENEZ/jimenez.831.5hz$OAtoOC_PM1_AMS
  
  # ------- append PI to colnames ----------
  cc = colnames(co2.831.5hz)
  cc[2] = paste(cc[2],'_DISKIN',sep='')
  colnames(co2.831.5hz) = cc 
  cc = colnames(co.ch4.831.5hz)
  cc[2:4] =  paste(cc[2:4],'_DISKIN',sep='')
  colnames(co.ch4.831.5hz) = cc
  cc=colnames(met.831.5hz)
  cc[2:36] = paste(cc[2:36],'_YANG',sep='')
  colnames(met.831.5hz) = cc
  cc = colnames(warneke.831.5hz)
  cc[2:43] = paste(cc[2:43],'_WARNEKE', sep='')
  colnames(warneke.831.5hz) = cc
  cc = colnames(rollinsno.831.5hz)
  cc[1] = "Time_Start" # really Time_Mid but need for merge
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsno.831.5hz) = cc
  cc = colnames(rollinsso2.831.5hz)
  cc[1] = "Time_Start" # really Time_Mid but need for merge
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsso2.831.5hz) = cc
  # make GTCIMS consistent with 1s merge
  cc=colnames(gtcims.831.5hz)
  cc[38:39] = c("APAN_GTCIMS_HUEY" , "PAN_GTCIMS_HUEY")
  colnames(gtcims.831.5hz)=cc
  # since ISAF, CIT, and GTCIMS came from merge tool, alread has PI's appended.
  # Titanic?
  # ------- get fuel moisture data --------
  if (doFM == 1){
    f1 = '/Users/ktravis1/Library/CloudStorage/Box-Box/FuelMoisture/fuel_moisture_content-20210715T1049Z/fmc_20190831_20Z.nc'
    fid = nc_open(f1)
    fuelMDead = ncvar_get(fid, varid = 'FMCG2D')
    fuelMLive = ncvar_get(fid, varid = 'FMCGLH2D')
    xlon = ncvar_get(fid, varid="XLONG_M")
    xlat = ncvar_get(fid, varid="XLAT_M")
    nc_close(fid)
  }  
  # ========================== JAWS  =====================
  fire="Jaws"; fuel="rice"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 45; stop = stopO -27.4; startB = startO -2; stopB = startO -1}
    if (i == 2){start = startO + 5.6; stop = stopO -62.6; startB = startO  -2; stopB = startO -1} # high background
    
    tmp = time_align(start,stop,co2.831.5hz,  co.ch4.831.5hz, 
                     warneke.831.5hz,  isaf.831.5hz,
                     rollinsno.831.5hz, rollinsso2.831.5hz, 
                     cit.831.5hz, gtcims.831.5hz,moore.831fast, jimenez.831.5hz,met.831.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1; startB=startB-1; stopB=stopB-1}
    if (i ==2){start=start-1; startB=startB-1; stopB=stopB-1}
    aug31st.fire= time_alignSLOWNOCANS(start,stop, co2.831.1hz, co.ch4.831.1hz, 
                                       warneke.831.1hz, isaf.831.1hz,rollinsno.831.1hz, rollinsso2.831.1hz, 
                                       cit.831.1hz,gtcims.831.1hz,ryerson.831.1hz , jimenez.831.1hz, 
                                       schwarz.831.1hz, freid.c2h6.831.1hz, freid.ch2o.831.1hz, 
                                       womack.831.1hz, stclair.831.1hz,veres.831.1hz, #wisthaler.831.1hz,
                                       moore.831, met.831.1hz)
    ind2 = which(aug31st.fire$Time_Start >= start & aug31st.fire$Time_Start  <= stop	) 
    ind2B = which(aug31st.fire$Time_Start >= startB & aug31st.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug31st.fire$fire = ''; aug31st.fire$fuel = ''
    tmp$fuel[ind] = 'rice' ; tmp$fire[ind] = 'Jaws'
    aug31st.fire$fuel[ind2] = 'rice' ; aug31st.fire$fire[ind2] = 'Jaws'
    # -------- Blake cans 
    indBLAKE = which(blake.831.merge$Time_Start >= (start-5) & blake.831.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.831.merge[indBLAKE,]
      blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.831.merge$Time_Start >= (start-0) & gilman.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.831.merge[indGILMAN,]
      GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])    
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.831.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.831.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.831.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.831.merge$Time_Start >= (start-0) & apel.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.831.merge[indAPEL,]
      APELBG = apel.831.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.831.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.831.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.831.merge[1,]
      APELBG = apel.831.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.831$Time_Start >= (start-0) & newTOGA.831$Time_Start <= stop & newTOGA.831$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.831[indBECKY,]
      BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.831$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.831[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.831[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug31st.fire[ind2,])
      plotpass5hz(aug31st.fire[ind2,])
      plotpass1hz(aug31st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      plotpass1hzJUSTNH3(aug31st.fire[ind2,])
      
      # plot background?
      plotpass5hzJUSTCO(aug31st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug31st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug31st.fire[ind2,],aug31st.fire[ind2B,],xspecies,'Jaws', 'rice',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug31st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug31st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug31st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Jaws', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop   
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Jaws.1hz.EF = tmpEF1 ;Jaws.5hz.EF = tmpEF5 }
    if (i > 1){Jaws.1hz.EF = rbind(Jaws.1hz.EF, tmpEF1) ; Jaws.5hz.EF = rbind(Jaws.5hz.EF, tmpEF5) }
    
  }
  Jaws.1hz.EF = Jaws.1hz.EF[order(Jaws.1hz.EF$variable),]
  Jaws.1hz.EF$transect_source_fire_ID = indA[1]
  Jaws.5hz.EF = Jaws.5hz.EF[order(Jaws.5hz.EF$variable),]
  Jaws.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== KINGPIN  =====================
  fire="Kingpin"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 51; stop = stopO -12; startB = startO -4; stopB = startO+0}
    if (i == 2){start = startO + 1.8; stop = stopO -27.6; startB = startO -3; stopB = startO-2}
    if (i == 3){;start = startO + 2; stop = stopO -28; startB = startO + 0; stopB = startO+1}
    tmp = time_align(start,stop,co2.831.5hz,  co.ch4.831.5hz, 
                     warneke.831.5hz,  isaf.831.5hz,
                     rollinsno.831.5hz, rollinsso2.831.5hz, 
                     cit.831.5hz, gtcims.831.5hz,moore.831fast, jimenez.831.5hz,met.831.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1; startB=startB-1; stopB=stopB-1}
    if (i ==2){start=start-1; startB=startB-1; stopB=stopB-1}
    if (i ==3){start=start-1; startB=startB-1; stopB=stopB-1}
    aug31st.fire= time_alignSLOWNOCANS(start,stop, co2.831.1hz, co.ch4.831.1hz, 
                                       warneke.831.1hz, isaf.831.1hz,rollinsno.831.1hz, rollinsso2.831.1hz, 
                                       cit.831.1hz,gtcims.831.1hz,ryerson.831.1hz , jimenez.831.1hz, 
                                       schwarz.831.1hz, freid.c2h6.831.1hz, freid.ch2o.831.1hz, 
                                       womack.831.1hz, stclair.831.1hz,veres.831.1hz, #wisthaler.831.1hz,
                                       moore.831, met.831.1hz)
    ind2 = which(aug31st.fire$Time_Start >= start & aug31st.fire$Time_Start  <= stop	) 
    ind2B = which(aug31st.fire$Time_Start >= startB & aug31st.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug31st.fire$fire = ''; aug31st.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Kingpin'
    aug31st.fire$fuel[ind2] = fuel ; aug31st.fire$fire[ind2] = 'Kingpin'
    # -------- Blake cans 
    indBLAKE = which(blake.831.merge$Time_Start >= (start-5) & blake.831.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.831.merge[indBLAKE,]
      blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.831.merge$Time_Start >= (start-0) & gilman.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.831.merge[indGILMAN,]
      GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])     
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.831.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.831.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.831.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.831.merge$Time_Start >= (start-0) & apel.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.831.merge[indAPEL,]
      APELBG = apel.831.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.831.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.831.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.831.merge[1,]
      APELBG = apel.831.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.831$Time_Start >= (start-0) & newTOGA.831$Time_Start <= stop & newTOGA.831$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.831[indBECKY,]
      BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.831$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.831[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.831[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug31st.fire[ind2,])
      plotpass5hz(aug31st.fire[ind2,])
      plotpass1hz(aug31st.fire[ind2,])
      plotpass1hzBLAKE(aug31st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCH4(tmp[ind,])
      plotpass5hz(tmp[ind,])
      plotslowstuff(tmp,startO,stopO,'Kingpin',i,blake,blakeBG,GILMAN,GILMANBG,BECKY,BECKYBG)
      # plot background?
      plotpass5hzJUSTCO(aug31st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    #plot(tmp$CO_DACOM_DISKIN[ind], tmp$Benzene_NOAAPTR_ppbv_WARNEKE[ind]*1E3, pch=19, ylab='Benzene, ppt', xlab='CO, ppb')
    #points(BECKY$CO_DACOM_DISKIN_BECKY-BECKYBG$CO_DACOM_DISKIN_BECKY, BECKY$Benzene_ppt- BECKYBG$Benzene_ppt, col='green', pch=19)
    #points(GILMAN$CO_DACOM_DISKIN_GILMAN-GILMANBG$CO_DACOM_DISKIN_GILMAN, (GILMAN$Benzene_NOAAiWAS_GILMAN- GILMANBG$Benzene_NOAAiWAS_GILMAN)*1E3, col='purple', pch=19)
    #points(GILMAN$CO_DACOM_DISKIN_GILMAN, (GILMAN$Benzene_NOAAiWAS_GILMAN)*1E3, col='purple', pch=2)
    #points(blake$CO_DACOM_DISKIN_BLAKE-blakeBG$CO_DACOM_DISKIN_BLAKE, blake$Benzene_WAS_BLAKE - blakeBG$Benzene_WAS_BLAKE, col='red', pch=19)
    #tt = lm( tmp$Benzene_NOAAPTR_ppbv_WARNEKE[ind]*1E3 ~ tmp$CO_DACOM_DISKIN[ind])
    #abline(tt)
    #abline(b=1.75,a=-60.6, lty=2)
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug31st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug31st.fire[ind2,],aug31st.fire[ind2B,],xspecies,'Kingpin', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug31st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug31st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug31st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Kingpin', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop 
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Kingpin.1hz.EF = tmpEF1 ;Kingpin.5hz.EF = tmpEF5 }
    if (i > 1){Kingpin.1hz.EF = rbind(Kingpin.1hz.EF, tmpEF1) ; Kingpin.5hz.EF = rbind(Kingpin.5hz.EF, tmpEF5) }
    
  }
  Kingpin.1hz.EF = Kingpin.1hz.EF[order(Kingpin.1hz.EF$variable),]
  Kingpin.1hz.EF$transect_source_fire_ID = indA[1]
  Kingpin.5hz.EF = Kingpin.5hz.EF[order(Kingpin.5hz.EF$variable),]
  Kingpin.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== LEON  =====================
  fire="Leon"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 8.6; stop = stopO -0; startB = startO -6; stopB = startO -4}
   # if (i == 1){start = startO + 0; stop = stopO -0; startB = startO -6; stopB = startO -4}
    if (i == 2){start = startO + 0.4; stop = stopO -5; startB = 68644; stopB = 68646} # split from above
    tmp = time_align(start,stop,co2.831.5hz,  co.ch4.831.5hz, 
                     warneke.831.5hz,  isaf.831.5hz,
                     rollinsno.831.5hz, rollinsso2.831.5hz, 
                     cit.831.5hz, gtcims.831.5hz,moore.831fast, jimenez.831.5hz,met.831.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1}
    if (i ==2){start=start-1}
    aug31st.fire= time_alignSLOWNOCANS(start,stop, co2.831.1hz, co.ch4.831.1hz, 
                                       warneke.831.1hz, isaf.831.1hz,rollinsno.831.1hz, rollinsso2.831.1hz, 
                                       cit.831.1hz,gtcims.831.1hz,ryerson.831.1hz , jimenez.831.1hz, 
                                       schwarz.831.1hz, freid.c2h6.831.1hz, freid.ch2o.831.1hz, 
                                       womack.831.1hz, stclair.831.1hz,veres.831.1hz, #wisthaler.831.1hz,
                                       moore.831, met.831.1hz)
    ind2 = which(aug31st.fire$Time_Start >= start & aug31st.fire$Time_Start  <= stop	) 
    ind2B = which(aug31st.fire$Time_Start >= startB & aug31st.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug31st.fire$fire = ''; aug31st.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Leon'
    aug31st.fire$fuel[ind2] = fuel ; aug31st.fire$fire[ind2] = 'Leon'
    # -------- Blake cans 
    indBLAKE = which(blake.831.merge$Time_Start >= (start-5) & blake.831.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.831.merge[indBLAKE,]
      blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.831.merge$Time_Start >= (start-0) & gilman.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.831.merge[indGILMAN,]
      GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])     
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.831.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.831.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.831.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.831.merge$Time_Start >= (start-0) & apel.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.831.merge[indAPEL,]
      APELBG = apel.831.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.831.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.831.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.831.merge[1,]
      APELBG = apel.831.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.831$Time_Start >= (start-0) & newTOGA.831$Time_Start <= stop & newTOGA.831$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.831[indBECKY,]
      BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.831$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.831[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.831[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug31st.fire[ind2,])
      plotpass5hz(aug31st.fire[ind2,])
      plotpass1hz(aug31st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug31st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug31st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug31st.fire[ind2,],aug31st.fire[ind2B,],xspecies,'Leon', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug31st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug31st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug31st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Leon', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop   
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Leon.1hz.EF = tmpEF1 ;Leon.5hz.EF = tmpEF5 }
    if (i > 1){Leon.1hz.EF = rbind(Leon.1hz.EF, tmpEF1) ; Leon.5hz.EF = rbind(Leon.5hz.EF, tmpEF5) }
    
  }
  Leon.1hz.EF = Leon.1hz.EF[order(Leon.1hz.EF$variable),]
  Leon.1hz.EF$transect_source_fire_ID = indA[1]
  Leon.5hz.EF = Leon.5hz.EF[order(Leon.5hz.EF$variable),]
  Leon.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== MEATBALLS  =====================
  fire="Meatballs"; fuel="soybean"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 0; stop = stopO -3; startB = startO -1; stopB = startO +.8}
    if (i == 2){start = startO + 2.8; stop = stopO -2; startB = startO + 0; stopB = startO+2}
    tmp = time_align(start,stop,co2.831.5hz,  co.ch4.831.5hz, 
                     warneke.831.5hz,  isaf.831.5hz,
                     rollinsno.831.5hz, rollinsso2.831.5hz, 
                     cit.831.5hz, gtcims.831.5hz,moore.831fast, jimenez.831.5hz,met.831.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i ==1){start=start-1}
    if (i ==2){start=start-1; stopB=stopB-1;startB=startB-1}
    aug31st.fire= time_alignSLOWNOCANS(start,stop, co2.831.1hz, co.ch4.831.1hz, 
                                       warneke.831.1hz, isaf.831.1hz,rollinsno.831.1hz, rollinsso2.831.1hz, 
                                       cit.831.1hz,gtcims.831.1hz,ryerson.831.1hz , jimenez.831.1hz, 
                                       schwarz.831.1hz, freid.c2h6.831.1hz, freid.ch2o.831.1hz, 
                                       womack.831.1hz, stclair.831.1hz,veres.831.1hz, #wisthaler.831.1hz,
                                       moore.831, met.831.1hz)
    ind2 = which(aug31st.fire$Time_Start >= start & aug31st.fire$Time_Start  <= stop	) 
    ind2B = which(aug31st.fire$Time_Start >= startB & aug31st.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug31st.fire$fire = ''; aug31st.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Meatballs'
    aug31st.fire$fuel[ind2] = fuel ; aug31st.fire$fire[ind2] = 'Meatballs'
    # -------- Blake cans 
    indBLAKE = which(blake.831.merge$Time_Start >= (start-5) & blake.831.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.831.merge[indBLAKE,]
      blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.831.merge$Time_Start >= (start-0) & gilman.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.831.merge[indGILMAN,]
      GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])    
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.831.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.831.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.831.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.831.merge$Time_Start >= (start-0) & apel.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.831.merge[indAPEL,]
      APELBG = apel.831.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.831.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.831.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.831.merge[1,]
      APELBG = apel.831.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.831$Time_Start >= (start-0) & newTOGA.831$Time_Start <= stop & newTOGA.831$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.831[indBECKY,]
      BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.831$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.831[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.831[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug31st.fire[ind2,])
      plotpass5hz(aug31st.fire[ind2,])
      plotpass1hz(aug31st.fire[ind2,])
      plotpass5hz(aug31st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug31st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug31st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug31st.fire[ind2,],aug31st.fire[ind2B,],xspecies,'Meatballs', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug31st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug31st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug31st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Meatballs', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Meatballs.1hz.EF = tmpEF1 ;Meatballs.5hz.EF = tmpEF5 }
    if (i > 1){Meatballs.1hz.EF = rbind(Meatballs.1hz.EF, tmpEF1) ; Meatballs.5hz.EF = rbind(Meatballs.5hz.EF, tmpEF5) }
    
  }
  Meatballs.1hz.EF = Meatballs.1hz.EF[order(Meatballs.1hz.EF$variable),]
  Meatballs.1hz.EF$transect_source_fire_ID = indA[1]
  Meatballs.5hz.EF = Meatballs.5hz.EF[order(Meatballs.5hz.EF$variable),]
  Meatballs.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== NIKITA  =====================
  fire="Nikita"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 6.6; stop = stopO -6; startB = startO + 0; stopB = startO + 2}
    if (i == 2){start = startO + 8.6; stop = stopO -30; startB = 69725; stopB = 69727} # split from above
    if (i == 3){start = startO + 1; stop = stopO -9; startB = startO -1; stopB = startO + 1}
    if (i == 4){;start = startO + 3; stop = stopO -67; startB = startO -1; stopB = startO + 1}
    if (i == 5){start = startO + 3; stop = stopO -3; startB = startO + 0; stopB = startO + 2}
    tmp = time_align(start,stop,co2.831.5hz,  co.ch4.831.5hz, 
                     warneke.831.5hz,  isaf.831.5hz,
                     rollinsno.831.5hz, rollinsso2.831.5hz, 
                     cit.831.5hz, gtcims.831.5hz,moore.831fast, jimenez.831.5hz,met.831.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i == 1){start=start-1; stopB=stopB-1; startB=startB-1}
    if (i == 2){start=start-1; stopB=stopB-1; startB=startB-1}
    if (i == 3){start=start-1; stopB=stopB-1; startB=startB-1}
    if (i == 4){start=start-1; stopB=stopB-1; startB=startB-1}
    if (i == 5){start=start-1; stopB=stopB-1; startB=startB-1}
    aug31st.fire= time_alignSLOWNOCANS(start,stop, co2.831.1hz, co.ch4.831.1hz, 
                                       warneke.831.1hz, isaf.831.1hz,rollinsno.831.1hz, rollinsso2.831.1hz, 
                                       cit.831.1hz,gtcims.831.1hz,ryerson.831.1hz , jimenez.831.1hz, 
                                       schwarz.831.1hz, freid.c2h6.831.1hz, freid.ch2o.831.1hz, 
                                       womack.831.1hz, stclair.831.1hz,veres.831.1hz, #wisthaler.831.1hz,
                                       moore.831, met.831.1hz)
    ind2 = which(aug31st.fire$Time_Start >= start & aug31st.fire$Time_Start  <= stop	) 
    ind2B = which(aug31st.fire$Time_Start >= startB & aug31st.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug31st.fire$fire = ''; aug31st.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Nikita'
    aug31st.fire$fuel[ind2] = fuel ; aug31st.fire$fire[ind2] = 'Nikita'
    # -------- Blake cans 
    indBLAKE = which(blake.831.merge$Time_Start >= (start-5) & blake.831.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.831.merge[indBLAKE,]
      blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.831.merge$Time_Start >= (start-0) & gilman.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.831.merge[indGILMAN,]
      GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])     
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.831.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.831.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.831.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.831.merge$Time_Start >= (start-0) & apel.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.831.merge[indAPEL,]
      APELBG = apel.831.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.831.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.831.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.831.merge[1,]
      APELBG = apel.831.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.831$Time_Start >= (start-0) & newTOGA.831$Time_Start <= stop & newTOGA.831$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.831[indBECKY,]
      BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.831$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.831[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.831[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug31st.fire[ind2,])
      plotpass5hz(aug31st.fire[ind2,])
      plotpass1hz(aug31st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug31st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug31st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug31st.fire[ind2,],aug31st.fire[ind2B,],xspecies,'Nikita', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug31st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug31st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug31st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Nikita', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop 
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Nikita.1hz.EF = tmpEF1 ;Nikita.5hz.EF = tmpEF5 }
    if (i > 1){Nikita.1hz.EF = rbind(Nikita.1hz.EF, tmpEF1) ; Nikita.5hz.EF = rbind(Nikita.5hz.EF, tmpEF5) }
    
  }
  Nikita.1hz.EF = Nikita.1hz.EF[order(Nikita.1hz.EF$variable),]
  Nikita.1hz.EF$transect_source_fire_ID = indA[1]
  Nikita.5hz.EF = Nikita.5hz.EF[order(Nikita.5hz.EF$variable),]
  Nikita.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== OBLIVION  =====================
  fire="Oblivion"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 2; stop = stopO -3; startB = startO-1; stopB = startO + 1}
    tmp = time_align(start,stop,co2.831.5hz,  co.ch4.831.5hz, 
                     warneke.831.5hz,  isaf.831.5hz,
                     rollinsno.831.5hz, rollinsso2.831.5hz, 
                     cit.831.5hz, gtcims.831.5hz,moore.831fast, jimenez.831.5hz,met.831.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start = start - 1; startB=startB-1;stopB=stopB-1}
    aug31st.fire= time_alignSLOWNOCANS(start,stop, co2.831.1hz, co.ch4.831.1hz, 
                                       warneke.831.1hz, isaf.831.1hz,rollinsno.831.1hz, rollinsso2.831.1hz, 
                                       cit.831.1hz,gtcims.831.1hz,ryerson.831.1hz , jimenez.831.1hz, 
                                       schwarz.831.1hz, freid.c2h6.831.1hz, freid.ch2o.831.1hz, 
                                       womack.831.1hz, stclair.831.1hz,veres.831.1hz, #wisthaler.831.1hz,
                                       moore.831, met.831.1hz)
    ind2 = which(aug31st.fire$Time_Start >= start & aug31st.fire$Time_Start  <= stop	) 
    ind2B = which(aug31st.fire$Time_Start >= startB & aug31st.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug31st.fire$fire = ''; aug31st.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Oblivion'
    aug31st.fire$fuel[ind2] = fuel ; aug31st.fire$fire[ind2] = 'Oblivion'
   
    # -------- Blake cans 
    indBLAKE = which(blake.831.merge$Time_Start >= (start-5) & blake.831.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.831.merge[indBLAKE,]
      blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.831.merge$Time_Start >= (start-0) & gilman.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.831.merge[indGILMAN,]
      GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])     
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.831.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.831.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.831.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.831.merge$Time_Start >= (start-0) & apel.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.831.merge[indAPEL,]
      APELBG = apel.831.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.831.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.831.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.831.merge[1,]
      APELBG = apel.831.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.831$Time_Start >= (start-0) & newTOGA.831$Time_Start <= stop & newTOGA.831$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.831[indBECKY,]
      BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.831$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.831[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.831[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug31st.fire[ind2,])
      plotpass5hz(aug31st.fire[ind2,])
      plotpass1hz(aug31st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug31st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug31st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug31st.fire[ind2,],aug31st.fire[ind2B,],xspecies,'Oblivion', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug31st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug31st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug31st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Oblivion', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop 
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Oblivion.1hz.EF = tmpEF1 ;Oblivion.5hz.EF = tmpEF5 }
    if (i > 1){Oblivion.1hz.EF = rbind(Oblivion.1hz.EF, tmpEF1) ; Oblivion.5hz.EF = rbind(Oblivion.5hz.EF, tmpEF5) }
  }
  Oblivion.1hz.EF = Oblivion.1hz.EF[order(Oblivion.1hz.EF$variable),]
  Oblivion.1hz.EF$transect_source_fire_ID = indA[1]
  Oblivion.5hz.EF = Oblivion.5hz.EF[order(Oblivion.5hz.EF$variable),]
  Oblivion.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== OBLIVION  Jr=====================
  fire="OblivionJr"; fuel="soybean"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    start = startO + 1.4; stop = stopO -2.6; startB = startO -1; stopB = startO +1
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    tmp = time_align(start,stop,co2.831.5hz,  co.ch4.831.5hz, 
                     warneke.831.5hz,  isaf.831.5hz,
                     rollinsno.831.5hz, rollinsso2.831.5hz, 
                     cit.831.5hz, gtcims.831.5hz,moore.831fast, jimenez.831.5hz,met.831.5hz)
    
    start = start - 2; stop=stop+1;startB=startB-1;stopB=stopB-1
    aug31st.fire= time_alignSLOWNOCANS(start,stop, co2.831.1hz, co.ch4.831.1hz, 
                                       warneke.831.1hz, isaf.831.1hz,rollinsno.831.1hz, rollinsso2.831.1hz, 
                                       cit.831.1hz,gtcims.831.1hz,ryerson.831.1hz , jimenez.831.1hz, 
                                       schwarz.831.1hz, freid.c2h6.831.1hz, freid.ch2o.831.1hz, 
                                       womack.831.1hz, stclair.831.1hz,veres.831.1hz, #wisthaler.831.1hz,
                                       moore.831, met.831.1hz)
    ind2 = which(aug31st.fire$Time_Start >= start & aug31st.fire$Time_Start  <= stop	) 
    ind2B = which(aug31st.fire$Time_Start >= startB & aug31st.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug31st.fire$fire = ''; aug31st.fire$fuel = ''
  
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'OblivionJr'
    aug31st.fire$fuel[ind2] = fuel ; aug31st.fire$fire[ind2] = 'OblivionJr' 
    # -------- Blake cans 
    indBLAKE = which(blake.831.merge$Time_Start >= (start-5) & blake.831.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.831.merge[indBLAKE,]
      blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.831.merge$Time_Start >= (start-0) & gilman.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.831.merge[indGILMAN,]
      GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.831.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.831.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.831.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.831.merge$Time_Start >= (start-0) & apel.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.831.merge[indAPEL,]
      APELBG = apel.831.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.831.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.831.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.831.merge[1,]
      APELBG = apel.831.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.831$Time_Start >= (start-0) & newTOGA.831$Time_Start <= stop & newTOGA.831$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.831[indBECKY,]
      BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.831$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.831[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.831[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug31st.fire[ind2,])
      plotpass5hz(aug31st.fire[ind2,])
      plotpass1hz(aug31st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug31st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug31st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug31st.fire[ind2,],aug31st.fire[ind2B,],xspecies,'OblivionJr', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug31st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug31st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug31st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'OblivionJr', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop 
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){OblivionJr.1hz.EF = tmpEF1 ;OblivionJr.5hz.EF = tmpEF5 }
    if (i > 1){OblivionJr.1hz.EF = rbind(OblivionJr.1hz.EF, tmpEF1) ; OblivionJr.5hz.EF = rbind(OblivionJr.5hz.EF, tmpEF5) }
  }
  OblivionJr.1hz.EF = OblivionJr.1hz.EF[order(OblivionJr.1hz.EF$variable),]
  OblivionJr.1hz.EF$transect_source_fire_ID = indA[1]
  OblivionJr.5hz.EF = OblivionJr.5hz.EF[order(OblivionJr.5hz.EF$variable),]
  OblivionJr.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== PSYCHO  =====================
  fire="Psycho"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # pass 4 is bad
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 2; stop = stopO -7; startB = startO + 0; stopB = startO + 2}
    if (i == 2){start = startO + 6; stop = stopO -5; startB = startO -1; stopB = startO + 1}
    if (i == 3){start = startO + 51; stop = stopO -20; startB = startO -1; stopB = startO+1}
    if (i == 4){start = startO + 8.4; stop = stopO -104; startB = startO + 0; stopB = startO+ 2}
    if (i == 5){;start = startO + 0.6; stop = stopO -12; startB = startO -1; stopB = startO + .6}
    if (i == 6){;start = startO + 3; stop = stopO -5; startB = startO-1; stopB = startO + 1}
    tmp = time_align(start,stop,co2.831.5hz,  co.ch4.831.5hz, 
                     warneke.831.5hz,  isaf.831.5hz,
                     rollinsno.831.5hz, rollinsso2.831.5hz, 
                     cit.831.5hz, gtcims.831.5hz,moore.831fast, jimenez.831.5hz,met.831.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1; stopB=stopB-1; startB=startB-1}
    if (i ==2){start=start-1; stopB=stopB-1; startB=startB-1}
    if (i ==3){ stopB=stopB-1; startB=startB-1}
    if (i ==5){start=start-1}
    if (i ==6){start=start-1; stopB=stopB-1; startB=startB-1}
    aug31st.fire= time_alignSLOWNOCANS(start,stop, co2.831.1hz, co.ch4.831.1hz, 
                                       warneke.831.1hz, isaf.831.1hz,rollinsno.831.1hz, rollinsso2.831.1hz, 
                                       cit.831.1hz,gtcims.831.1hz,ryerson.831.1hz , jimenez.831.1hz, 
                                       schwarz.831.1hz, freid.c2h6.831.1hz, freid.ch2o.831.1hz, 
                                       womack.831.1hz, stclair.831.1hz,veres.831.1hz, #wisthaler.831.1hz,
                                       moore.831, met.831.1hz)
    ind2 = which(aug31st.fire$Time_Start >= start & aug31st.fire$Time_Start  <= stop	) 
    ind2B = which(aug31st.fire$Time_Start >= startB & aug31st.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug31st.fire$fire = ''; aug31st.fire$fuel = ''
    
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Psycho'
    aug31st.fire$fuel[ind2] = fuel ; aug31st.fire$fire[ind2] = 'Psycho'
    # -------- Blake cans 
    indBLAKE = which(blake.831.merge$Time_Start >= (start-5) & blake.831.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.831.merge[indBLAKE,]
      blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.831.merge$Time_Start >= (start-0) & gilman.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.831.merge[indGILMAN,]
      GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])    
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.831.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.831.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.831.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.831.merge$Time_Start >= (start-0) & apel.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.831.merge[indAPEL,]
      APELBG = apel.831.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.831.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.831.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.831.merge[1,]
      APELBG = apel.831.merge[1,]
    }
    # Becky
    # Filter for > 1ppm per Becky
    indBECKY = which(newTOGA.831$Time_Start >= (start-0) & newTOGA.831$Time_Start <= stop & newTOGA.831$Percent_plume > 0.5
                     & newTOGA.831$CO_DACOM_DISKIN_BECKY > 1E3) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.831[indBECKY,]
      BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.831$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.831[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.831[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug31st.fire[ind2,])
      plotpass5hz(aug31st.fire[ind2,])
      plotpass1hz(aug31st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug31st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug31st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug31st.fire[ind2,],aug31st.fire[ind2B,],xspecies,'Psycho', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug31st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug31st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug31st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Psycho', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop   
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Psycho.1hz.EF = tmpEF1 ;Psycho.5hz.EF = tmpEF5 }
    if (i > 1){Psycho.1hz.EF = rbind(Psycho.1hz.EF, tmpEF1) ; Psycho.5hz.EF = rbind(Psycho.5hz.EF, tmpEF5) }
    
  }
  Psycho.1hz.EF = Psycho.1hz.EF[order(Psycho.1hz.EF$variable),]
  Psycho.1hz.EF$transect_source_fire_ID = indA[1]
  Psycho.5hz.EF = Psycho.5hz.EF[order(Psycho.5hz.EF$variable),]
  Psycho.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== QUARANTINE  =====================
  fire="Quarantine"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 0; stop = stopO -5; startB = startO -1; stopB = startO +1}
    if (i == 2){start = startO + 1; stop = stopO -16; startB = startO-1; stopB = startO +1}
    tmp = time_align(start,stop,co2.831.5hz,  co.ch4.831.5hz, 
                     warneke.831.5hz,  isaf.831.5hz,
                     rollinsno.831.5hz, rollinsso2.831.5hz, 
                     cit.831.5hz, gtcims.831.5hz,moore.831fast, jimenez.831.5hz,met.831.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1; startB=startB-1;stopB=stopB-1}
    if (i ==2){start=start-1; startB=startB-1;stopB=stopB-1}
    aug31st.fire= time_alignSLOWNOCANS(start,stop, co2.831.1hz, co.ch4.831.1hz, 
                                       warneke.831.1hz, isaf.831.1hz,rollinsno.831.1hz, rollinsso2.831.1hz, 
                                       cit.831.1hz,gtcims.831.1hz,ryerson.831.1hz , jimenez.831.1hz, 
                                       schwarz.831.1hz, freid.c2h6.831.1hz, freid.ch2o.831.1hz, 
                                       womack.831.1hz, stclair.831.1hz,veres.831.1hz, #wisthaler.831.1hz,
                                       moore.831, met.831.1hz)
    ind2 = which(aug31st.fire$Time_Start >= start & aug31st.fire$Time_Start  <= stop	) 
    ind2B = which(aug31st.fire$Time_Start >= startB & aug31st.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug31st.fire$fire = ''; aug31st.fire$fuel = ''
    
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Quarantine'
    aug31st.fire$fuel[ind2] = fuel ; aug31st.fire$fire[ind2] = 'Quarantine'
    # -------- Blake cans 
    indBLAKE = which(blake.831.merge$Time_Start >= (start-5) & blake.831.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.831.merge[indBLAKE,]
      blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.831.merge$Time_Start >= (start-0) & gilman.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.831.merge[indGILMAN,]
      GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])    
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.831.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.831.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.831.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.831.merge$Time_Start >= (start-0) & apel.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.831.merge[indAPEL,]
      APELBG = apel.831.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.831.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.831.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.831.merge[1,]
      APELBG = apel.831.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.831$Time_Start >= (start-0) & newTOGA.831$Time_Start <= stop & newTOGA.831$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.831[indBECKY,]
      BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.831$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.831[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.831[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug31st.fire[ind2,])
      plotpass5hz(aug31st.fire[ind2,])
      plotpass1hz(aug31st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug31st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug31st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug31st.fire[ind2,],aug31st.fire[ind2B,],xspecies,'Quarantine', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug31st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug31st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug31st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Quarantine', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop  
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Quarantine.1hz.EF = tmpEF1 ;Quarantine.5hz.EF = tmpEF5 }
    if (i > 1){Quarantine.1hz.EF = rbind(Quarantine.1hz.EF, tmpEF1) ; Quarantine.5hz.EF = rbind(Quarantine.5hz.EF, tmpEF5) }
    
  }
  Quarantine.1hz.EF = Quarantine.1hz.EF[order(Quarantine.1hz.EF$variable),]
  Quarantine.1hz.EF$transect_source_fire_ID = indA[1]
  Quarantine.5hz.EF = Quarantine.5hz.EF[order(Quarantine.5hz.EF$variable),]
  Quarantine.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== RATATOUILLE  =====================
  fire="Ratatouille"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 1; stop = stopO -3; startB = startO-1; stopB = startO +.7}
    if (i == 2){;start = startO + 1.2; stop = stopO -2; startB = startO -1; stopB = startO +1}
    tmp = time_align(start,stop,co2.831.5hz,  co.ch4.831.5hz, 
                     warneke.831.5hz,  isaf.831.5hz,
                     rollinsno.831.5hz, rollinsso2.831.5hz, 
                     cit.831.5hz, gtcims.831.5hz,moore.831fast, jimenez.831.5hz,met.831.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1; stop=stop+1; startB=startB-1; stopB=stopB-1}
    if (i ==2){start=start-1; stop=stop+1; startB=startB-1; stopB=stopB-1}
    aug31st.fire= time_alignSLOWNOCANS(start,stop, co2.831.1hz, co.ch4.831.1hz, 
                                       warneke.831.1hz, isaf.831.1hz,rollinsno.831.1hz, rollinsso2.831.1hz, 
                                       cit.831.1hz,gtcims.831.1hz,ryerson.831.1hz , jimenez.831.1hz, 
                                       schwarz.831.1hz, freid.c2h6.831.1hz, freid.ch2o.831.1hz, 
                                       womack.831.1hz, stclair.831.1hz,veres.831.1hz, #wisthaler.831.1hz,
                                       moore.831, met.831.1hz)
    ind2 = which(aug31st.fire$Time_Start >= start & aug31st.fire$Time_Start  <= stop	) 
    ind2B = which(aug31st.fire$Time_Start >= startB & aug31st.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug31st.fire$fire = ''; aug31st.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Ratatouille'
    aug31st.fire$fuel[ind2] = fuel ; aug31st.fire$fire[ind2] = 'Ratatouille'
    # -------- Blake cans 
    indBLAKE = which(blake.831.merge$Time_Start >= (start-5) & blake.831.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.831.merge[indBLAKE,]
      blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.831.merge$Time_Start >= (start-0) & gilman.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.831.merge[indGILMAN,]
      GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])    
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.831.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.831.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.831.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.831.merge$Time_Start >= (start-0) & apel.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.831.merge[indAPEL,]
      APELBG = apel.831.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.831.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.831.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.831.merge[1,]
      APELBG = apel.831.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.831$Time_Start >= (start-0) & newTOGA.831$Time_Start <= stop & newTOGA.831$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.831[indBECKY,]
      BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.831$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.831[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.831[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug31st.fire[ind2,])
      plotpass5hz(aug31st.fire[ind2,])
      plotpass1hz(aug31st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug31st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug31st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug31st.fire[ind2,],aug31st.fire[ind2B,],xspecies,'Ratatouille', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug31st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug31st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug31st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Ratatouille', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop 
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Ratatouille.1hz.EF = tmpEF1 ;Ratatouille.5hz.EF = tmpEF5 }
    if (i > 1){Ratatouille.1hz.EF = rbind(Ratatouille.1hz.EF, tmpEF1) ; Ratatouille.5hz.EF = rbind(Ratatouille.5hz.EF, tmpEF5) }
    
  }
  Ratatouille.1hz.EF = Ratatouille.1hz.EF[order(Ratatouille.1hz.EF$variable),]
  Ratatouille.1hz.EF$transect_source_fire_ID = indA[1]
  Ratatouille.5hz.EF = Ratatouille.5hz.EF[order(Ratatouille.5hz.EF$variable),]
  Ratatouille.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== SPACEBALLS  =====================
  fire="Spaceballs"; fuel="rice"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 2; stop = stopO -0; startB = startO + 0; stopB = startO+2}
    if (i == 2){start = startO + 0; stop = stopO -55; startB = 75283; stopB = 75285}# split from above
    if (i == 3){start = startO + 2; stop = stopO -3; startB = startO + 0; stopB = startO+2}
    tmp = time_align(start,stop,co2.831.5hz,  co.ch4.831.5hz, 
                     warneke.831.5hz,  isaf.831.5hz,
                     rollinsno.831.5hz, rollinsso2.831.5hz, 
                     cit.831.5hz, gtcims.831.5hz,moore.831fast, jimenez.831.5hz,met.831.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start = start-1; stopB=stopB-1; startB=startB-1}
    if (i ==2){start = start-1; stopB=stopB-1; startB=startB-1}
    if (i ==3){start = start-1; stopB=stopB-1; startB=startB-1}
    aug31st.fire= time_alignSLOWNOCANS(start,stop, co2.831.1hz, co.ch4.831.1hz, 
                                       warneke.831.1hz, isaf.831.1hz,rollinsno.831.1hz, rollinsso2.831.1hz, 
                                       cit.831.1hz,gtcims.831.1hz,ryerson.831.1hz , jimenez.831.1hz, 
                                       schwarz.831.1hz, freid.c2h6.831.1hz, freid.ch2o.831.1hz, 
                                       womack.831.1hz, stclair.831.1hz,veres.831.1hz, #wisthaler.831.1hz,
                                       moore.831, met.831.1hz)
    ind2 = which(aug31st.fire$Time_Start >= start & aug31st.fire$Time_Start  <= stop	) 
    ind2B = which(aug31st.fire$Time_Start >= startB & aug31st.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug31st.fire$fire = ''; aug31st.fire$fuel = ''
    tmp$fuel[ind] = 'rice' ; tmp$fire[ind] = 'Spaceballs'
    aug31st.fire$fuel[ind2] = 'rice' ; aug31st.fire$fire[ind2] = 'Spaceballs'
    # -------- Blake cans 
    indBLAKE = which(blake.831.merge$Time_Start >= (start-5) & blake.831.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.831.merge[indBLAKE,]
      blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.831.merge$Time_Start >= (start-0) & gilman.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.831.merge[indGILMAN,]
      GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])   
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.831.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.831.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.831.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.831.merge$Time_Start >= (start-0) & apel.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.831.merge[indAPEL,]
      APELBG = apel.831.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.831.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.831.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.831.merge[1,]
      APELBG = apel.831.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.831$Time_Start >= (start-0) & newTOGA.831$Time_Start <= stop & newTOGA.831$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.831[indBECKY,]
      BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.831$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.831[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.831[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug31st.fire[ind2,])
      plotpass5hz(aug31st.fire[ind2,])
      plotpass1hz(aug31st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug31st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug31st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug31st.fire[ind2,],aug31st.fire[ind2B,],xspecies,'Spaceballs', 'rice',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug31st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug31st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug31st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Spaceballs', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop  
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Spaceballs.1hz.EF = tmpEF1 ;Spaceballs.5hz.EF = tmpEF5 }
    if (i > 1){Spaceballs.1hz.EF = rbind(Spaceballs.1hz.EF, tmpEF1) ; Spaceballs.5hz.EF = rbind(Spaceballs.5hz.EF, tmpEF5) }
    
  }
  Spaceballs.1hz.EF = Spaceballs.1hz.EF[order(Spaceballs.1hz.EF$variable),]
  Spaceballs.1hz.EF$transect_source_fire_ID = indA[1]
  Spaceballs.5hz.EF = Spaceballs.5hz.EF[order(Spaceballs.5hz.EF$variable),]
  Spaceballs.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== TREMORS  =====================
  fire="Tremors"; fuel="rice"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){;start = startO + 1; stop = stopO -55; startB = startO -1; stopB = startO + 1}
    if (i == 2){start = startO + 1; stop = stopO -4; startB = startO -1; stopB = startO + .8}
    if (i == 3){start = startO + 4; stop = stopO -4; startB = startO + 0; stopB = startO + 2}
    if (i == 4){start = startO + 6.4; stop = stopO -9.6; startB = startO -1; stopB = startO + 1}
    if (i == 5){start = startO + 3; stop = stopO -3; startB = startO + 0; stopB = startO + 2}
    tmp = time_align(start,stop,co2.831.5hz,  co.ch4.831.5hz, 
                     warneke.831.5hz,  isaf.831.5hz,
                     rollinsno.831.5hz, rollinsso2.831.5hz, 
                     cit.831.5hz, gtcims.831.5hz,moore.831fast, jimenez.831.5hz,met.831.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start = start-1; stopB=stopB-1; startB=startB-1}
    if (i ==2){start = start-1; stopB=stopB-1; startB=startB-1}
    if (i ==3){start = start-1; stopB=stopB-1; startB=startB-1}
    if (i ==4){start = start-2; stop=stop+1;stopB=stopB-1; startB=startB-1}
    aug31st.fire= time_alignSLOWNOCANS(start,stop, co2.831.1hz, co.ch4.831.1hz, 
                                       warneke.831.1hz, isaf.831.1hz,rollinsno.831.1hz, rollinsso2.831.1hz, 
                                       cit.831.1hz,gtcims.831.1hz,ryerson.831.1hz , jimenez.831.1hz, 
                                       schwarz.831.1hz, freid.c2h6.831.1hz, freid.ch2o.831.1hz, 
                                       womack.831.1hz, stclair.831.1hz,veres.831.1hz, #wisthaler.831.1hz,
                                       moore.831, met.831.1hz)
    ind2 = which(aug31st.fire$Time_Start >= start & aug31st.fire$Time_Start  <= stop	) 
    ind2B = which(aug31st.fire$Time_Start >= startB & aug31st.fire$Time_Start <= stopB	) 

    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug31st.fire$fire = ''; aug31st.fire$fuel = ''
    tmp$fuel[ind] = 'rice' ; tmp$fire[ind] = 'Tremors'
    aug31st.fire$fuel[ind2] = 'rice' ; aug31st.fire$fire[ind2] = 'Tremors'
    # -------- Blake cans 
    indBLAKE = which(blake.831.merge$Time_Start >= (start-5) & blake.831.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.831.merge[indBLAKE,]
      blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.831.merge$Time_Start >= (start-0) & gilman.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.831.merge[indGILMAN,]
      GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])   
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.831.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.831.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.831.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.831.merge$Time_Start >= (start-0) & apel.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.831.merge[indAPEL,]
      APELBG = apel.831.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.831.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.831.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.831.merge[1,]
      APELBG = apel.831.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.831$Time_Start >= (start-0) & newTOGA.831$Time_Start <= stop & newTOGA.831$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.831[indBECKY,]
      BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.831$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.831[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.831[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug31st.fire[ind2,])
      plotpass5hz(aug31st.fire[ind2,])
      plotpass1hz(aug31st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug31st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug31st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug31st.fire[ind2,],aug31st.fire[ind2B,],xspecies,'Tremors', 'rice',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug31st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug31st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug31st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Tremors', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop 
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Tremors.1hz.EF = tmpEF1 ;Tremors.5hz.EF = tmpEF5 }
    if (i > 1){Tremors.1hz.EF = rbind(Tremors.1hz.EF, tmpEF1) ; Tremors.5hz.EF = rbind(Tremors.5hz.EF, tmpEF5) }
    
  }
  Tremors.1hz.EF = Tremors.1hz.EF[order(Tremors.1hz.EF$variable),]
  Tremors.1hz.EF$transect_source_fire_ID = indA[1]
  Tremors.5hz.EF = Tremors.5hz.EF[order(Tremors.5hz.EF$variable),]
  Tremors.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== UP  =====================

  fire="Up"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 5.4; stop = stopO -21; startB = startO -1; stopB = startO +.6}
    if (i == 2){start = startO + 2.6; stop = stopO -4.6; startB = startO -1; stopB = startO+.8}
    tmp = time_align(start,stop,co2.831.5hz,  co.ch4.831.5hz, 
                     warneke.831.5hz,  isaf.831.5hz,
                     rollinsno.831.5hz, rollinsso2.831.5hz, 
                     cit.831.5hz, gtcims.831.5hz,moore.831fast, jimenez.831.5hz,met.831.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	)
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start = start-2}
    if (i ==2){start = start-1; startB=startB-1;stopB=stopB-1}
    aug31st.fire= time_alignSLOWNOCANS(start,stop, co2.831.1hz, co.ch4.831.1hz, 
                                       warneke.831.1hz, isaf.831.1hz,rollinsno.831.1hz, rollinsso2.831.1hz, 
                                       cit.831.1hz,gtcims.831.1hz,ryerson.831.1hz , jimenez.831.1hz, 
                                       schwarz.831.1hz, freid.c2h6.831.1hz, freid.ch2o.831.1hz, 
                                       womack.831.1hz, stclair.831.1hz,veres.831.1hz, #wisthaler.831.1hz,
                                       moore.831, met.831.1hz)
    ind2 = which(aug31st.fire$Time_Start >= start & aug31st.fire$Time_Start  <= stop	) 
    ind2B = which(aug31st.fire$Time_Start >= startB & aug31st.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug31st.fire$fire = ''; aug31st.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = fire
    aug31st.fire$fuel[ind2] = fuel ; aug31st.fire$fire[ind2] = fire
    # -------- Blake cans 
    indBLAKE = which(blake.831.merge$Time_Start >= (start-5) & blake.831.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.831.merge[indBLAKE,]
      blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.831.merge$Time_Start >= (start-0) & gilman.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.831.merge[indGILMAN,]
      GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])    
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.831.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.831.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.831.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.831.merge$Time_Start >= (start-0) & apel.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.831.merge[indAPEL,]
      APELBG = apel.831.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.831.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.831.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.831.merge[1,]
      APELBG = apel.831.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.831$Time_Start >= (start-0) & newTOGA.831$Time_Start <= stop & newTOGA.831$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.831[indBECKY,]
      BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.831$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.831[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.831[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug31st.fire[ind2,])
      plotpass5hz(aug31st.fire[ind2,])
      plotpass1hz(aug31st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
       # plot background?
      plotpass5hzJUSTCO(aug31st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug31st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug31st.fire[ind2,],aug31st.fire[ind2B,],xspecies,'Up', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug31st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug31st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug31st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Up', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop  
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Up.1hz.EF = tmpEF1 ;Up.5hz.EF = tmpEF5 }
    if (i > 1){Up.1hz.EF = rbind(Up.1hz.EF, tmpEF1) ; Up.5hz.EF = rbind(Up.5hz.EF, tmpEF5) }
    
  }
  Up.1hz.EF = Up.1hz.EF[order(Up.1hz.EF$variable),]
  Up.1hz.EF$transect_source_fire_ID = indA[1]
  Up.5hz.EF = Up.5hz.EF[order(Up.5hz.EF$variable),]
  Up.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== VERTIGO  =====================

  fire="Vertigo"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO +5; stop = stopO -0; startB = startO -1; stopB = startO+.8}
    if (i == 2){start = startO + 0.6; stop = stopO -5; startB =  77555; stopB =  77556.8}# split from above
    if (i == 3){start = startO + 2.6; stop = stopO -48; startB = startO -1; stopB = startO+1}
    tmp = time_align(start,stop,co2.831.5hz,  co.ch4.831.5hz, 
                     warneke.831.5hz,  isaf.831.5hz,
                     rollinsno.831.5hz, rollinsso2.831.5hz, 
                     cit.831.5hz, gtcims.831.5hz,moore.831fast, jimenez.831.5hz,met.831.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start = start-1}
    if (i ==2){start = start-1}
    if (i ==3){start = start-1; startB=startB-1; stopB=stopB-1}
    aug31st.fire= time_alignSLOWNOCANS(start,stop, co2.831.1hz, co.ch4.831.1hz, 
                                       warneke.831.1hz, isaf.831.1hz,rollinsno.831.1hz, rollinsso2.831.1hz, 
                                       cit.831.1hz,gtcims.831.1hz,ryerson.831.1hz , jimenez.831.1hz, 
                                       schwarz.831.1hz, freid.c2h6.831.1hz, freid.ch2o.831.1hz, 
                                       womack.831.1hz, stclair.831.1hz,veres.831.1hz, #wisthaler.831.1hz,
                                       moore.831, met.831.1hz)
    ind2 = which(aug31st.fire$Time_Start >= start & aug31st.fire$Time_Start  <= stop	) 
    ind2B = which(aug31st.fire$Time_Start >= startB & aug31st.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug31st.fire$fire = ''; aug31st.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Vertigo'
    aug31st.fire$fuel[ind2] = fuel ; aug31st.fire$fire[ind2] = 'Vertigo'
    # -------- Blake cans 
    indBLAKE = which(blake.831.merge$Time_Start >= (start-5) & blake.831.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.831.merge[indBLAKE,]
      blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.831.merge$Time_Start >= (start-0) & gilman.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.831.merge[indGILMAN,]
      GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]) 
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.831.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.831.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.831.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.831.merge$Time_Start >= (start-0) & apel.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.831.merge[indAPEL,]
      APELBG = apel.831.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.831.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.831.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.831.merge[1,]
      APELBG = apel.831.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.831$Time_Start >= (start-0) & newTOGA.831$Time_Start <= stop & newTOGA.831$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.831[indBECKY,]
      BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.831$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.831[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.831[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug31st.fire[ind2,])
      plotpass5hz(aug31st.fire[ind2,])
      plotpass1hz(aug31st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug31st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug31st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug31st.fire[ind2,],aug31st.fire[ind2B,],xspecies,'Vertigo', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug31st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug31st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug31st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Vertigo', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop 
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Vertigo.1hz.EF = tmpEF1 ;Vertigo.5hz.EF = tmpEF5 }
    if (i > 1){Vertigo.1hz.EF = rbind(Vertigo.1hz.EF, tmpEF1) ; Vertigo.5hz.EF = rbind(Vertigo.5hz.EF, tmpEF5) }
    
  }
  Vertigo.1hz.EF = Vertigo.1hz.EF[order(Vertigo.1hz.EF$variable),]
  Vertigo.1hz.EF$transect_source_fire_ID = indA[1]
  Vertigo.5hz.EF = Vertigo.5hz.EF[order(Vertigo.5hz.EF$variable),]
  Vertigo.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== WILLOW  =====================
  fire="Willow"; fuel="slash"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){ # 
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){;start = startO + 2; stop = stopO -.8; startB = startO + 0; stopB = startO+2}
    if (i == 2){start = startO + 0; stop = stopO -.2; startB = 81522; stopB = 81524} # split from above
    if (i == 3){start = startO + 2.4; stop = stopO -3.6; startB = 81522; stopB = 81524} # split from above
    if (i == 4){start = startO + 0; stop = stopO -91; startB = 81522; stopB = 81524} # split from above
    
    if (i == 5){start = startO + 1; stop = stopO -2; startB = startO -1; stopB = startO+1}
    if (i == 6){start = startO + 1; stop = stopO -25; startB = 81968; stopB = 81970} # split from above
    
    if (i == 7){start = startO + 2; stop = stopO -4; startB = startO + 0; stopB = startO+2}
    if (i == 8){start = startO + 0; stop = stopO -0; startB = startO -1; stopB = startO+1}
    if (i == 9){start = startO + 29; stop = stopO -5; startB = 82410; stopB =  82412}# split from above
    
    if (i == 10){start = startO + 2; stop = stopO -5; startB = startO + 0; stopB = startO+2}
    
    if (i == 11){start = startO + 6; stop = stopO -23; startB = startO + 0; stopB = startO+2} # 
    if (i == 12){start = startO + 30.2; stop = stopO -2; startB =  82896; stopB =  82898} # split from above
    
    if (i == 13){start = startO + 1; stop = stopO -1; startB = startO -1; stopB = startO+1}
    if (i == 14){start = startO + 5; stop = stopO -90; startB = 83249; stopB = 83251} # split from above
    
    if (i == 15){start = startO + 2; stop = stopO -5; startB = startO + 0; stopB = startO+2}
    tmp = time_align(start,stop,co2.831.5hz,  co.ch4.831.5hz, 
                     warneke.831.5hz,  isaf.831.5hz,
                     rollinsno.831.5hz, rollinsso2.831.5hz, 
                     cit.831.5hz, gtcims.831.5hz,moore.831fast, jimenez.831.5hz,met.831.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1;stopB=stopB-1;startB=startB-1}
    if (i ==2){start=start-1;stopB=stopB-1;startB=startB-1}
    if (i ==3){start=start-1;stopB=stopB-1;startB=startB-1}
    if (i ==4){start=start-1;stopB=stopB-1;startB=startB-1}
    if (i ==5){start=start-1;stopB=stopB-1;startB=startB-1}
    if (i ==6){start=start-1;stopB=stopB-1;startB=startB-1}
    if (i ==7){start=start-1;stopB=stopB-1;startB=startB-1}
    if (i ==10){start=start-1;stopB=stopB-1;startB=startB-1}
    if (i ==12){start=start-1;stopB=stopB-1;startB=startB-1}
    if (i ==13){start=start-1;stop=stop+2;stopB=stopB-1;startB=startB-1}
    if (i ==14){start=start-1;stop=stop+2;stopB=stopB-1;startB=startB-1}
    if (i ==15){start=start-1;stopB=stopB-1;startB=startB-1}
    aug31st.fire= time_alignSLOWNOCANS(start,stop, co2.831.1hz, co.ch4.831.1hz, 
                                       warneke.831.1hz, isaf.831.1hz,rollinsno.831.1hz, rollinsso2.831.1hz, 
                                       cit.831.1hz,gtcims.831.1hz,ryerson.831.1hz , jimenez.831.1hz, 
                                       schwarz.831.1hz, freid.c2h6.831.1hz, freid.ch2o.831.1hz, 
                                       womack.831.1hz, stclair.831.1hz,veres.831.1hz, #wisthaler.831.1hz,
                                       moore.831, met.831.1hz)
    ind2 = which(aug31st.fire$Time_Start >= start & aug31st.fire$Time_Start  <= stop	) 
    ind2B = which(aug31st.fire$Time_Start >= startB & aug31st.fire$Time_Start <= stopB	) 
  
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    aug31st.fire$fire = ''; aug31st.fire$fuel = ''
    tmp$fuel[ind] = 'slash' ; tmp$fire[ind] = 'Willow'
    aug31st.fire$fuel[ind2] = 'slash' ; aug31st.fire$fire[ind2] = 'Willow'
    # -------- Blake cans 
    indBLAKE = which(blake.831.merge$Time_Start >= (start-5) & blake.831.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.831.merge[indBLAKE,]
      blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.831.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.831.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.831.merge$Time_Start >= (start-0) & gilman.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.831.merge[indGILMAN,]
      GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]) 
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.831.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.831.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.831.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.831.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.831.merge$Time_Start >= (start-0) & apel.831.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.831.merge[indAPEL,]
      APELBG = apel.831.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.831.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.831.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.831.merge[1,]
      APELBG = apel.831.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.831$Time_Start >= (start-0) & newTOGA.831$Time_Start <= stop & newTOGA.831$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.831[indBECKY,]
      BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.831$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.831[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.831[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.831[1,]
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(aug31st.fire[ind2,])
      plotpass5hz(aug31st.fire[ind2,])
      plotpass1hz(aug31st.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(aug31st.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    #par(mfrow=c(1,2),cex=1.5)
    #plot(aug31st.fire$Time_Start[ind2], aug31st.fire$C2H6_CAMS_pptv_FRIED[ind2], ylab='C2H6, ppt', xlab='Time', pch=19, main='Willow (slash fire), pass 7')
    #points(aug31st.fire$Time_Start[ind2], aug31st.fire$Ethane_NOAAiWAS_GILMAN[ind2]*1E3, col='blue', pch=19)
    #points(aug31st.fire$Time_Start[ind2], aug31st.fire$Ethane_WAS_BLAKE[ind2], col='red', pch=19)
    #legend("topleft", c("CAMS","Gilman","Blake"), col=c("black","blue","red"), pch=19)
    
    #plot(aug31st.fire$Time_Start[ind2], aug31st.fire$AcetonePropanal_NOAAPTR_ppbv_WARNEKE[ind2],
    #     ylab='Acetone, ppb', xlab='Time', pch=19, main='Willow (slash fire), pass 7')
    #points(aug31st.fire$Time_Start[ind2], aug31st.fire$Acetone_NOAAiWAS_GILMAN[ind2], col='blue', pch=19)
    #points(aug31st.fire$Time_Start[ind2], aug31st.fire$AcetonePropanal_WAS_BLAKE[ind2]/1E3, col='red', pch=19)
    #legend("topleft", c("CAMS","Gilman","Blake"), col=c("black","blue","red"), pch=19)
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, aug31st.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(aug31st.fire[ind2,],aug31st.fire[ind2B,],xspecies,'Willow', 'slash',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(aug31st.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(aug31st.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(aug31st.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Willow', 'slash',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(aug31st.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(aug31st.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop 
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Willow.1hz.EF = tmpEF1 ;Willow.5hz.EF = tmpEF5 }
    if (i > 1){Willow.1hz.EF = rbind(Willow.1hz.EF, tmpEF1) ; Willow.5hz.EF = rbind(Willow.5hz.EF, tmpEF5) }
    
  }
  Willow.1hz.EF = Willow.1hz.EF[order(Willow.1hz.EF$variable),]
  Willow.1hz.EF$transect_source_fire_ID = indA[1]
  Willow.5hz.EF = Willow.5hz.EF[order(Willow.5hz.EF$variable),]
  Willow.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ----- [[[[[[[[[[[[[[[[[[[[ Sep 3rd ]]]]]]]]]]]]]]]]]]]]] -------
  # ------------- Get 1 Hz Data -------------
  # --------#######--------- Get 1 Hz Data individual 8/21------#######---------
  # plume tags
  tags = getICARTTdataSIMPLE('InputFiles/firexaq-fire-Flags-1HZ_DC8_20190903_R9.ict') ;tags$Time_Start = tags$TIME_START
  # MET DATA
  met.903.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-MetNav_DC8_20190903_R1.ict')
  met.903.1hz = merge(met.903.1hz, tags, by='Time_Start')
  # CO2
  co2.903.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CO2-7000_DC8_20190903_R2.ict')
  # -------- DISKIN -----CO, CH4
  co.ch4.903.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-DACOM_DC8_20190903_R1.ict')
  # --------- WARNEKE ----  VOCs
  warneke.903.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-NOAAPTR-VOCs-1Hz_DC8_20190903_R3.ict')
  # ------ HANISCO - ISAF HCHO - merged to 5Hz from the online merge
  isaf.903.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-ISAF-CH2O-1Hz_DC8_20190903_R0.ict')
  #  ------- ROLLINS - SO2 and NO
  rollinsno.903.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-NO_DC8_20190903_R1.ict')
  rollinsno.903.1hz$Time_Start = rollinsno.903.1hz$time_mid
  rollinsso2.903.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-SO2_DC8_20190903_R1.ict')
  rollinsso2.903.1hz$Time_Start = rollinsso2.903.1hz$time_mid
  #  ----- WENNBERG - CIT VOCs - 
  cit.903.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190903_R0_CIT.ict')
  # ------ HUEY - GTCIMS PANs - not sure how to match up peaks here
  gtcims.903.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190903_R0_Huey.ict')
  # ------ RYERSON
  ryerson.A = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NO_DC8_20190903_R1.ict')
  ryerson.B = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NO2_DC8_20190903_R1.ict')
  ryerson.C = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-NOy_DC8_20190903_R1.ict')
  ryerson.D = getICARTTdataSIMPLE('InputFiles/firexaq-NOyO3-O3_DC8_20190903_R1.ict')
  ryerson.903.1hz = cbind(ryerson.A,ryerson.B,ryerson.C,ryerson.D) ; ryerson.903.1hz$Time_Start = ryerson.903.1hz$Time_start
  # ----- JIMENEZ ---
  jimenez.903.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-EESI_DC8_20190903_R1.ict')
  
  # average to 1s
  #jimenez.903.1hz = aggregate(jimenez.903.1hz, by=list(round(jimenez.903.1hz$Time_Start)), FUN='mean', na.rm=TRUE)
  # ----- SCHWARZ ---
  schwarz.903.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-SP2-BC-1HZ_DC8_20190903_R2.ict')
  # ----- FREID ---
  freid.c2h6.903.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-C2H6_DC8_20190903_R3.ict')
  freid.ch2o.903.1hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CH2O_DC8_20190903_R3.ict')
  # ------ WOMACK ---
  womackA = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-CH3COCHO_DC8_20190903_R1.ict')
  womackB = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-CHOCHO_DC8_20190903_R1.ict')
  womackC = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-HNO2_DC8_20190903_R1.ict')
  womackD = getICARTTdataSIMPLE('InputFiles/firexaq-ACES-NO2_DC8_20190903_R1.ict')
  womack.903.1hz = cbind(womackA, womackB, womackC, womackD)
  # -------St Clair
  stclair.903.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-CANOE-NO2_DC8_20190903_R0.ict')
  stclair.903.1hz$Time_Start = stclair.903.1hz$Time_start
  # ------- VERES
  veres.A = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-ClNO2_DC8_20190903_R0.ict')
  veres.B = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HCOOH_DC8_20190903_R1.ict')
  veres.C = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HNO2_DC8_20190903_R1.ict')
  veres.D = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-N2O5_DC8_20190903_R0.ict')
  veres.E = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HPMTF_DC8_20190903_R0.ict')
  veres.F = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-CH3COOCl_DC8_20190903_R0.ict')
  veres.G = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-Cl2_DC8_20190903_R0.ict')
  veres.H = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrCl_DC8_20190903_R0.ict')
  veres.I = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrCN_DC8_20190903_R0.ict')
  veres.J = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-BrO_DC8_20190903_R0.ict')
  veres.K = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HCN_DC8_20190903_R0.ict')
  veres.L = getICARTTdataSIMPLE('InputFiles/FIREXAQ-NOAACIMS-HNCO_DC8_20190903_R0.ict')
  veres.903.1hz = cbind(veres.A,veres.B,veres.C,veres.D,veres.E,veres.F,veres.G,veres.H,veres.I,veres.J,veres.K,veres.L)
  
  # --- WISTHALER
  #wisthaler.903.1hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg1_dc8_20190903_R0_Wisthaler.ict')

  # ---- BLAKE
  blake.903.1hz = getICARTTdataSIMPLE('InputFiles/WAS-MERGE/firexaq-mrgWAS-dc8_merge_20190903_R1.ict')
  cc = colnames(blake.903.1hz)
  blake.903.merge = blake.903.1hz[,c(1,2,96:225)]
  blake.903.merge$CO_DACOM_DISKIN_BLAKE = blake.903.1hz$CO_DACOM_DISKIN
  blake.903.merge$CO2_7000_ppm_DISKIN_BLAKE = blake.903.1hz$CO2_7000_ppm_DISKIN
  
  # ------ APEL
  apel.903.1hz = getICARTTdataSIMPLE('InputFiles/TOGA-MERGE/firexaq-mrgTOGA-dc8_merge_20190903_R1.ict')
  cc = colnames(apel.903.1hz)
  apel.903.merge = apel.903.1hz[,c(1,2,226:315)]
  apel.903.merge$CO_DACOM_DISKIN_APEL = apel.903.1hz$CO_DACOM_DISKIN
  apel.903.merge$CO2_7000_ppm_DISKIN_APEL =apel.903.1hz$CO2_7000_ppm_DISKIN
  
  # Becky's better merge
  file = 'InputFiles/Hornbrook/FIREX-AQ weighted TOGA merge 2022-01-24_0903.xlsx'
  newTOGA.903 = readxl::read_xlsx(file); newTOGA.903[newTOGA.903==-999] = NaN; newTOGA.903[newTOGA.903==-888] = NaN
  newTOGA.903$CO_DACOM_DISKIN_BECKY = newTOGA.903$CO_DACOM_DISKIN
  newTOGA.903$CO2_7000_ppm_DISKIN_BECKY = NaN
  newTOGA.903$Time_Start=newTOGA.903$Time_Start...4
  
  # ----GILMAN
  gilman.903.1hz = getICARTTdataSIMPLE('InputFiles/iWAS-MERGE/firexaq-mrgiWAS-dc8_merge_20190903_R1.ict')
  cc = colnames(gilman.903.1hz)
  gilman.903.merge = gilman.903.1hz[,c(1,2,316:361)]
  gilman.903.merge$CO_DACOM_DISKIN_GILMAN = gilman.903.1hz$CO_DACOM_DISKIN
  gilman.903.merge$CO2_7000_ppm_DISKIN_GILMAN = gilman.903.1hz$CO2_7000_ppm_DISKIN

  # ------ Moore 
  moore.903fast = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190903_R0_MOORE.ict')
  
  moore.903p1 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-AerosolCloudConc_DC8_20190903_R0.ict')
  moore.903p2 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-LAScold_DC8_20190903_R0.ict')
  moore.903p3 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-LAShot_DC8_20190903_R0.ict')
  moore.903p4 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-CPSPD_DC8_20190903_R0.ict')
  moore.903p5 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-CDP_DC8_20190903_R0.ict')
  moore.903p6 = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LARGE-SMPS_DC8_20190903_R0.ict')
  moore.903 =merge(moore.903p1, moore.903p2, by='Time_mid', all = TRUE, incomparables = NA)
  moore.903 =merge(moore.903, moore.903p3, by='Time_mid', all = TRUE, incomparables = NA)
  moore.903 =merge(moore.903, moore.903p4, by='Time_mid', all = TRUE, incomparables = NA)
  moore.903 =merge(moore.903, moore.903p5, by='Time_mid', all = TRUE, incomparables = NA)
  moore.903 =merge(moore.903, moore.903p6, by='Time_mid', all = TRUE, incomparables = NA)
  
  # ------- append PI to colnames 1hz ----------
  cc = colnames(co2.903.1hz)
  cc[2] = paste(cc[2],'_DISKIN',sep='')
  colnames(co2.903.1hz) = cc 
  cc = colnames(co.ch4.903.1hz)
  cc[2:4] =  paste(cc[2:4],'_DISKIN',sep='')
  colnames(co.ch4.903.1hz) = cc
  cc=colnames(met.903.1hz)
  
  colnames(isaf.903.1hz) = c("Time_Start"  , "CH2O_ISAF_HANISCO" ,"CH2O_ISAF_precision_HANISCO")
  cc[2:36] = paste(cc[2:36],'_YANG',sep='')
  colnames(met.903.1hz) = cc
  cc = colnames(warneke.903.1hz)
  cc[2:43] = paste(cc[2:43],'_WARNEKE', sep='')
  colnames(warneke.903.1hz) = cc
  cc = colnames(rollinsno.903.1hz)
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsno.903.1hz) = cc
  cc = colnames(rollinsso2.903.1hz)
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsso2.903.1hz) = cc
  # make GTCIMS consistent with 1s merge
  cc=colnames(gtcims.903.1hz)
  cc[4:7] = c("PAN_GTCIMS_HUEY" , "PPN_GTCIMS_HUEY"  ,"APAN_GTCIMS_HUEY" ,"PBN_GTCIMS_HUEY")
  colnames(gtcims.903.1hz)=cc
  # since ISAF, CIT, and GTCIMS came from merge tool, alread has PI's appended.
  ryerson.903.1hz <- ryerson.903.1hz[, !duplicated(colnames(ryerson.903.1hz))]
  
  cc = colnames(ryerson.903.1hz)
  cc[2:9] = paste(cc[2:9],'_RYERSON',sep='')
  colnames(ryerson.903.1hz)=cc
  
  cc = colnames(schwarz.903.1hz)
  cc[2:3] = paste(cc[2:3],'_SCHWARZ', sep='')
  colnames(schwarz.903.1hz) =cc
  
  cc = colnames(freid.c2h6.903.1hz)
  cc[4:7] = paste(cc[4:7], '_FRIED',sep='')
  colnames(freid.c2h6.903.1hz) = cc
  
  cc = colnames(freid.ch2o.903.1hz)
  cc[4:7] = paste(cc[4:7], '_FRIED',sep='')
  colnames(freid.ch2o.903.1hz) = cc
  
  colnames(freid.ch2o.903.1hz) = cc
  womack.903.1hz <- womack.903.1hz[, !duplicated(colnames(womack.903.1hz))]
  cc = colnames(womack.903.1hz)
  cc[2:5] = paste(cc[2:5], '_WOMACK',sep='')
  colnames(womack.903.1hz) = cc
  
  veres.903.1hz <- veres.903.1hz[, !duplicated(colnames(veres.903.1hz))]
  cc = colnames(veres.903.1hz)
  cc[2:13] = paste(cc[2:13], '_VERES',sep='')
  colnames(veres.903.1hz) = cc
  
  #cc = colnames(wisthaler.903.1hz)
  #cc[4:5] = paste(cc[4:5], '_WISTHALER',sep='')
  #colnames(wisthaler.903.1hz) = cc
  cc = colnames(jimenez.903.1hz)
  cc[2:8] = paste(cc[2:8], '_JIMENEZ',sep='')
  colnames(jimenez.903.1hz) = cc
  # --------#######--------- Get 5 or 10 Hz Data 9/03 ------#######---------
  # MET DATA - BUI + YANG + DLH
  met.903.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190903_R0_met.ict')
  #CO, CH4
  co.ch4.903.5hz = getICARTTdataSIMPLE('InputFiles//FIREXAQ-DACOM-5Hz_DC8_20190903_R1.ict')
  # CO2
  co2.903.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-CO2-7000-5Hz_DC8_20190903_R1.ict')
  # WARNEKE VOCs
  warneke.903.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-NOAAPTR-VOCs-5Hz_DC8_20190903_R3.ict')
  # ISAF HCHO - merged to 5Hz from the online merge
  isaf.903.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190903_R0_ISAF.ict')
  # ROLLINS SO2 and NO
  rollinsno.903.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-NO-5Hz_DC8_20190903_R0.ict')
  rollinsso2.903.5hz = getICARTTdataSIMPLE('InputFiles/FIREXAQ-LIF-SO2-5Hz_DC8_20190903_R1.ict')
  # CIT VOCs - 
  cit.903.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190903_R0_CIT.ict')
  # GTCIMS PANs - not sure how to match up peaks here
  gtcims.903.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_dc8_20190903_R0_huey.ict')
  # ----- Jimenez ---
  jimenez.903.5hz = getICARTTdataSIMPLE('InputFiles/firexaq-mrg5hz_AMS_20190903_R0_20230314T134139.ict')
  jimenez.903.5hz$OC_PM1_AMS_JIMENEZ = jimenez.903.5hz$OA_PM1_AMS_JIMENEZ/jimenez.903.5hz$OAtoOC_PM1_AMS
  
  # ------- append PI to colnames ----------
  cc = colnames(co2.903.5hz)
  cc[2] = paste(cc[2],'_DISKIN',sep='')
  colnames(co2.903.5hz) = cc 
  cc = colnames(co.ch4.903.5hz)
  cc[2:4] =  paste(cc[2:4],'_DISKIN',sep='')
  colnames(co.ch4.903.5hz) = cc
  cc=colnames(met.903.5hz)
  cc[2:36] = paste(cc[2:36],'_YANG',sep='')
  colnames(met.903.5hz) = cc
  cc = colnames(warneke.903.5hz)
  cc[2:43] = paste(cc[2:43],'_WARNEKE', sep='')
  colnames(warneke.903.5hz) = cc
  cc = colnames(rollinsno.903.5hz)
  cc[1] = "Time_Start" # really Time_Mid but need for merge
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsno.903.5hz) = cc
  cc = colnames(rollinsso2.903.5hz)
  cc[1] = "Time_Start" # really Time_Mid but need for merge
  cc[2] = paste(cc[2],'_ROLLINS',sep='')
  colnames(rollinsso2.903.5hz) = cc
  # make GTCIMS consistent with 1s merge
  cc=colnames(gtcims.903.5hz)
  cc[38:39] = c("APAN_GTCIMS_HUEY" , "PAN_GTCIMS_HUEY")
  colnames(gtcims.903.5hz)=cc
  # since ISAF, CIT, and GTCIMS came from merge tool, alread has PI's appended.
  # ------- get fuel moisture data --------
  if (doFM == 1){
    f1 = '/Users/ktravis1/Library/CloudStorage/Box-Box/FuelMoisture/fuel_moisture_content-20210715T1049Z/fmc_20190903_20Z.nc'
    fid = nc_open(f1)
    fuelMDead = ncvar_get(fid, varid = 'FMCG2D')
    fuelMLive = ncvar_get(fid, varid = 'FMCGLH2D')
    xlon = ncvar_get(fid, varid="XLONG_M")
    xlat = ncvar_get(fid, varid="XLAT_M")
    nc_close(fid)
  }  
  # ========================== ASTERISK  * need to check=====================
  fire="Asterisk"; fuel="?"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 166; stop = stopO -5; startB = startO + 0; stopB = startO + 2} # maybe cut down more or split?
   
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = '?' ; tmp$fire[ind] = 'Asterisk'
    sep03rd.fire$fuel[ind2] = '?' ; sep03rd.fire$fire[ind2] = 'Asterisk'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])    
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,'Asterisk', '?',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Asterisk', '?',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Asterisk.1hz.EF = tmpEF1 ;Asterisk.5hz.EF = tmpEF5 }
    if (i > 1){Asterisk.1hz.EF = rbind(Asterisk.1hz.EF, tmpEF1) ; Asterisk.5hz.EF = rbind(Asterisk.5hz.EF, tmpEF5) }
    
  }
  Asterisk.1hz.EF = Asterisk.1hz.EF[order(Asterisk.1hz.EF$variable),]
  Asterisk.1hz.EF$transect_source_fire_ID = indA[1]
  Asterisk.5hz.EF = Asterisk.5hz.EF[order(Asterisk.5hz.EF$variable),]
  Asterisk.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== BUGS BUNNY =====================
  # use ethane for Blake and Gilman
  #plot(sep03rd.fire$Ethane_WAS_BLAKE, sep03rd.fire$C2H6_CAMS_pptv_FRIED)
  #cor.test(sep03rd.fire$Ethane_WAS_BLAKE, sep03rd.fire$C2H6_CAMS_pptv_FRIED)
  #plot(sep03rd.fire$Ethane_NOAAiWAS_GILMAN, sep03rd.fire$C2H6_CAMS_pptv_FRIED)
  #cor.test(sep03rd.fire$Ethane_NOAAiWAS_GILMAN, sep03rd.fire$C2H6_CAMS_pptv_FRIED)
  # use benzene for Apel
  #plot(sep03rd.fire$Benzene_TOGA_APEL, sep03rd.fire$Benzene_NOAAPTR_ppbv_WARNEKE)
  #cor.test(sep03rd.fire$Benzene_TOGA_APEL, sep03rd.fire$Benzene_NOAAPTR_ppbv_WARNEKE)
  fire="Bugs Bunny"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 256; stop = stopO -10.6; startB = startO + 0; stopB = startO + 2}
    if (i == 2){start = startO + 4.6; stop = stopO -7.6; startB = startO + 0; stopB = startO + 2}
    if (i == 3){start = startO + 1; stop = stopO -71;   startB = startO -2; stopB = startO +0}
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i ==1){start=start-1}
    if (i ==2){start=start-1}
    if (i ==3){start=start-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'BugsBunny'
    sep03rd.fire$fuel[ind2] = fuel ; sep03rd.fire$fire[ind2] = 'BugsBunny'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])    
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,'BugsBunny', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'BugsBunny', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){BugsBunny.1hz.EF = tmpEF1 ;BugsBunny.5hz.EF = tmpEF5 }
    if (i > 1){BugsBunny.1hz.EF = rbind(BugsBunny.1hz.EF, tmpEF1) ; BugsBunny.5hz.EF = rbind(BugsBunny.5hz.EF, tmpEF5) }
    
  }
  BugsBunny.1hz.EF = BugsBunny.1hz.EF[order(BugsBunny.1hz.EF$variable),]
  BugsBunny.1hz.EF$transect_source_fire_ID = indA[1]
  BugsBunny.5hz.EF = BugsBunny.5hz.EF[order(BugsBunny.5hz.EF$variable),]
  BugsBunny.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== CHARLIE BROWN =====================
  fire="Charlie Brown"; fuel="shrub"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 1; stop = stopO -5.6; startB = startO -1; stopB = startO + .8}
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1; startB=startB-1;stopB=stopB-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = 'shrub' ; tmp$fire[ind] = 'CharlieBrown'
    sep03rd.fire$fuel[ind2] = 'shrub' ; sep03rd.fire$fire[ind2] = 'CharlieBrown'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])  
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,'CharlieBrown', 'grass',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'CharlieBrown', 'grass',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop    
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){CharlieBrown.1hz.EF = tmpEF1 ;CharlieBrown.5hz.EF = tmpEF5 }
    if (i > 1){CharlieBrown.1hz.EF = rbind(CharlieBrown.1hz.EF, tmpEF1) ; CharlieBrown.5hz.EF = rbind(CharlieBrown.5hz.EF, tmpEF5) }
    
  }
  CharlieBrown.1hz.EF = CharlieBrown.1hz.EF[order(CharlieBrown.1hz.EF$variable),]
  CharlieBrown.1hz.EF$transect_source_fire_ID = indA[1]
  CharlieBrown.5hz.EF = CharlieBrown.5hz.EF[order(CharlieBrown.5hz.EF$variable),]
  CharlieBrown.5hz.EF$transect_source_fire_ID = indA[1]
  # ========================== SHAWNEE =====================
  fire="Shawnee"; fuel="rice"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 5; stop = stopO -46; startB = startO -1; stopB = startO + 0.4} # split?
    if (i == 2){start = startO + 1; stop = stopO -6; startB = startO -1; stopB = startO + 1} 
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i==2){start=start-1; stopB=stopB-1;startB=startB-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Shawnee'
    sep03rd.fire$fuel[ind2] = fuel ; sep03rd.fire$fire[ind2] = 'Shawnee'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])    
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,'Shawnee', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Shawnee', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop   
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Shawnee.1hz.EF = tmpEF1 ;Shawnee.5hz.EF = tmpEF5 }
    if (i > 1){Shawnee.1hz.EF = rbind(Shawnee.1hz.EF, tmpEF1) ; Shawnee.5hz.EF = rbind(Shawnee.5hz.EF, tmpEF5) }
    
  }
  Shawnee.1hz.EF = Shawnee.1hz.EF[order(Shawnee.1hz.EF$variable),]
  Shawnee.1hz.EF$transect_source_fire_ID = indA[1]
  Shawnee.5hz.EF = Shawnee.5hz.EF[order(Shawnee.5hz.EF$variable),]
  Shawnee.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== DAFFY DUCK =====================
  fire="Daffy Duck"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 1; stop = stopO -10; startB = startO -1; stopB = startO + 1}
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i==1){start=start-1; stopB=stopB-1;startB=startB-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'DaffyDuck'
    sep03rd.fire$fuel[ind2] = fuel ; sep03rd.fire$fire[ind2] = 'DaffyDuck'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])    
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,'DaffyDuck', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'DaffyDuck', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop   
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){DaffyDuck.1hz.EF = tmpEF1 ;DaffyDuck.5hz.EF = tmpEF5 }
    if (i > 1){DaffyDuck.1hz.EF = rbind(DaffyDuck.1hz.EF, tmpEF1) ; DaffyDuck.5hz.EF = rbind(DaffyDuck.5hz.EF, tmpEF5) }
    
  }
  DaffyDuck.1hz.EF = DaffyDuck.1hz.EF[order(DaffyDuck.1hz.EF$variable),]
  DaffyDuck.1hz.EF$transect_source_fire_ID = indA[1]
  DaffyDuck.5hz.EF = DaffyDuck.5hz.EF[order(DaffyDuck.5hz.EF$variable),]
  DaffyDuck.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== EEYORE =====================
  fire="Eeyore"; fuel="rice"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 0.4; stop = stopO -0; startB = startO -1; stopB = startO+1}
    if (i == 2){start = startO + 0.4; stop = stopO -5.8; startB = 71089; stopB = 71091} # split from above
    if (i == 3){start = startO + 5.8; stop = stopO -4.2; startB = 71089; stopB = 71091} # split from above
    if (i == 4){start = startO + 2.2; stop = stopO -9; startB = startO -1; stopB = startO+1}
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i==1){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i==2){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i==3){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i==4){start=start-2;startB=startB-1;stopB=stopB-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = 'rice' ; tmp$fire[ind] = 'Eeyore'
    sep03rd.fire$fuel[ind2] = 'rice' ; sep03rd.fire$fire[ind2] = 'Eeyore'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])   
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])   
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,'Eeyore', 'rice',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Eeyore', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)

    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop  
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Eeyore.1hz.EF = tmpEF1 ;Eeyore.5hz.EF = tmpEF5 }
    if (i > 1){Eeyore.1hz.EF = rbind(Eeyore.1hz.EF, tmpEF1) ; Eeyore.5hz.EF = rbind(Eeyore.5hz.EF, tmpEF5) }
    
  }
  Eeyore.1hz.EF = Eeyore.1hz.EF[order(Eeyore.1hz.EF$variable),]
  Eeyore.1hz.EF$transect_source_fire_ID = indA[1]
  Eeyore.5hz.EF = Eeyore.5hz.EF[order(Eeyore.5hz.EF$variable),]
  Eeyore.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== FAT ALBERT =====================
  fire="Fat Albert"; fuel="rice"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 0; stop = stopO -23; startB = startO -1; stopB = startO + .6}
    if (i == 2){start = startO + 10; stop = stopO -5; startB = startO + 0; stopB = startO + 2} # split?
    if (i == 3){start = startO + 0; stop = stopO -20; startB = startO -1; stopB = startO + 1}
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i==1){startB=startB-1;stopB=stopB-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = 'rice' ; tmp$fire[ind] = 'FatAlbert'
    sep03rd.fire$fuel[ind2] = 'rice' ; sep03rd.fire$fire[ind2] = 'FatAlbert'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]) 
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 300 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,'FatAlbert', 'rice',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'FatAlbert', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){FatAlbert.1hz.EF = tmpEF1 ;FatAlbert.5hz.EF = tmpEF5 }
    if (i > 1){FatAlbert.1hz.EF = rbind(FatAlbert.1hz.EF, tmpEF1) ; FatAlbert.5hz.EF = rbind(FatAlbert.5hz.EF, tmpEF5) }
    
  }
  FatAlbert.1hz.EF = FatAlbert.1hz.EF[order(FatAlbert.1hz.EF$variable),]
  FatAlbert.1hz.EF$transect_source_fire_ID = indA[1]
  FatAlbert.5hz.EF = FatAlbert.5hz.EF[order(FatAlbert.5hz.EF$variable),]
  FatAlbert.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== HOBBES =====================
  fire="Hobbes"; fuel="pile"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i))
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 11.8; stop = stopO -3.4; startB = startO -1; stopB = startO + 1}
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-1; stopB=stopB-1; startB=startB-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Hobbes'
    sep03rd.fire$fuel[ind2] = fuel ; sep03rd.fire$fire[ind2] = 'Hobbes'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])    
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,'Hobbes', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Hobbes', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
   tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop    
   tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
   #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
   #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
   
    if (i == 1){Hobbes.1hz.EF = tmpEF1 ;Hobbes.5hz.EF = tmpEF5 }
    if (i > 1){Hobbes.1hz.EF = rbind(Hobbes.1hz.EF, tmpEF1) ; Hobbes.5hz.EF = rbind(Hobbes.5hz.EF, tmpEF5) }
    
  }
  Hobbes.1hz.EF = Hobbes.1hz.EF[order(Hobbes.1hz.EF$variable),]
  Hobbes.1hz.EF$transect_source_fire_ID = indA[1]
  Hobbes.5hz.EF = Hobbes.5hz.EF[order(Hobbes.5hz.EF$variable),]
  Hobbes.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== GRINCH =====================
  fire="Grinch"; fuel="rice"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO +1.8; stop = stopO -5; startB = startO -1; stopB = startO +.8}
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i ==1){start=start-2; startB=startB-1; stopB=stopB-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = 'rice' ; tmp$fire[ind] = 'Grinch'
    sep03rd.fire$fuel[ind2] = 'rice' ; sep03rd.fire$fire[ind2] = 'Grinch'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,fire,fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire,fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop   
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Grinch.1hz.EF = tmpEF1 ;Grinch.5hz.EF = tmpEF5 }
    if (i > 1){Grinch.1hz.EF = rbind(Grinch.1hz.EF, tmpEF1) ; Grinch.5hz.EF = rbind(Grinch.5hz.EF, tmpEF5) }
    
  }
  Grinch.1hz.EF = Grinch.1hz.EF[order(Grinch.1hz.EF$variable),]
  Grinch.1hz.EF$transect_source_fire_ID = indA[1]
  Grinch.5hz.EF = Grinch.5hz.EF[order(Grinch.5hz.EF$variable),]
  Grinch.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== IAGO =====================
  fire="Iago"; fuel="pile"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 38; stop = stopO -11.4; startB = startO-1; stopB = startO+1}
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i==1){start=start-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Iago'
    sep03rd.fire$fuel[ind2] = fuel ; sep03rd.fire$fire[ind2] = 'Iago'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
     indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])      
     if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,'Iago', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Iago', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop    
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Iago.1hz.EF = tmpEF1 ;Iago.5hz.EF = tmpEF5 }
    if (i > 1){Iago.1hz.EF = rbind(Iago.1hz.EF, tmpEF1) ; Iago.5hz.EF = rbind(Iago.5hz.EF, tmpEF5) }
    
  }
  Iago.1hz.EF = Iago.1hz.EF[order(Iago.1hz.EF$variable),]
  Iago.1hz.EF$transect_source_fire_ID = indA[1]
  Iago.5hz.EF = Iago.5hz.EF[order(Iago.5hz.EF$variable),]
  Iago.5hz.EF$transect_source_fire_ID = indA[1]
  # ========================== JANE AND JUDY JETSON =====================
  fire="Jane and Judy Jetson"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 1; stop = stopO -41.4; startB = startO -1; stopB = startO + 1}
    if (i == 2){start = startO + 8.2; stop = stopO -1.8; startB = startO -1; stopB = startO +.6}
    if (i == 3){;start = startO + 1; stop = stopO -0; startB = startO -1; stopB = startO + 1} 
    if (i == 4){start = startO + 1; stop = stopO -1; startB = 73979; stopB = 73981} #split from above
    if (i == 5){start = startO + 0; stop = stopO -29; startB = 73979; stopB = 73981} #split from above
    if (i == 6){start =start+32; stop = stop - 15; startB = stopO-1.4; stopB = stopO} # missing data at the beginning
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i==1){start=start-1}
    if (i==2){start=start-2}
    if (i==3){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i==4){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i==5){start=start-1;startB=startB-1;stopB=stopB-1}
    if (i==6){start=start-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Jetson'
    sep03rd.fire$fuel[ind2] = fuel ; sep03rd.fire$fire[ind2] = 'Jetson'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])      
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])     
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,'Jetson', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'Jetson', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop   
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Jetson.1hz.EF = tmpEF1 ;Jetson.5hz.EF = tmpEF5 }
    if (i > 1){Jetson.1hz.EF = rbind(Jetson.1hz.EF, tmpEF1) ; Jetson.5hz.EF = rbind(Jetson.5hz.EF, tmpEF5) }
    
  }
  Jetson.1hz.EF = Jetson.1hz.EF[order(Jetson.1hz.EF$variable),]
  Jetson.1hz.EF$transect_source_fire_ID = indA[1]
  Jetson.5hz.EF = Jetson.5hz.EF[order(Jetson.5hz.EF$variable),]
  Jetson.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== KIM POSSIBLE =====================
  fire="Kim Possible"; fuel ="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){   
    print(c(fire,i))
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 4; stop = stopO -1; startB = stopO -1; stopB = stopO+1} #missing data at beginning
    if (i == 2){start = startO + 5; stop = stopO -62; startB = startO -1; stopB = startO+1} 
    if (i == 3){start = startO + 3; stop = stopO -6; startB = startO -1; stopB = startO+1}
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i==1){start=start-1}
    if (i==2){start=start-1}
    if (i==3){start=start-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'KimPossible'
    sep03rd.fire$fuel[ind2] = fuel ; sep03rd.fire$fire[ind2] = 'KimPossible'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1;peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      if (is.nan(GILMANBG$CO_DACOM_DISKIN_GILMAN[1]) | is.na(GILMANBG$CO_DACOM_DISKIN_GILMAN[1])){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,'KimPossible', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'KimPossible', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop  
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){KimPossible.1hz.EF = tmpEF1 ;KimPossible.5hz.EF = tmpEF5 }
    if (i > 1){KimPossible.1hz.EF = rbind(KimPossible.1hz.EF, tmpEF1) ; KimPossible.5hz.EF = rbind(KimPossible.5hz.EF, tmpEF5) }
    
  }
  KimPossible.1hz.EF = KimPossible.1hz.EF[order(KimPossible.1hz.EF$variable),]
  KimPossible.1hz.EF$transect_source_fire_ID = indA[1]
  KimPossible.5hz.EF = KimPossible.5hz.EF[order(KimPossible.5hz.EF$variable),]
  KimPossible.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== LISA SIMPSON =====================
  fire="Lisa Simpson"; fuel="rice"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 0; stop = stopO -4; startB = startO -1; stopB = startO + .6}
    if (i == 2){start = startO + 76.6; stop = stopO -6; startB = startO + 0; stopB =startO + 2}
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i==1){start=start-1; startB=startB-1; stopB=stopB-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = 'rice' ; tmp$fire[ind] = 'LisaSimpson'
    sep03rd.fire$fuel[ind2] = 'rice' ; sep03rd.fire$fire[ind2] = 'LisaSimpson'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]) 
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      while (is.nan(GILMANBG$CO_DACOM_DISKIN_GILMAN[1]) | is.na(GILMANBG$CO_DACOM_DISKIN_GILMAN[1])){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        print("HERE")
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (i == 1 | i == 2){indGILMANBACKGROUND = indGILMANBACKGROUND - 2} #kludge
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,]) 
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,'LisaSimpson', fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'LisaSimpson', fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop   
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){LisaSimpson.1hz.EF = tmpEF1 ;LisaSimpson.5hz.EF = tmpEF5 }
    if (i > 1){LisaSimpson.1hz.EF = rbind(LisaSimpson.1hz.EF, tmpEF1) ; LisaSimpson.5hz.EF = rbind(LisaSimpson.5hz.EF, tmpEF5) }
    
  }
  LisaSimpson.1hz.EF = LisaSimpson.1hz.EF[order(LisaSimpson.1hz.EF$variable),]
  LisaSimpson.1hz.EF$transect_source_fire_ID = indA[1]
  LisaSimpson.5hz.EF = LisaSimpson.5hz.EF[order(LisaSimpson.5hz.EF$variable),]
  LisaSimpson.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== UNKOWN KIM =====================
  fire="unknown Kim"; fuel="?"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 7; stop = stopO -29; startB = startO + 0; stopB =startO + 2}
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i==1){start=start-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = '?' ; tmp$fire[ind] = 'Unknown Kim'
    sep03rd.fire$fuel[ind2] = '?' ; sep03rd.fire$fire[ind2] = 'Unknown Kim'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])   
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      # go backwards till we get a background if necessary
      if (is.nan(GILMANBG$CO_DACOM_DISKIN_GILMAN[1]) | is.na(GILMANBG$CO_DACOM_DISKIN_GILMAN[1])){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,'unknownKim', '?',1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,'unknownKim', '?',5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop    
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){unknownKim.1hz.EF = tmpEF1 ;unknownKim.5hz.EF = tmpEF5 }
    if (i > 1){unknownKim.1hz.EF = rbind(unknownKim.1hz.EF, tmpEF1) ; unknownKim.5hz.EF = rbind(unknownKim.5hz.EF, tmpEF5) }
    
  }
  unknownKim.1hz.EF = unknownKim.1hz.EF[order(unknownKim.1hz.EF$variable),]
  unknownKim.1hz.EF$transect_source_fire_ID = indA[1]
  unknownKim.5hz.EF = unknownKim.5hz.EF[order(unknownKim.5hz.EF$variable),]
  unknownKim.5hz.EF$transect_source_fire_ID = indA[1]
  # ========================== MARGE =====================
  fire="Marge"; fuel = "rice"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i))
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 1; stop = stopO -5; startB = startO -1; stopB = startO +1}
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i==1){start=start-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = 'rice' ; tmp$fire[ind] = 'Marge'
    sep03rd.fire$fuel[ind2] = 'rice' ; sep03rd.fire$fire[ind2] = 'Marge'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      if (is.nan(GILMANBG$CO_DACOM_DISKIN_GILMAN[1]) | is.na(GILMANBG$CO_DACOM_DISKIN_GILMAN[1])){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }      
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,fire,fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire,fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop   
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Marge.1hz.EF = tmpEF1 ;Marge.5hz.EF = tmpEF5 }
    if (i > 1){Marge.1hz.EF = rbind(Marge.1hz.EF, tmpEF1) ; Marge.5hz.EF = rbind(Marge.5hz.EF, tmpEF5) }
    
  }
  Marge.1hz.EF = Marge.1hz.EF[order(Marge.1hz.EF$variable),]
  Marge.1hz.EF$transect_source_fire_ID = indA[1]
  Marge.5hz.EF = Marge.5hz.EF[order(Marge.5hz.EF$variable),]
  Marge.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== NEMO  =====================
  
  fire="Nemo"; fuel = "corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 2; stop = stopO -4; startB = startO +0; stopB = startO +2}
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i==1){start=start-1; startB=startB-1; stopB=stopB-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Nemo'
    sep03rd.fire$fuel[ind2] = fuel ; sep03rd.fire$fire[ind2] = 'Nemo'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])    
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      if (is.nan(GILMANBG$CO_DACOM_DISKIN_GILMAN[1]) | is.na(GILMANBG$CO_DACOM_DISKIN_GILMAN[1])){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] =fuel ; tmp$fire[ind] = fire
    sep03rd.fire$fuel[ind2] = fuel ; sep03rd.fire$fire[ind2] = fire
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,fire,fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire,fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop    
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Nemo.1hz.EF = tmpEF1 ;Nemo.5hz.EF = tmpEF5 }
    if (i > 1){Nemo.1hz.EF = rbind(Nemo.1hz.EF, tmpEF1) ; Nemo.5hz.EF = rbind(Nemo.5hz.EF, tmpEF5) }
    
  }
  Nemo.1hz.EF = Nemo.1hz.EF[order(Nemo.1hz.EF$variable),]
  Nemo.1hz.EF$transect_source_fire_ID = indA[1]
  Nemo.5hz.EF = Nemo.5hz.EF[order(Nemo.5hz.EF$variable),]
  Nemo.5hz.EF$transect_source_fire_ID = indA[1]
  # ========================== OBELIX =====================
  fire="Obelix"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 1; stop = stopO -0; startB = startO +0; stopB = startO +1.6}
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i==1){start=start-1; startB=startB-1; stopB=stopB-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = fuel; tmp$fire[ind] = 'Obelix'
    sep03rd.fire$fuel[ind2] =fuel ; sep03rd.fire$fire[ind2] = 'Obelix'   
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])    
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      if (is.nan(GILMANBG$CO_DACOM_DISKIN_GILMAN[1]) | is.na(GILMANBG$CO_DACOM_DISKIN_GILMAN[1])){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,fire,fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire,fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop    
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Obelix.1hz.EF = tmpEF1 ;Obelix.5hz.EF = tmpEF5 }
    if (i > 1){Obelix.1hz.EF = rbind(Obelix.1hz.EF, tmpEF1) ; Obelix.5hz.EF = rbind(Obelix.5hz.EF, tmpEF5) }
    
  }
  Obelix.1hz.EF = Obelix.1hz.EF[order(Obelix.1hz.EF$variable),]
  Obelix.1hz.EF$transect_source_fire_ID = indA[1]
  Obelix.5hz.EF = Obelix.5hz.EF[order(Obelix.5hz.EF$variable),]
  Obelix.5hz.EF$transect_source_fire_ID = indA[1]
  # ========================== POPEYE  =====================
  fire="Popeye"; fuel='rice'
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 2; stop = stopO -3; startB = startO +0; stopB = startO +2}
    if (i == 2){start = startO + 5; stop = stopO -0; startB = 76376; stopB = 76378}#split from above
    if (i == 3){start = startO + 2; stop = stopO -6.8; startB = 76376; stopB = 76378}#split from above
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    
    if (i==1){start=start-1; startB=startB-1; stopB=stopB-1}
    if (i==2){start=start-1; startB=startB-1; stopB=stopB-1}
    if (i==3){start=start-1; startB=startB-1; stopB=stopB-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] =fuel; tmp$fire[ind] = 'Popeye'
    sep03rd.fire$fuel[ind2] =fuel ; sep03rd.fire$fire[ind2] = 'Popeye'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])  
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      if (is.nan(GILMANBG$CO_DACOM_DISKIN_GILMAN[1]) | is.na(GILMANBG$CO_DACOM_DISKIN_GILMAN[1])){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,fire,fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire,fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop    
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Popeye.1hz.EF = tmpEF1 ;Popeye.5hz.EF = tmpEF5 }
    if (i > 1){Popeye.1hz.EF = rbind(Popeye.1hz.EF, tmpEF1) ; Popeye.5hz.EF = rbind(Popeye.5hz.EF, tmpEF5) }
    
  }
  Popeye.1hz.EF = Popeye.1hz.EF[order(Popeye.1hz.EF$variable),]
  Popeye.1hz.EF$transect_source_fire_ID = indA[1]
  Popeye.5hz.EF = Popeye.5hz.EF[order(Popeye.5hz.EF$variable),]
  Popeye.5hz.EF$transect_source_fire_ID = indA[1]
  # ========================== ROADRUNNER  =====================
  fire="Roadrunner"; fuel="corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 0.6; stop = stopO -15; startB = startO -1; stopB = startO + .6}
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i==1){start=start-2; stop=stop+1;startB=startB-1; stopB=stopB-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Roadrunner'
    sep03rd.fire$fuel[ind2] = fuel ; sep03rd.fire$fire[ind2] = 'Roadrunner'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]]) 
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      if (is.nan(GILMANBG$CO_DACOM_DISKIN_GILMAN[1]) | is.na(GILMANBG$CO_DACOM_DISKIN_GILMAN[1])){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hzJUSTCN(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,fire,fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire,fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop 
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Roadrunner.1hz.EF = tmpEF1 ;Roadrunner.5hz.EF = tmpEF5 }
    if (i > 1){Roadrunner.1hz.EF = rbind(Roadrunner.1hz.EF, tmpEF1) ; Roadrunner.5hz.EF = rbind(Roadrunner.5hz.EF, tmpEF5) }
    
  }
  Roadrunner.1hz.EF = Roadrunner.1hz.EF[order(Roadrunner.1hz.EF$variable),]
  Roadrunner.1hz.EF$transect_source_fire_ID = indA[1]
  Roadrunner.5hz.EF = Roadrunner.5hz.EF[order(Roadrunner.5hz.EF$variable),]
  Roadrunner.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ========================== ROADRUNNER SPONGEBOB  =====================

  fire="Spongebob"; fuel = "corn"
  indA = which(flags$transec_source_fire_namestr == fire)
  for (i in 1:length(indA)){
    print(c(fire,i) )
    startO = flags$`transect_start_time (UTC s from midnight)`[indA[i]] ; start = startO 
    stopO  = flags$`transect_end_time  (UTC s from midnight)`[indA[i]]; stop = stopO
    if (i == 1){start = startO + 6; stop = stopO -1; startB = startO + 1; stopB = startO + 2}
    if (i == 2){start = startO + 0.8; stop = stopO -74.8; startB = 76591; stopB =76592} #split from above
    tmp = time_align(start,stop,co2.903.5hz,  co.ch4.903.5hz, 
                     warneke.903.5hz,  isaf.903.5hz,
                     rollinsno.903.5hz, rollinsso2.903.5hz, 
                     cit.903.5hz, gtcims.903.5hz,moore.903fast, jimenez.903.5hz,met.903.5hz)
    ind = which(tmp$Time_Start >= start & tmp$Time_Start  <= stop	) 
    indB = which(tmp$Time_Start >= startB & tmp$Time_Start  <= stopB	) 
    if (i==1){start=start-1}
    if (i==2){start=start-1}
    sep03rd.fire= time_alignSLOWNOCANS(start,stop, co2.903.1hz, co.ch4.903.1hz, 
                                       warneke.903.1hz, isaf.903.1hz,rollinsno.903.1hz, rollinsso2.903.1hz, 
                                       cit.903.1hz,gtcims.903.1hz,ryerson.903.1hz , jimenez.903.1hz, 
                                       schwarz.903.1hz, freid.c2h6.903.1hz, freid.ch2o.903.1hz, 
                                       womack.903.1hz, stclair.903.1hz,veres.903.1hz, #wisthaler.903.1hz,
                                       moore.903, met.903.1hz)
    ind2 = which(sep03rd.fire$Time_Start >= start & sep03rd.fire$Time_Start  <= stop	) 
    ind2B = which(sep03rd.fire$Time_Start >= startB & sep03rd.fire$Time_Start <= stopB	) 
    # collect all ag plumes together at 5hz and 1hz
    tmp$fire = ''; tmp$fuel = ''
    sep03rd.fire$fire = '';sep03rd.fire$fuel = ''
    tmp$fuel[ind] = fuel ; tmp$fire[ind] = 'Spongebob'
    sep03rd.fire$fuel[ind2] = fuel ; sep03rd.fire$fire[ind2] = 'Spongebob'
    # -------- Blake cans 
    indBLAKE = which(blake.903.merge$Time_Start >= (start-5) & blake.903.merge$Time_Start <= stop) # check at least 5 sec before
    if (length(indBLAKE) > 0){
      doBlake = 1
      indBLAKEBACKGROUND = indBLAKE[1] - 1
      blake = blake.903.merge[indBLAKE,]
      blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      # go backwards till we get a background if necessary
      while (blakeBG$CO_DACOM_DISKIN_BLAKE[1] > 200 ){
        indBLAKEBACKGROUND = indBLAKEBACKGROUND - 1
        blakeBG = blake.903.merge[indBLAKEBACKGROUND,]
      }
    } else{
      doBlake = 0
      # just dummy variables so things dont break
      blake = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
      blakeBG = blake.903.merge[1,]; blake$CO_DACOM_DISKIN_BLAKE=NaN
    }
    
    # Gilman cans
    indGILMAN = which(gilman.903.merge$Time_Start >= (start-0) & gilman.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indGILMAN) > 0){
      doGILMAN = 1
      peakCO = which(tmp$CO_DACOM_DISKIN[ind] == max(tmp$CO_DACOM_DISKIN[ind], na.rm=TRUE))
      indGILMANBACKGROUND = indGILMAN[1] - 1
      GILMAN = gilman.903.merge[indGILMAN,]
      GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      indGG = which(GILMAN$Time_Start > tmp$Time_Start[ind[peakCO]])   
      if (length(indGG) == 0){ doGILMAN = 0} # don't use GILMAN if didn't capture peak       
      if (is.nan(GILMANBG$CO_DACOM_DISKIN_GILMAN[1]) | is.na(GILMANBG$CO_DACOM_DISKIN_GILMAN[1])){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
      # go backwards till we get a background if necessary
      while (GILMANBG$CO_DACOM_DISKIN_GILMAN[1] > 200 ){
        indGILMANBACKGROUND = indGILMANBACKGROUND - 1
        if (indGILMANBACKGROUND < 1){ # start from the end
          indGILMANBACKGROUND = length(gilman.903.merge$CO_DACOM_DISKIN_GILMAN)
        }
        GILMANBG = gilman.903.merge[indGILMANBACKGROUND,]
      }
    } else{
      doGILMAN = 0
      # just dummy variables so things dont break
      GILMAN = gilman.903.merge[1,]; GILMAN$CO_DACOM_DISKIN_GILMAN=  NaN
      GILMANBG = gilman.903.merge[1,]
    }
    
    # Apel  
    indAPEL = which(apel.903.merge$Time_Start >= (start-0) & apel.903.merge$Time_Start <= stop) # check at least 5 sec before?
    if (length(indAPEL) > 0){
      doAPEL = 1
      indAPELBACKGROUND = indAPEL[1] - 1
      APEL = apel.903.merge[indAPEL,]
      APELBG = apel.903.merge[indAPELBACKGROUND,]
      # go backwards till we get a background if necessary
      while (APELBG$CO_DACOM_DISKIN_APEL[1] > 200 ){
        indAPELBACKGROUND = indAPELBACKGROUND - 1
        if (indAPELBACKGROUND < 1){ # start from the end
          indAPELBACKGROUND = length(apel.903.merge$CO_DACOM_DISKIN_APEL)
        }
        APELBG = apel.903.merge[indAPELBACKGROUND,]
      }
    } else{
      doAPEL = 0
      # just dummy variables so things dont break
      APEL = apel.903.merge[1,]
      APELBG = apel.903.merge[1,]
    }
    # Becky
    indBECKY = which(newTOGA.903$Time_Start >= (start-0) & newTOGA.903$Time_Start <= stop & newTOGA.903$Percent_plume > 0.5) # check at least 5 sec before?
    if (length(indBECKY) > 0){
      doBECKY = 1
      indBECKYBACKGROUND = indBECKY[1] - 1
      BECKY = newTOGA.903[indBECKY,]
      BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      # go backwards till we get a background if necessary
      while (BECKYBG$CO_DACOM_DISKIN_BECKY[1] > 200 ){
        indBECKYBACKGROUND = indBECKYBACKGROUND - 1
        if (indBECKYBACKGROUND < 1){ # start from the end
          indBECKYBACKGROUND = length(newTOGA.903$CO_DACOM_DISKIN_BECKY)
        }
        BECKYBG = newTOGA.903[indBECKYBACKGROUND,]
      }
    } else{
      doBECKY = 0
      # just dummy variables so things dont break
      BECKY = newTOGA.903[1,]; BECKY$CO_DACOM_DISKIN_BECKY = NaN
      BECKYBG = newTOGA.903[1,]; BECKYBG$CO_DACOM_DISKIN_BECKY=NaN
    }
    if (doplot == 1){
      plotpass5hzJUSTCO(sep03rd.fire[ind2,])
      plotpass5hz(sep03rd.fire[ind2,])
      plotpass1hz(sep03rd.fire[ind2,])
      plotpass5hzJUSTCO(tmp[ind,])
      plotpass5hz(tmp[ind,])
      # plot background?
      plotpass5hzJUSTCO(sep03rd.fire[ind2B,])
      plotpass5hzJUSTCO(tmp[indB,])
    }
    allfires.5hz = rbind.fill(allfires.5hz, tmp[ind,]); BECKY$fire = fire;BECKY$fuel=fuel;BECKY$pass = i; GILMAN$fire = fire;GILMAN$fuel=fuel;GILMAN$pass = i;blake$fire = fire;blake$fuel=fuel;blake$pass = i     
    allfires.1hz = rbind.fill(allfires.1hz, sep03rd.fire[ind2,]); if (doBECKY == 1){toga.all=rbind(toga.all,BECKY)};  if (doBlake == 1){was.all = rbind(was.all, blake)};  if (doGILMAN == 1){iwas.all = rbind(iwas.all, GILMAN)}
    
    # ------------------ 1Hz
    tmpEF1 = ERsEFsALLhzv2(sep03rd.fire[ind2,],sep03rd.fire[ind2B,],xspecies,fire,fuel,1,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF1$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE); tmpEF1$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    tmpEF1$transect_fuel_class = median(sep03rd.fire$transect_fuel_class[ind2], na.rm=TRUE);  tmpEF1$transect_dominant_fuel = median(sep03rd.fire$transect_dominant_fuel[ind2], na.rm=TRUE); tmpEF1$transect_fuel_confidence = median(sep03rd.fire$transect_fuel_confidence[ind2], na.rm=TRUE)
    # ------------------ 5Hz
    tmpEF5  = ERsEFsALLhzv2(tmp[ind,],tmp[indB,],xspecies,fire,fuel,5,SLOW,blake, blakeBG,doBlake,GILMAN,GILMANBG,doGILMAN,BECKY, BECKYBG, doBECKY)
    tmpEF5$age = median(sep03rd.fire$transect_smoke_age[ind2], na.rm=TRUE) ; tmpEF5$transect_type = median(sep03rd.fire$transect_type[ind2], na.rm=TRUE)
    
    tmpEF1$pass = i ; tmpEF1$StartO = startO; tmpEF1$Start = start; tmpEF1$StopO = stopO; tmpEF1$Stop = stop 
    tmpEF5$pass = i ; tmpEF5$StartO = startO; tmpEF5$Start = start; tmpEF5$StopO = stopO; tmpEF5$Stop = stop     
    #tmpEF5$fuelMDead = fuelMDead[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    #tmpEF5$fuelMLive = fuelMLive[as.numeric(flags$II[indA[i]]),as.numeric(flags$JJ[indA[i]])]
    
    if (i == 1){Spongebob.1hz.EF = tmpEF1 ;Spongebob.5hz.EF = tmpEF5 }
    if (i > 1){Spongebob.1hz.EF = rbind(Spongebob.1hz.EF, tmpEF1) ; Spongebob.5hz.EF = rbind(Spongebob.5hz.EF, tmpEF5) }
    
  }
  Spongebob.1hz.EF = Spongebob.1hz.EF[order(Spongebob.1hz.EF$variable),]
  Spongebob.1hz.EF$transect_source_fire_ID = indA[1]
  Spongebob.5hz.EF = Spongebob.5hz.EF[order(Spongebob.5hz.EF$variable),]
  Spongebob.5hz.EF$transect_source_fire_ID = indA[1]
  
  # ---- Save data for later analysis --------
  save.image(file='AgFires.RData')
  
 } else( load(file='AgFires.RData'))
