# -------- ### All  EFs ### -------- 
doclear=1
if (doclear == 1){ rm(list=ls()) } # clear all analysis
load('AgFires.RData')
 source('speciateSpecies.R')
require(dplyr); require(plyr); require(GMCM)

R2filter = 0.70; R2filterCO = 0.90 # stricter criteria for CO as this defines the plume
R2Bot = 0.5 ; COcutoff = 400 # ppb 
doprocess=0;doprocessSTEP2 =0
OHval = 5E6
lifetimecutoff = 6 #hours for short vs. long-lived
# ----- Some ggplot settings ----
fuelshapes = c(19,15,2,18,7,8,4,2,16,12,13)
fuellimits = c("corn","soybean","rice","winter wheat","grass","pile","slash","shrub","forest")

cbp1a <- c( "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", 
            "#CC79A7","#000000","#999999","#CC0000")
# colorblind safe option
cbp1b = c('#377eb8','#e41a1c','#4daf4a','#f781bf','#a65628','#984ea3','#ff7f00','#ffff33')
cbp1c=c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9',
        '#74add1','#4575b4','#313695')
cbp1d <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbp1 = cbp1b
# ----- +++++++++++ ANALYSIS +++++++++++++ ------
# How do emission factors change at a single fire over time, if at all?
# How do emission factors vary across various fires sampled?
#  Are there metrics (i.e. MCE etc) that provide some explanatory fpower for the variation that may be present?
require(reshape); require(ggmap) ; require(OrgMassSpecR); library(readxl) ;require(plyr)# ; library(ggbrace)
require(dplyr); require(ggpubr); require(ncdf4)
# ---- EFs from Xiaoxi Liu (2016 ACP SEAC4RS rice straw)
source("xioaxi.R")
xiaoxi.avg=as.data.frame(colMeans(xiaoxi, na.rm=TRUE))
xiaoxi.sd=apply(xiaoxi, 2, sd)

xiaoxi.avg$var =rownames(xiaoxi.avg)
colnames(xiaoxi.avg)  = c('mean','var')
xiaoxi.avg$sd = xiaoxi.sd
# --------- Get Andreae emission factors ------------
f2 = 'InputFiles/OtherStudies/Andreae-BB-EMFactors-14Apr2021_justtable1.csv'
andreae = read.csv(f2)
# --------- Get Akagi emission factors ---------
akagi=readxl::read_xlsx('InputFiles/OtherStudies/Akagi_acp-11-4039-2011-supplement/Tables 1-5_4.27.11.xlsx')
akagi$CropEF = as.numeric(akagi$CropEF)
akagi$SavannahEF= as.numeric(akagi$Savannah)
akagi$TemperateEF= as.numeric(akagi$Temperate)

if (doprocess == 1){
  # ------ combine 5hz & 1Hz fires ----
  all5hz = rbind.fill(CopperBreaks.5hz.EF,Vivian.5hz.EF,Halfpint.5hz.EF,Loretta.5hz.EF,LilDebbie.5hz.EF,Ricearoni.5hz.EF,CrawbabyHouse.5hz.EF,Crawdaddy.5hz.EF, Gumbo.5hz.EF,Jumbalaya.5hz.EF,JumbalayaJr.1hz.EF,Crouton.5hz.EF,PoBoy.5hz.EF, # aug 21
                      Ant.5hz.EF, Blanket.5hz.EF, Chips.5hz.EF, Dip.5hz.EF, Escargot.5hz.EF,Frisbee.5hz.EF, Guac.5hz.EF, Hamburger.5hz.EF, IPA.5hz.EF, Jello.5hz.EF, Kebab.5hz.EF, 
                      Limoncello.5hz.EF, Mustard.5hz.EF, # aug 23
                      DaysLater.5hz.EF,Alien.5hz.EF,Bambi.5hz.EF,BambiJr.5hz.EF,Deadpool.5hz.EF,Elf.5hz.EF,Fargo.5hz.EF,Hellboy.5hz.EF,Invictus.5hz.EF,InvictusU.5hz.EF,# aug 26
                      HickoryRidge.5hz.EF, Tallgrass.5hz.EF,Boxer.5hz.EF,Chessie.5hz.EF, Dingo.5hz.EF, Elkhound.5hz.EF, # aug 29th,
                      Blackwaterriver.5hz.EF, WIGGINS.5hz.EF, WIGGINSNEIGHBORS.5hz.EF, ZZTop.5hz.EF, YoungMC.5hz.EF, XTC.5hz.EF, Weezer.5hz.EF, ViolentFemmes.5hz.EF, U2.5hz.EF, Toto.5hz.EF, Supertramp.5hz.EF, # Aug 30
                      Jaws.5hz.EF, Kingpin.5hz.EF, Leon.5hz.EF, Meatballs.5hz.EF, Nikita.5hz.EF, Oblivion.5hz.EF, OblivionJr.5hz.EF, Psycho.5hz.EF, Quarantine.5hz.EF, Ratatouille.5hz.EF, Spaceballs.5hz.EF,Tremors.5hz.EF, Up.5hz.EF, Vertigo.5hz.EF, Willow.5hz.EF,# aug 31st
                      Asterisk.5hz.EF, BugsBunny.5hz.EF, CharlieBrown.5hz.EF, Shawnee.5hz.EF, DaffyDuck.5hz.EF, Eeyore.5hz.EF, FatAlbert.5hz.EF, Grinch.5hz.EF, Hobbes.5hz.EF, Iago.5hz.EF, Jetson.5hz.EF, KimPossible.5hz.EF, LisaSimpson.5hz.EF, unknownKim.5hz.EF, Marge.5hz.EF, Nemo.5hz.EF,Obelix.5hz.EF,Popeye.5hz.EF,Roadrunner.5hz.EF,Spongebob.5hz.EF)# Aug 31 Akito was not a good fire, W-V & XTC-WEEZER-UNKNOWN bad fire
  all5hz = all5hz[order(all5hz$variable),]
  all5hz$uniqueid =(all5hz$pass/10+all5hz$transect_source_fire_ID)
  #ind = which(all5hz$variable == 'CO_DACOM_DISKIN')
  #tnd2 = which(duplicated(all5hz$uniqueid[ind]))
  #all5hz$uniqueid[ind[tnd2]]
  
  # ---- make 5hz VOC EF ------- this is including duplicates, so only for comparison to MCE, shouldn't be used as a total
  #ind = which(all5hz$Category == 1)
  #all5hz.VOC = all5hz[ind,]
  #test=aggregate(all5hz.VOC$EF1, by=list(all5hz.VOC$uniqueid), FUN=sum, na.rm=TRUE)
  #ind = which(test$x == 0) ; test$x[ind] = NaN # if there are zero VOCs for a fire 
  #test2=aggregate(all5hz.VOC, by=list(all5hz.VOC$uniqueid), FUN='mean', na.rm=TRUE)
  ## don't want a zero EF
  #cc = colnames(all5hz)
  #fullline = data.frame(matrix(vector(), length(test$Group.1), length(cc)))
  #colnames(fullline) = cc
  #fullline$uniqueid = test$Group.1
  #fullline$EF1 = test$x
  #fullline$variable = 'All5HzVOC'
  #fullline$mce = test2$mce
  #fullline$Category=1
  #fullline$age = test2$age
  #fullline$R2toX = test2$R2toX
  #for (i in 1:length(fullline$fuel)){
  #  ind = which(fullline$uniqueid[i] == all5hz.VOC$uniqueid)
  #  fullline$fuel[i] = all5hz.VOC$fuel[ind[1]]
  #  fullline$fire[i] = all5hz.VOC$fire[ind[1]]
  #}
  #all5hz = rbind(all5hz, fullline)
  # ----------------------------------------
  all1hz = rbind.fill(CopperBreaks.1hz.EF,Vivian.1hz.EF,Halfpint.1hz.EF,Loretta.1hz.EF,LilDebbie.1hz.EF,Ricearoni.1hz.EF,CrawbabyHouse.1hz.EF,Crawdaddy.1hz.EF, Gumbo.1hz.EF,
                      Jumbalaya.1hz.EF,JumbalayaJr.1hz.EF,Crouton.1hz.EF,PoBoy.1hz.EF, # aug 21
                      Ant.1hz.EF, Blanket.1hz.EF, Chips.1hz.EF, Dip.1hz.EF, Escargot.1hz.EF,Frisbee.1hz.EF, Guac.1hz.EF, Hamburger.1hz.EF, IPA.1hz.EF, Jello.1hz.EF, Kebab.1hz.EF,
                      Limoncello.1hz.EF, Mustard.1hz.EF, # aug 23
                      DaysLater.1hz.EF,Alien.1hz.EF,Bambi.1hz.EF,BambiJr.1hz.EF,Deadpool.1hz.EF,Elf.1hz.EF,Fargo.1hz.EF,Hellboy.1hz.EF,Invictus.1hz.EF,InvictusU.1hz.EF,# aug 26
                      HickoryRidge.1hz.EF, Tallgrass.1hz.EF,Boxer.1hz.EF,Chessie.1hz.EF, Dingo.1hz.EF, Elkhound.1hz.EF, # aug 29th,
                      Blackwaterriver.1hz.EF, WIGGINS.1hz.EF, WIGGINSNEIGHBORS.1hz.EF, ZZTop.1hz.EF, YoungMC.1hz.EF, XTC.1hz.EF, Weezer.1hz.EF, ViolentFemmes.1hz.EF, U2.1hz.EF, Toto.1hz.EF, Supertramp.1hz.EF, # Aug 30
                      Jaws.1hz.EF, Kingpin.1hz.EF, Leon.1hz.EF, Meatballs.1hz.EF, Nikita.1hz.EF, Oblivion.1hz.EF, OblivionJr.1hz.EF,Psycho.1hz.EF, Quarantine.1hz.EF, Ratatouille.1hz.EF, Spaceballs.1hz.EF, Tremors.1hz.EF,Up.1hz.EF, Vertigo.1hz.EF, Willow.1hz.EF, # Aug 31st
                      Asterisk.1hz.EF, BugsBunny.1hz.EF, CharlieBrown.1hz.EF, Shawnee.1hz.EF, DaffyDuck.1hz.EF, Eeyore.1hz.EF, FatAlbert.1hz.EF, Grinch.1hz.EF, Hobbes.1hz.EF, Iago.1hz.EF, Jetson.1hz.EF, KimPossible.1hz.EF, LisaSimpson.1hz.EF, unknownKim.1hz.EF, Marge.1hz.EF,  Nemo.1hz.EF,Obelix.1hz.EF,Popeye.1hz.EF, Roadrunner.1hz.EF, Spongebob.1hz.EF)# Aug 31 Akito was not a good fire, W-V & XTC-WEEZER-UNKNOWN bad fire
  all1hz = all1hz[order(all1hz$variable),]
  
  
  # -----  get rid of negatives ------ eventually might want to try and fix it, this is TOGA, GILMAN, BLAKE
  ind = which(all1hz$ERtoX < 0)
  all1hz$EF1[ind] = NaN
  all1hz$ERtoX[ind] = NaN
  ind = which(all1hz$ERtoCO < 0)
  toinvestigate = all1hz[ind,]
  all1hz$ERtoCO[ind] = NaN
  all1hz$EF1CO[ind] = NaN
  
  all1hz$uniqueid =(all1hz$pass/10+all1hz$transect_source_fire_ID)
  
  # remove -888, Inf, negatives???
  LOD = which(all1hz$maxval == -888 | all1hz$maxval == -0.888 | all1hz$maxval == -Inf | all1hz$maxval < 0)
  all1hz$maxval[LOD] = NaN
  tt = unique(all5hz$fuel)

  # --- Plot fires analyzed ------
  # Investigate O3/CO
  # ind = which(allfires.1hz$fuel2 != '?' &allfires.1hz$fuel2 != '' & 
  #               allfires.1hz$fuel2 != 'forest'  & 
  #               allfires.1hz$fuel2 != 'savannah'  & 
  #               allfires.1hz$fuel2 != 'timber'  & 
  #               allfires.1hz$fuel2 != 'Understory mixed, shrub,rice' &
  #               allfires.1hz$fuel2 != 'coniferous/decidous')
  # O31=ggplot(allfires.1hz[ind,])+geom_point(aes(x=CO_DACOM_DISKIN, y=O3_CL_RYERSON, col=fuel2))+theme_classic()
  # ind = which(allfires.1hz$fuel2 != '?' &allfires.1hz$fuel2 != '' & 
  #               allfires.1hz$fuel2 != 'forest'  & 
  #               allfires.1hz$fuel2 != 'savannah'  & 
  #               allfires.1hz$fuel2 != 'timber'  & 
  #               allfires.1hz$fuel2 != 'Understory mixed, shrub,rice' &
  #               allfires.1hz$fuel2 != 'coniferous/decidous' & is.finite(allfires.1hz$NO_LIF_ROLLINS) &
  #               allfires.1hz$CO_DACOM_DISKIN < 1E3)
  # O32=ggplot(allfires.1hz[ind,])+geom_point(aes(x=CO_DACOM_DISKIN, y=O3_CL_RYERSON, col=fire))+theme_classic()
  # O32=ggplot(allfires.1hz[ind,])+geom_point(aes(x=CO_DACOM_DISKIN,y=O3_CL_RYERSON, col=NO_LIF_ROLLINS/1E3))+theme_classic()+
  #   scale_color_viridis_c(trans='log10')
  
  # merge names, flags, and data so I can plot the ag fires I actually analyzed
  # should probably cut though to criteria for both

  # ------- How many fuels? --------
  #fuels = unique(all5hz.map$fuel)
  #length(unique(all5hz.map$fire))
  
  ind = which(all1hz$Category == 1 & is.finite(as.numeric(all1hz$ERtoCO)) & as.numeric(all1hz$ERtoCO) > 0
              & as.numeric(all1hz$R2toCO) > 0.75)
  tmp = unique(all1hz$variable[ind])
  newlist = c()
  for (i in 1:length(tmp)){
    tt = strsplit(tmp[i], "_")
    tt = tt[[1]][1]
    newlist = c(newlist,tt)
  }
  uniqueLIST = unique(newlist)
  unique(all1hz$variable[ind])
  
  # ------ Merge 1hz and 5hz together ------
  allBOTH = merge(all5hz, all1hz, by=c('variable','fire','fuel', 'transect_source_fire_ID','mWs',
                                       'nCs',"Start","Stop","StartO","StopO","pass","uniqueid"),
                  all = TRUE, suffixes = c(".5hz", ".1hz"))#, incomparables = NA) # x = 1Hz, y=5hz
  allBOTH$lifetime_5hz_hr = 1/(OHval*as.numeric(allBOTH$OHrate.5hz))/60/60
  allBOTH$lifetime_1hz_hr = 1/(OHval*as.numeric(allBOTH$OHrate.1hz))/60/60
  
  allBOTH$PI = ''
  for (i in 1:length(allBOTH$PI)){
    tmp = strsplit(allBOTH$variable[i], '_')
    tt = tmp[[1]]
    allBOTH$PI[i] = tt[length(tt)]
  }
  ind = which(allBOTH$PI != "APEL")
  allBOTH = allBOTH[ind,]
  # ----------- Provide the 5hz mce, MAtoF, and 5hz CO+CH4+CO2 ERtoCO  for the 1hz data --------------------------
  ff = unique(allBOTH$uniqueid)
  allBOTH$MCE = NaN
  for (i in 1:length(ff)){
    ind = which(allBOTH$uniqueid == ff[i] )
    allBOTH$MCE[ind] = max(as.numeric(allBOTH$mce.5hz[ind]), na.rm=TRUE) # they should all be the same, just fill in
    allBOTH$MAtoF.5hz[ind] = max(as.numeric(allBOTH$MAtoF.5hz[ind]), na.rm=TRUE) # they should all be the same, just fill in
    allBOTH$TC1.5hz[ind] = max(as.numeric(allBOTH$TC1.5hz[ind]), na.rm=TRUE) # they should all be the same, just fill in)
    allBOTH$TC1CO.5hz[ind] = max(as.numeric(allBOTH$TC1CO.5hz[ind]), na.rm=TRUE) # they should all be the same, just fill in)
  }
  # Do my 1hz EFs correlate better with 5Hz if I use 5hz TC?
  #ind = which(allBOTH$variable == 'Benzene_NOAAPTR_ppbv_WARNEKE')
  #plot(allBOTH$EF1.5hz[ind],allBOTH$EF1.1hz[ind])
  #cContentCorn = 500 # Assuming 50 % C, 500 g/kg
  #allBOTH$EF1CO.1hz.5hz = NaN
  #allBOTH$EF1.1hz.5hz = NaN
  #allBOTH$EF1CO.1hz.5hz =  cContentCorn * allBOTH$mWs/12 * allBOTH$ERtoCO.1hz /allBOTH$TC1CO.5hz
  #allBOTH$EF1.1hz.5hz =  cContentCorn * allBOTH$mWs/12 * allBOTH$ERtoX.1hz /allBOTH$TC1.5hz
  #ind = which(allBOTH$variable == 'Benzene_NOAAPTR_ppbv_WARNEKE')
  #points(allBOTH$EF1.5hz[ind],allBOTH$EF1.1hz.5hz[ind], col='red' ) # barely improves#

  # ---- Cut data to > 400 ppb where CO to CO2 R2 > 0.9
  dofilter=1
  if (dofilter==1){
    ind = which(allBOTH$variable == 'CO_DACOM_DISKIN' & 
                  as.numeric(allBOTH$R2toX.5hz) >= R2filterCO & as.numeric(allBOTH$maxval.5hz) > COcutoff) 
    goodpasses = allBOTH$uniqueid[ind]
    
    cc = c()
    for (i in 1:length(goodpasses)){
      ind = which(allBOTH$uniqueid == goodpasses[i])
      cc = c(cc, ind)
    }
    allBOTH.filter = allBOTH[cc,]
  
  }  
  allBOTH.filter$OHrate.1hz=as.numeric(allBOTH.filter$OHrate.1hz)
  allBOTH.filter$OHrate.5hz=as.numeric(allBOTH.filter$OHrate.5hz)
  allBOTH.filter$maxval.5hz=as.numeric(allBOTH.filter$maxval.5hz)
  allBOTH.filter$maxval.1hz=as.numeric(allBOTH.filter$maxval.1hz)
  
  allBOTH.filter$Category.5hz= as.numeric(allBOTH.filter$Category.5hz)
  allBOTH.filter$Category.1hz= as.numeric(allBOTH.filter$Category.1hz)
  allBOTH.filter$BGSpecies.5hz= as.numeric(allBOTH.filter$BGSpecies.5hz)
  allBOTH.filter$BGSpecies.1hz= as.numeric(allBOTH.filter$BGSpecies.1hz)
  allBOTH.filter$BGX.1hz = as.numeric(allBOTH.filter$BGX.1hz)
  allBOTH.filter$BGX.5hz = as.numeric(allBOTH.filter$BGX.5hz)
  allBOTH.filter$BGCO.1hz = as.numeric(allBOTH.filter$BGCO.1hz)
  allBOTH.filter$BGCO.5hz = as.numeric(allBOTH.filter$BGCO.5hz)
  allBOTH.filter$intercept.1hz = as.numeric(allBOTH.filter$intercept.1hz)
  allBOTH.filter$intercept.5hz = as.numeric(allBOTH.filter$intercept.5hz)
  
  # Provide same MCE to all passes from 5hz
  ff = unique(allBOTH.filter$uniqueid) # need to redo this
  allBOTH.filter$MCE = NaN
  for (i in 1:length(ff)){
    ind = which(allBOTH.filter$uniqueid == ff[i] )
    allBOTH.filter$MCE[ind] = max(as.numeric(allBOTH.filter$mce.5hz[ind]), na.rm=TRUE) # they should all be the same, just fill in
  }
  
  # Remove EF and ER below R2 filter
  ind = which(round(allBOTH.filter$R2toCO.5hz, digits = 2) < R2filter)
  allBOTH.filter$EF1CO.5hz[ind] = NaN
  allBOTH.filter$ERtoCO.5hz[ind] = NaN
  allBOTH.filter$MCE[ind] = NaN
  
  ind = which(round(allBOTH.filter$R2toCO.1hz, digits=2) < R2filter)
  allBOTH.filter$EF1CO.1hz[ind] = NaN
  allBOTH.filter$ERtoCO.1hz[ind] = NaN
  allBOTH.filter$MCE[ind] = NaN
  # start with 5hz
  allBOTH.filter$FinalEF = allBOTH.filter$EF1CO.5hz
  allBOTH.filter$FinalERtoCO = allBOTH.filter$ERtoCO.5hz
  allBOTH.filter$FinalR2 = allBOTH.filter$R2toCO.5hz
  
  ind = which(!is.finite(allBOTH.filter$FinalEF) & is.finite(allBOTH.filter$EF1CO.1hz))
  allBOTH.filter$FinalEF[ind] = allBOTH.filter$EF1CO.1hz[ind]
  ind = which(!is.finite(allBOTH.filter$FinalERtoCO) & is.finite(allBOTH.filter$ERtoCO.1hz))
  allBOTH.filter$FinalERtoCO[ind] = allBOTH.filter$ERtoCO.1hz[ind]
  allBOTH.filter$FinalR2[ind] = allBOTH.filter$R2toCO.1hz[ind]
  
  # If R2 is between 0.5 and 0.75, assume we can use the integration method
  # actually don't
  doINT = 0
  if (doINT == 1){
    ind = which(!is.finite(allBOTH.filter$FinalEF) & allBOTH.filter$R2toCO.5hz > R2Bot  &
                  allBOTH.filter$R2toCO.5hz < R2filter & is.finite(allBOTH.filter$EF1COintfill.5hz))
    allBOTH.filter$FinalEF[ind] = allBOTH.filter$EF1COintfill.5hz[ind]
    ind = which(!is.finite(allBOTH.filter$FinalERtoCO) & allBOTH.filter$R2toCO.5hz > R2Bot  &
                  allBOTH.filter$R2toCO.5hz < R2filter & is.finite(allBOTH.filter$ERtoCOintfill.5hz))
    allBOTH.filter$FinalERtoCO[ind] = allBOTH.filter$ERtoCOintfill.5hz[ind]
   
    ind = which(!is.finite(allBOTH.filter$FinalEF) & allBOTH.filter$R2toCO.1hz > R2Bot  &
                  allBOTH.filter$R2toCO.1hz < R2filter & is.finite(allBOTH.filter$EF1COintfill.1hz))
    allBOTH.filter$FinalEF[ind] = allBOTH.filter$EF1COintfill.1hz[ind]
    ind = which(!is.finite(allBOTH.filter$FinalERtoCO) & allBOTH.filter$R2toCO.1hz > R2Bot  &
                  allBOTH.filter$R2toCO.1hz < R2filter & is.finite(allBOTH.filter$ERtoCOintfill.1hz))
    allBOTH.filter$FinalERtoCO[ind] = allBOTH.filter$ERtoCOintfill.1hz[ind]
  }   
  # ----------- Now, choose EFs that had the best correlation to either CO2 or CO --------
  #allBOTH.filter$ChosenEF.5hz = NaN; allBOTH.filter$ChosenEF.R2.5hz = NaN
  #allBOTH.filter$ChosenEF.1hz = NaN; allBOTH.filter$ChosenEF.R2.1hz = NaN
  # Pick highest correlation for either the CO or CO2 EF
  #for (i in 1:length(allBOTH.filter$variable)){
    # Pick for 5hz data
  #  tmpR = c(as.numeric(allBOTH.filter$R2toX.5hz[i]),as.numeric(allBOTH.filter$R2toCO.5hz[i]))
  #  tmpEF= c(as.numeric(allBOTH.filter$EF1.5hz[i]),as.numeric(allBOTH.filter$EF1CO.5hz[i]))
  #  # Ok, actually I just want to use CO
  #  tmpR = c(as.numeric(allBOTH.filter$R2toCO.5hz[i]))
  #  tmpEF= c(as.numeric(allBOTH.filter$EF1CO.5hz[i]))
  #  ind = which(is.finite(tmpR))
  #  if (length(ind) > 0){
  #    ind2 = which(tmpR == max(tmpR, na.rm=TRUE))
  #    if (length(ind2) ==1){
  #      allBOTH.filter$ChosenEF.5hz[i] = tmpEF[ind2]
  #      allBOTH.filter$ChosenEF.R2.5hz[i] = tmpR[ind2]
  #    } else{
  #      allBOTH.filter$ChosenEF.5hz[i] = mean(tmpEF[ind2], na.rm=TRUE)
  #      allBOTH.filter$ChosenEF.R2.5hz[i] = mean(tmpR[ind2], na.rm=TRUE)
  #    }
  #  }
  #  
  #  # Pick for 1hz data
  #  tmpR = c(as.numeric(allBOTH.filter$R2toX.1hz[i]),as.numeric(allBOTH.filter$R2toCO.1hz[i]))
  #  tmpEF= c(as.numeric(allBOTH.filter$EF1.1hz[i]),as.numeric(allBOTH.filter$EF1CO.1hz[i]))
    # Ok, actually I just want to use CO
  #  tmpR = c(as.numeric(allBOTH.filter$R2toCO.1hz[i]))
  #  tmpEF= c(as.numeric(allBOTH.filter$EF1CO.1hz[i]))
  #  ind = which(is.finite(tmpR))
  #  if (length(ind) > 0){
  #    ind2 = which(tmpR == max(tmpR, na.rm=TRUE))
  #    if (length(ind2) ==1){
  #      allBOTH.filter$ChosenEF.1hz[i] = tmpEF[ind2]
  #      allBOTH.filter$ChosenEF.R2.1hz[i] = tmpR[ind2]
  #    } else{
        #if (max(tmpEF[ind2],na.rm=TRUE) < 0){print(c("i",i,tmpEF[ind2]))}
        # for cans, R2 is both 1.  But the ratio to CO seems to work better
  #      if (tmpR[1] == 1 & tmpR[2] == 1 & allBOTH.filter$PI == 'APEL'){
  #          allBOTH.filter$ChosenEF.1hz[i] = as.numeric(allBOTH.filter$EF1CO.1hz[i])
  #        allBOTH.filter$ChosenEF.R2.1hz[i] = as.numeric(allBOTH.filter$R2toCO.1hz[i])
  #      }
  #      if (tmpR[1] == 1 & tmpR[2] == 1 & allBOTH.filter$PI == 'ppt'){
  #        allBOTH.filter$ChosenEF.1hz[i] = as.numeric(allBOTH.filter$EF1CO.1hz[i])
  #        allBOTH.filter$ChosenEF.R2.1hz[i] = as.numeric(allBOTH.filter$R2toCO.1hz[i])
  #      }
  #      if (tmpR[1] == 1 & tmpR[2] == 1 & allBOTH.filter$PI == 'BLAKE'){
  #        allBOTH.filter$ChosenEF.1hz[i] = as.numeric(allBOTH.filter$EF1CO.1hz[i])
  #        allBOTH.filter$ChosenEF.R2.1hz[i] = as.numeric(allBOTH.filter$R2toCO.1hz[i])
  #      }
  #      if (tmpR[1] == 1 & tmpR[2] == 1 & allBOTH.filter$PI == 'GILMAN'){
  #        allBOTH.filter$ChosenEF.1hz[i] = as.numeric(allBOTH.filter$EF1CO.1hz[i])
  #        allBOTH.filter$ChosenEF.R2.1hz[i] = as.numeric(allBOTH.filter$R2toCO.1hz[i])
  #      }
  #    }
  #  }
  #  if (allBOTH.filter$variable[i] == 'NO2_ACES_WOMACK'){
  #    print(c(i,allBOTH.filter$variable[i],allBOTH.filter$ChosenEF.5hz[i],allBOTH.filter$ChosenEF.1hz[i]))
  #  }
  #}
  
  ff = unique(allBOTH.filter$uniqueid)
  gg = unique(allBOTH.filter$variable)
  for (i in 1:length(ff)){
    for (j in 1:length(gg)){
      ind = which(allBOTH.filter$uniqueid == ff[i] & allBOTH.filter$variable == gg[j])
      allBOTH.filter$Category.5hz[ind] = max(c(as.numeric(allBOTH.filter$Category.1hz[ind]),as.numeric(allBOTH.filter$Category.5hz[ind])), na.rm=TRUE) # they should all be the same, just fill in
    }
  }
  
  # I think I need to get rid of plumes at 1hz < 400 ppb CO too for WAS
  ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' &  as.numeric(allBOTH.filter$maxval.1hz) < COcutoff) 
  badpasses= allBOTH.filter$uniqueid[ind]
  ind = which(allBOTH.filter$PI == 'BLAKE' & allBOTH.filter$uniqueid %in% badpasses) # WAS
  allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
  ind = which(allBOTH.filter$PI == 'ppt' & allBOTH.filter$uniqueid %in% badpasses) # TOGA
  allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
  ind = which(allBOTH.filter$PI == 'GILMAN' & allBOTH.filter$uniqueid %in% badpasses) # iWAS
  allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
  
  save(allBOTH.filter, file='AllBOTH.filter.RData')
} else{
  load('AllBOTH.filter.RData')
}# end do process


# ---- Remove very strange outliers -----
# TOGA struggles with Wiggins-neighbor pass 2 & Limoncello 3
ind = which(allBOTH.filter$fire == 'Wiggins-neighbor' & allBOTH.filter$pass == 2 & allBOTH.filter$PI == 'ppt')
allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
ind = which(allBOTH.filter$fire == 'Limoncello' & allBOTH.filter$pass == 3 & allBOTH.filter$PI == 'ppt')
allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
# iWAS struggles with Ratatouille pass 2
ind = which(allBOTH.filter$fire == 'Ratatouille' & allBOTH.filter$pass == 2 & allBOTH.filter$PI == 'GILMAN')
allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN

if (doprocessSTEP2 == 1){
  
  # ----------- Get rid of emission factors at 1hz and 5hz that had poor correlation ----------
  #ind = which(allBOTH.filter$ChosenEF.R2.1hz < R2filter)
  #allBOTH.filter$ChosenEF.1hz[ind] = NaN
  #ind = which(allBOTH.filter$ChosenEF.R2.5hz < R2filter)
  #allBOTH.filter$ChosenEF.5hz[ind] = NaN
  
  # ------ Make a column just of Final EF where first we use 5hz
  #allBOTH.filter$FinalEF = allBOTH.filter$ChosenEF.5hz
  #allBOTH.filter$FinalR2 = allBOTH.filter$ChosenEF.R2.5hz
  #ind = which(is.na(allBOTH.filter$FinalEF) & is.finite(allBOTH.filter$ChosenEF.1hz))
  #allBOTH.filter$FinalEF[ind] = allBOTH.filter$ChosenEF.1hz[ind]
  #allBOTH.filter$FinalR2[ind] = allBOTH.filter$ChosenEF.R2.1hz[ind]
  # ------ Make a column just of Final ERtoCO where first we use 5hz
  #allBOTH.filter$FinalERtoCO = allBOTH.filter$ERtoCO.5hz
  #ind = which(is.na(allBOTH.filter$FinalERtoCO) & is.finite(allBOTH.filter$ERtoCO.1hz))
  #allBOTH.filter$FinalERtoCO[ind] = allBOTH.filter$ERtoCO.1hz[ind]
  # ---- Can have finite ER to CO but NaN EF if CH4 or CO2 are missing
  # ---------------remove these for consistency
  #ind = which(!is.finite(allBOTH.filter$FinalEF))
  #allBOTH.filter$FinalERtoCO[ind] = NaN
  
  allBOTH.filter$kind = allBOTH.filter$kind.5hz
  allBOTH.filter$formula = allBOTH.filter$formula.5hz
  allBOTH.filter$names = allBOTH.filter$names.5hz
  allBOTH.filter$lifetime = allBOTH.filter$lifetime_5hz_hr
  ind = which(is.na(allBOTH.filter$lifetime))
  allBOTH.filter$lifetime[ind]= allBOTH.filter$lifetime_1hz_hr[ind]
  ind = which(is.na(allBOTH.filter$kind))
  allBOTH.filter$kind[ind] = allBOTH.filter$kind.1hz[ind]
  ind = which(is.na(allBOTH.filter$formula))
  allBOTH.filter$formula[ind] = allBOTH.filter$formula.1hz[ind]
  ind = which(is.na(allBOTH.filter$names))
  allBOTH.filter$names[ind] = allBOTH.filter$names.1hz[ind]
  
  # get rid of Inf
  ind = which(allBOTH.filter$maxval.5hz == -Inf)
  allBOTH.filter$maxval.5hz[ind] = NaN
  
  
  # Probs should get rid of negative emission factors, however need to go back and check these
  ind = which(allBOTH.filter$ERtoCO.1hz < 0 & allBOTH.filter$Category.1hz == 1)
  allBOTH.filter$EF1CO.1hz[ind] = NaN
  allBOTH.filter$ERtoCO.1hz[ind] = NaN
  
  allBOTH.filter.allfuels = allBOTH.filter
  # provide the dominant fuel class for the 1hz data
 # ff = unique(allBOTH.filter.allfuels$uniqueid)
#  for (i in 1:length(ff)){
#    ind = which(allBOTH.filter.allfuels$uniqueid == ff[i] )
#    allBOTH.filter.allfuels$transect_dominant_fuel.5hz[ind] = max(allBOTH.filter.allfuels$transect_dominant_fuel.1hz[ind], na.rm=TRUE) # they should all be the same, just fill in
#    allBOTH.filter.allfuels$transect_fuel_class.5hz[ind] = max(allBOTH.filter.allfuels$transect_fuel_class.1hz[ind], na.rm=TRUE) # they should all be the same, just fill in
#    allBOTH.filter.allfuels$transect_fuel_confidence.5hz[ind] = max(allBOTH.filter.allfuels$transect_fuel_confidence.1hz[ind], na.rm=TRUE) # they should all be the same, just fill in
 # }
  # Dominant fuel class
  #1	Forest
  #2	Savanna
  #3	Shrubland
  #4	Grassland
  #5	Cropland
  #6	Pile
  #7	Slash
  #8	Understory
  #9	Urban/Barren
  # ----- For all EFs that are NaN, set MCE to NaN -----
  ind = which(is.na(allBOTH.filter.allfuels$FinalEF))
  allBOTH.filter.allfuels$MCE[ind] = NaN

  # ------- get rid of non-ag fuels ------------
  # Use the new flags from Amber
  #ind = which(allBOTH.filter.allfuels$transect_fuel_class.5hz >= 3 &
  #              allBOTH.filter.allfuels$transect_fuel_class.5hz <= 7)
  ind = which(allBOTH.filter.allfuels$fuel != '?' & allBOTH.filter.allfuels$fuel != 'house')
  # lets keep all the fuels but '?' for now, and filter after calculating more stuff
  #ind = which(allBOTH.filter.allfuels$fuel == 'slash' | allBOTH.filter.allfuels$fuel == 'pile' | 
  #              allBOTH.filter.allfuels$fuel == 'shrub' |
  #              allBOTH.filter.allfuels$fuel == 'grass' |
  #              allBOTH.filter.allfuels$fuel == 'rice'  | allBOTH.filter.allfuels$fuel == 'corn' |
  #              allBOTH.filter.allfuels$fuel == 'soybean' | 
  #              allBOTH.filter.allfuels$fuel == 'winter wheat' )
  allBOTH.filter = allBOTH.filter.allfuels[ind,]

  # check for zeros that will mess up averages
  ind = which(allBOTH.filter$ERtoCO.1hz == 0)
  allBOTH.filter$ERtoCO.1hz[ind] = NaN
  allBOTH.filter$EF1CO.1hz[ind] = NaN

  # ---- SET USEME ----
  # 1: include in the table and total VOC EF
  # 2: include in the table not total EF
  # -1 : Dont report at all, not an emission
  # 0: Dont include in the table but include for analysis of MCE dependence etc.
  allBOTH.filter$USEME = 1  
  # No longer reported by CALTECH
  ind1 = which(allBOTH.filter$names == 'C4 Hydroxyperoxide')
  allBOTH.filter$USEME[ind1] = 0
  # ------ These PM1 species don't correlate with CO ------
  ind1 = which(allBOTH.filter$variable =="Iodine_JIMENEZ")
  allBOTH.filter$USEME[ind1] = -1
  ind1 = which(allBOTH.filter$variable == "ClO4_JIMENEZ")
  allBOTH.filter$USEME[ind1] = -1
  ind1 = which(allBOTH.filter$variable =="Bromine_JIMENEZ")
  allBOTH.filter$USEME[ind1] = -1
  ind1 = which(allBOTH.filter$variable == "Seasalt_JIMENEZ")
  allBOTH.filter$USEME[ind1] = -1
  ind1 = which(allBOTH.filter$variable == "MSA_JIMENEZ")
  allBOTH.filter$USEME[ind1] = -1
  # ---------------- Get PM1 EF --------------------------------
  ind1 = which(allBOTH.filter$variable == 'OC_JIMENEZ')
  oc = allBOTH.filter[ind1,]
  oa = oc
  oa$FinalEF = oa$FinalEF*oa$OAtoOC.5hz
  oa$FinalERtoCO = oa$FinalERtoCO*oa$OAtoOC.5hz
  oa$names = 'Organic aerosol'
  oa$variable = 'OA_JIMENEZ'
  
  ind1 = which(allBOTH.filter$variable == 'BC_SCHWARZ')
  bc = allBOTH.filter[ind1,]
  ind1 = which(allBOTH.filter$variable == 'Ammonium_JIMENEZ')
  ammonium= allBOTH.filter[ind1,]
  ind1 = which(allBOTH.filter$variable == "Sulfate_JIMENEZ")
  sulf = allBOTH.filter[ind1,]
  ind1 = which(allBOTH.filter$variable == "Nitrate_JIMENEZ")
  nit = allBOTH.filter[ind1,]
  ind1 = which(allBOTH.filter$variable == "NR_Chloride_JIMENEZ")
  nrcl = allBOTH.filter[ind1,]
  ind1 = which(allBOTH.filter$variable == "Potassium_JIMENEZ")
  pot = allBOTH.filter[ind1,]
  for (i in 1:length(oc$variable)){
    newline = oc[i,]
    vars = c(oa$FinalEF[i],bc$FinalEF[i],ammonium$FinalEF[i],
             sulf$FinalEF[i],nit$FinalEF[i],nrcl$FinalEF[i], pot$FinalEF[i])
    varsER = c(oc$FinalERtoCO[i],bc$FinalERtoCO[i],ammonium$FinalERtoCO[i],
             sulf$FinalERtoCO[i],nit$FinalERtoCO[i],nrcl$FinalERtoCO[i],pot$FinalERtoCO[i])
    # only sum PM1 if we have OC
    newline$FinalEF = NaN
    if (is.finite(vars[1])){newline$FinalEF = sum(vars, na.rm=TRUE)}
    if (is.finite(vars[1])){newline$FinalERtoCO = sum(varsER, na.rm=TRUE)}
    newline$variable = 'PM1'
    newline$names = 'PM1'
    newline$formula = 'PM1'
    # give it a big MW so its in the right order
    newline$mWs= 500
    allBOTH.filter = rbind(allBOTH.filter,newline)
  }
  # also append OA
  allBOTH.filter = rbind(allBOTH.filter,oa)
  

  # --- actually don't do this
  dobeckyoutliers = 0
  if (dobeckyoutliers == 1){
    # ---- Get rid of outlier points from Becky ------
    ind = which(allBOTH.filter$fire == "U2" & allBOTH.filter$variable == 'nButane_ppt')
    allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    ind = which(allBOTH.filter$fire == "U2" & allBOTH.filter$variable == 'nPentane_ppt')
    allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    
    ind = which(allBOTH.filter$fire == "Supertramp" & allBOTH.filter$variable == 'nButane_ppt')
    allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    ind = which(allBOTH.filter$fire == "Supertramp" & allBOTH.filter$variable == 'nPentane_ppt')
    allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    
    #ind = which(allBOTH.filter$fire == "Blackwater" & allBOTH.filter$variable == 'nButane_ppt')
    #allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    #ind = which(allBOTH.filter$fire == "Blackwater" & allBOTH.filter$variable == 'nPentane_ppt')
    #allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    #ind = which(allBOTH.filter$fire == "Blackwater" & allBOTH.filter$variable == 'iButane_ppt')
    #allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    #ind = which(allBOTH.filter$fire == "Blackwater" & allBOTH.filter$variable == 'iButene1Butene_ppt')
    #allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    
    ind = which(allBOTH.filter$fire == "Willow" & allBOTH.filter$names == 'n-Butane')
    allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    ind = which(allBOTH.filter$fire == "Willow" & allBOTH.filter$variable == 'n-Pentane')
    allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    ind = which(allBOTH.filter$fire == "Willow" & allBOTH.filter$variable == 'Isobutane')
    allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    ind = which(allBOTH.filter$fire == "Willow" & allBOTH.filter$variable == 'i-Butene/1Butene')
    allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    
    ind = which(allBOTH.filter$fire == "BugsBunny" & allBOTH.filter$variable == 'n-Butane')
    allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    ind = which(allBOTH.filter$fire == "BugsBunny" & allBOTH.filter$variable == 'n-Pentane')
    allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    ind = which(allBOTH.filter$fire == "BugsBunny" & allBOTH.filter$variable == 'Isobutane')
    allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    ind = which(allBOTH.filter$fire == "BugsBunny" & allBOTH.filter$variable == 'i-Butene/1Butene')
    allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    ind = which(allBOTH.filter$fire == "BugsBunny" & allBOTH.filter$variable == '2-Methylpentane')
    allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    ind = which(allBOTH.filter$fire == "BugsBunny" & allBOTH.filter$variable == '3-Methylpentane')
    allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
    
    ind = which(allBOTH.filter$fire == "FatAlbert" & allBOTH.filter$variable == 'n-Pentane')
    allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN
  
    # just for plotting
    ind = which(was.all$fire == "U2" )
    was.all$nButane_WAS_BLAKE[ind] = NaN 
    was.all$nPentane_WAS_BLAKE[ind] = NaN 
  
    ind = which(was.all$fire == "Supertramp" )
    was.all$nButane_WAS_BLAKE[ind] = NaN 
    was.all$nPentane_WAS_BLAKE[ind] = NaN
    
    ind = which(was.all$fire == "Willow" )
    was.all$nButane_WAS_BLAKE[ind] = NaN 
    was.all$nPentane_WAS_BLAKE[ind] = NaN
    was.all$iButane_WAS_BLAKE[ind] = NaN
    was.all$iButene_WAS_BLAKE[ind] = NaN
    
    ind = which(was.all$fire == "BugsBunny" )
    was.all$nButane_WAS_BLAKE[ind] = NaN 
    was.all$nPentane_WAS_BLAKE[ind] = NaN
    was.all$iButane_WAS_BLAKE[ind] = NaN
    was.all$iButene_WAS_BLAKE[ind] = NaN
    was.all$x2MePentane_WAS_BLAKE[ind] = NaN
    was.all$x3MePentane_WAS_BLAKE[ind] = NaN
    
    ind = which(was.all$fire == "FatAlbert" )
    was.all$nPentane_WAS_BLAKE[ind] = NaN
  }
    
  # ------- Removing these measurements entirely
  # No emission factors come out for these species
  ind = which(allBOTH.filter$formula == 'C2H4O3S')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$formula == 'MSA')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$formula == 'Na')
  allBOTH.filter$USEME[ind] = -1
  # Don't correlate with CO
  ind = which(allBOTH.filter$formula == 'BrO')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$formula == 'BrCl')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$formula == 'BrCN')
  allBOTH.filter$USEME[ind] = -1
  # ---------- These TOGA species don't correlate with CO
  ind = which(allBOTH.filter$variable  == 'iPropONO2_ppt')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CHBr3_ppt')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable =='CHBrCl2_ppt')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CH3CCl3_ppt')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CHCl3_ppt')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'HFC134a_ppt')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'HCFC142b_ppt')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'HCFC141b_ppt')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CHBr2Cl_ppt')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CH2ClI_ppt')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'LimoneneD3Carene_ppt')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'Propane_ppt')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'Propene_ppt')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'MBO_ppt')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'iButONO2and2ButONO2_ppt')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'C2H5OH_ppt')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CH2ClCH2Cl_ppt' ) 
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'HCFC22_ppt' ) 
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CH2Cl2_ppt' ) # TOGA doesn't correlate with CO but WAS does
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'C2Cl4_ppt' ) # 
  allBOTH.filter$USEME[ind] = -1
  # These iWAS species dont correlate with CO
  ind = which(allBOTH.filter$variable == 'CycHexane_NOAAiWAS_GILMAN')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable  == 'x3MePentane_NOAAiWAS_GILMAN')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable  == 'x224TriMePentane_NOAAiWAS_GILMAN')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable  == 'x22DiMeButane_NOAAiWAS_GILMAN')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable  == 'CHCl3_NOAAiWAS_GILMAN')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'C2Cl4_NOAAiWAS_GILMAN')
  allBOTH.filter$USEME[ind] = -1
  
  #These WAS species dont correlate with CO
  ind = which(allBOTH.filter$variable == 'C2Cl4_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CHBrCl2_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'C2HCl3_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CHCl3_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'Limonene_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'ClBenzene_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'H1211_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CFC11_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CFC12_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'x234TrimePentane_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CH2ClCH2Cl_WAS_BLAKE' ) 
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CCl4_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CH3CCl3_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CFC12_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'HFC134a_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'HFC152a_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'HCFC142b_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'HCFC141b_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'HFC365mfc_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'x2MePentane_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CFC114_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'x3MePentane_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'H1301_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'H2402_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CHBr2Cl_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CFC113_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable =='HCFC22_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'x23Dimebutane_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'CHBr3_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable =='x2PentONO2_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable =='x2ButONO2_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable =='x3PentONO2_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'iPropONO2_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable == 'x3Me2ButONO2_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = -1
#  ind = which(allBOTH.filter$variable =='CycPentane_WAS_BLAKE')
 # allBOTH.filter$USEME[ind] = 0
  
  # For these species, set USEME == 2 to use in table, not total VOC
  ind2 = which(allBOTH.filter$variable == 'MVK_ppt')
  allBOTH.filter$USEME[ind2] = 2
  ind3 = which(allBOTH.filter$variable== 'MAC_ppt')
  allBOTH.filter$USEME[ind3] = 2
  ind4 = which(allBOTH.filter$variable== 'x2Butenals_ppt')
  allBOTH.filter$USEME[ind4] = 2
  
  ind2 = which(allBOTH.filter$variable== 'Propanal_ppt')
  allBOTH.filter$USEME[ind2] = 2
  ind3 = which(allBOTH.filter$variable== 'Acetone_ppt')
  allBOTH.filter$USEME[ind3] =2
  
  ind3 = which(allBOTH.filter$variable== 'CRESOL_WENNBERG')
  allBOTH.filter$USEME[ind3] =2
  
  # ----- These species just really don't correlate in my opinion!
  ind = which(allBOTH.filter$variable =='Cl2_NOAACIMS_VERES')
  allBOTH.filter$USEME[ind] = -1
  ind = which(allBOTH.filter$variable =='ISOPN_WENNBERG')
  allBOTH.filter$USEME[ind] = -1
  
  # Speciate C9 aromatics with Blake (C9H12) - maybe not enough, so just dont include the speciated in the total VOC
  ind = which(allBOTH.filter$variable== 'C9Aromatics_NOAAPTR_ppbv_WARNEKE')
  allBOTH.filter$USEME[ind] = 1
  ind2 = which(allBOTH.filter$names== '1,2,4-Trimethylbenzene')
  allBOTH.filter$USEME[ind2] = 2
  ind3 = which(allBOTH.filter$names== '1,3,5-trimethylbenzene' | allBOTH.filter$names== '1,3,5-Trimethylbenzene' )
  allBOTH.filter$USEME[ind3] = 2
  ind4 = which(allBOTH.filter$names== 'i-Propylbenzene')
  allBOTH.filter$USEME[ind4] = 2
  ind5 = which(allBOTH.filter$names== 'n-Propylbenzene')
  allBOTH.filter$USEME[ind5] = 2
  ind6 = which(allBOTH.filter$names== '2-Ethyltoluene')
  allBOTH.filter$USEME[ind6] = 2
  ind7 = which(allBOTH.filter$names== '3-Ethyltoluene')
  allBOTH.filter$USEME[ind7] = 2
  ind8 = which(allBOTH.filter$names== '4-Ethyltoluene')
  allBOTH.filter$USEME[ind8] = 2
  
  # Toga CH2Br2 seems weird - maybe use blake?
  ind = which(allBOTH.filter$variable == 'CH2Br2_ppt')
  allBOTH.filter$USEME[ind] =-1
  # -- Actually just use TOGA furan, methyl furan, and furfural per GIGI's paper 
  ind = which(allBOTH.filter$names == ' Furan and fragments' & allBOTH.filter$PI == 'WARNEKE'  ) 
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$names == 'sum of 2-methylfuran 3-methylfuran and fragments' & allBOTH.filter$PI == 'WARNEKE'  ) 
  allBOTH.filter$USEME[ind] = 0
  # -- Warneke has the fragments so keep it for total VOC for  and furfural
  ind = which(allBOTH.filter$formula == 'C5H4O2' & allBOTH.filter$PI == 'WARNEKE'  ) 
  allBOTH.filter$USEME[ind] = 0
  
  # don't use BLAKE MVK and MACR at all
  ind = which(allBOTH.filter$variable == 'MVK_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'MAC_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  
  # -- Warneke - keep for total VOC
  ind = which(allBOTH.filter$names == 'Acetone')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$names == 'Propanal')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$names == 'Acetone/Propanal' & allBOTH.filter$PI != 'WARNEKE')
  allBOTH.filter$USEME[ind] = 0
  
  ind = which(allBOTH.filter$PI == 'STCLAIR')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$names == 'NOy')
  allBOTH.filter$USEME[ind] = 0
  # all measurements agree so just keep WARNEKE
  ind = which(allBOTH.filter$names == 'Benzene' & allBOTH.filter$PI != 'WARNEKE') 
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$names == 'Toluene' & allBOTH.filter$PI != 'WARNEKE') 
  allBOTH.filter$USEME[ind] = 0
  # ---- Weird toluene outlier
  ind = which(allBOTH.filter$variable == 'Toluene_NOAAPTR_ppbv_WARNEKE' & allBOTH.filter$FinalEF > 2)
  allBOTH.filter$FinalEF[ind] = allBOTH.filter$EF1CO.5hz[ind]
  allBOTH.filter$FinalERtoCO[ind] = allBOTH.filter$ERtoCO.5hz[ind]
  
  ind = which(allBOTH.filter$names == 'Acetaldehyde' & allBOTH.filter$PI != 'WARNEKE') 
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$names == 'Acrolein' & allBOTH.filter$PI != 'WARNEKE') #APEL and WARNEKE agree, BLAKE is low
  allBOTH.filter$USEME[ind] = 0
  # Don't use WARNEKE or TOGA CH2O
  ind = which(allBOTH.filter$names == 'Formaldehyde' & allBOTH.filter$PI == 'ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$names == 'Formaldehyde' & allBOTH.filter$PI == 'WARNEKE')
  allBOTH.filter$USEME[ind] = 0
  
  # All but WAS agree so just keep warneke for MEK, but average iWAS and TOGA for table
  #allBOTH.filter = mergelines(allBOTH.filter,  'MEK_ppt','MEK_NOAAiWAS_GILMAN')
  
  ind = which(allBOTH.filter$names == 'MEK/ 2-methyl propanal' & allBOTH.filter$PI == 'WARNEKE')
  allBOTH.filter$USEME[ind] = 1
  ind = which(allBOTH.filter$names == 'Methyl ethyl ketone' & allBOTH.filter$PI == 'GILMAN')#
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$names == 'Methyl ethyl ketone' & allBOTH.filter$PI == 'BLAKE' )
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$names == 'Methyl ethyl ketone' &  allBOTH.filter$PI == 'ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$names == 'Methyl ethyl ketone' & allBOTH.filter$PI != 'WARNEKE' & allBOTH.filter$PI != 'GILMAN' &
                allBOTH.filter$PI != 'ppt' & allBOTH.filter$PI != 'BLAKE')
  allBOTH.filter$USEME[ind] = 2
  # Stategy - maybe average the 'short' measurements?
  # ------ Make averages of specific species, set USME for individuals == 0 -----------
  # --- Average Fried and Hanisco HCHO
  allBOTH.filter = mergelines(allBOTH.filter,  'CH2O_CAMS_pptv_FRIED','CH2O_ISAF_HANISCO')
  # --- Average Ryerson and Rollins NO
  allBOTH.filter = mergelines(allBOTH.filter,  'NO_LIF_ROLLINS','NO_RYERSON')
  # --- Average Ryerson and Womack NO2
  allBOTH.filter = mergelines(allBOTH.filter,  'NO2_ACES_WOMACK','NO2_RYERSON')
  # --- Average VERES and WOMACK HONO
  allBOTH.filter = mergelines(allBOTH.filter,  'HNO2_ACES_WOMACK','HNO2_NOAACIMS_VERES')
  # --- Average GILMAN AND BLAKE ethene
  allBOTH.filter = mergelines(allBOTH.filter, 'Ethene_NOAAiWAS_GILMAN','Ethene_WAS_BLAKE')
  # --- Average GILMAN AND BLAKE Propene
  allBOTH.filter = mergelines(allBOTH.filter, 'Propene_NOAAiWAS_GILMAN','Propene_WAS_BLAKE')
  # --- Average GILMAN AND BLAKE i-Butene
  allBOTH.filter = mergelines(allBOTH.filter, 'iButene_NOAAiWAS_GILMAN','iButene_WAS_BLAKE')
  # --- Average GILMAN AND BLAKE   trans-2-Butene
  allBOTH.filter = mergelines(allBOTH.filter, 't2butene_NOAAiWAS_GILMAN','t2Butene_WAS_BLAKE')
  # --- Average GILMAN AND BLAKE  cis-2-Butene
  allBOTH.filter = mergelines(allBOTH.filter, 'c2Butene_NOAAiWAS_GILMAN','c2Butene_WAS_BLAKE')
  # --- Average GILMAN AND BLAKE 1-Butene
  allBOTH.filter = mergelines(allBOTH.filter, 'x1Butene_NOAAiWAS_GILMAN','x1Butene_WAS_BLAKE')
  # --- Average GILMAN AND BLAKE 1-Pentene
  allBOTH.filter = mergelines(allBOTH.filter, 'x1Pentene_NOAAiWAS_GILMAN','x1Pentene_WAS_BLAKE')
  # --- Average GILMAN AND BLAKE t2-Pentene
  allBOTH.filter = mergelines(allBOTH.filter, 't2Pentene_NOAAiWAS_GILMAN','t2Pentene_WAS_BLAKE')
  # --- Average GILMAN AND BLAKE c2-Pentene
  allBOTH.filter = mergelines(allBOTH.filter, 'c2Pentene_NOAAiWAS_GILMAN','c2Pentene_WAS_BLAKE')
  # --- Average GILMAN AND BLAKE 1,3-pentadiene
  allBOTH.filter = mergelines(allBOTH.filter, 't13Pentadiene_NOAAiWAS_GILMAN','x13Pentadienes_WAS_BLAKE')
    # --- Average GILMAN AND BLAKE 2-methyl-1-butene
  allBOTH.filter = mergelines(allBOTH.filter, 'X2Me1Butene_WAS_BLAKE','x2Me1Butene_NOAAiWAS_GILMAN')
  # --- Average GILMAN AND BLAKE 3-methyl-1-butene
  allBOTH.filter = mergelines(allBOTH.filter, 'X3Me1Butene_WAS_BLAKE','x3Me1Butene_NOAAiWAS_GILMAN')
  # --- Average GILMAN AND BLAKE Methylcyclopentane
  allBOTH.filter = mergelines(allBOTH.filter, 'MeCycPentane_NOAAiWAS_GILMAN','MeCycPentane_WAS_BLAKE')
  # --- Average GILMAN AND BLAKE Methylcyclohexane
  allBOTH.filter = mergelines(allBOTH.filter, 'MeCycHexane_NOAAiWAS_GILMAN','MeCycHexane_WAS_BLAKE')
  # --- Average Toga AND BLAKE Camphene
  allBOTH.filter = mergelines(allBOTH.filter, 'Camphene_WAS_BLAKE','Camphene_ppt')
  # --- Average Toga AND BLAKE Isobutanal # actually just use TOGA
  ind = which(allBOTH.filter$variable == 'iButanal_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
#  allBOTH.filter = mergelines(allBOTH.filter, 'iButanal_WAS_BLAKE','iButanal_ppt')
  # --- Average Toga AND BLAKE Butanal  # actually just use TOGA
  #allBOTH.filter = mergelines(allBOTH.filter, 'Butanal_WAS_BLAKE','Butanal_ppt')
  ind = which(allBOTH.filter$variable == 'Butanal_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  # --- Average Toga AND BLAKE CH3I
  allBOTH.filter = mergelines(allBOTH.filter, 'CH3I_WAS_BLAKE','CH3I_ppt')
  # --- Average Toga AND BLAKE CH2Cl2 - does seem emitted from corn
  allBOTH.filter = mergelines(allBOTH.filter, 'CH2Cl2_WAS_BLAKE','CH2Cl2_ppt')
  # --- Average Toga AND BLAKE ethono2
  allBOTH.filter = mergelines(allBOTH.filter, 'EthONO2_WAS_BLAKE','EthONO2_ppt')
  # --- Average Toga AND BLAKE meono2
  allBOTH.filter = mergelines(allBOTH.filter, 'MeONO2_WAS_BLAKE','MeONO2_ppt')
  # --- Average GILMAN AND BLAKE Cyclohexane
  # actually, iwas cyclohexane is bad, just use blake
  ind = which(allBOTH.filter$variable == 'CycHexane_NOAAiWAS_GILMAN')
  allBOTH.filter$USEME[ind] = 0
  #allBOTH.filter = mergelines(allBOTH.filter, 'CycHexane_WAS_BLAKE','CycHexane_NOAAiWAS_GILMAN')
  # --- Average GILMAN AND BLAKE n-Decane
  allBOTH.filter = mergelines(allBOTH.filter, 'nDecane_NOAAiWAS_GILMAN','nDecane_WAS_BLAKE')
  # --- Average GILMAN AND BLAKE ethyne
  allBOTH.filter = mergelines(allBOTH.filter, 'Ethyne_NOAAiWAS_GILMAN','Ethyne_WAS_BLAKE')
  # --- Average GILMAN, BLAKE, APEL 3-methylpentane
  # GILMAN is bad, dont use
  ind = which(allBOTH.filter$variable == 'x3MePentane_NOAAiWAS_GILMAN')
  allBOTH.filter$USEME[ind] = 0
  # Blake shows a negative dependence on CO...
  ind = which(allBOTH.filter$variable == 'x3MePentane_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
#  allBOTH.filter = mergelines3(allBOTH.filter, 'x3MePentane_ppt',','')
  # --- Average GILMAN, BLAKE, APEL 2-methylpentane
  allBOTH.filter = mergelines3(allBOTH.filter, 'x2MePentane_ppt','x2MePentane_WAS_BLAKE','x2MePentane_NOAAiWAS_GILMAN')
  # --- Average GILMAN, BLAKE, APEL n-hexane
  allBOTH.filter = mergelines3(allBOTH.filter, 'nHexane_ppt','nHexane_WAS_BLAKE','nHexane_NOAAiWAS_GILMAN')
  # --- Average  BLAKE, APEL n-heptane
  allBOTH.filter = mergelines(allBOTH.filter, 'nHeptane_ppt','nHeptane_WAS_BLAKE')
  # --- Average GILMAN, BLAKE, APEL m,p-xylene
  allBOTH.filter = mergelines3(allBOTH.filter, 'mpXylene_ppt','mpXylene_WAS_BLAKE','mpXylene_NOAAiWAS_GILMAN')
  # --- Average GILMAN, BLAKE, APEL o-xylene
  allBOTH.filter = mergelines3(allBOTH.filter, 'oXylene_ppt','oXylene_WAS_BLAKE','oXylene_NOAAiWAS_GILMAN')
  # --- Average GILMAN, BLAKE, APEL ethylbenzene
  allBOTH.filter = mergelines3(allBOTH.filter, 'EthBenzene_ppt','EthBenzene_WAS_BLAKE','EthBenzene_NOAAiWAS_GILMAN')
   # --- Average WARNEKE, BLAKE, APEL Nitromethane
 # allBOTH.filter = mergelines3(allBOTH.filter, 'CH3NO2_NOAAPTR_ppbv_WARNEKE','Nitromethane_WAS_BLAKE','Nitromethane_ppt')
  # TOGA nitromethane is way low compaared to WAS, which is within 30% of WARNEKE. USE Warneke
  ind = which(allBOTH.filter$variable == 'Nitromethane_WAS_BLAKE' | allBOTH.filter$variable == 'Nitromethane_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CH3NO2_NOAAPTR_ppbv_WARNEKE')
  allBOTH.filter$USEME[ind] = 1
  
    # --- Average GBlake, TOGA, CH3Br
  allBOTH.filter = mergelines(allBOTH.filter, 'CH3Br_ppt','CH3Br_WAS_BLAKE')
  # --- Average GILMAN BLAKE, APEL 2,2,4-Trimethylpentane
  allBOTH.filter = mergelines3(allBOTH.filter, 'x224TrimePentane_ppt','x224TrimePentane_WAS_BLAKE','x224TriMePentane_NOAAiWAS_GILMAN')
  # --- Average GILMAN, BLAKE, APEL octane
  allBOTH.filter = mergelines3(allBOTH.filter, 'nOctane_ppt','nOctane_WAS_BLAKE','nOctane_NOAAiWAS_GILMAN')
  # --- Average GILMAN AND BLAKE   nonane
  allBOTH.filter = mergelines(allBOTH.filter, 'nNonane_NOAAiWAS_GILMAN','nNonane_WAS_BLAKE')
  # --- Average GILMAN, BLAKE, APEL alphapinene
  allBOTH.filter = mergelines3(allBOTH.filter, 'aPinene_ppt','aPinene_WAS_BLAKE','aPinene_NOAAiWAS_GILMAN')
  # --- Average BLAKE, APEL isopropanol
  allBOTH.filter = mergelines(allBOTH.filter, 'iPropanol_ppt','iPropanol_WAS_BLAKE')
  # --- Average BLAKE, APEL methyl acetate - actually just use TOGA
  ind = which(allBOTH.filter$variable == 'MeAcetate_WAS_BLAKE')
  allBOTH.filter$USEME == 0
  #  allBOTH.filter = mergelines(allBOTH.filter, 'MeAcetate_ppt','MeAcetate_WAS_BLAKE')
  # Use warneke for the total VOC, but keep for table
  ind = which(allBOTH.filter$variable == 'MeAcetate_ppt')
  allBOTH.filter$USEME[ind] == 2
  # For betapinene/myrcene, lets add was together, then average with TOGA
  ind = which(allBOTH.filter$variable == 'Myrcene_WAS_BLAKE')
  ind2 = which(allBOTH.filter$variable == 'bPinene_WAS_BLAKE')
  tmp1 = allBOTH.filter[ind,]; tmp2 = allBOTH.filter[ind2,]

  tmp1$FinalEF = rowSums(cbind(tmp1$FinalEF,tmp2$FinalEF), na.rm=TRUE)
  tmp1$FinalERtoCO = rowSums(cbind(tmp1$FinalERtoCO,tmp2$FinalERtoCO), na.rm=TRUE)
  tmp1$variable = 'bPinene/Myrcene_WAS_BLAKE'
  tmp1$names = 'beta-Pinene/Myrcene'
  allBOTH.filter = rbind(allBOTH.filter,tmp1)
  # have zeros though now instead of NaNs
  ind = which(allBOTH.filter$variable == 'bPinene/Myrcene_WAS_BLAKE' & allBOTH.filter$FinalEF == 0.0)
  allBOTH.filter$FinalEF[ind] = NaN
  allBOTH.filter = mergelines(allBOTH.filter, 'bPinene/Myrcene_WAS_BLAKE',
                              'bPineneMyrcene_ppt')
  ind = which(allBOTH.filter$variable == 'bPinene/Myrcene_WAS_BLAKE' |
                allBOTH.filter$variable =='bPineneMyrcene_ppt')
  allBOTH.filter$USEME[ind] = 0
  # --- Average VERES and WARNEKE HCOOH
  allBOTH.filter = mergelines(allBOTH.filter, 'HCOOH_NOAACIMS_VERES','HCOOH_NOAAPTR_ppbv_WARNEKE')
  # --- Average GILMAN, BLAKE, APEL nbutane
  allBOTH.filter = mergelines3(allBOTH.filter, 'nButane_ppt','nButane_WAS_BLAKE','nButane_NOAAiWAS_GILMAN')
  # --- Average GILMAN, BLAKE, APEL ibutane
  allBOTH.filter = mergelines3(allBOTH.filter, 'iButane_ppt','iButane_WAS_BLAKE','iButane_NOAAiWAS_GILMAN')
  # --- Average APEL and GILMAN methylformate
  allBOTH.filter = mergelines(allBOTH.filter, 'MeFormate_NOAAiWAS_GILMAN','MeFormate_ppt')
  # --- Average GILMAN, BLAKE, APEL ipentane
  allBOTH.filter = mergelines3(allBOTH.filter, 'iPentane_ppt','iPentane_WAS_BLAKE','iPentane_NOAAiWAS_GILMAN')
  # --- Average GILMAN, BLAKE, APEL npentane
  allBOTH.filter = mergelines3(allBOTH.filter, 'nPentane_ppt','nPentane_WAS_BLAKE','nPentane_NOAAiWAS_GILMAN')
  # --- Average GILMAN AND BLAKE   22-dimethylbutane
  allBOTH.filter = mergelines(allBOTH.filter, 'x22DiMeButane_NOAAiWAS_GILMAN','x22Dimebutane_WAS_BLAKE')
  # --- Average Apel AND BLAKE  ethynlbenzene
  allBOTH.filter = mergelines(allBOTH.filter, 'EthynylBenzene_ppt','EthynylBenzene_WAS_BLAKE')
  # --- Average Apel AND BLAKE  tricyclene
  allBOTH.filter = mergelines(allBOTH.filter, 'Tricyclene_ppt','Tricyclene_WAS_BLAKE')
  # --- Average warneke, veres, HNCO
  allBOTH.filter = mergelines(allBOTH.filter, 'HNCO_NOAACIMS_VERES','HNCO_NOAAPTR_ppbv_WARNEKE')
  # --- Average wennberg, gilman, toga acrylonitrile
  allBOTH.filter = mergelines4(allBOTH.filter,'Acrylonitrile_NOAAPTR_ppbv_WARNEKE', 'Acrylonitrile_ppt','Acrylonitrile_NOAAiWAS_GILMAN', 'Acrylonitrile_WAS_BLAKE')
  # --- Average wennberg, veres, warneke, and apel hydrogen cyanide - Ok, per Lu Xu, don't use WARNEKE as it may have a water interference
  # water intereference is now fixed
  allBOTH.filter = mergelines4(allBOTH.filter, 'HCN_NOAACIMS_VERES','HCN_WENNBERG','HCN_ppt','HCN_NOAAPTR_ppbv_WARNEKE')
  # --- Average Apel AND BLAKE propionitrile
  allBOTH.filter = mergelines(allBOTH.filter, 'PropNitrile_ppt','PropNitrile_WAS_BLAKE')
  # --- Average Apel AND BLAKE DMS
  allBOTH.filter = mergelines(allBOTH.filter, 'PropNitrile_ppt','PropNitrile_WAS_BLAKE')
  # --- Average BLAKE and APEL 2-methylfuran ***actually just use TOGA
  ind = which(allBOTH.filter$variable == 'x2MeFuran_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'x3MeFuran_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  
#  allBOTH.filter = mergelines(allBOTH.filter, 'x2MeFuran_WAS_BLAKE','x2MeFuran_ppt')
  # --- Average BLAKE and APEL 3-methylfuran
#  allBOTH.filter = mergelines(allBOTH.filter, 'x3MeFuran_WAS_BLAKE','x3MeFuran_ppt')
  # keep for the table, use WARNEKE for the total
  # actually use TOGA per Gigi's table
  #ind = which(allBOTH.filter$names== '2-Methylfuran' & allBOTH.filter$USEME == 1)
  #allBOTH.filter$USEME[ind] =2
  #ind = which(allBOTH.filter$names== '3-Methylfuran' & allBOTH.filter$USEME == 1)
  #allBOTH.filter$USEME[ind] =2
  # --- Average BLAKE and APEL propane
  allBOTH.filter = mergelines(allBOTH.filter, 'Propane_WAS_BLAKE','Propane_NOAAiWAS_GILMAN')
  # --- *****  Average GILMAN, BLAKE, APEL Furan - actually just use TOGA *********
  ind = which(allBOTH.filter$variable == 'Furan_WAS_BLAKE' | allBOTH.filter$variable == 'Furan_NOAAiWAS_GILMAN')
  allBOTH.filter$USEME[ind] = 0
#  allBOTH.filter = mergelines3(allBOTH.filter, 'Furan_ppt','Furan_WAS_BLAKE','Furan_NOAAiWAS_GILMAN')
  # For furan, use TOGA per Gigi's paper
  #ind2 = which(allBOTH.filter$names == 'Furan' & allBOTH.filter$USEME == 1 )
  #allBOTH.filter$USEME[ind] = 2 # keep for the table, not for the total
  # ----- For whatever reason, PTRMS phenol is factor of 4 higher than CIT-CIMS.  Using CIT-CIMS for now
  ind = which(allBOTH.filter$variable == 'Phenol_NOAAPTR_ppbv_WARNEKE')
  allBOTH.filter$USEME[ind] = 2
  # - For acetonitrile, just use WARNEKE
  ind = which(allBOTH.filter$formula == 'CH3CN' & allBOTH.filter$PI != 'WARNEKE')
  allBOTH.filter$USEME[ind] = 0
  # ----- No data on HPMTF
  ind = which(allBOTH.filter$formula == 'C2H4O3S')
  allBOTH.filter$USEME[ind] = 0
  
  # -------- For benzofuran, just use WARNEKE ----
  ind = which(allBOTH.filter$names == 'Benzofuran' & allBOTH.filter$PI != 'WARNEKE')
  allBOTH.filter$USEME[ind] = 0
  # ---- Don't have corrected TOGA methanol yet so just use WARNEKE
  ind = which(allBOTH.filter$variable == 'CH3OH_ppt')
  allBOTH.filter$USEME[ind] = 0
  # ---- For i/Butene/1butene, Gilman + Blake add up to Apel, so keep speciated from Gilman + Blake
  ind = which(allBOTH.filter$variable == 'iButene1Butene_ppt')
  allBOTH.filter$USEME[ind] = 0
  # ---- For styrene just use Warneke -----
  ind = which(allBOTH.filter$variable == 'Styrene_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'Styrene_ppt')
  allBOTH.filter$USEME[ind] = 0
  # ----- For isoprene, dont use  Warneke - possible interferences
  ind = which(allBOTH.filter$variable == 'Isoprene_NOAAPTR_ppbv_WARNEKE')
  allBOTH.filter$USEME[ind] = 0
  allBOTH.filter = mergelines3(allBOTH.filter, 'Isoprene_WAS_BLAKE','Isoprene_ppt','Isoprene_NOAAiWAS_GILMAN')
  # But, there are some outliers that look odd for the above average.  If it is 3x higher than PTRMS, remove.
  #ttt = unique(allBOTH.filter$uniqueid)
  #for ( i in 1:length(ttt)){
  #  ind = which(allBOTH.filter$uniqueid == ttt[i] & allBOTH.filter$variable == 'Isoprene_NOAAPTR_ppbv_WARNEKE')
  #  try1 =allBOTH.filter[ind,]
  #  ind2 = which(allBOTH.filter$uniqueid == ttt[i] & allBOTH.filter$variable == 'Isoprene_BLAKE_ppt_')
  #  try2= allBOTH.filter[ind2,]
  #  #if (length(try1$FinalEF) > 1){print(c(i,try1$FinalEF,try2$FinalEF))}
  #  if (is.finite(try2$FinalEF) & is.finite(try1$FinalEF)){
  #    if (try2$FinalEF > (try1$FinalEF*3)){
  #      print(c(try2$fire,try2$pass,try2$FinalEF, try1$FinalEF))
  #      allBOTH.filter$FinalEF[ind2] = NaN
  #      allBOTH.filter$FinalERtoCO[ind2] = NaN
  #    }
  #  }
  #}
  
  # ---- For DMS just use Warneke -----
#  ind = which(allBOTH.filter$variable == 'DMS_WAS_BLAKE')
#  allBOTH.filter$USEME[ind] = 0
#  ind = which(allBOTH.filter$variable == 'DMS_ppt')
#  allBOTH.filter$USEME[ind] = 0
  allBOTH.filter = mergelines3(allBOTH.filter, 'DMS_ppt','DMS_WAS_BLAKE','DMS_NOAAPTR_ppbv_WARNEKE')
  
  allBOTH.filter = mergelines3(allBOTH.filter, 'Ethane_WAS_BLAKE','Ethane_NOAAiWAS_GILMAN','C2H6_CAMS_pptv_FRIED')
  # Oh actually want to report NOy ER
  ind = which(allBOTH.filter$formula == 'NOy')
  allBOTH.filter$USEME[ind] = 1
  # ------ Lifetime category ---------
  # ------- Add jvalues --------------
  # jvalues
  hall1 = getICARTTdataSIMPLE('InputFiles/firexaq-HALL_dc8_20190821_R0_20221004T143858.ict')
  hall2 = getICARTTdataSIMPLE('InputFiles/firexaq-HALL_dc8_20190823_R0_20221004T143901.ict')
  hall3 = getICARTTdataSIMPLE('InputFiles/firexaq-HALL_dc8_20190826_R0_20221004T143904.ict')
  hall4 = getICARTTdataSIMPLE('InputFiles/firexaq-HALL_dc8_20190829_R0_20221004T143907.ict')
  hall5 = getICARTTdataSIMPLE('InputFiles/firexaq-HALL_dc8_20190830_R0_20221004T143909.ict')
  hall6 = getICARTTdataSIMPLE('InputFiles/firexaq-HALL_dc8_20190831_R0_20221004T143913.ict')
  hall7 = getICARTTdataSIMPLE('InputFiles/firexaq-HALL_dc8_20190903_R0_20221004T143916.ict')
  hall = rbind(hall1,hall2,hall3,hall4,hall5,hall6,hall7)
  hall$LST = hall$Time_Start/60/60-5
  allBOTH.filter$lifetime_jval = NaN
  ind2 = which(hall$LST >= 10 & hall$LST <= 17)
  
  ind = which(allBOTH.filter$names == 'Formaldehyde')
  allBOTH.filter$lifetime_jval[ind] =1/ mean(hall$jCH2O_H2_CO_CAFS_HALL[ind2]+hall$jCH2O_H_HCO_CAFS_HALL[ind2], na.rm=TRUE)/60/60
  ind = which(allBOTH.filter$names == 'Acetaldehyde')
  allBOTH.filter$lifetime_jval[ind] = 1/mean(hall$jCH3CHO_CH3_HCO_CAFS_HALL[ind2], na.rm=TRUE)/60/60
  ind = which(allBOTH.filter$names == 'Propanal')
  allBOTH.filter$lifetime_jval[ind] = 1/mean(hall$jCH3CHO_CH3_HCO_CAFS_HALL[ind2], na.rm=TRUE)/60/60
  ind = which(allBOTH.filter$names == 'Methyl Ethyl Ketone' | allBOTH.filter$variable == 'C4Carbonyls_NOAAPTR_ppbv_WARNEKE')
  allBOTH.filter$lifetime_jval[ind] = 1/mean(hall$jMEK_CH3CO_CH2CH3_CAFS_HALL[ind2], na.rm=TRUE)/60/60
  ind = which(allBOTH.filter$variable == 'x23Butanedione_NOAAPTR_ppbv_WARNEKE')
  allBOTH.filter$lifetime_jval[ind] = 1/mean(hall$j23Butanedione_NoProductsSpecified_CAFS_HALL[ind2], na.rm=TRUE)/60/60
  ind = which(allBOTH.filter$names == 'Acetone')
  allBOTH.filter$lifetime_jval[ind] = 1/mean(hall$jAcetone_CH3CO_CH3_CAFS_HALL[ind2], na.rm=TRUE)/60/60
  ind = which(allBOTH.filter$names == 'Methylglyoxal')
  allBOTH.filter$lifetime_jval[ind] = 1/mean(hall$jCH3COCHO_CH3CO_HCO_CAFS_HALL[ind2], na.rm=TRUE)/60/60
  ind = which(allBOTH.filter$variable == 'C3H6O2_NOAAPTR_ppbv_WARNEKE')
  allBOTH.filter$lifetime_jval[ind] = 1/mean(hall$jHydroxyacetone_CH3CO_CH3O_CAFS_HALL[ind2]+
                                               hall$jHydroxyacetone_CH3COO_CH3_CAFS_HALL[ind2], na.rm=TRUE)/60/60
  ind = which(allBOTH.filter$names == 'Glyoxal')
  allBOTH.filter$lifetime_jval[ind] = 1/mean(hall$jCHOCHO_CH2O_CO_CAFS_HALL[ind2]+
                                               hall$jCHOCHO_H2_2CO_CAFS_HALL[ind2]+ 
                                               hall$jCHOCHO_HCO_HCO_CAFS_HALL[ind2], na.rm=TRUE)/60/60
  ind = which(allBOTH.filter$names == 'Methyl vinyl ketone')
  allBOTH.filter$lifetime_jval[ind] = 1/mean(hall$jMVK_NoProductsSpecified_CAFS_HALL[ind2], na.rm=TRUE)/60/60
  ind = which(allBOTH.filter$names == 'Methacrolein')
  allBOTH.filter$lifetime_jval[ind] = 1/mean(hall$jMAC_NoProductsSpecified_CAFS_HALL[ind2], na.rm=TRUE)/60/60
  ind = which(allBOTH.filter$names == 'Nitrous acid')
  allBOTH.filter$lifetime_jval[ind] = 1/mean(hall$jHNO2_OH_NO_CAFS_HALL[ind2], na.rm=TRUE)/60/60
  ind = which(allBOTH.filter$names == 'Nitric acid')
  allBOTH.filter$lifetime_jval[ind] = 1/mean(hall$jHNO3_OH_NO2_CAFS_HALL[ind2], na.rm=TRUE)/60/60
  ind = which(allBOTH.filter$names == 'Nitrogen dioxide')
  allBOTH.filter$lifetime_jval[ind] = 1/mean(hall$jNO2_NO_O3P_CAFS_HALL[ind2], na.rm=TRUE)/60/60
  ind = which(allBOTH.filter$names == 'Nitryl chloride')
  allBOTH.filter$lifetime_jval[ind] = 1/mean(hall$jClNO2_Cl_NO2_CAFS_HALL[ind2][ind2], na.rm=TRUE)/60/60
  
  # add lifetimes in parallel
  allBOTH.filter$lifetime_oh_1hz_hr = allBOTH.filter$lifetime_1hz_hr # don't overwrite this
  allBOTH.filter$lifetime_oh_5hz_hr = allBOTH.filter$lifetime_5hz_hr # don't overwrite this
  ind = which(is.finite(allBOTH.filter$lifetime_jval))
  allBOTH.filter$lifetime_1hz_hr[ind] = 1/( 1/allBOTH.filter$lifetime_oh_1hz_hr[ind] + 1/allBOTH.filter$lifetime_jval[ind])
  allBOTH.filter$lifetime_5hz_hr[ind] = 1/( 1/allBOTH.filter$lifetime_oh_5hz_hr[ind] + 1/allBOTH.filter$lifetime_jval[ind])
  
  allBOTH.filter$LifetimeCat = 1 # assume fast if I haven't found an OH rate yet
  ind = which(allBOTH.filter$lifetime_1hz_hr <= lifetimecutoff | allBOTH.filter$lifetime_5hz_hr <= lifetimecutoff)
  allBOTH.filter$LifetimeCat[ind] = 1
  ind = which(allBOTH.filter$lifetime_1hz_hr > lifetimecutoff | allBOTH.filter$lifetime_5hz_hr > lifetimecutoff)
  allBOTH.filter$LifetimeCat[ind] = 2
  
  # Need category to be the same
  allBOTH.filter$Category = allBOTH.filter$Category.5hz
  ind = which(!is.finite(allBOTH.filter$Category))
  allBOTH.filter$Category[ind] = allBOTH.filter$Category.1hz[ind]
  
  # For monoterpenes, make everything USEME = 2 except for WARNEKE
  ind = which(allBOTH.filter$names == 'Camphene' | allBOTH.filter$names == 'beta-Pinene/Myrcene' | 
                allBOTH.filter$names == 'beta-Pinene'| allBOTH.filter$names == 'Myrcene' | 
                allBOTH.filter$names == 'alpha-Pinene' | allBOTH.filter$names == 'Tricyclene')
  ind2 = which(allBOTH.filter$USEME[ind] == 1)
  allBOTH.filter$USEME[ind[ind2]] =2
  # ------------  Make a total VOC emission factor -------------
  ind1 = which(allBOTH.filter$variable =="Isoprene_NOAAPTR_ppbv_WARNEKE")
  isop= allBOTH.filter[ind1,]
  for (i in 1:length(isop$variable)){
    newline = isop[i,]
    newline2 = isop[i,]
    newline3 = isop[i,]
    
    # find all the VOCs for this fire/pass
    ind = which(allBOTH.filter$USEME == 1 & allBOTH.filter$Category == 1 &
                  allBOTH.filter$uniqueid == newline$uniqueid & allBOTH.filter$formula != 'N/A')
    ind2 = which(allBOTH.filter$USEME == 1 & allBOTH.filter$Category == 1 &
                  allBOTH.filter$uniqueid == newline$uniqueid& allBOTH.filter$LifetimeCat == 1 & allBOTH.filter$formula != 'N/A')
    ind3 = which(allBOTH.filter$USEME == 1 & allBOTH.filter$Category == 1 &
                   allBOTH.filter$uniqueid == newline$uniqueid& allBOTH.filter$LifetimeCat == 2 & allBOTH.filter$formula != 'N/A')
    
    # kludge for C8 aromatics (short-lived), MEK (long-lived), Methyl acetate (long-lived)
    ind1B = which(allBOTH.filter$names[ind] != 'Ethylbenzene' & allBOTH.filter$names[ind] != 'o-Xylene' & 
                    allBOTH.filter$names[ind] != 'm,p-Xylene' &
#                    allBOTH.filter$names[ind] != 'Methyl Ethyl Ketone' & 
                    allBOTH.filter$names[ind] != 'Methyl acetate' &
                    allBOTH.filter$names[ind] != ' Furan and fragments')
    ind2B = which(allBOTH.filter$names[ind2] != 'Ethylbenzene' & allBOTH.filter$names[ind2] != 'o-Xylene' & 
                    allBOTH.filter$names[ind2] != 'm,p-Xylene' & allBOTH.filter$names[ind2] != ' Furan and fragments')
    ind3B = which(allBOTH.filter$names[ind3] != 'Methyl acetate')
    
    tN1 = allBOTH.filter$FinalEF[ind[ind1B]]
    tN1ER = allBOTH.filter$FinalERtoCO[ind[ind1B]]
    
    tN1names = allBOTH.filter$names[ind[ind1B]]
    tNQ = as.data.frame(rbind(tN1))
    colnames(tNQ) = tN1names
      
    if (is.finite(tNQ$Formaldehyde) & is.finite(tNQ$Acetaldehyde) & 
        is.finite(tNQ$Methylglyoxal) &
        is.finite(tNQ$`2,3-Butanedione/2-Oxobutanal/1,4-Butanedial`) & 
        # long-lived
        is.finite(tNQ$`Methyl acetate/Ethyl formate/Hydroxyacetone`) & 
        is.finite(tNQ$`Acetone/Propanal`) & is.finite(tNQ$`Acetic acid/Glycolaldehyde`) & is.finite(tNQ$Methanol)){
      sFVOC = rowSums(tNQ, na.rm=TRUE)
      sRVOC = sum(tN1ER, na.rm=TRUE)
    } else{ sFVOC = NaN; sRVOC = NaN}
     # All
    allnames = allBOTH.filter$names[ind[ind1B]]
    allPI = allBOTH.filter$PI[ind[ind1B]]
    NMVOC.table = as.data.frame(cbind(allnames, allPI))
    write.csv(NMVOC.table, 'NMVOCtable.csv')
    newline$FinalEF = sFVOC
    newline$FinalERtoCO = sRVOC
    newline$variable = 'VOC'
    newline$names = 'VOC'
    newline$formula = 'N/A'
    newline$USEME = 2
    newline$FinalR2 =1 
    newline$LifetimeCat = 2 # to put at the end
    newline$mWs = 1000
    allBOTH.filter = rbind.fill(allBOTH.filter,newline)
    # ----- Short-lived
    tN1 = allBOTH.filter$FinalEF[ind2[ind2B]]
    tN1ER = allBOTH.filter$FinalERtoCO[ind2[ind2B]]
    
    tN1names = allBOTH.filter$names[ind2[ind2B]]
    tNQ = as.data.frame(rbind(tN1))
    colnames(tNQ) = tN1names
    if (is.finite(tNQ$Formaldehyde) & is.finite(tNQ$Acetaldehyde) & 
        is.finite(tNQ$Methylglyoxal) &
        is.finite(tNQ$`2,3-Butanedione/2-Oxobutanal/1,4-Butanedial`)){
      sFVOC = rowSums(tNQ, na.rm=TRUE)
      sRVOC = sum(tN1ER, na.rm=TRUE)
    } else{ sFVOC = NaN; sRVOC = NaN}
    
    newline2$FinalEF = sFVOC
    newline2$FinalERtoCO = sRVOC
    newline2$variable = 'Short-lived VOC'
    newline2$names = 'Short-lived VOC'
    newline2$formula = 'N/A'
    newline2$USEME = 2
    newline2$LifetimeCat = 2
    newline2$FinalR2 =1 
    #newline2$PI = ''
    newline2$mWs = 500
    allBOTH.filter = rbind.fill(allBOTH.filter,newline2)
    # Make a total VOC emission factor - long-lived
    tN1 = allBOTH.filter$FinalEF[ind3[ind3B]]
    tN1ER = allBOTH.filter$FinalERtoCO[ind3[ind3B]]
    
    tN1names = allBOTH.filter$names[ind3[ind3B]]
    tNQ = as.data.frame(rbind(tN1))
    colnames(tNQ) = tN1names
    if (is.finite(tNQ$`Methyl acetate/Ethyl formate/Hydroxyacetone`) & 
        is.finite(tNQ$`Acetone/Propanal`) & is.finite(tNQ$`Acetic acid/Glycolaldehyde`) & is.finite(tNQ$Methanol)){
      sFVOC = rowSums(tNQ, na.rm=TRUE)
      sRVOC = sum(tN1ER, na.rm=TRUE)
    } else{ sFVOC = NaN; sRVOC = NaN}
    
    newline3$FinalEF = sFVOC
    newline3$FinalERtoCO = sRVOC
    newline3$variable = 'Long-lived VOC'
    newline3$names = 'Long-lived VOC'
    newline3$formula = 'N/A'
    #newline3$PI = ''
    newline3$mWs = 800
    newline3$FinalR2 =1 
    newline3$USEME = 2
    newline3$LifetimeCat = 2
    allBOTH.filter = rbind.fill(allBOTH.filter,newline3)
  }
    
  # ----- save for processing
  save(allBOTH.filter, file='AllBOTH.filterP2.RData')
} else (load('AllBOTH.filterP2.RData'))
# For OA, fix variable
ind = which(allBOTH.filter$names == 'Organic aerosol')
if (length(ind) > 0){ allBOTH.filter$variable[ind] = 'OA_JIMENEZ'}
# ---- Give all passes the 5hz MCE ----
ff = unique(allBOTH.filter$uniqueid) 
allBOTH.filter$MCE = NaN
for (i in 1:length(ff)){
  ind = which(allBOTH.filter$uniqueid == ff[i] )
  allBOTH.filter$MCE[ind] = max(as.numeric(allBOTH.filter$mce.5hz[ind]), na.rm=TRUE) # they should all be the same, just fill in
}

# Get rid of any zero emission factors
ind = which(allBOTH.filter$FinalEF == 0)
allBOTH.filter$FinalEF[ind] = NaN

# Also need to check that we have a good CH4 EF for each
ind = which(allBOTH.filter$variable == 'CH4_DACOM_DISKIN' & !is.finite(allBOTH.filter$FinalEF)) 
badCH4passes = allBOTH.filter$uniqueid[ind]

cc = c()
for (i in 1:length(badCH4passes)){
  ind = which(allBOTH.filter$uniqueid == badCH4passes[i])
  allBOTH.filter$FinalEF[ind] = NaN
  allBOTH.filter$FinalERtoCO[ind] = NaN
}

# ------- Make a NOx as NO -------
ind = which(allBOTH.filter$formula == 'NO' & allBOTH.filter$USEME == 1)
tmpNO = allBOTH.filter[ind,]

ind = which(allBOTH.filter$formula == 'NO2' & allBOTH.filter$USEME == 1)
tmpNO2 = allBOTH.filter[ind,]

for (i in 1:length(tmpNO$variable)){
  ind = which(tmpNO2$uniqueid == tmpNO$uniqueid[i])
  tmpNO$FinalEF[i] = tmpNO$FinalEF[i] + tmpNO2$FinalEF[ind] * mWNO/mWNO2
  tmpNO$FinalERtoCO[i] = tmpNO$FinalERtoCO[i] + tmpNO2$FinalERtoCO[ind] 
  tmpNO$variable[i] = 'NOx (as NO)'
  tmpNO$names[i] = 'NOx (as NO)'
  tmpNO$formula[i] = 'NOx (as NO)'
  
}
allBOTH.filter = rbind(allBOTH.filter,tmpNO)

# ----- Make a pNO3 + HNO3 as HNO3 -------
ind = which(allBOTH.filter$formula == 'NO3' & allBOTH.filter$USEME == 1)# & allBOTH.filter$fuel != 'forest')
tmpNO3 = allBOTH.filter[ind,]

ind = which(allBOTH.filter$formula == 'HNO3' & allBOTH.filter$USEME == 1)#& allBOTH.filter$fuel != 'forest')
tmpHNO3 = allBOTH.filter[ind,]

for (i in 1:length(tmpNO3$variable)){
  ind = which(tmpHNO3$uniqueid == tmpNO3$uniqueid[i])
  tmpNO3$FinalEF[i] = tmpHNO3$FinalEF[i] + tmpNO3$FinalEF[ind] * mWHNO3/mWNO3
  tmpNO3$FinalERtoCO[i] = tmpHNO3$FinalERtoCO[i] + tmpNO3$FinalERtoCO[ind] 
  tmpNO3$variable[i] = 'pN (as HNO3)'
  tmpNO3$names[i] = 'pN (as HNO3)'
  tmpNO3$formula[i] = 'pN (as HNO3)'
  
}
allBOTH.filter = rbind(allBOTH.filter,tmpNO3)

# ---- Total carbon with CO+CO2+CH4 vs. everything -----
#plot(allBOTH.filter$TC2CO.1hz, allBOTH.filter$TC1CO.1hz, xlim=c(0,45), ylim=c(0,45))
#abline(lm(allBOTH.filter$TC1CO.1hz~allBOTH.filter$TC2CO.1hz+0))
# -------- min vs. max MCE across each fire -----
ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' &
              allBOTH.filter$fire != "Copper Breaks" &
              allBOTH.filter$fire != "Vivian" &
              allBOTH.filter$fire != "Invictus" &
             allBOTH.filter$fuel != 'coniferous/decidous')
fireMCEmin = aggregate(allBOTH.filter$MCE[ind], by=list(allBOTH.filter$fire[ind], allBOTH.filter$fuel[ind]), FUN='min', na.rm=TRUE)
fireMCEmax = aggregate(allBOTH.filter$MCE[ind], by=list(allBOTH.filter$fire[ind], allBOTH.filter$fuel[ind]), FUN='max', na.rm=TRUE)
fireMCEmin = as.data.frame(fireMCEmin)
colnames(fireMCEmin)=c('Fire','Fuel','MCEmin')
fireMCEmin$MCEmax = fireMCEmax$x
fireMCEmin$CT = NaN
for (i in 1:length(fireMCEmin$Fire)){
  ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' & is.finite(allBOTH.filter$FinalEF)&
                allBOTH.filter$fire == fireMCEmin$Fire[i])
  fireMCEmin$CT[i] = length(ind)
}
fireMCEmin$Diff = (fireMCEmin$MCEmax-fireMCEmin$MCEmin)*100/fireMCEmin$MCEmin
fireMCEmin = fireMCEmin[order(fireMCEmin$Diff,decreasing = TRUE),]

# Hamburger example
doHamburger = 0
if (doHamburger ==1){
  ind = which(allfires.5hz$fire == 'Hamburger' & allfires.5hz$Time_Start >= 72647)
  pass = allfires.5hz[ind,]
  cz = 1.2
  par(mfrow=c(1,1),mar = c(3, 4, 2, 4), cex=cz)
  yy=max(pass$CO2_7000_ppm, na.rm=TRUE)
  tt=3
  pass$Time_Start2 = seq(1,length(pass$Time_Start))
  plot(pass$Time_Start, pass$CO2_7000_ppm_DISKIN, ylab = "CO2, ppm",type='o',lwd=tt,
       main = "", xlab = "Time", ylim=c(min(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE), yy))
  par(new = TRUE, cex=cz)
  plot(pass$Time_Start, pass$CO_DACOM_DISKIN, xaxt = "n", yaxt = "n",type='o',lwd=tt,pch=16,cex=0.5,
       ylab = "", xlab = "", col = "red", lty = 2, ylim=c(0,max(pass$CO_DACOM_DISKIN, na.rm=TRUE)))
  axis(side = 4, cex=cz)
  mtext("CO, ppb", side = 4, line = 2, cex=cz)
  legend("topright", c("CO2", "CO"), col = c("black", "red"), bty='n',lty = c(1, 2))
       
  pass$Pass = NaN
  ind = which(pass$Time_Start >= 72647& pass$Time_Start < 72650.0 )
  pass$Pass[ind] = 1
  ind = which(pass$Time_Start >= 72651 & pass$Time_Start < 72655.0 )
  pass$Pass[ind] = 2
  ind = which(pass$Time_Start >= 72655 & pass$Time_Start < 72659.8 )
  pass$Pass[ind] = 3
  ind = which(pass$Time_Start >= 72870 & pass$Time_Start < 72879.0 )
  pass$Pass[ind] = 4
  ind = which(pass$Time_Start >= 72880 & pass$Time_Start < 72885.0 )
  pass$Pass[ind] = 5
  ind = which(pass$Time_Start >= 73305 & pass$Time_Start < 73335.0 )
  pass$Pass[ind] = 6
  ind = which(pass$Time_Start >= 73336 & pass$Time_Start < 73341.0 )
  pass$Pass[ind] = 7
  ind = which(pass$Time_Start >=  73341 & pass$Time_Start < 73347.0 )
  pass$Pass[ind] = 8
  ind = which(pass$Time_Start >=  73567 & pass$Time_Start < 73578.0)
  pass$Pass[ind] = 9
  ggplot(pass) + geom_point(aes(x=CO2_7000_ppm_DISKIN, y=CO_DACOM_DISKIN, col=factor(Pass)))+theme_classic()
  ind = which(is.nan(pass$Pass))
  pass$Pass[ind] = 1
  # Value used to transform the data
  coeff <- 15
  a.diff = max(pass$CO_DACOM_DISKIN, na.rm=TRUE) - min(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  b.diff = max(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE) - min(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE)
  a.min = min(pass$CO_DACOM_DISKIN, na.rm=TRUE)
  b.min= min(pass$CO2_7000_ppm_DISKIN, na.rm=TRUE)
  b = pass$CO2_7000_ppm_DISKIN
  ggplot(pass, aes(x=Time_Start2,y=CO_DACOM_DISKIN)) +
    geom_point( aes(col=factor(Pass)), size=2) + 
    geom_point(aes(y = (b - b.min) / b.diff * a.diff + a.min)) +
    #geom_line( aes(y=CO2_7000_ppm_DISKIN * coeff), size=2) +
    #scale_y_continuous(    name = "CO, ppb",
    #  sec.axis = sec_axis(~./coeff, name="CO2, ppm", limits=c(400,560))) + 
    scale_y_continuous(name='CO, ppb',sec.axis = sec_axis(trans = ~((. -a.min) * b.diff / a.diff) + b.min,
                                                          name = "CO2, ppm"))+
    theme_classic() +  theme(axis.title.y = element_text( size=13),
      axis.title.y.right = element_text( size=13)  ) + labs(col="Plume #")
}  
# Really want max CO for each plume
ind = which(allBOTH.filter$names == 'Carbon Monoxide')
maxCO.5hz = allBOTH.filter$maxval.5hz[ind]
maxCO.1hz = allBOTH.filter$maxval.1[ind]
passes = allBOTH.filter$uniqueid[ind]
allBOTH.filter$maxCO.5hz = NaN;allBOTH.filter$maxCO.1hz = NaN
for (i in 1:length(passes)){
  ind = which(allBOTH.filter$uniqueid == passes[i])
  allBOTH.filter$maxCO.5hz[ind] = maxCO.5hz[i]
  allBOTH.filter$maxCO.1hz[ind] = maxCO.1hz[i]
}

# Do we still have negative EFs????
ind = which(allBOTH.filter$FinalEF < 0 & allBOTH.filter$Category != 5 & is.finite(allBOTH.filter$Category))
allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN; allBOTH.filter$MCE[ind] = NaN
# --Remove aged all plumes not just blackwater
ind = which(allBOTH.filter$MAtoF_rq.5hz < 0.7 | is.nan(allBOTH.filter$MAtoF.5hz))
allBOTH.filter$MAtoF.5hz[ind] = .2 # kludge to not get rid of poorly calculated MAtoF
ind = which( allBOTH.filter$MAtoF.5hz > 0.2)
write.csv(allBOTH.filter[ind,],file='removedforMAtoF.csv')
allBOTH.filter$FinalEF[ind] = NaN; allBOTH.filter$FinalERtoCO[ind] = NaN; allBOTH.filter$MCE[ind] = NaN

dokludge = 1
if (dokludge == 1){
  # ------ kludge for particle number since I forgot to precalculate it -----
  indA = which(allBOTH.filter$variable == 'CNgt6nm' & allBOTH.filter$fire != 'Invictus')
  indB= which( allBOTH.filter$variable == 'Ngt100nm_LAS_stdPT' & allBOTH.filter$fire != 'Invictus')
  indC = which(allBOTH.filter$names == 'CCN_034' & allBOTH.filter$fire != 'Invictus')
  indD = which(allBOTH.filter$variable == 'CNgt20nm' & allBOTH.filter$fire != 'Invictus')
  indE = which(allBOTH.filter$variable == 'CNgt3nm'& allBOTH.filter$fire != 'Invictus')
  
  ind1 = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN'& allBOTH.filter$fire != 'Invictus')
  fuel = allBOTH.filter$fuel[ind1]
  tmpCContent = allBOTH.filter$Start*0 + 500 # Assuming 50 % C, 500 g/kg
  ind = which(fuel == 'corn' | fuel == 'soybean' | fuel == 'winter wheat' | fuel == 'rice')
  tmpCContent[ind] = 415.1 
  ind =which(fuel == 'grass' | fuel == 'shrub')
  tmpCContent[ind] = 462.7
  ind = which(fuel == 'pile' | fuel == 'slash') 
  tmpCContent[ind] =  511.1
  
  ind2 = which(allBOTH.filter$variable == 'CH4_DACOM_DISKIN'& allBOTH.filter$fire != 'Invictus')
  ind3 = which(allBOTH.filter$variable == 'CO2_ppb'& allBOTH.filter$fire != 'Invictus')
  c1 = allBOTH.filter$FinalERtoCO[indA]/(allBOTH.filter$FinalERtoCO[ind1]*allBOTH.filter$nCs[ind1]+
                                           allBOTH.filter$FinalERtoCO[ind2]*allBOTH.filter$nCs[ind2]+  # just CO, CO2, CH4
                                           allBOTH.filter$FinalERtoCO[ind3] *allBOTH.filter$nCs[ind3]) # particles/cm3/ppbC
  allBOTH.filter$FinalEF[indA] =  (c1*6.022E23/2.6E10/12)*tmpCContent[indA] #conversion at 273K and 1atm
  allBOTH.filter$PI[indA] = 'MOORE'
  
  c1 = allBOTH.filter$FinalERtoCO[indB]/(allBOTH.filter$FinalERtoCO[ind1]*allBOTH.filter$nCs[ind1]+
                                           allBOTH.filter$FinalERtoCO[ind2]*allBOTH.filter$nCs[ind2]+  # just CO, CO2, CH4
                                           allBOTH.filter$FinalERtoCO[ind3] *allBOTH.filter$nCs[ind3]) # particles/cm3/ppbC
  allBOTH.filter$FinalEF[indB] =  (c1*6.022E23/2.6E10/12)*tmpCContent[indB] #conversion at 273K and 1atm
  allBOTH.filter$PI[indB] = 'MOORE'
  
  c1 = allBOTH.filter$FinalERtoCO[indC]/(allBOTH.filter$FinalERtoCO[ind1]*allBOTH.filter$nCs[ind1]+
                                           allBOTH.filter$FinalERtoCO[ind2]*allBOTH.filter$nCs[ind2]+  # just CO, CO2, CH4
                                           allBOTH.filter$FinalERtoCO[ind3] *allBOTH.filter$nCs[ind3]) # particles/cm3/ppbC
  allBOTH.filter$FinalEF[indC] =  (c1*6.022E23/2.6E10/12)*tmpCContent[indC] #conversion at 273K and 1atm
  allBOTH.filter$PI[indC] = 'MOORE'
  
  c1 = allBOTH.filter$FinalERtoCO[indD]/(allBOTH.filter$FinalERtoCO[ind1]*allBOTH.filter$nCs[ind1]+
                                           allBOTH.filter$FinalERtoCO[ind2]*allBOTH.filter$nCs[ind2]+  # just CO, CO2, CH4
                                           allBOTH.filter$FinalERtoCO[ind3] *allBOTH.filter$nCs[ind3]) # particles/cm3/ppbC
  allBOTH.filter$FinalEF[indD] =  (c1*6.022E23/2.6E10/12)*tmpCContent[indD] #conversion at 273K and 1atm
  allBOTH.filter$PI[indD] = 'MOORE'
  
  c1 = allBOTH.filter$FinalERtoCO[indE]/(allBOTH.filter$FinalERtoCO[ind1]*allBOTH.filter$nCs[ind1]+
                                           allBOTH.filter$FinalERtoCO[ind2]*allBOTH.filter$nCs[ind2]+  # just CO, CO2, CH4
                                           allBOTH.filter$FinalERtoCO[ind3] *allBOTH.filter$nCs[ind3]) # particles/cm3/ppbC
  allBOTH.filter$FinalEF[indE] =  (c1*6.022E23/2.6E10/12)*tmpCContent[indE] #conversion at 273K and 1atm
  allBOTH.filter$PI[indE] = 'MOORE'
}

# --------- Supplemental ----------
# MEK vs. MEK/2-propanal
plot( toga.all$C4Carbonyls_NOAAPTR_TM,toga.all$MEK_ppt, pch=19, 
      xlim=c(0,7.5E3), ylim=c(0, 7.5E3),ylab='MEK, ppt', xlab='C4Carbonyls (NOAA PTRMS), ppt')
yy=as.numeric(toga.all$MEK_ppt)
xx=as.numeric(toga.all$C4Carbonyls_NOAAPTR_TM)
tt = lm(yy~xx)
abline(tt, col='black')
text(2000,6000,paste("TOGA slope = ", round(tt$coefficients[2], digits = 2)))

points(toga.all$C4Carbonyls_NOAAPTR_TM, toga.all$MEK_NOAAiWAS_TM, pch=19, col='red')
yy=as.numeric(toga.all$MEK_NOAAiWAS_TM)
tt2 = lm(yy~xx)
abline(tt2, col='red')
text(2000,5700,paste("iWAS slope = ", round(tt2$coefficients[2], digits = 2)), col='red')

points(toga.all$C4Carbonyls_NOAAPTR_TM, toga.all$MEK_WAS_TM, pch=19, col='blue')
yy=as.numeric(toga.all$MEK_WAS_TM)
tt2 = lm(yy~xx)
abline(tt2, col='red')
text(2000,5400,paste("WAS slope = ", round(tt2$coefficients[2], digits = 2)), col='blue')

# ------ Range of Dates ------
require(chron)
dS1= chron('2019-01-01', '00:00:00', format=c('y-m-d','h:m:s'))-1

ind = which(allfires.1hz$fuel == 'corn')
dS1+unique(allfires.1hz$Julian_date[ind], na.rm=TRUE)
ind = which(allfires.1hz$fuel == 'soybean')
dS1+unique(allfires.1hz$Julian_date[ind], na.rm=TRUE)
ind = which(allfires.1hz$fuel == 'winter wheat')
dS1+unique(allfires.1hz$Julian_date[ind], na.rm=TRUE)
ind = which(allfires.1hz$fuel == 'rice')
dS1+unique(allfires.1hz$Julian_date[ind], na.rm=TRUE)
ind = which(allfires.1hz$fuel == 'slash')
dS1+unique(allfires.1hz$Julian_date[ind], na.rm=TRUE)
ind = which(allfires.1hz$fuel == 'pile')
dS1+unique(allfires.1hz$Julian_date[ind], na.rm=TRUE)
ind = which(allfires.1hz$fuel == 'grass')
dS1+unique(allfires.1hz$Julian_date[ind], na.rm=TRUE)
ind = which(allfires.1hz$fuel == 'shrub')
dS1+unique(allfires.1hz$Julian_date[ind], na.rm=TRUE)


