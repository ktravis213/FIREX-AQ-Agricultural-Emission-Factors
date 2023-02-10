# -------- ### All  EFs ### -------- 
doclear=1
if (doclear == 1){ rm(list=ls()) } # clear all analysis
load('AgFires.RData')
source('getERsv2.R') ; require(GMCM)
source('makeplots.R'); source('speciateSpecies.R')
require(dplyr); require(plyr)
source('/Users/ktravis1/OneDrive - NASA/FIREX/plotSpeciesMCE.R')

R2filter = 0.75; R2filterCO = 0.90 # stricter criteria for CO as this defines the plume
R2Bot = 0.5 ; COcutoff = 400 # ppb 
doprocess=0;doprocessSTEP2 =0

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
# ---------- Get McCarty et al., 2011 EFs -------
fuel = c('corn',"rice",'soybean','winter wheat','grass','slash','pile')
co2 = c(1520,1520,1520,1630,1550,NaN,NaN)
co2SD = c(NaN,NaN,NaN,NaN,50,NaN,NaN)
co = c(53,52,69,55,91,NaN,NaN)
coSD = c(24,28,25,22,43,NaN,NaN)
ch4 = c(2.24,2.09,3.15,2.12,5.11,NaN,NaN)
ch4SD = c(0.49,0.94,1.00,1.20,4.32,NaN,NaN)
mccarty = as.data.frame(cbind(co2,co2SD,co,coSD,ch4,ch4SD))
mccarty$mce = (mccarty$co2/mWCO2)/(mccarty$co2/mWCO2 + mccarty$co/mWCO)
mccarty$fuel = fuel
# --------- Get Andreae emission factors ------------
f2 = '/Users/ktravis1/OneDrive - NASA/FIREX/Andreae-BB-EMFactors-14Apr2021_justtable1.csv'
andreae = read.csv(f2)
# --------- Get Akagi emission factors ---------
akagi=readxl::read_xlsx('Akagi_acp-11-4039-2011-supplement/Tables 1-5_4.27.11.xlsx')
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
  allBOTH$lifetime_5hz_hr = 1/(5E6*as.numeric(allBOTH$OHrate.5hz))/60/60
  allBOTH$lifetime_1hz_hr = 1/(5E6*as.numeric(allBOTH$OHrate.1hz))/60/60
  
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
  ind = which(round(allBOTH.filter$R2toCO.5hz, digits = 1) < R2filter)
  allBOTH.filter$EF1CO.5hz[ind] = NaN
  allBOTH.filter$ERtoCO.5hz[ind] = NaN
  allBOTH.filter$MCE[ind] = NaN
  
  ind = which(round(allBOTH.filter$R2toCO.1hz, digits=1) < R2filter)
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
      allBOTH.filter$Category.5hz[ind] = max(as.numeric(allBOTH.filter$Category.1hz[ind]), na.rm=TRUE) # they should all be the same, just fill in
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
  ff = unique(allBOTH.filter.allfuels$uniqueid)
  for (i in 1:length(ff)){
    ind = which(allBOTH.filter.allfuels$uniqueid == ff[i] )
    allBOTH.filter.allfuels$transect_dominant_fuel.5hz[ind] = max(allBOTH.filter.allfuels$transect_dominant_fuel.1hz[ind], na.rm=TRUE) # they should all be the same, just fill in
    allBOTH.filter.allfuels$transect_fuel_class.5hz[ind] = max(allBOTH.filter.allfuels$transect_fuel_class.1hz[ind], na.rm=TRUE) # they should all be the same, just fill in
    allBOTH.filter.allfuels$transect_fuel_confidence.5hz[ind] = max(allBOTH.filter.allfuels$transect_fuel_confidence.1hz[ind], na.rm=TRUE) # they should all be the same, just fill in
  }
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
  ind = which(allBOTH.filter.allfuels$transect_fuel_class.5hz >= 3 &
                allBOTH.filter.allfuels$transect_fuel_class.5hz <= 7)
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
    vars = c(oc$FinalEF[i]*oc$OAtoOC.1hz[i],bc$FinalEF[i],ammonium$FinalEF[i],
             sulf$FinalEF[i],nit$FinalEF[i],nrcl$FinalEF[i], pot$FinalEF[i])
    varsER = c(oc$FinalERtoCO[i]*oc$OAtoOC.1hz[i],bc$FinalERtoCO[i],ammonium$FinalERtoCO[i],
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
  ind3 = which(allBOTH.filter$names== '1,3,5-trimethylbenzene')
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
  ind = which(allBOTH.filter$names == 'Furan' & allBOTH.filter$PI == 'WARNEKE'  ) 
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
  
  # All but WAS agree so just keep warneke
  ind = which(allBOTH.filter$names == 'Methyl Ethyl Ketone' & allBOTH.filter$PI != 'WARNEKE')
  allBOTH.filter$USEME[ind] = 0
  
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
  # --- Average Toga AND BLAKE Isobutanal
  allBOTH.filter = mergelines(allBOTH.filter, 'iButanal_WAS_BLAKE','iButanal_ppt')
  # --- Average Toga AND BLAKE Butanal
  allBOTH.filter = mergelines(allBOTH.filter, 'Butanal_WAS_BLAKE','Butanal_ppt')
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
  # --- Average BLAKE, APEL methyl acetate
  allBOTH.filter = mergelines(allBOTH.filter, 'MeAcetate_ppt','MeAcetate_WAS_BLAKE')
  # Use warneke for the total VOC, but keep for table
  ind = which(allBOTH.filter$names == 'Methyl acetate')
  allBOTH.filter$USEME[ind] == 2
  # For betapinene/myrcene, lets add was together, then average with TOGA
  ind = which(allBOTH.filter$variable == 'Myrcene_WAS_BLAKE')
  ind2 = which(allBOTH.filter$variable == 'bPinene_WAS_BLAKE')
  allBOTH.filter$USEME[ind2] = 0
  allBOTH.filter$FinalEF[ind] = rowSums(cbind(allBOTH.filter$FinalEF[ind],allBOTH.filter$FinalEF[ind2]), na.rm=TRUE)
  allBOTH.filter$FinalERtoCO[ind] = rowSums(cbind(allBOTH.filter$FinalERtoCO[ind],allBOTH.filter$FinalERtoCO[ind2]), na.rm=TRUE)
  allBOTH.filter$variable[ind] = 'bPinene/Myrcene_WAS_BLAKE'
  allBOTH.filter$names[ind] = 'beta-Pinene/Myrcene'
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
  ind = which(allBOTH.filter$variable == 'HCN_NOAAPTR_ppbv_WARNEKE')
  allBOTH.filter$USEME[ind] =0
  allBOTH.filter = mergelines3(allBOTH.filter, 'HCN_NOAACIMS_VERES','HCN_WENNBERG','HCN_ppt')
  # --- Average Apel AND BLAKE propionitrile
  allBOTH.filter = mergelines(allBOTH.filter, 'PropNitrile_ppt','PropNitrile_WAS_BLAKE')
  # --- Average Apel AND BLAKE DMS
  allBOTH.filter = mergelines(allBOTH.filter, 'PropNitrile_ppt','PropNitrile_WAS_BLAKE')
  # --- Average BLAKE and APEL 2-methylfuran
  allBOTH.filter = mergelines(allBOTH.filter, 'x2MeFuran_WAS_BLAKE','x2MeFuran_ppt')
  # --- Average BLAKE and APEL 3-methylfuran
  allBOTH.filter = mergelines(allBOTH.filter, 'x3MeFuran_WAS_BLAKE','x3MeFuran_ppt')
  # keep for the table, use WARNEKE for the total
  # actually use TOGA per Gigi's table
  #ind = which(allBOTH.filter$names== '2-Methylfuran' & allBOTH.filter$USEME == 1)
  #allBOTH.filter$USEME[ind] =2
  #ind = which(allBOTH.filter$names== '3-Methylfuran' & allBOTH.filter$USEME == 1)
  #allBOTH.filter$USEME[ind] =2
  # --- Average BLAKE and APEL propane
  allBOTH.filter = mergelines(allBOTH.filter, 'Propane_WAS_BLAKE','Propane_NOAAiWAS_GILMAN')
  
  # --- Average GILMAN, BLAKE, APEL Furan
  allBOTH.filter = mergelines3(allBOTH.filter, 'Furan_ppt','Furan_WAS_BLAKE','Furan_NOAAiWAS_GILMAN')
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
  # ------ Lifetime category ---------
  ind = which(allBOTH.filter$variable == 'x23Butanedione_NOAAPTR_ppbv_WARNEKE') # fast photolysis
  allBOTH.filter$lifetime_5hz_hr[ind] = 4 # 
  allBOTH.filter$lifetime_1hz_hr[ind] = 4 # 1
  allBOTH.filter$lifetime[ind] = 4 # 1
  
  allBOTH.filter$LifetimeCat = 1 # assume fast if I haven't found an OH rate yet
  ind = which(allBOTH.filter$lifetime_1hz_hr <= 12 | allBOTH.filter$lifetime_5hz_hr <= 12)
  allBOTH.filter$LifetimeCat[ind] = 1
  ind = which(allBOTH.filter$lifetime_1hz_hr > 12 | allBOTH.filter$lifetime_5hz_hr > 12)
  allBOTH.filter$LifetimeCat[ind] = 2
  
  # Need category to be the same
  allBOTH.filter$Category = allBOTH.filter$Category.5hz
  ind = which(!is.finite(allBOTH.filter$Category))
  allBOTH.filter$Category[ind] = allBOTH.filter$Category.1hz[ind]
  
  # For monoterpenes, make everything USEME = 2 except for WARNEKE
  ind = which(allBOTH.filter$names == 'Camphene' | allBOTH.filter$names == 'beta-Pinene/Myrcene' | allBOTH.filter$names == 'alpha-Pinene' | allBOTH.filter$names == 'Tricyclene')
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

# Make a NOx as NO
ind = which(allBOTH.filter$formula == 'NO' & allBOTH.filter$USEME == 1)
tmpNO = allBOTH.filter[ind,]

ind = which(allBOTH.filter$formula == 'NO2' & allBOTH.filter$USEME == 1)
tmpNO2 = allBOTH.filter[ind,]

for (i in 1:length(tmpNO$variable)){
  ind = which(tmpNO2$uniqueid == tmpNO$uniqueid[i])
  tmpNO$FinalEF[i] = tmpNO$FinalEF[i] + tmpNO2$FinalEF[ind] * 30/46
  tmpNO$FinalERtoCO[i] = tmpNO$FinalERtoCO[i] + tmpNO2$FinalERtoCO[ind] 
  tmpNO$variable[i] = 'NOx (as NO)'
  tmpNO$names[i] = 'NOx (as NO)'
  tmpNO$formula[i] = 'NOx (as NO)'
  
}
allBOTH.filter = rbind(allBOTH.filter,tmpNO)
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
# --Remove aged blackwater
ind = which(allBOTH.filter$fire == 'BlackwaterRiver' & allBOTH.filter$MAtoF.5hz > 0.2)
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
# ----------------- AVERAGES ------------------------------
# Average by fire
#allBOTH.filter.avg.fire = aggregate(allBOTH.filter, by=list(allBOTH.filter$fire,allBOTH.filter$fuel, allBOTH.filter$variable), FUN='mean', na.rm=TRUE)
#allBOTH.filter.sd.fire = aggregate(allBOTH.filter, by=list(allBOTH.filter$fire, allBOTH.filter$fuel,allBOTH.filter$variable), FUN='sd', na.rm=TRUE)
# 
   # Average by  by individual fuel and species
#allBOTH.filter.avg.fire.fuel = aggregate(allBOTH.filter.avg.fire, by=list(allBOTH.filter.avg.fire$Group.2, allBOTH.filter.avg.fire$Group.3), FUN='mean', na.rm=TRUE)
#allBOTH.filter.sd.fire.fuel = aggregate(allBOTH.filter.avg.fire, by=list(allBOTH.filter.avg.fire$Group.2, allBOTH.filter.avg.fire$Group.3), FUN='sd', na.rm=TRUE)
# ----- Anything set to USEME = 0 here is already in the VOC EF, so this is double counted. Need to go back and fix at the end probably.
#       
#ind = which(allBOTH.filter$variable == 'PM1' & allBOTH.filter$names == 'Organic Carbon')
#allBOTH.filter <- allBOTH.filter[-c(ind), ]
q1 = 0.25; q2=0.75

# Just for ag + prescribed
ind = which(allBOTH.filter$fuel == 'corn' | allBOTH.filter$fuel == 'soybean' |  allBOTH.filter$fuel == 'rice' |
              allBOTH.filter$fuel == 'pile' | allBOTH.filter$fuel == 'slash' )
quantile(allBOTH.filter$MCE[ind], na.rm=TRUE)
allBOTH.filter.median = aggregate(allBOTH.filter[ind,], by=list(allBOTH.filter$variable[ind]), FUN='median', na.rm=TRUE)
allBOTH.filter.25 = aggregate(allBOTH.filter$FinalEF[ind], by=list(allBOTH.filter$variable[ind]), 'quantile',probs=c(q1),na.rm=TRUE)
allBOTH.filter.75 = aggregate(allBOTH.filter$FinalEF[ind], by=list(allBOTH.filter$variable[ind]), 'quantile',probs=c(q2),na.rm=TRUE)
allBOTH.filter.25ER = aggregate(allBOTH.filter$FinalERtoCO[ind], by=list(allBOTH.filter$variable[ind]), 'quantile',probs=c(q1),na.rm=TRUE)
allBOTH.filter.75ER = aggregate(allBOTH.filter$FinalERtoCO[ind], by=list(allBOTH.filter$variable[ind]), 'quantile',probs=c(q2),na.rm=TRUE)
allBOTH.filter.median$FinalEF_25 = allBOTH.filter.25$x
allBOTH.filter.median$FinalEF_75 = allBOTH.filter.75$x
allBOTH.filter.median$FinalERtoCO_25 = allBOTH.filter.25ER$x
allBOTH.filter.median$FinalERtoCO_75 = allBOTH.filter.75ER$x

ind = which(allBOTH.filter$fuel == 'corn')
allBOTH.filter.corn.median = aggregate(allBOTH.filter[ind,],
    by=list(allBOTH.filter$variable[ind]), FUN='median', na.rm=TRUE)
# PM1 breakdown
ind = which(allBOTH.filter.median$Category == 4)
tmp = allBOTH.filter.median[ind,]
# Average by individual fuel and species
# medians
allBOTH.filter.median.fuel = aggregate(allBOTH.filter, by=list(allBOTH.filter$fuel,allBOTH.filter$variable), FUN='median', na.rm=TRUE)
allBOTH.filter.mean.fuel = aggregate(allBOTH.filter, by=list(allBOTH.filter$fuel,allBOTH.filter$variable),   FUN='mean', na.rm=TRUE)
allBOTH.filter.sd.fuel = aggregate(allBOTH.filter, by=list(allBOTH.filter$fuel,allBOTH.filter$variable), FUN='sd', na.rm=TRUE)

allBOTH.filter.median.fuel$FinalEF_mean = allBOTH.filter.mean.fuel$FinalEF
allBOTH.filter.median.fuel$FinalEF_sd = allBOTH.filter.sd.fuel$FinalEF

allBOTH.filter.25.fuel = aggregate(allBOTH.filter$FinalEF, by=list(allBOTH.filter$fuel,allBOTH.filter$variable), 'quantile',probs=c(q1),na.rm=TRUE)
allBOTH.filter.75.fuel = aggregate(allBOTH.filter$FinalEF, by=list(allBOTH.filter$fuel,allBOTH.filter$variable), 'quantile',probs=c(q2),na.rm=TRUE)
allBOTH.filter.25.fuelER = aggregate(allBOTH.filter$FinalERtoCO, by=list(allBOTH.filter$fuel,allBOTH.filter$variable), 'quantile',probs=c(q1),na.rm=TRUE)
allBOTH.filter.75.fuelER = aggregate(allBOTH.filter$FinalERtoCO, by=list(allBOTH.filter$fuel,allBOTH.filter$variable), 'quantile',probs=c(q2),na.rm=TRUE)

allBOTH.filter.median.fuel$FinalEF_25 = allBOTH.filter.25.fuel$x
allBOTH.filter.median.fuel$FinalEF_75 = allBOTH.filter.75.fuel$x
allBOTH.filter.median.fuel$FinalERtoCO_25 = allBOTH.filter.25.fuelER$x
allBOTH.filter.median.fuel$FinalERtoCO_75 = allBOTH.filter.75.fuelER$x

# need to recover fire, kind, formula, and names
allBOTH.filter.median.fuel$variable = allBOTH.filter.median.fuel$Group.2
for (i in 1:length(allBOTH.filter.median.fuel$kind)){
  ind = which(allBOTH.filter$variable == allBOTH.filter.median.fuel$Group.2[i])
  allBOTH.filter.median.fuel$kind[i] = allBOTH.filter$kind[ind[1]]
  allBOTH.filter.median.fuel$formula[i] = allBOTH.filter$formula[ind[1]]
  allBOTH.filter.median.fuel$names[i] = allBOTH.filter$names[ind[1]]
  allBOTH.filter.median.fuel$PI[i] = allBOTH.filter$PI[ind[1]]
 # print(c(allBOTH.filter$variable[ind[1]], allBOTH.filter.median.fuel$Group.2[i]))
}
for (i in 1:length(allBOTH.filter.median$kind)){
  ind = which(allBOTH.filter$variable == allBOTH.filter.median$Group.1[i])
  allBOTH.filter.median$kind[i] = allBOTH.filter$kind[ind[1]]
  allBOTH.filter.median$formula[i] = allBOTH.filter$formula[ind[1]]
  allBOTH.filter.median$names[i] = allBOTH.filter$names[ind[1]]
  allBOTH.filter.median$PI[i] = allBOTH.filter$PI[ind[1]]
  # print(c(allBOTH.filter$variable[ind[1]], allBOTH.filter.median.fuel$Group.2[i]))
}
for (i in 1:length(allBOTH.filter.mean.fuel$kind)){
  ind = which(allBOTH.filter$variable == allBOTH.filter.mean.fuel$Group.2[i])
  allBOTH.filter.mean.fuel$kind[i] = allBOTH.filter$kind[ind[1]]
  allBOTH.filter.mean.fuel$formula[i] = allBOTH.filter$formula[ind[1]]
  allBOTH.filter.mean.fuel$names[i] = allBOTH.filter$names[ind[1]]
  allBOTH.filter.mean.fuel$PI[i] = allBOTH.filter$PI[ind[1]]
}


# --------- Get Plume Counts by FUEL ---------------------------------------------
# REDO
ff = unique(allBOTH.filter$uniqueid) # need to redo this
allBOTH.filter$MCE = NaN
for (i in 1:length(ff)){
  ind = which(allBOTH.filter$uniqueid == ff[i] )
  allBOTH.filter$MCE[ind] = max(as.numeric(allBOTH.filter$mce.5hz[ind]), na.rm=TRUE) # they should all be the same, just fill in
}
allBOTH.filter.median.fuel= getplumesANDmcebyfuel(allBOTH.filter.median.fuel, allBOTH.filter )

# --------------------------- Speciate these NOAA PTRMS species based on TOGA -----------------------------------------
# --------- Sum of m-xylene p-xylene o-xylene and ethyl benzene -------
# -------- AcetonePropanal --------------
# -------- MVK/MACR/2Butenals -------
allBOTH.filter.median.fuel = speciateSpecies(allBOTH.filter.median.fuel)
allBOTH.filter.median = speciateSpecies(allBOTH.filter.median)
 # C4Carbonyls_NOAAPTR_ppbv
ind = which(allBOTH.filter.median.fuel$variable == 'C4Carbonyls_NOAAPTR_ppbv_WARNEKE')
ind = which(allBOTH.filter.median.fuel$names == 'MEK')

# Make larger fuel categories
allBOTH.filter$fuelORIG = allBOTH.filter$fuel
allBOTH.filter$fuel2 = allBOTH.filter$fuel
ind = which(allBOTH.filter$fuel2 == 'corn' | allBOTH.filter$fuel2 == 'soybean' | allBOTH.filter$fuel2 == 'rice' |
              allBOTH.filter$fuel2 == 'winter wheat')
allBOTH.filter$fuel2[ind] = 'agriculture'
ind = which(allBOTH.filter$fuel2 == 'pile' | allBOTH.filter$fuel2 == 'slash' | allBOTH.filter$fuel2 == 'shrub')
allBOTH.filter$fuel2[ind] = 'prescribed'
ind = which(allBOTH.filter$fire == 'BlackwaterRiver')
allBOTH.filter$fuel2[ind] = 'Blackwater'

# MCE hist
ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN')
tmpCO = allBOTH.filter[ind,]
ind = which(tmpCO$fuel2 == 'Blackwater' & tmpCO$MAtoF.5hz > 0.2)
tmpCO$FinalEF[ind] = NaN
ind = which(tmpCO$fuel2 != 'forest' & tmpCO$fuel2 != 'coniferous/decidous' & is.finite(tmpCO$FinalEF))
min(tmpCO$MCE[ind])
# ------- names --------
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
# --------- which TOGA species dont correlate with CO? ---------------
cors.toga = c(); cors.toga.soybean = c(); pval.toga=c(); counts = c()
cc = colnames(toga.all)
for (i in 1:length(anames)){
  ind = which(cc == anames[i])
  yy=unlist(toga.all[,ind])
  ind2 = which(toga.all$fuel == 'corn')
  ind3 = which(toga.all$fuel == 'rice')
  ind4 = which(toga.all$fuel == 'soybean')
  ind5 = which(toga.all$fuel == 'grass')
  ind6 = which(toga.all$fuel == 'pile')
  ind7 = which(toga.all$fuel == 'slash')
  
  xx = as.numeric(toga.all$CO_DACOM_DISKIN_BECKY)
  xx = xx[ind2]; yyCORN = yy[ind2] # cut to just corn
  yyRICE = yy[ind2] ;yySOYBEAN= yy[ind2] ;yyGRASS = yy[ind2] ;yyPILE = yy[ind2] ;yySLASH= yy[ind2] 
  tt = which(is.finite(yyCORN)); tt2 = which(is.finite(yyRICE)); tt3 = which(is.finite(yySOYBEAN))
  tt4 = which(is.finite(yyGRASS)); tt5 = which(is.finite(yyPILE)); tt6 = which(is.finite(yySLASH))
  if (length(tt) > 2){
    cors = cor.test(xx,yyCORN)
    cors.toga = c(cors.toga, cors$estimate)
    pval.toga = c(pval.toga, cors$p.value)
    counts = c(counts, length(tt))
  } else{ cors.toga = c(cors.toga,NaN); pval.toga = c(pval.toga, NaN);    counts = c(counts, length(tt))}
  if (length(tt2) > 2){
    cors = cor.test(xx,yySOYBEAN)
    cors.toga = c(cors.toga, cors$estimate)
    pval.toga = c(pval.toga, cors$p.value)
    counts = c(counts, length(tt))
  } else{ cors.toga = c(cors.toga,NaN); pval.toga = c(pval.toga, NaN);    counts = c(counts, length(tt))}
 
  if (length(tt3) > 2){
    cors = cor.test(xx,yyCORN)
    cors.toga = c(cors.toga, cors$estimate)
    pval.toga = c(pval.toga, cors$p.value)
    counts = c(counts, length(tt))
  } else{ cors.toga = c(cors.toga,NaN); pval.toga = c(pval.toga, NaN);    counts = c(counts, length(tt))}
 
  if (length(tt4) > 2){
    cors = cor.test(xx,yyCORN)
    cors.toga = c(cors.toga, cors$estimate)
    pval.toga = c(pval.toga, cors$p.value)
    counts = c(counts, length(tt))
  } else{ cors.toga = c(cors.toga,NaN); pval.toga = c(pval.toga, NaN);    counts = c(counts, length(tt))}
 
  if (length(tt5) > 2){
    cors = cor.test(xx,yyCORN)
    cors.toga = c(cors.toga, cors$estimate)
    pval.toga = c(pval.toga, cors$p.value)
    counts = c(counts, length(tt))
  } else{ cors.toga = c(cors.toga,NaN); pval.toga = c(pval.toga, NaN);    counts = c(counts, length(tt))}
 
  if (length(tt6) > 2){
    cors = cor.test(xx,yyCORN)
    cors.toga = c(cors.toga, cors$estimate)
    pval.toga = c(pval.toga, cors$p.value)
    counts = c(counts, length(tt))
  } else{ cors.toga = c(cors.toga,NaN); pval.toga = c(pval.toga, NaN);    counts = c(counts, length(tt))}
 
   if (length(tt) > 2){
    par(mfrow=c(1,2))
    plot(xx,yyCORN, pch=19, main=anames[i], xlab='CO, ppb')
    text(550,max(yyCORN, na.rm=TRUE),paste("R=",round(cors$estimate,2)))
  }
}
cors.toga = as.data.frame(cors.toga)
cors.toga$name = anames
cors.toga$pval = pval.toga
cors.toga$counts = counts
write.csv(cors.toga, 'cors.togaCORN.csv')
# Don't report negative TOGA correlations
ind = which(is.finite(cors.toga$cors.toga) & cors.toga$cors.toga < 0)

# ----- which iWAS  species dont correlate with CO? ------
dev.off()
cors.iwas = c(); pval.iwas=c(); counts = c()
cc = colnames(iwas.all)
for (i in 1:length(gnames)){
  ind = which(cc == gnames[i])
  yy=unlist(iwas.all[,ind])
  tt = which(is.finite(yy))
  if (length(tt) > 0){
    par(mfrow=c(1,2))
    plot(as.numeric(iwas.all$CO_DACOM_DISKIN_GILMAN),yy, pch=19, main=gnames[i], xlab='CO, ppb')
    cors = cor.test(as.numeric(iwas.all$CO_DACOM_DISKIN_GILMAN), unlist(iwas.all[,ind]))

    cors.iwas = c(cors.iwas, cors$estimate)
    pval.iwas = c(pval.iwas, cors$p.value)
    counts = c(counts, length(tt))
    text(550,max(yy, na.rm=TRUE),paste("R=",round(cors$estimate,2)))
  } else{
    cors.iwas = c(cors.iwas,NaN)
    pval.iwas = c(pval.iwas, NaN)
    counts = c(counts, length(tt))
  }
}

cors.iwas = as.data.frame(cors.iwas)
cors.iwas$name = gnames
cors.iwas$pval = pval.iwas
cors.iwas$counts = counts

# ----- which WAS  species dont correlate with CO? ------
dev.off()
cors.was.ag = c(); pval.was.ag=c(); counts.ag = c()
cors.was.pb = c(); pval.was.pb=c(); counts.pb = c()
cc = colnames(was.all)
par(mfrow=c(2,2))
for (i in 1:length(bnames)){
  ind = which(cc == bnames[i])
  yy=unlist(was.all[,ind])
  ff = was.all$fuel
  ind1 = which(ff == 'slash' | ff == 'pile' | ff == 'grass')
  ff[ind1] = 'prescribed'
  ind2 = which(ff == 'corn' | ff == 'rice' | ff == 'soybean'| ff == 'wheat')
  ff[ind2] = 'agriculture'
  tt = which(is.finite(yy) & ff == 'agriculture')
  if (length(tt) > 2){
    plot(as.numeric(was.all$CO_DACOM_DISKIN_BLAKE[tt]),yy[tt], pch=19, main='Agriculture',ylab=bnames[i], xlab='CO, ppb')
    cors = cor.test(as.numeric(was.all$CO_DACOM_DISKIN_BLAKE[tt]), unlist(was.all[tt,ind]))
    
    cors.was.ag = c(cors.was.ag, cors$estimate)
    pval.was.ag = c(pval.was.ag, cors$p.value)
    counts.ag = c(counts.ag, length(tt))
    text(550,max(yy[tt], na.rm=TRUE),paste("R=",round(cors$estimate,2)))
  } else{
    cors.was.ag = c(cors.was.ag,NaN)
    pval.was.ag = c(pval.was.ag, NaN)
    counts.ag = c(counts.ag, length(tt))
  }
  tt = which(is.finite(yy) & ff == 'prescribed')
  if (length(tt) > 2){
    plot(as.numeric(was.all$CO_DACOM_DISKIN_BLAKE[tt]),yy[tt], pch=19,main='Prescribed', ylab=bnames[i], xlab='CO, ppb')
    cors = cor.test(as.numeric(was.all$CO_DACOM_DISKIN_BLAKE[tt]), unlist(was.all[tt,ind]))
    
    cors.was.pb = c(cors.was.pb, cors$estimate)
    pval.was.pb = c(pval.was.pb, cors$p.value)
    counts.pb = c(counts.pb, length(tt))
    text(550,max(yy[tt], na.rm=TRUE),paste("R=",round(cors$estimate,2)))
  } else{
    cors.was.pb = c(cors.was.pb,NaN)
    pval.was.pb = c(pval.was.pb, NaN)
    counts.pb = c(counts.pb, length(tt))
  }
  
}
cors.was.ag = as.data.frame(cors.was.ag)
cors.was.ag$name = bnames
cors.was.ag$pval = pval.was.ag
cors.was.ag$counts = counts.ag

cors.was.pb = as.data.frame(cors.was.pb)
cors.was.pb$name = bnames
cors.was.pb$pval = pval.was.pb
cors.was.pb$counts = counts.pb

# ------- What is unique from Blake
ind = which(allBOTH.filter.median.fuel$Group.1 == 'corn')
blake = allBOTH.filter.median.fuel[ind,]
blakenames = unique(blake$names)
ind = which(allBOTH.filter.median.fuel$Group.1 == 'corn' & allBOTH.filter.median.fuel$PI != 'BLAKE')
notblake = allBOTH.filter.median.fuel[ind,]
notblakenames = unique(notblake$names)

# --------- Get Andreae emission factors ------------
allBOTH.filter.median.fuel$AndreaeEF = NaN
allBOTH.filter.median.fuel$AndreaeEFsd = NaN
allBOTH.filter.median.fuel$AndreaeName = NaN
allBOTH.filter.median.fuel$AndreaeNN = NaN
fix=c()
for (i in 1:length(allBOTH.filter.median.fuel$AndreaeEF)){
  tt = strsplit(allBOTH.filter.median.fuel$Group.2[i], '_')
  ind = which(tt[[1]][1] == andreae$Katie)
  ind2 = which(allBOTH.filter.median.fuel$names[i] == andreae$Katie)
  if (length(ind) == 0 & length(ind2) == 0){fix=c(fix,allBOTH.filter.median.fuel$names[i])}
  print(c(tt[[1]][1], andreae$Katie[ind]))
  
  if (length(ind2) ==1){
    allBOTH.filter.median.fuel$AndreaeEF[i] = andreae$average[ind2]
    allBOTH.filter.median.fuel$AndreaeEFsd[i] = andreae$std.dev.[ind2]
    allBOTH.filter.median.fuel$AndreaeName[i] = andreae$Species[ind2]
    allBOTH.filter.median.fuel$AndreaeNN[i] = andreae$N[ind2]
  } else if (length(ind) ==1){
    allBOTH.filter.median.fuel$AndreaeEF[i] = andreae$average[ind]
    allBOTH.filter.median.fuel$AndreaeEFsd[i] = andreae$std.dev.[ind]
    allBOTH.filter.median.fuel$AndreaeName[i] = andreae$Species[ind]
    allBOTH.filter.median.fuel$AndreaeNN[i] = andreae$N[ind]
  } else if (length(ind) ==2){
    allBOTH.filter.median.fuel$AndreaeEF[i] = sum(andreae$average[ind])
    allBOTH.filter.median.fuel$AndreaeEFsd[i] = mean(andreae$std.dev.[ind], na.rm=TRUE)
    allBOTH.filter.median.fuel$AndreaeName[i] = paste(andreae$Species[ind[1]],andreae$Species[ind[2]])
    allBOTH.filter.median.fuel$AndreaeNN[i] = sum(andreae$N[ind])
  }
  if (length(ind) >2){print(">2")}
}
print(c("Done Andreae"))
# ---- ---
# # --------- Get Akagi emission factors ------------
# allBOTH.filter.median.fuel$AkagiEF = NaN
# allBOTH.filter.median.fuel$AkagiEFsd = NaN
# fix=c()
# for (i in 1:length(allBOTH.filter.median.fuel$AkagiEF)){
#   ind = which(allBOTH.filter.median.fuel$names[i] == akagi$Species)
#   if (length(ind) == 0 ){fix=c(fix,allBOTH.filter.median.fuel$names[i])}
#   print(c(allBOTH.filter.median.fuel$names[i], akagi$Species[ind]))
#   
#   if (length(ind) ==1){
#     allBOTH.filter.median.fuel$AkagiEF[i] = akagi$CropEF[ind]
#     allBOTH.filter.median.fuel$AkagiEFsd[i] = akagi$CropSD[ind]
#   }
#   if (length(ind) >2){print(">2")}
# }
# 
# ------ Want R2 to CO from the data------
# allBOTH.filter.median.fuel$RtoCOoverall = NaN
# for ( i in 1:length(allBOTH.filter.median.fuel$Group.2)){
#   cc1 = colnames(allfires.1hz); cc2 = colnames(allfires.5hz)
#   cc3 = colnames(toga.all); cc4 = colnames(was.all); cc5 = colnames(iwas.all)
#   # just do corn to remove any fuel dependence
#   ind1 = which(cc1 == allBOTH.filter.median.fuel$variable[i])
#   ind1B = which(allfires.1hz$fuel == 'corn')
#   ind2 = which(cc2 == allBOTH.filter.median.fuel$variable[i])
#   ind2B = which(allfires.5hz$fuel == 'corn')
#   ind3 = which(cc3 == allBOTH.filter.median.fuel$variable[i])
#   ind3B = which(toga.all$fuel == 'corn')
#   ind4 = which(cc4 == allBOTH.filter.median.fuel$variable[i])
#   ind4B = which(was.all$fuel == 'corn')
#   ind5 = which(cc5 == allBOTH.filter.median.fuel$variable[i])
#   ind5B = which(iwas.all$fuel == 'corn')
#   if (length(ind1) > 0){
#     tmp = as.numeric(allfires.1hz[ind1B,ind1])
#     iq = which(is.finite(tmp))
#     if (length(iq) > 2){tt = cor.test(as.numeric(allfires.1hz$CO_DACOM_DISKIN[ind1B]), tmp)}
#     allBOTH.filter.median.fuel$RtoCOoverall[i] = tt$estimate 
#   }
#   if (length(ind2) > 0){
#     tmp = as.numeric(allfires.5hz[ind2B,ind2])
#     iq = which(is.finite(tmp))
#     if (length(iq) > 2){tt = cor.test(as.numeric(allfires.5hz$CO_DACOM_DISKIN[ind2B]), tmp)}
#     allBOTH.filter.median.fuel$RtoCOoverall[i] = tt$estimate 
#   }
#   if (length(ind3) > 0){
#     tmp = (toga.all[ind3B,ind3])
#     tmp = unlist(tmp)
#     tmp = as.numeric(tmp)
#     iq = which(is.finite(tmp))
#     if (length(iq) > 2){tt = cor.test(as.numeric(toga.all$CO_DACOM_DISKIN_BECKY[ind3B]), tmp)}
#     allBOTH.filter.median.fuel$RtoCOoverall[i] = tt$estimate 
#   }
#   if (length(ind4) > 0){
#     tmp = as.numeric(was.all[ind4B,ind4])
#     iq = which(is.finite(tmp))
#     if (length(iq) > 2){tt = cor.test(as.numeric(was.all$CO_DACOM_DISKIN_BLAKE[ind4B]), tmp)}
#     allBOTH.filter.median.fuel$RtoCOoverall[i] = tt$estimate 
#   }
#   if (length(ind5) > 0){
#     tmp = as.numeric(as.numeric(iwas.all[ind5B,ind5]))
#     iq = which(is.finite(tmp))
#     if (length(iq) > 2){tt = cor.test(as.numeric(iwas.all$CO_DACOM_DISKIN_GILMAN[ind5B]), tmp)}
#     allBOTH.filter.median.fuel$RtoCOoverall[i] = tt$estimate 
#   }
# }

# ----------------------------------------------------------------------------
# --------- here run TestTableJustEF.R -------
# --------- then run SupplementXioaxi.R -------
# --------- then run SlopeDifferences.R -------
# --------- then run UniqueSpecies.R -------

# Get average size distribution by species

# ----------------------------------------------------------------------------

dothis=0
if (dothis == 1){
    ind = which(allBOTH.filter$fuel == 'corn')
    allBOTH.filter.avg.corn = aggregate(allBOTH.filter[ind,], by=list(allBOTH.filter$kind[ind]), FUN='mean', na.rm=TRUE)
    
    # Average by species
    allBOTH.filter.avg = aggregate(allBOTH.filter, by=list(allBOTH.filter$variable), FUN='mean', na.rm=TRUE)
    allBOTH.filter.sd = aggregate(allBOTH.filter, by=list(allBOTH.filter$variable), FUN='sd', na.rm=TRUE)
    allBOTH.filter.avg$EF1.5hz.sd = allBOTH.filter.sd$EF1.5hz
    allBOTH.filter.avg$EF1.1hz.sd = allBOTH.filter.sd$EF1.1hz
    
    # to compare with Xiaoxi's table
    #ind = which(allBOTH.filter.avg.fuel$Group.1 == 'corn' |
    #              allBOTH.filter.avg.fuel$Group.1 == 'rice' |
    #              allBOTH.filter.avg.fuel$Group.1 == 'soybean'  )
    #cornricesoybean = allBOTH.filter.avg.fuel[ind,]
    #ind = which(cornricesoybean$Group.2 == 'NO2_ACES_WOMACK') # why is Womack NO2 this brokenCH3OH_NOAAPTR_ppbv_WARNEKE
    #cornricesoybean$Group.1[ind]
    #cornricesoybean$mce.5hz[ind]
    #cornricesoybean$FinalEF[ind]
    #cornricesoybean$FinalEF.sd[ind]
    
    # Average by agriculture/land clearing and species
    # ---- Select agricultural fuels vs. silviculture -----
    ind = which(allBOTH.filter$fuel == 'corn' | allBOTH.filter$fuel == 'soybean' | allBOTH.filter$fuel == 'rice' | allBOTH.filter$fuel == 'winter wheat' )
    allBOTH.filter.ag = allBOTH.filter[ind,]
    ind = which(  allBOTH.filter$fuel == 'slash' | allBOTH.filter$fuel == 'pile' | allBOTH.filter$fuel == 'grass' )
    allBOTH.filter.sc = allBOTH.filter[ind,]
    # -------- agriculture
    allBOTH.filter.ag.avg = aggregate(allBOTH.filter.ag, by=list(allBOTH.filter.ag$variable), FUN='mean', na.rm=TRUE)
    allBOTH.filter.ag.sd = aggregate(allBOTH.filter.ag, by=list(allBOTH.filter.ag$variable), FUN='sd', na.rm=TRUE)
    allBOTH.filter.ag.avg$EF1.5hz.ag.sd = allBOTH.filter.ag.sd$EF1.5hz
    allBOTH.filter.ag.avg$EF1.1hz.ag.sd = allBOTH.filter.ag.sd$EF1.1hz
    allBOTH.filter.ag.avg$FinalEF.ag.sd = allBOTH.filter.ag.sd$FinalEF
    allBOTH.filter.ag.avg$mce.ag.sd = allBOTH.filter.ag.sd$mce.5hz
    # -------- land-clearing
    allBOTH.filter.sc.avg = aggregate(allBOTH.filter.sc, by=list(allBOTH.filter.sc$variable), FUN='mean', na.rm=TRUE)
    allBOTH.filter.sc.sd = aggregate(allBOTH.filter.sc, by=list(allBOTH.filter.sc$variable), FUN='sd', na.rm=TRUE)
    allBOTH.filter.sc.avg$EF1.5hz.sc.sd = allBOTH.filter.sc.sd$EF1.5hz
    allBOTH.filter.sc.avg$EF1.1hz.sc.sd = allBOTH.filter.sc.sd$EF1.1hz
    allBOTH.filter.sc.avg$FinalEF.sc.sd = allBOTH.filter.sc.sd$FinalEF
    allBOTH.filter.sc.avg$mce.sc.sd = allBOTH.filter.sc.sd$mce.5hz
    
    allBOTH.filter.ag.avg  = getplumesANDmce(allBOTH.filter.ag.avg, allBOTH.filter.ag )
    allBOTH.filter.sc.avg  = getplumesANDmce(allBOTH.filter.sc.avg, allBOTH.filter.sc )
    allBOTH.filter.avg     = getplumesANDmce(allBOTH.filter.avg, allBOTH.filter )
    
    
    # ------- Is there Distinct MCE ------- 
    # corn vs. others
    cornvsfuel= fuelvsfuel(allBOTH.filter,"corn","CO_DACOM_DISKIN")
    ricevsfuel= fuelvsfuel(allBOTH.filter,"rice","CO_DACOM_DISKIN")
    soyvsfuel= fuelvsfuel(allBOTH.filter,"soybean","CO_DACOM_DISKIN")
    pilevsfuel= fuelvsfuel(allBOTH.filter,"pile","CO_DACOM_DISKIN")
    slashvsfuel= fuelvsfuel(allBOTH.filter,"slash","CO_DACOM_DISKIN")
    grassvsfuel= fuelvsfuel(allBOTH.filter,"grass","CO_DACOM_DISKIN")
    
    distinctMCE = cbind(cornvsfuel, ricevsfuel, soyvsfuel,pilevsfuel, slashvsfuel, grassvsfuel)
    runplots =1
    if (runplots == 1){
      

    # all variables
      
      # Andreae comparison
      
      vars = unique(allBOTH.filter$variable)
      for (i in 1:length(vars)){
        vv = vars[i]
        # need at least 25% of the data?
        ind = which(allBOTH.filter$variable == vv & is.finite(allBOTH.filter$FinalEF) & allBOTH.filter$USEME == 1)
        ind3 = which(allBOTH.filter$variable == vv )
        nn = allBOTH.filter$names[ind3]
        if (length(ind)/length(ind3) > 0.25){
          ff=plotSpeciesMCE(allBOTH.filter, vv,nn[1],vv,vv)
          print(ff)
          ggsave(filename=paste(vv, '.ps', sep=''),ff, path='/Users/ktravis1/OneDrive - NASA/FIREX/FiguresMCE/')

        }
      }
      # ----------- Average agricultural table ---------- -----
    
    
    # Make a histogram plot or something of the correlation with MCE
    ind = which(is.finite(allBOTH.filter.avg$EF1.1hz) | is.finite(allBOTH.filter.avg$EF1.5hz))
    allBOTH.filter.avg = allBOTH.filter.avg[ind,]
    allBOTH.filter.sd = allBOTH.filter.sd[ind,]
    library(tidyverse)
    # --------- Get Andreae emission factors ------------
    allBOTH.filter.avg$AndreaeEF = NaN
    allBOTH.filter.avg$AndreaeEFsd = NaN
    allBOTH.filter.avg$AndreaeName = NaN
    allBOTH.filter.avg$AndreaeNN = NaN
    for (i in 1:length(allBOTH.filter.avg$AndreaeEF)){
      tt = strsplit(allBOTH.filter.avg$Group.1[i], '_')
      ind = which(tt[[1]][1] == andreae$Katie)
      print(c(tt[[1]][1], andreae$Katie[ind]))
      if (length(ind) ==1){
        allBOTH.filter.avg$AndreaeEF[i] = andreae$average[ind]
        allBOTH.filter.avg$AndreaeEFsd[i] = andreae$std.dev.[ind]
        allBOTH.filter.avg$AndreaeName[i] = andreae$Species[ind]
        allBOTH.filter.avg$AndreaeNN[i] = andreae$N[ind]
      }
      if (length(ind) ==2){
        allBOTH.filter.avg$AndreaeEF[i] = sum(andreae$average[ind])
        allBOTH.filter.avg$AndreaeEFsd[i] = mean(andreae$std.dev.[ind], na.rm=TRUE)
        allBOTH.filter.avg$AndreaeName[i] = paste(andreae$Species[ind[1]],andreae$Species[ind[2]])
        allBOTH.filter.avg$AndreaeNN[i] = sum(andreae$N[ind])
      }
      if (length(ind) >2){print(">2")}
    }
    dd = as.data.frame(cbind(allBOTH.filter.avg$lifetime_5hz_hr, allBOTH.filter.avg$lifetime_1hz_hr))
    allBOTH.filter.avg$Lifetime=  rowMeans(dd, na.rm=TRUE)
    dd2 = as.data.frame(cbind(allBOTH.filter.avg$OHrate.5hz, allBOTH.filter.avg$OHrate.1hz))
    allBOTH.filter.avg$OHrate=  rowMeans(dd2, na.rm=TRUE)
    
    tt = cbind(allBOTH.filter.avg$EF1.1hz, allBOTH.filter.avg$EF1.5hz)
    allBOTH.filter.avg$maxEF = rowMeans(tt, na.rm=TRUE)
    
    allBOTH.filter.avg$EF1.5hz.carbon = allBOTH.filter.avg$EF1.5hz*
      (allBOTH.filter.avg$nCs*12/allBOTH.filter.avg$mWs)
    # ------- Average Table ----------
    AverageTable <- tibble(
      Species = paste0(allBOTH.filter.avg$Group.1," "),
      #EF5Hz = round(allBOTH.filter.avg$EF1.5hz, digits=2),
      EF5Hz = round(allBOTH.filter.avg$ChosenEF.5hz, digits=2),
      SD5 = paste0("(",round(allBOTH.filter.sd$EF1.5hz, digits=2), ") "),
      n5=  paste0(round(allBOTH.filter.avg$COUNT_EF5hz, digits=2), " "),
      rMCE=  paste0(round(allBOTH.filter.avg$corMCE_5hz, digits=2), " "),
      
      #EF1Hz = round(allBOTH.filter.avg$EF1.1hz, digits=2),
      EF1Hz = round(allBOTH.filter.avg$ChosenEF.1hz, digits=2),
      SD = paste0("(",round(allBOTH.filter.sd$EF1.1hz, digits=2), ") "),
      n=  paste0(round(allBOTH.filter.avg$COUNT_EF1hz, digits=2), " "),
      rMCE1hz=  paste0(round(allBOTH.filter.avg$corMCE_1hz, digits=2), " "),
      OHrate = paste0(allBOTH.filter.avg$OHrate, " "),
      Lifetime = paste0(allBOTH.filter.avg$Lifetime, " "),
      LifetimeCat = allBOTH.filter.avg$LifetimeCat,
      Category = rowMeans(cbind(allBOTH.filter.avg$Category.1hz,allBOTH.filter.avg$Category.5hz), na.rm=TRUE),
      AndreaeEF = paste0(" ",allBOTH.filter.avg$AndreaeEF, " "),
      AndreaeEFsd = paste0(allBOTH.filter.avg$AndreaeEFsd, " "),
      AndreaeName = paste0(allBOTH.filter.avg$AndreaeName, ""),
      AndreaeNN = allBOTH.filter.avg$AndreaeNN,
      maxEF = allBOTH.filter.avg$maxEF,
      maxVal.1hz = allBOTH.filter.avg$maxval.1hz,
      maxVal.5hz = allBOTH.filter.avg$maxval.5hz,
      ERtoCO.5Hz = allBOTH.filter.avg$ERtoCO.5hz,
      ERtoCO.51z = allBOTH.filter.avg$ERtoCO.1hz,
      distinctCorn1 = allBOTH.filter.avg$CornvsFuel1hz,
      distinctCorn5 = allBOTH.filter.avg$CornvsFuel5hz,
      mw = allBOTH.filter.avg$mWs,
      nC = allBOTH.filter.avg$nCs)
    
    AverageTable$SpeciesShort = ""
    for (i in 1:length(AverageTable$SpeciesShort)){
      tt = strsplit(AverageTable$Species[i], "_")
      tt = tt[[1]][1]
      AverageTable$SpeciesShort[i]= tt
    }
    AverageTable = arrange(AverageTable, (Category), desc(LifetimeCat), (SpeciesShort),desc(maxEF))
    #AverageTable = getAverageTable(allBOTH.filter.avg,)
    write.table(AverageTable, file='AverageTable.csv')
    
    # ----- Want to compare average ag to akagi
    
    # plume length
    ind = which(all5hz$variable == 'CO_DACOM_DISKIN')
    plumelength = all5hz$Stop[ind] - all5hz$Start[ind] 
    
    # ---------- Average Table Ag -------------
    AverageTableAG = getAverageTable(allBOTH.filter.ag.avg, allBOTH.filter.ag.sd)
    ind = which(allBOTH.filter.ag.avg$Group.1 == "CO_DACOM_DISKIN")
    allBOTH.filter.ag.avg$mce.5hz[ind]
    allBOTH.filter.ag.sd$mce.5hz[ind]
    # --------- Get Akagi emission factors ------------
    allBOTH.filter.ag.avg$AkagiEF = NaN
    allBOTH.filter.ag.avg$AkagiEFsd = NaN
    allBOTH.filter.ag.avg$AkagiName = NaN
    allBOTH.filter.ag.avg$AkagiNN = NaN
    for (i in 1:length(allBOTH.filter.ag.avg$AkagiEF)){
      tt = strsplit(allBOTH.filter.ag.avg$Group.1[i], '_')
      ind = which(tt[[1]][1] == Akagi$Katie)
      print(c(tt[[1]][1], Akagi$Katie[ind]))
      if (length(ind) ==1){
        allBOTH.filter.ag.avg$AkagiEF[i] = Akagi$average[ind]
        allBOTH.filter.ag.avg$AkagiEFsd[i] = Akagi$std.dev.[ind]
        allBOTH.filter.ag.avg$AkagiName[i] = Akagi$Species[ind]
        allBOTH.filter.ag.avg$AkagiNN[i] = Akagi$N[ind]
      }
      if (length(ind) ==2){
        allBOTH.filter.ag.avg$AkagiEF[i] = sum(Akagi$average[ind])
        allBOTH.filter.ag.avg$AkagiEFsd[i] = mean(Akagi$std.dev.[ind], na.rm=TRUE)
        allBOTH.filter.ag.avg$AkagiName[i] = paste(Akagi$Species[ind[1]],Akagi$Species[ind[2]])
        allBOTH.filter.ag.avg$AkagiNN[i] = sum(Akagi$N[ind])
      }
      if (length(ind) >2){print(">2")}
    }
    dd = as.data.frame(cbind(allBOTH.filter.ag.avg$lifetime_5hz_hr, allBOTH.filter.ag.avg$lifetime_1hz_hr))
    allBOTH.filter.ag.avg$Lifetime=  rowMeans(dd, na.rm=TRUE)
    dd2 = as.data.frame(cbind(allBOTH.filter.ag.avg$OHrate.5hz, allBOTH.filter.ag.avg$OHrate.1hz))
    allBOTH.filter.ag.avg$OHrate=  rowMeans(dd2, na.rm=TRUE)
    
    tt = cbind(allBOTH.filter.ag.avg$EF1.1hz, allBOTH.filter.ag.avg$EF1.5hz)
    allBOTH.filter.ag.avg$maxEF = rowMeans(tt, na.rm=TRUE)
    
    allBOTH.filter.ag.avg$EF1.5hz.carbon = allBOTH.filter.ag.avg$EF1.5hz*
      (allBOTH.filter.ag.avg$nCs*12/allBOTH.filter.ag.avg$mWs)
    
    # ---------- Average Table Land clearing -------------
    AverageTableSC = getAverageTable(allBOTH.filter.sc.avg, allBOTH.filter.sc.sd)
    ind = which(allBOTH.filter.sc.avg$Group.1 == "CO_DACOM_DISKIN")
    allBOTH.filter.sc.avg$mce.5hz[ind]
    allBOTH.filter.sc.sd$mce.5hz[ind]
    
  
    # -----------
    
    # ---- Compare to Andreae, 2019 ------
    # ---- recreate Gigi's plot
    # make shorter names
    allBOTH.filter.avg$shortname <- 
      gsub(x = allBOTH.filter.avg$Group.1, pattern = "_NOAAPTR_ppbv_WARNEKE", replacement = "")  
    allBOTH.filter.avg$shortname <- 
      gsub(x = allBOTH.filter.avg$shortname, pattern = "_DACOM_DISKIN", replacement = "")  
    allBOTH.filter.avg$shortname <- 
      gsub(x = allBOTH.filter.avg$shortname, pattern = "_WENNBERG", replacement = "")  
    allBOTH.filter.avg$shortname <- 
      gsub(x = allBOTH.filter.avg$shortname, pattern = "_LIF_ROLLINS", replacement = "")  
    allBOTH.filter.avg$shortname <- 
      gsub(x = allBOTH.filter.avg$shortname, pattern = "_HUEY", replacement = "") 
    allBOTH.filter.avg$shortname <- 
      gsub(x = allBOTH.filter.avg$shortname, pattern = "_ISAF_HANISCO", replacement = "") 
    allBOTH.filter.avg$shortname <- 
      gsub(x = allBOTH.filter.avg$shortname, pattern = "_JIMENEZ", replacement = "") 
    allBOTH.filter.avg$shortname <- 
      gsub(x = allBOTH.filter.avg$shortname, pattern = "_NOAACIMS_VERES", replacement = "") 
    allBOTH.filter.avg$shortname <- 
      gsub(x = allBOTH.filter.avg$shortname, pattern = "_CAMS_pptv_FRIED", replacement = "") 
    allBOTH.filter.avg$shortname <- 
      gsub(x = allBOTH.filter.avg$shortname, pattern = "_ACES_WOMACK", replacement = "") 
    
    library(viridis)
    
    # ----- plot ----
    allBOTH.filter.avg= subset(allBOTH.filter.avg, select = -c( Group.1, Group.2))
    allBOTH.filter.avg$variableNUM = seq(1,length(allBOTH.filter.avg$variable))
    
    # ---- sort from highest to lowest EF
    allBOTH.filter.avg =allBOTH.filter.avg[order(allBOTH.filter.avg$maxEF,decreasing = TRUE),]
    ind = which(is.finite(allBOTH.filter.avg$maxEF) & allBOTH.filter.avg$variableNUM > 2
                & is.finite(allBOTH.filter.avg$AndreaeEF))
    tmp = allBOTH.filter.avg[ind,]
    tmp$ChosenEF = tmp$EF1.5hz
    ind = which(!is.finite(tmp$EF1.5hz))
    tmp$ChosenEF[ind] = tmp$EF1.1hz[ind]
    tmp$ChosenEF.sd = tmp$EF1.5hz.sd
    tmp$ChosenEF.sd[ind] = tmp$EF1.1hz.sd[ind]
    
    ggplot(tmp, aes(fill=factor(nCs),y=ChosenEF, x=variableNUM)) + 
      geom_bar(position="dodge", stat="identity") + theme_classic()+ # ylim(-3,7)+
      geom_errorbar(data=tmp,aes(x=variableNUM,ymin=ChosenEF-ChosenEF.sd, ymax=ChosenEF+ChosenEF.sd), width=.2,
                    position=position_dodge(.9), col='grey') +
      geom_point(data=tmp,aes(x=variableNUM, y=AndreaeEF), show.legend=FALSE)+
      geom_errorbar(data=tmp,aes(x=variableNUM,ymin=AndreaeEF-AndreaeEFsd, ymax=AndreaeEF+AndreaeEFsd), width=.2,
                    position=position_dodge(.9)) +
      scale_color_manual(values = c("purple", "blue", "green","yellow","red"))+
      guides(fill=guide_legend(title="nCs"))+
      theme(text = element_text(size=20))+
      geom_text(aes(x =variableNUM,y = 0,label = shortname), size=5,
                vjust = 0,  hjust = 1,angle = 90, nudge_y = -.01, nudge_x = 0.2)+ylab('EF, g/kg') 
    
    # samples per fire
    fires = as.data.frame(unique(allBOTH.filter.fire$Group.1))
    colnames(fires) = "fire"
    fires$npass = NaN; fires$minEF = NaN; fires$maxEF = NaN
    for (i in 1:length(fires$fire)){
      ind = which(allBOTH.filter$fire == fires$fire[i] & allBOTH.filter$variable == 'CO_DACOM_DISKIN' &
                    is.finite(allBOTH.filter$EF1.5hz))
      fires$npass[i] = length(allBOTH.filter$pass[ind])
      fires$minEF[i] = min(allBOTH.filter$EF1.5hz[ind])
      fires$maxEF[i] = max(allBOTH.filter$EF1.5hz[ind])
      
    }
    
    # can the EFs vary more across a fire then by fuel?
    ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN')
    ggplot(allBOTH.filter[ind,]) + 
      theme_classic()+geom_label(aes(x=mce.5hz, y=EF1.5hz, label=fire, col=fuel))+
      geom_point(aes(x=mce.5hz, y=EF1.5hz, col=fuel))
    ind = which(allBOTH.filter$fuel == 'corn' & allBOTH.filter$variable == 'CO_DACOM_DISKIN')
    plot(allBOTH.filter$mce.5hz[ind], allBOTH.filter$EF1.5hz[ind], pch=19, main='All Corn EFs')
    # then plot all 

    # Plot sampled fires
    co.ch4.821.5hz = getICARTTdataSIMPLE('Aircraft/FASTDATA/FIREXAQ-DACOM-5Hz_DC8_20190821_R1.ict')
    co.ch4.823.5hz = getICARTTdataSIMPLE('Aircraft/FASTDATA/FIREXAQ-DACOM-5Hz_DC8_20190823_R1.ict')
    co.ch4.826.5hz = getICARTTdataSIMPLE('Aircraft/FASTDATA/FIREXAQ-DACOM-5Hz_DC8_20190826_R1.ict')
    co.ch4.829.5hz = getICARTTdataSIMPLE('Aircraft/FASTDATA/FIREXAQ-DACOM-5Hz_DC8_20190829_R1.ict')
    co.ch4.830.5hz = getICARTTdataSIMPLE('Aircraft/FASTDATA/FIREXAQ-DACOM-5Hz_DC8_20190830_R1.ict')
    co.ch4.831.5hz = getICARTTdataSIMPLE('Aircraft/FASTDATA/FIREXAQ-DACOM-5Hz_DC8_20190831_R1.ict')
    co.ch4.903.5hz = getICARTTdataSIMPLE('Aircraft/FASTDATA/FIREXAQ-DACOM-5Hz_DC8_20190903_R1.ict')
    
    ind = which(is.na(allfires.5hz$Day_Of_Year_YANG ))
    ggplot(co.ch4.821.5hz)+geom_line(aes(x=Time_Start,y=CO_DACOM)) + scale_y_log10()+theme_classic()+
      geom_point(data=allfires.5hz[ind,], aes(x=Time_Start, y=CO_DACOM_DISKIN), col='red', cex=0.5)+
      geom_hline(aes(yintercept=400))+ggtitle("8/21")
    ind = which(allfires.5hz$Day_Of_Year_YANG == 235)
    ggplot(co.ch4.823.5hz)+geom_line(aes(x=Time_Start,y=CO_DACOM)) + scale_y_log10()+theme_classic()+
      geom_point(data=allfires.5hz[ind,], aes(x=Time_Start, y=CO_DACOM_DISKIN), col='red', cex=0.5)+
      geom_hline(aes(yintercept=400))+ggtitle("8/23")
    ind = which(allfires.5hz$Day_Of_Year_YANG == 238)
    ggplot(co.ch4.826.5hz)+geom_line(aes(x=Time_Start,y=CO_DACOM)) + scale_y_log10()+theme_classic()+
      geom_point(data=allfires.5hz[ind,], aes(x=Time_Start, y=CO_DACOM_DISKIN), col='red', cex=0.5)+
      geom_hline(aes(yintercept=400))+ggtitle("8/26")
    ind = which(allfires.5hz$Day_Of_Year_YANG == 241)
    ggplot(co.ch4.829.5hz)+geom_line(aes(x=Time_Start,y=CO_DACOM)) + scale_y_log10()+theme_classic()+
      geom_point(data=allfires.5hz[ind,], aes(x=Time_Start, y=CO_DACOM_DISKIN), col='red', cex=0.5)+
      geom_hline(aes(yintercept=400))+ggtitle("8/29")
    ind = which(allfires.5hz$Day_Of_Year_YANG == 242)
    ggplot(co.ch4.830.5hz)+geom_line(aes(x=Time_Start,y=CO_DACOM)) + scale_y_log10()+theme_classic()+
      geom_point(data=allfires.5hz[ind,], aes(x=Time_Start, y=CO_DACOM_DISKIN), col='red', cex=0.5)+
      geom_hline(aes(yintercept=400))+ggtitle("8/30")
    ind = which(allfires.5hz$Day_Of_Year_YANG == 243)
    ggplot(co.ch4.831.5hz)+geom_line(aes(x=Time_Start,y=CO_DACOM)) + scale_y_log10()+theme_classic()+
      geom_point(data=allfires.5hz[ind,], aes(x=Time_Start, y=CO_DACOM_DISKIN), col='red', cex=0.5)+
      geom_hline(aes(yintercept=400))+ggtitle("8/31")
    ind = which(allfires.5hz$Day_Of_Year_YANG == 246)
    ggplot(co.ch4.903.5hz)+geom_line(aes(x=Time_Start,y=CO_DACOM)) + scale_y_log10()+theme_classic()+
      geom_point(data=allfires.5hz[ind,], aes(x=Time_Start, y=CO_DACOM_DISKIN), col='red', cex=0.5)+
      geom_hline(aes(yintercept=400))+ggtitle("9/3")
    
    doblackwater=0
    # -------------------Blackwater--------------------------------------------
    if (doblackwater == 1){
      # Just Blackwater
      blackwater1 = Blackwaterriver.5hz.EF[order(Blackwaterriver.5hz.EF$variable),]
      blackwater2 = Blackwaterriver.1hz.EF[order(Blackwaterriver.1hz.EF$variable),]
      
      allBOTH.blackwater = merge(blackwater1, blackwater2, by=c('variable','fire','fuel', 'transect_source_fire_ID','mWs',
                                                                'nCs',"Start","Stop","StartO","StopO","pass","uniqueid"),
                                 all = TRUE, suffixes = c(".5hz", ".1hz"))
      # ----------- Provide the 5hz mce, MAtoF, and TC1.5hz for the 1hz data
      ff = unique(allBOTH.blackwater$uniqueid)
      for (i in 1:length(ff)){
        ind = which(allBOTH.blackwater$uniqueid == ff[i] )
        allBOTH.blackwater$mce.1hz[ind] = max(allBOTH.blackwater$mce.5hz[ind], na.rm=TRUE) # they should all be the same, just fill in
        allBOTH.blackwater$MAtoF.1hz[ind] = max(allBOTH.blackwater$MAtoF.5hz[ind], na.rm=TRUE) # they should all be the same, just fill in
        allBOTH.blackwater$TC1.1hz[ind] = max(allBOTH.blackwater$TC1.5hz[ind], na.rm=TRUE) # they should all be the same, just fill in
      }
      # Make EFs with the 5hz TC1.5hz from the 1hz ERtoCO?
      allBOTH.blackwater$EF1CO.1hz.5hz = NaN
      allBOTH.blackwater$EF1.1hz.5hz = NaN
      # lifetime category
      allBOTH.blackwater$LifetimeCat = 1 # assume fast if I haven't found an OH rate yet
      ind = which(allBOTH.blackwater$lifetime_1hz_hr <= 12 | allBOTH.blackwater$lifetime_5hz_hr <= 12)
      allBOTH.blackwater$LifetimeCat[ind] = 1
      ind = which(allBOTH.blackwater$lifetime_1hz_hr > 12 | allBOTH.blackwater$lifetime_5hz_hr > 12)
      allBOTH.blackwater$LifetimeCat[ind] = 2
      dofilter=1
      if (dofilter==1){
        # blackwater
        ind = which(allBOTH.blackwater$variable == 'CO_DACOM_DISKIN' & 
                      allBOTH.blackwater$R2toX.5hz >= R2filterCO & allBOTH.blackwater$maxval.5hz > 400 ) 
        goodpasses = allBOTH.blackwater$uniqueid[ind]
        
        cc = c()
        for (i in 1:length(goodpasses)){
          ind = which(allBOTH.blackwater$uniqueid == goodpasses[i])
          cc = c(cc, ind)
        }
        allBOTH.blackwater.filter = allBOTH.blackwater[cc,]
        
      }  
      # ----------- Blackwater Now, choose EFs that had the best correlation to either CO2 or CO
      allBOTH.blackwater.filter$ChosenEF.5hz = NaN; allBOTH.blackwater.filter$ChosenEF.R2.5hz = NaN
      allBOTH.blackwater.filter$ChosenEF.1hz = NaN; allBOTH.blackwater.filter$ChosenEF.1hz.5hz = NaN; allBOTH.blackwater.filter$ChosenEF.R2.1hz = NaN
      # Pick highest correlation for either the CO or CO2 EF
      for (i in 1:length(allBOTH.blackwater.filter$variable)){
        # Pick for 5hz data
        tmpR = c(allBOTH.blackwater.filter$R2toX.5hz[i],allBOTH.blackwater.filter$R2toCO.5hz[i])
        tmpEF= c(allBOTH.blackwater.filter$EF1.5hz[i],allBOTH.blackwater.filter$EF1CO.5hz[i])
        MAtoF = c(allBOTH.blackwater.filter$MAtoF.5hz[i],allBOTH.blackwater.filter$MAtoF.1hz[i])
        ind = which(is.finite(tmpR))
        if (length(ind) > 0){
          ind2 = which(tmpR == max(tmpR, na.rm=TRUE))
          if (length(ind2) ==1){
            allBOTH.blackwater.filter$ChosenEF.5hz[i] = tmpEF[ind2]
            allBOTH.blackwater.filter$ChosenEF.R2.5hz[i] = tmpR[ind2]
          } else{
            allBOTH.blackwater.filter$ChosenEF.5hz[i] = mean(tmpEF[ind2], na.rm=TRUE)
            allBOTH.blackwater.filter$ChosenEF.R2.5hz[i] = mean(tmpR[ind2], na.rm=TRUE)
          }
        }
        # Pick for 1hz data
        tmpR = c(allBOTH.blackwater.filter$R2toX.1hz[i],allBOTH.blackwater.filter$R2toCO.1hz[i])
        tmpEF= c(allBOTH.blackwater.filter$EF1.1hz[i],allBOTH.blackwater.filter$EF1CO.1hz[i])
        tmpEF2=c(allBOTH.blackwater.filter$EF1CO.1hz.5hz[i],allBOTH.blackwater.filter$EF1.1hz.5hz[i])
        ind = which(is.finite(tmpR))
        if (length(ind) > 0){
          ind2 = which(tmpR == max(tmpR, na.rm=TRUE))
          if (length(ind2) ==1){
            allBOTH.blackwater.filter$ChosenEF.1hz[i] = tmpEF[ind2]
            allBOTH.blackwater.filter$ChosenEF.1hz.5hz[i] = tmpEF[ind2]
            allBOTH.blackwater.filter$ChosenEF.R2.1hz[i] = tmpR[ind2]
            
          } else{
            allBOTH.blackwater.filter$ChosenEF.1hz[i] = mean(tmpEF[ind2], na.rm=TRUE)
            allBOTH.blackwater.filter$ChosenEF.1hz.5hz[i] = mean(tmpEF2[ind2], na.rm=TRUE)
            allBOTH.blackwater.filter$ChosenEF.R2.1hz[i] = mean(tmpR[ind2], na.rm=TRUE)
          }
        }
      }
      
      ind = which(allBOTH.blackwater.filter$transect_type.1hz  == 1 | allBOTH.blackwater.filter$transect_type.5hz == 1 )
      allBOTH.blackwater.filter=allBOTH.blackwater.filter[ind,]
      
      # ----------- filter
      ind = which(allBOTH.blackwater.filter$ChosenEF.R2.1hz < R2filter & allBOTH.blackwater.filter$MAtoF.5hz <= 0.2)
      allBOTH.blackwater.filter$ChosenEF.1hz[ind] = NaN
      ind = which(allBOTH.blackwater.filter$ChosenEF.R2.5hz < R2filter & allBOTH.blackwater.filter$MAtoF.5hz <= 0.2)
      allBOTH.blackwater.filter$ChosenEF.5hz[ind] = NaN
      allBOTH.blackwater.filter$FinalEF = allBOTH.blackwater.filter$ChosenEF.5hz
      allBOTH.blackwater.filter$FinalR2 = allBOTH.blackwater.filter$ChosenEF.R2.5hz
      ind = which(is.na(allBOTH.blackwater.filter$FinalEF) & is.finite(allBOTH.blackwater.filter$ChosenEF.1hz))
      allBOTH.blackwater.filter$FinalEF[ind] = allBOTH.blackwater.filter$ChosenEF.1hz[ind]
      allBOTH.blackwater.filter$FinalERtoCO = allBOTH.blackwater.filter$ERtoCO.5hz
      ind = which(is.na(allBOTH.blackwater.filter$FinalERtoCO) & is.finite(allBOTH.blackwater.filter$ERtoCO.1hz))
      allBOTH.blackwater.filter$FinalERtoCO[ind] = allBOTH.blackwater.filter$ERtoCO.1hz[ind]
      ind = which(allBOTH.blackwater.filter$maxval.5hz == -Inf)
      allBOTH.blackwater.filter$maxval.5hz[ind] = NaN
      write.csv(allBOTH.blackwater.filter,'Allboth.blackwater.filter.csv')
      ind = which(allBOTH.filter.allfuels$variable == 'OC_JIMENEZ' & allBOTH.filter.allfuels$fire != 'Blackwater')
      ind2 = which(allBOTH.blackwater.filter$variable == 'OC_JIMENEZ' ) ; allBOTH.blackwater.filter$fuel = 'Blackwater'
      
      ggplot(allBOTH.filter.allfuels[ind,])+geom_point(aes(x=mce.5hz,y=FinalEF, col=fuel), size=3)+theme_classic()+
        geom_point(data=allBOTH.blackwater.filter[ind2,],aes(x=mce.5hz,y=FinalEF), col='black')+xlab("MCE")+ylab("EF OC, g/kg")
      
      # average blackwater?
      #cut out aged
      ind = which(allBOTH.blackwater.filter$MAtoF.5hz <= 0.2) # remove aged
      allBOTH.blackwater.filter = allBOTH.blackwater.filter[ind,]
      allBOTH.blackwater.filter.avg =  aggregate(allBOTH.blackwater.filter, by=list(allBOTH.blackwater.filter$fire,allBOTH.blackwater.filter$fuel,
                                                                                    allBOTH.blackwater.filter$variable), FUN='mean', na.rm=TRUE)
      allBOTH.blackwater.filter.sd = aggregate(allBOTH.blackwater.filter, by=list(allBOTH.blackwater.filter$fire, allBOTH.blackwater.filter$fuel,
                                                                                  allBOTH.blackwater.filter$variable), FUN='sd', na.rm=TRUE)
      
      allBOTH.blackwater.filter.avg$Group.1 = allBOTH.blackwater.filter.avg$Group.3
      allBOTH.blackwater.filter.avg   = getplumesANDmce(allBOTH.blackwater.filter.avg, allBOTH.blackwater.filter )
      
      # ---------- Average Table Blackwater 
      AverageTableBlackwater = getAverageTable(allBOTH.blackwater.filter.avg, allBOTH.blackwater.filter.sd)
      ind = which(allBOTH.blackwater.filter.avg$Group.1 == "CO_DACOM_DISKIN")
      allBOTH.blackwater.filter.avg$mce.5hz[ind]
      allBOTH.blackwater.filter.sd$mce.5hz[ind]
      write.csv(AverageTableBlackwater, file='AverageTableBlackwater.csv')
    }
}

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

ind = which(is.finite(as.numeric(toga.all$C4Carbonyls_NOAAPTR_TM))& is.finite(as.numeric(toga.all$MEK_NOAAiWAS_TM)) &
             is.finite(as.numeric(toga.all$MEK_ppt)) & is.finite(as.numeric(toga.all$MEK_WAS_TM)))




  # ------
ind2 = which(allBOTH.filter$formula == 'HNCO'  & allBOTH.filter$USEME == 1)
ind1 = which(allBOTH.filter$formula == 'HNO2' & allBOTH.filter$USEME == 1 )
fuelshapes = c(19,15,17,18,7,8,4,2,16,12)
cbp1 <- c( "#E69F00", "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00",
           "#CC79A7","#000000","#999999","#CC0000")
fuellimits = c("corn","soybean","rice","winter wheat","grass","pile","slash","shrub")
ind = which(allfires.1hz$fuel != '?' & allfires.1hz$fuel != 'forest' & allfires.1hz$fuel != 'coniferous/decidous' &
              allfires.1hz$fuel != '' & allfires.1hz$CO_DACOM_DISKIN > 1500)
tmp = allfires.1hz[ind,]
tmp$fuel <- factor(tmp$fuel,
                            levels = c("winter wheat", "soybean","rice","corn","grass","slash","pile","shrub"))
ggplot(tmp)+geom_point(aes(x=HNCO_NOAACIMS_VERES, y=HNO2_NOAACIMS_VERES, col=fuel), size=3) + theme_classic()+
  scale_color_manual(values = c(cbp1 ))+
  scale_shape_manual(values =fuelshapes)

# HNCO/HNCO + HCN
ind = which(allBOTH.filter$formula == 'HCN' & allBOTH.filter$USEME == 1)
HCN = allBOTH.filter[ind,]
ind = which(allBOTH.filter$formula == 'HNCO' & allBOTH.filter$USEME == 1)
HNCO= allBOTH.filter[ind,]
ind = which(allBOTH.filter$formula == 'NH3' & allBOTH.filter$USEME == 1)
NH3= allBOTH.filter[ind,]
ind = which(allBOTH.filter$formula == 'NH4' & allBOTH.filter$USEME == 1)
NH4 =allBOTH.filter[ind,]
ind = which(allBOTH.filter$formula == 'NO2' & allBOTH.filter$USEME == 1& allBOTH.filter$fuel != 'forest' & allBOTH.filter$fuel != 'coniferous/decidous')
NO2 =allBOTH.filter[ind,]
ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' & allBOTH.filter$USEME == 1 & allBOTH.filter$fuel != 'forest'& allBOTH.filter$fuel != 'coniferous/decidous')
CO=allBOTH.filter[ind,]
ggplot()+geom_point(aes(x= HCN$MCE, y=HNCO$FinalERtoCO/(HCN$FinalERtoCO), col=HNCO$fuel), size=3)

ggplot()+geom_point(aes(x= HCN$MCE, y=HCN$FinalERtoCO/(NH3$FinalERtoCO), col=HNCO$fuel), size=3)
CO$fuel <- factor(CO$fuel,levels = c("winter wheat", "soybean","rice","corn","grass","slash","pile","shrub"))

p1=ggplot()+geom_point(aes(x= CO$MCE, y=NO2$FinalERtoCO/(CO$FinalERtoCO), col=CO$fuel), size=3)+
  theme_classic()+xlab("MCE")+ylab(expression(paste("NO"[2],"/CO")))+labs(col="",shape="")+
  theme(legend.background=element_blank())+
  theme(text = element_text(size = 20)) +
  scale_color_manual(values = c(cbp1 ))

p2=ggplot()+geom_point(aes(x= CO$MCE, y=NO2$FinalERtoCO, col=CO$fuel), size=3)+
  theme_classic()+xlab("MCE")+ylab(expression(paste("NO"[2])))+labs(col="",shape="")+
  theme(legend.background=element_blank())+
  theme(text = element_text(size = 20)) +
  scale_color_manual(values = c(cbp1 ))

tmp = merge(CO, NO2, by="uniqueid", suffixes = c("CO","NO2"))
FinalERtoCOCO = as.numeric(tmp$FinalERtoCOCO)
FinalERtoCONO2 = as.numeric(tmp$FinalERtoCONO2)
MCECO = as.numeric(tmp$MCECO)
tt = as.data.frame(cbind(FinalERtoCOCO, FinalERtoCONO2, MCECO))
tt$fuel = tmp$fuelCO
tt.med = aggregate(tt, by=list(tmp$fuelCO), FUN='median', na.rm=TRUE)
f1 = 'Aircraft/1s_MERGES/firexaq-mrg01-dc8_merge_20190830_RL.ict'
aug30merge = read.csv(f1, skip=674, header=TRUE)

ind = which(allfires.1hz$fuel != 'forest' & allfires.1hz$fuel != '?' & allfires.1hz$fuel != 'coniferous/decidous' & allfires.1hz$fuel != 'savannah' & allfires.1hz$fuel != '')
ggplot(allfires.1hz[ind,])+geom_point(aes(x=NO2_CL_RYERSON, y=CH4_DACOM_DISKIN, col=fuel), size=3)+xlab("NO2, ppb")+ylab("CH4, ppb")+theme_classic()

ind = which(allBOTH.filter$names == 'NH3' & allBOTH.filter$USEME == 1)
NH3= allBOTH.filter[ind,]
ind = which(allBOTH.filter$formula == 'NH4' & allBOTH.filter$USEME == 1)
NH4 =allBOTH.filter[ind,]

ind = which(allBOTH.filter$names == 'Furan' & allBOTH.filter$USEME == 1)
FURAN= allBOTH.filter[ind,]
ind = which(allBOTH.filter$formula == 'C2H2' & allBOTH.filter$USEME == 1)
C2H2 =allBOTH.filter[ind,]
FURANC2H2 = merge(FURAN, C2H2, suffixes = c("FURAN","C2H2"), by="uniqueid")
#FURANC2H2$fuelFURAN <- factor(FURANC2H2$fuelFURAN,levels = c("winter wheat", "soybean","rice","corn","grass","slash","pile","shrub"))
ind =which(FURANC2H2$fuelC2H2 != 'forest')
ggplot(FURANC2H2[ind,]) + geom_point(aes(x=MCEFURAN,y=(FinalERtoCOC2H2/0.0393)/(FinalERtoCOFURAN/0.0159),  col=fuelFURAN))+theme_classic()+
  scale_color_manual(values = c(cbp1 ))
tmp = aggregate(FURANC2H2, by=list(FURANC2H2$fuelC2H2), FUN='median', na.rm=TRUE)
tmp$TFACT = (tmp$FinalERtoCOC2H2/0.0393)/(tmp$FinalERtoCOFURAN/0.0159)

tmp$BOTH = tmp$FinalERtoCOFURAN+tmp$FinalERtoCOC2H2
newDATA = as.data.frame(rbind( tmp$MCEC2H2, tmp$FinalERtoCOC2H2/0.0393, tmp$FinalERtoCOFURAN/0.0159))
colnames(newDATA) = tmp$Group.1
rownames(newDATA) = c("MCE","C2H2","Furan")
library(ggplot2); library(reshape)
test =  melt(newDATA)
test$CHAR = rep(c("MCE","C2H2","Furan"),9)
test$CHAR2 = rep(c("MCE","var","var"),9)

# Basic barplot
ind = which(test$CHAR == 'MCE')
ind2 = which(test$CHAR == 'C2H2')
ind3 = which(test$CHAR == 'Furan')

p<-ggplot() +
  geom_bar(data=test[ind,], aes(x=test$value[ind], y=test$value[ind2]),stat="identity")+
  geom_bar(data=test[ind,], aes(x=test$value[ind], y=test$value[ind3]),stat="identity")

ind = which(test$CHAR != 'MCE' & test$variable != 'forest' & test$variable != 'winter wheat' & test$variable != 'shrub')
ind2 = which(test$CHAR == 'MCE'& test$variable != 'forest' & test$variable != 'winter wheat' & test$variable != 'shrub')

p=ggplot(data=test[ind,], aes(x=variable, y=value, fill=CHAR, label=c(test$value[c(ind2,ind2)]))) +
  geom_bar(stat="identity")
p+theme_classic()+ylab('ER, ppt/ppb')# + coord_flip()


# --------- EC/OC
ind = which(allBOTH.filter$variable == 'OC_JIMENEZ' & allBOTH.filter$fuel != 'forest' )#& allBOTH.filter$fuel == 'corn')
OC = allBOTH.filter$FinalEF[ind]
ind = which(allBOTH.filter$variable == 'BC_SCHWARZ' & allBOTH.filter$fuel != 'forest' )#& allBOTH.filter$fuel == 'corn')
EC = allBOTH.filter$FinalEF[ind]
fuel = allBOTH.filter$fuel[ind]
MCE = allBOTH.filter$MCE[ind]
tmp =as.data.frame( cbind(OC,EC,fuel,MCE))
tmp$fuel <- factor(tmp$fuel,
                   levels = c("winter wheat", "soybean","rice","corn","grass","slash","pile","shrub"))
ggplot(tmp)+geom_point(aes(x=as.numeric(MCE), y=as.numeric(EC)/as.numeric(OC), col=fuel), size=3)+theme_classic()+
  geom_point(data=xiaoxi,aes(x=xiaoxi$MCE, xiaoxi$BC/xiaoxi$OC), col='green', pch=0, size=3, lwd=3)+
  scale_color_manual(values = c(cbp1 ))+ylab("EC/OC")+xlab("MCE")
cor.test(as.numeric(tmp$EC)/as.numeric(tmp$OC), as.numeric(tmp$MCE))
# ------------------------
ggplot(tmp)+geom_point(aes(x=HNCO_NOAACIMS_VERES, y=HNO2_NOAACIMS_VERES, col=fuel), size=3) + theme_classic()+
  scale_color_manual(values = c(cbp1 ))+
  scale_shape_manual(values =fuelshapes)

# Range of Dates
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

ind = which(allBOTH.filter$fire == 'Ant' |
              allBOTH.filter$fire == 'Blanket' | allBOTH.filter$fire == 'Chips' |
              allBOTH.filter$fire == 'Dip' | allBOTH.filter$fire == 'ChipsDip' |
              allBOTH.filter$fire == 'Escargot' | allBOTH.filter$fire == 'Frisbee' |
              allBOTH.filter$fire == 'Guac' | allBOTH.filter$fire == 'Hamburger' |
              allBOTH.filter$fire == 'IPA' | allBOTH.filter$fire == 'Jello' |
              allBOTH.filter$fire == 'Mustard' | allBOTH.filter$fire == 'Pushmataha' |
              allBOTH.filter$fire == 'Kebab' | allBOTH.filter$fire == 'Limoncello')

#aug23 = allBOTH.filter[ind,]
#ind = which(aug23$formula == 'CO' | aug23$formula == 'CO2' | aug23$formula == 'CH4' | aug23$formula == 'EC' |
#              aug23$PI == 'MOORE' | aug23$formula == 'OC')
#aug23 = aug23[ind,]
#write.csv(aug23,'aug23.csv')
}

