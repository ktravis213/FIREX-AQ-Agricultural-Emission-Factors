# -------- ### All  EFs ### -------- 
doclear=1
if (doclear == 1){ rm(list=ls()) } # clear all analysis
load("AgFires.RData")
source('getERs.R') ; require(GMCM)
source('makeplots.R'); source('speciateSpecies.R')
require(dplyr); require(plyr)
source('/Users/ktravis1/OneDrive - NASA/FIREX/plotSpeciesMCE.R')

R2filter = 0.75; R2filterCO = 0.90 # stricter criteria for CO as this defines the plume
R2Bot = 0.5 ; COcutoff = 400 # ppb 
doprocess=0;doprocessSTEP2 =0

fuelshapes = c(19,15,17,18,7,8,4,2,16,12)
fuellimits = c("corn","soybean","rice","winter wheat","grass","pile","slash","shrub")

cbp1 <- c( "#E69F00", "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", 
           "#CC79A7","#000000","#999999","#CC0000")
# ----- +++++++++++ ANALYSIS +++++++++++++ ------
# How do emission factors change at a single fire over time, if at all?
# How do emission factors vary across various fires sampled?
#  Are there metrics (i.e. MCE etc) that provide some explanatory fpower for the variation that may be present?
require(reshape); require(ggmap) ; require(OrgMassSpecR); library(readxl) ;require(plyr) ; library(ggbrace)
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
akagi$PastureEF = as.numeric(akagi$PastureEF)
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
  require(ggmap)
  ind = which(all5hz$fuel != "?" & all5hz$fuel != "forest" & all5hz$fuel != "savannah?" &
                all5hz$fuel != "coniferous/decidous" & all5hz$fuel != "house" & 
                all5hz$fuel != 'Understory mixed, shrub,rice' &
                #all5hz$fuel != 'shrub' &
                all5hz$fuel != 'Post-blackjack oak forest' &
                all5hz$variable == 'CO_DACOM_DISKIN' &
                as.numeric(all5hz$R2toX) >= R2filterCO & as.numeric(all5hz$maxval) > COcutoff ) 
  
  all5hz.map = all5hz[ind,]
  all5hz.map$maxval = as.numeric(all5hz.map$maxval)
  all5hz.map$catCO = NaN
  ind = which(all5hz.map$maxval < 1E3)
  all5hz.map$catCO[ind] = 1 # 52/237
  ind = which(all5hz.map$maxval >= 1E3 & all5hz.map$maxval < 5E3)
  all5hz.map$catCO[ind] = 5#123
  ind = which(all5hz.map$maxval >= 5E3 & all5hz.map$maxval < 10E3)
  all5hz.map$catCO[ind] = 10 # 42
  ind = which(all5hz.map$maxval >= 10E3 & all5hz.map$maxval < 20E3)
  all5hz.map$catCO[ind] = 30 # 14
  ind = which(all5hz.map$maxval >= 20E3 & all5hz.map$maxval < 40E3)
  all5hz.map$catCO[ind] = 40 # 6
  all5hz.map =all5hz.map[order(all5hz.map$maxval,decreasing = TRUE),]
  
  domap=0
  if (domap == 1){
    map.us <- get_map(c(-98.5,30,-82,41))
    mapfire = ggmap(map.us) +
      geom_point(data = all5hz.map, aes(x = lon,y = lat,colour =fuel, size=catCO))+ 
      ggtitle("Ag Fires, 8/21-9/03 ")+  scale_color_manual(values = c(cbp1 ), 
                                                           limits=fuellimits)
    ind = which(all5hz.map$fuel == 'corn' | all5hz.map$fuel == 'rice' | 
                  all5hz.map$fuel == 'soybean' | all5hz.map$fuel == 'winter wheat')   
    ind = which(all5hz.map$fuel == 'grass' | all5hz.map$fuel == 'shrub')   
    ind = which(all5hz.map$fuel == 'slash' | all5hz.map$fuel == 'pile')   
    
    ggmap(map.us) +
      geom_point(data = all5hz.map[ind,], aes(x = lon,y = lat,colour =fuel, size=catCO))
    
    # ---- MCE  histogram -----
    ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' & allBOTH.filter$fuel != 'coniferous/decidous' & allBOTH.filter$fuel != 'forest')
    MCEhist = ggplot(allBOTH.filter[ind,], aes(x=MCE, col=fuel, fill=fuel)) + geom_histogram()+theme_classic()+
      scale_color_manual(values = c(cbp1 ),limits=fuellimits) + scale_fill_manual(values = c(cbp1 ),limits=fuellimits) +
      theme(legend.position = "top")+labs(fill="")
    ggsave(filename = 'MCEhist.pdf',MCEhist,width = 6,height = 6,units = 'in')
    
    ggsave(filename = 'Figure1_FireMap.pdf',mapfire,width = 6,height = 6,units = 'in')
  }
  # ------- How many fuels? --------
  fuels = unique(all5hz.map$fuel)
  length(unique(all5hz.map$fire))
  
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
  
  ind = which(allBOTH.filter$R2toCO.5hz < R2filter)
  allBOTH.filter$EF1CO.5hz[ind] = NaN
  allBOTH.filter$ERtoCO.5hz[ind] = NaN
  allBOTH.filter$MCE[ind] = NaN
  
  ind = which(allBOTH.filter$R2toCO.1hz < R2filter)
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
  allBOTH.filter$USEME = 1  # 1: include in the table and total VOC EF, 2: include in the table not total EF, 0: dont include at all
  # ------ These PM1 species don't correlate with CO ------
  ind1 = which(allBOTH.filter$variable =="Iodine_JIMENEZ")
  allBOTH.filter$USEME[ind1] = 0
  ind1 = which(allBOTH.filter$variable == "ClO4_JIMENEZ")
  allBOTH.filter$USEME[ind1] = 0
  ind1 = which(allBOTH.filter$variable =="Bromine_JIMENEZ")
  allBOTH.filter$USEME[ind1] = 0
  ind1 = which(allBOTH.filter$variable == "Seasalt_JIMENEZ")
  allBOTH.filter$USEME[ind1] = 0
  ind1 = which(allBOTH.filter$variable == "MSA_JIMENEZ")
  allBOTH.filter$USEME[ind1] = 0
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
  
  # ------- Removing these measurements entirely
  # No emission factors come out for these species
  ind = which(allBOTH.filter$formula == 'C2H4O3S')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$formula == 'MSA')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$formula == 'Na')
  allBOTH.filter$USEME[ind] = 0
  # Don't correlate with CO
  ind = which(allBOTH.filter$formula == 'BrO')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$formula == 'BrCl')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$formula == 'BrCN')
  allBOTH.filter$USEME[ind] = 0
  # ---------- These TOGA species don't correlate with CO
  ind = which(allBOTH.filter$variable  == 'iPropONO2_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CHBr3_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable =='CHBrCl2_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CH3CCl3_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CHCl3_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'HFC134a_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'HCFC142b_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'HCFC141b_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CHBr2Cl_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CH2ClI_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'LimoneneD3Carene_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'Propane_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'Propene_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'MBO_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'iButONO2and2ButONO2_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'C2H5OH_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CH2ClCH2Cl_ppt' ) 
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'HCFC22_ppt' ) 
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CH2Cl2_ppt' ) # TOGA doesn't correlate with CO but WAS does
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'C2Cl4_ppt' ) # 
  allBOTH.filter$USEME[ind] = 0
  # These iWAS species dont correlate with CO
  ind = which(allBOTH.filter$variable == 'CycHexane_NOAAiWAS_GILMAN')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable  == 'x3MePentane_NOAAiWAS_GILMAN')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable  == 'x224TriMePentane_NOAAiWAS_GILMAN')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable  == 'x22DiMeButane_NOAAiWAS_GILMAN')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable  == 'CHCl3_NOAAiWAS_GILMAN')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'C2Cl4_NOAAiWAS_GILMAN')
  allBOTH.filter$USEME[ind] = 0
  
  #These WAS species dont correlate with CO
  ind = which(allBOTH.filter$variable == 'C2Cl4_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CHBrCl2_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'C2HCl3_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CHCl3_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'Limonene_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'ClBenzene_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'H1211_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CFC11_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CFC12_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'x234TrimePentane_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CH2ClCH2Cl_WAS_BLAKE' ) 
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CCl4_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CH3CCl3_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CFC12_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'HFC134a_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'HFC152a_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'HCFC142b_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'HCFC141b_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'HFC365mfc_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'x2MePentane_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CFC114_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'x3MePentane_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'H1301_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'H2402_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CHBr2Cl_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CFC113_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable =='HCFC22_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'x23Dimebutane_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'CHBr3_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable =='x2PentONO2_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable =='x2ButONO2_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable =='x3PentONO2_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'iPropONO2_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'x3Me2ButONO2_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
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
  allBOTH.filter$USEME[ind] = 0
  
  ind = which(allBOTH.filter$variable =='ISOPN_WENNBERG')
  allBOTH.filter$USEME[ind] = 0
  
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
  allBOTH.filter$USEME[ind] =0
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
  ind = which(allBOTH.filter$names == 'Acetaldehyde' & allBOTH.filter$PI != 'WARNEKE') 
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$names == 'Acrolein' & allBOTH.filter$PI != 'WARNEKE') #APEL and WARNEKE agree, BLAKE is low
  allBOTH.filter$USEME[ind] = 0
  # Don't use WARNEKE or TOGA CH2O
  ind = which(allBOTH.filter$names == 'Formaldehyde' & allBOTH.filter$PI == 'ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$names == 'Formaldehyde' & allBOTH.filter$PI == 'WARNEKE')
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
  # --- Average GILMAN, BLAKE, APEL MEK
  allBOTH.filter = mergelines3(allBOTH.filter, 'MEK_ppt','MEK_WAS_BLAKE','MEK_NOAAiWAS_GILMAN')
  # --- Average WARNEKE, BLAKE, APEL Nitromethane
  allBOTH.filter = mergelines3(allBOTH.filter, 'CH3NO2_NOAAPTR_ppbv_WARNEKE','Nitromethane_WAS_BLAKE','Nitromethane_ppt')
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
  ind = which(allBOTH.filter$variable == 'Isoprene_WAS_BLAKE')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'Isoprene_ppt')
  allBOTH.filter$USEME[ind] = 0
  ind = which(allBOTH.filter$variable == 'Isoprene_NOAAiWAS_GILMAN')
  allBOTH.filter$USEME[ind] = 0
  allBOTH.filter = mergelines3(allBOTH.filter, 'Isoprene_WAS_BLAKE','Isoprene_ppt','Isoprene_NOAAiWAS_GILMAN')
  
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
    ind1B = which(allBOTH.filter$names[ind] != 'Ethylbenzene' & allBOTH.filter$names[ind] != 'o-Xylene' & allBOTH.filter$names[ind] != 'm,p-Xylene' &
                    allBOTH.filter$names[ind] != 'Methyl Ethyl Ketone' & allBOTH.filter$names[ind] != 'Methyl acetate')
    ind2B = which(allBOTH.filter$names[ind2] != 'Ethylbenzene' & allBOTH.filter$names[ind2] != 'o-Xylene' & allBOTH.filter$names[ind2] != 'm,p-Xylene')
    ind3B = which(allBOTH.filter$names[ind3] != 'Methyl Ethyl Ketone' & allBOTH.filter$names[ind3] != 'Methyl acetate')
    
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
    #allnames = allBOTH.filter$names[ind[ind1B]]
    #allPI = allBOTH.filter$PI[ind[ind1B]]
    #write.csv(NMVOC.table, 'NMVOCtable.csv')
    newline$FinalEF = sFVOC
    newline$FinalERtoCO = sRVOC
    newline$variable = 'NMVOC'
    newline$names = 'NMVOC'
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
    newline2$variable = 'Short-lived NMVOC'
    newline2$names = 'Short-lived NMVOC'
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
    newline3$variable = 'Long-lived NMVOC'
    newline3$names = 'Long-lived NMVOC'
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
# ----------------- AVERAGES ------------------------------
# Average by fire
#allBOTH.filter.avg.fire = aggregate(allBOTH.filter, by=list(allBOTH.filter$fire,allBOTH.filter$fuel, allBOTH.filter$variable), FUN='mean', na.rm=TRUE)
#allBOTH.filter.sd.fire = aggregate(allBOTH.filter, by=list(allBOTH.filter$fire, allBOTH.filter$fuel,allBOTH.filter$variable), FUN='sd', na.rm=TRUE)
# 
   # Average by  by individual fuel and species
#allBOTH.filter.avg.fire.fuel = aggregate(allBOTH.filter.avg.fire, by=list(allBOTH.filter.avg.fire$Group.2, allBOTH.filter.avg.fire$Group.3), FUN='mean', na.rm=TRUE)
#allBOTH.filter.sd.fire.fuel = aggregate(allBOTH.filter.avg.fire, by=list(allBOTH.filter.avg.fire$Group.2, allBOTH.filter.avg.fire$Group.3), FUN='sd', na.rm=TRUE)
# ----- Anything set to USEME = 0 here is already in the NMVOC EF, so this is double counted. Need to go back and fix at the end probably.
#       
#ind = which(allBOTH.filter$variable == 'PM1' & allBOTH.filter$names == 'Organic Carbon')
#allBOTH.filter <- allBOTH.filter[-c(ind), ]

allBOTH.filter.median = aggregate(allBOTH.filter, by=list(allBOTH.filter$variable), FUN='median', na.rm=TRUE)
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
q1 = 0.25; q2=0.75
allBOTH.filter.25.fuel = aggregate(allBOTH.filter$FinalEF, by=list(allBOTH.filter$fuel,allBOTH.filter$variable), 'quantile',probs=c(q1),na.rm=TRUE)
allBOTH.filter.75.fuel = aggregate(allBOTH.filter$FinalEF, by=list(allBOTH.filter$fuel,allBOTH.filter$variable), 'quantile',probs=c(q2),na.rm=TRUE)
allBOTH.filter.25.fuelER = aggregate(allBOTH.filter$FinalERtoCO, by=list(allBOTH.filter$fuel,allBOTH.filter$variable), 'quantile',probs=c(q1),na.rm=TRUE)
allBOTH.filter.75.fuelER = aggregate(allBOTH.filter$FinalERtoCO, by=list(allBOTH.filter$fuel,allBOTH.filter$variable), 'quantile',probs=c(q2),na.rm=TRUE)

allBOTH.filter.median.fuel$FinalEF_25 = allBOTH.filter.25.fuel$x
allBOTH.filter.median.fuel$FinalEF_75 = allBOTH.filter.75.fuel$x
allBOTH.filter.median.fuel$FinalERtoCO_25 = allBOTH.filter.25.fuelER$x
allBOTH.filter.median.fuel$FinalERtoCO_75 = allBOTH.filter.75.fuelER$x

# need to recover kind, formula, and names
for (i in 1:length(allBOTH.filter.median.fuel$kind)){
  ind = which(allBOTH.filter$variable == allBOTH.filter.median.fuel$Group.2[i])
  allBOTH.filter.median.fuel$kind[i] = allBOTH.filter$kind[ind[1]]
  allBOTH.filter.median.fuel$formula[i] = allBOTH.filter$formula[ind[1]]
  allBOTH.filter.median.fuel$names[i] = allBOTH.filter$names[ind[1]]
  allBOTH.filter.median.fuel$PI[i] = allBOTH.filter$PI[ind[1]]
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
cors.was.ag = c(); pval.was.ag=c(); counts.ag = c()
cors.was.pb = c(); pval.was.pb=c(); counts.pb = c()
cc = colnames(was.all)
par(mfrow=c(5,4))
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
allBOTH.filter.median.fuel$RtoCOoverall = NaN
for ( i in 1:length(allBOTH.filter.median.fuel$Group.2)){
  cc1 = colnames(allfires.1hz); cc2 = colnames(allfires.5hz)
  cc3 = colnames(toga.all); cc4 = colnames(was.all); cc5 = colnames(iwas.all)
  # just do corn to remove any fuel dependence
  ind1 = which(cc1 == allBOTH.filter.median.fuel$variable[i])
  ind1B = which(allfires.1hz$fuel == 'corn')
  ind2 = which(cc2 == allBOTH.filter.median.fuel$variable[i])
  ind2B = which(allfires.5hz$fuel == 'corn')
  ind3 = which(cc3 == allBOTH.filter.median.fuel$variable[i])
  ind3B = which(toga.all$fuel == 'corn')
  ind4 = which(cc4 == allBOTH.filter.median.fuel$variable[i])
  ind4B = which(was.all$fuel == 'corn')
  ind5 = which(cc5 == allBOTH.filter.median.fuel$variable[i])
  ind5B = which(iwas.all$fuel == 'corn')
  if (length(ind1) > 0){
    tmp = as.numeric(allfires.1hz[ind1B,ind1])
    iq = which(is.finite(tmp))
    if (length(iq) > 2){tt = cor.test(as.numeric(allfires.1hz$CO_DACOM_DISKIN[ind1B]), tmp)}
    allBOTH.filter.median.fuel$RtoCOoverall[i] = tt$estimate 
  }
  if (length(ind2) > 0){
    tmp = as.numeric(allfires.5hz[ind2B,ind2])
    iq = which(is.finite(tmp))
    if (length(iq) > 2){tt = cor.test(as.numeric(allfires.5hz$CO_DACOM_DISKIN[ind2B]), tmp)}
    allBOTH.filter.median.fuel$RtoCOoverall[i] = tt$estimate 
  }
  if (length(ind3) > 0){
    tmp = (toga.all[ind3B,ind3])
    tmp = unlist(tmp)
    tmp = as.numeric(tmp)
    iq = which(is.finite(tmp))
    if (length(iq) > 2){tt = cor.test(as.numeric(toga.all$CO_DACOM_DISKIN_BECKY[ind3B]), tmp)}
    allBOTH.filter.median.fuel$RtoCOoverall[i] = tt$estimate 
  }
  if (length(ind4) > 0){
    tmp = as.numeric(was.all[ind4B,ind4])
    iq = which(is.finite(tmp))
    if (length(iq) > 2){tt = cor.test(as.numeric(was.all$CO_DACOM_DISKIN_BLAKE[ind4B]), tmp)}
    allBOTH.filter.median.fuel$RtoCOoverall[i] = tt$estimate 
  }
  if (length(ind5) > 0){
    tmp = as.numeric(as.numeric(iwas.all[ind5B,ind5]))
    iq = which(is.finite(tmp))
    if (length(iq) > 2){tt = cor.test(as.numeric(iwas.all$CO_DACOM_DISKIN_GILMAN[ind5B]), tmp)}
    allBOTH.filter.median.fuel$RtoCOoverall[i] = tt$estimate 
  }
}

# --- MCE analysis Fig X ----- 
ind = which(allBOTH.filter.median.fuel$Group.1 == 'corn' & allBOTH.filter.median.fuel$USEME == 1 &
              is.finite(allBOTH.filter.median.fuel$corMCE_FINAL) & allBOTH.filter.median.fuel$variable != 'CO_DACOM_DISKIN' &
              allBOTH.filter.median.fuel$variable != 'CO2_ppb' & is.finite(allBOTH.filter.median.fuel$Category) & allBOTH.filter.median.fuel$Category != 5 &
              allBOTH.filter.median.fuel$formula != 'N/A')
uVOC = allBOTH.filter.median.fuel[ind,]
ind = which(allBOTH.filter.median.fuel$Group.1 == 'slash'  & allBOTH.filter.median.fuel$USEME == 1)
uVOC.slash = allBOTH.filter.median.fuel[ind,]

# does corn have a statistically different slope from slash? 
uVOC$DiffSlashCorn = NaN
for (i in 1:length(uVOC$Group.1)){
  ind = which(uVOC$names[i] == uVOC.slash$names & uVOC.slash$USEME == 1)
  if (length(ind) == 1){
    slope1 = c(uVOC$slopeMCE[i] - uVOC$slopeError[i], uVOC$slopeMCE[i] + uVOC$slopeError[i])
    slope2 = c(uVOC.slash$slopeMCE[ind] - uVOC.slash$slopeError[ind], uVOC.slash$slopeMCE[ind] + uVOC.slash$slopeError[ind])
    diff = (uVOC.slash$slopeMCE[ind] - uVOC$slopeMCE[i])*100/uVOC$slopeMCE[i]
    uVOC$DiffSlashCorn[i] =diff
  }
}
ind = which(uVOC$DiffSlashCorn > 50)
uVOC$DiffSlashCorn[ind] = median(uVOC$DiffSlashCorn[ind]) # 119
ind = which(uVOC$DiffSlashCorn > -25 & uVOC$DiffSlashCorn < 25)
uVOC$DiffSlashCorn[ind] = median(uVOC$DiffSlashCorn[ind]) # -7
ind = which(uVOC$DiffSlashCorn < -50)
uVOC$DiffSlashCorn[ind] = median(uVOC$DiffSlashCorn[ind]) # -67
ind = which(uVOC$DiffSlashCorn < -25 & uVOC$DiffSlashCorn > -50)
uVOC$DiffSlashCorn[ind] = median(uVOC$DiffSlashCorn[ind]) # -33

uVOC = uVOC[order(uVOC$corMCE_FINAL),]
ind = which(uVOC$variable != 'PM1')
uVOC = uVOC[ind,]
uVOC$NUM = seq(1,length(uVOC$Group.1))
# color this plot by whether corn is unique
ind = which(uVOC$kind == 'nitrogen' | uVOC$kind == 'nitrate' | uVOC$kind == 'alkyl nitrate')
uVOC$kind[ind] = 'NOy'
ind = which(uVOC$kind == 'CH2O')
uVOC$kind[ind] = 'oVOC'

ind = which(uVOC$kind == 'CH4')
uVOC$kind[ind] = 'alkane'

# Try different categories
uVOC$kind2 = uVOC$kind
uVOC$kind2[3] = "aldehyde"
uVOC$kind2[4] = "furans"
uVOC$kind2[6] = "aldehyde"
uVOC$kind2[7] = "aldehyde"
uVOC$kind2[8] = "ketone"
uVOC$kind2[9] = "aldehyde"
uVOC$kind2[10] = "aldehyde"
uVOC$kind2[11] = "ketone"
uVOC$kind2[12] = "aldehyde"
uVOC$kind2[13] = "phenolics"
uVOC$kind2[14] = "phenolics"
uVOC$kind2[15] = "carboxylic acid" # acetic acid is majority
uVOC$kind2[17] = "phenolics" #o-cresol?
uVOC$kind2[20] = "allene"
uVOC$kind2[24] = "amine"
uVOC$kind2[25] = "furans"
uVOC$kind2[26] = "nitrile"
uVOC$kind2[27] = "furans"
uVOC$kind2[28] = "phenolics" 
uVOC$kind2[29] = "alcohol" 
uVOC$kind2[30] = "furans"
uVOC$kind2[31] = "ketone"
uVOC$kind2[32] = "furans"
uVOC$kind2[33] = "ketone"
uVOC$kind2[36] = "ketone"
uVOC$kind2[37] = "nitrile"
uVOC$kind2[40] = "ketone" # metyl glyoxal, also an aldehyde
uVOC$kind2[47] = "allene"
uVOC$kind2[48] = "aldehyde"
uVOC$kind2[50] = "peroxide"
uVOC$kind2[59] = "peroxide"
uVOC$kind2[65] = "phenolics" 
uVOC$kind2[67] = "terpenes" 
uVOC$kind2[70] = "aldehyde"
uVOC$kind2[72] = "quinone"
uVOC$kind2[73] = "thiol"
uVOC$kind2[74] = "peroxide"
uVOC$kind2[78] = "furans"
uVOC$kind2[83] = "ester"
uVOC$kind2[84] = "nitrile"
uVOC$kind2[88] = "furans"
uVOC$kind2[97] = "phenolics" 
uVOC$kind2[99] = "ketone"
uVOC$kind2[102] = "aldehyde"
uVOC$kind2[103] = "nitrile"
uVOC$kind2[105] = "furans"
uVOC$kind2[107] = "ester"
uVOC$kind2[108] = "alcohol"
uVOC$kind2[111] = "nitrile"
uVOC$kind2[112] = "furans" # furans?

ind = which(uVOC$names != ' Furan and fragments')
uVOC = uVOC[ind,]
uVOC$corMCE_FINAL = round(uVOC$corMCE_FINAL, digits = 2)
uVOC = uVOC[order(uVOC$corMCE_FINAL),]
uVOC$NUM = seq(1,length(uVOC$Group.1))

library(pals)
cornEFMCE = ggplot(uVOC , aes(y=corMCE_FINAL, x=NUM, fill=kind2)) + #ylim(0,-1)+
  geom_bar(position="dodge", stat="identity") + theme_classic()+ 
  geom_text(aes(x =NUM,y = 0,label = names), size=5,
            vjust = 0,  hjust = 0,angle = 90, nudge_y = 0.01, nudge_x = 0.2)+
  ylab('R EF MCE') + xlab(c(""))+ylab("Correlation with MCE") +
  scale_fill_manual(values=as.vector(polychrome(26)))+ theme(legend.position = "top")+labs(fill="")

# ----------------------------------------------------------------------------
# --------- here run TestTable.R -------
# ----------------------------------------------------------------------------
# ---------- Table 1 ------------
ind = which(allBOTH.filter$fuel == 'corn' & allBOTH.filter$variable == 'CO_DACOM_DISKIN')
cornplumes = length(ind)
cornfires = length(unique(allBOTH.filter$fire[ind]))
cornmce = quantile(allBOTH.filter$mce.5hz[ind], na.rm=TRUE)
ind = which(allBOTH.filter$fuel == 'soybean' & allBOTH.filter$variable == 'CO_DACOM_DISKIN')
soyplumes = length(ind)
soyfires = length(unique(allBOTH.filter$fire[ind]))
soymce = quantile(allBOTH.filter$mce.5hz[ind], na.rm=TRUE)
ind = which(allBOTH.filter$fuel == 'rice' & allBOTH.filter$variable == 'CO_DACOM_DISKIN')
riceplumes = length(ind)
ricefires = length(unique(allBOTH.filter$fire[ind]))
ricemce = quantile(allBOTH.filter$mce.5hz[ind], na.rm=TRUE)
ind = which(allBOTH.filter$fuel == 'winter wheat' & allBOTH.filter$variable == 'CO_DACOM_DISKIN')
wheatplumes = length(ind)
wheatfires = length(unique(allBOTH.filter$fire[ind]))
wheatmce = quantile(allBOTH.filter$mce.5hz[ind], na.rm=TRUE)
ind = which(allBOTH.filter$fuel == 'grass' & allBOTH.filter$variable == 'CO_DACOM_DISKIN')
grassplumes = length(ind)
grassfires = length(unique(allBOTH.filter$fire[ind]))
grassmce = quantile(allBOTH.filter$mce.5hz[ind], na.rm=TRUE)
ind = which(allBOTH.filter$fuel == 'slash' & allBOTH.filter$variable == 'CO_DACOM_DISKIN')
slashplumes = length(ind)
slashfires = length(unique(allBOTH.filter$fire[ind]))
slashmce = quantile(allBOTH.filter$mce.5hz[ind], na.rm=TRUE)
ind = which(allBOTH.filter$fuel == 'pile' & allBOTH.filter$variable == 'CO_DACOM_DISKIN')
pileplumes = length(ind)
pilefires = length(unique(allBOTH.filter$fire[ind]))
pilemce = quantile(allBOTH.filter$mce.5hz[ind], na.rm=TRUE)

ind = which(allBOTH.filter$fuel == 'shrub' & allBOTH.filter$variable == 'CO_DACOM_DISKIN')
shrubplumes = length(ind)
shrubfires = length(unique(allBOTH.filter$fire[ind]))
shrubmce = quantile(allBOTH.filter$mce.5hz[ind], na.rm=TRUE)
fires = c(cornfires,ricefires,soyfires, wheatfires, grassfires, slashfires, pilefires, shrubfires  )
plumes = c(cornplumes,riceplumes, soyplumes, wheatplumes, grassplumes, slashplumes, pileplumes, shrubplumes  )
mce = round(c(cornmce[3],ricemce[3],soymce[3],wheatmce[3],grassmce[3],slashmce[3],pilemce[3],shrubmce[3]),2)
mce25 = round(c(cornmce[2],ricemce[2],soymce[2],wheatmce[2],grassmce[2],slashmce[2],pilemce[2],shrubmce[2]),2)
mce75 = round(c(cornmce[4],ricemce[4],soymce[4],wheatmce[4],grassmce[4],slashmce[4],pilemce[4],shrubmce[4]),2)
table1 = cbind(fires, plumes,mce, mce25,mce75)
rownames(table1) = c("corn", "rice","soy","wheat","grass","slash","pile","shrub")
table1 = as.data.frame(table1)
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
    runplots =0
    if (runplots == 1){
        
      ind2 = which(allBOTH.filter$variable == 'HCN_NOAAPTR_ppbv_WARNEKE' & allBOTH.filter$fuel == 'corn')
      ind2 = which(allBOTH.filter$variable == 'CH3CN_NOAAPTR_ppbv_WARNEKE' & allBOTH.filter$fuel == 'corn')
      ind2 = which(allBOTH.filter$variable == 'NO2_ACES_WOMACK' )#& allBOTH.filter$fuel == 'corn')
      ind = which(allBOTH.filter$variable == 'HNO2_ACES_WOMACK' )#& allBOTH.filter$fuel == 'corn')
      
      ggplot(allBOTH.filter[ind,])+geom_point(aes(x=allBOTH.filter$FinalEF[ind], y=allBOTH.filter$FinalEF[ind2], col=fuel, shape=fuel),stroke=3, size=4)+
        ylab("NO2 (ACES)")+xlab("HONO (ACES)")+theme_classic()+
        #scale_shape_manual(values=c(19,15,17,18,7,8,4))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+labs(col="",shape="")+
        theme(legend.background=element_blank())+theme(legend.position = "none")+ theme(text = element_text(size = 20)) 
      # ---- Figure1: Plot CH4 EF vs. MCE for separated fuel types -----
      CH4vsMCE = plotSpeciesMCE(allBOTH.filter,'CH4_DACOM_DISKIN','Methane','CH4','CH4')
      # Figure CH4
      ggsave('CH4vsMCE.ps',CH4vsMCE,width = 7*1.25*1.25, height=7*1.25)
      # ------- Figure SO2 -------------
      SO2vsMCE = plotSpeciesMCE(allBOTH.filter,'SO2_LIF_ROLLINS','SO2','SO2','SO2')
      ggsave('SO2vsMCE.ps',SO2vsMCE,width = 7*1.25*1.25, height=7*1.25)
      
      # ------- Figure aerosols ----------------
      PM1vsMCE = plotSpeciesMCE(allBOTH.filter,'PM1','PM1','SO2','SO2')
      ECvsMCE = plotSpeciesMCE(allBOTH.filter,'BC_SCHWARZ','BlackCarbon','Black Carbon','BC')
      xiaoxi$OC = xiaoxi$OA/2
      OCvsMCE = plotSpeciesMCE(allBOTH.filter,'OC_JIMENEZ','OrganicCarbon','Organic Carbon','OC')
      LGvsMCE = plotSpeciesMCE(allBOTH.filter,'C6H10O5_JIMENEZ','Levoglucosan','Levoglucosan','Levoglucosan')
      
      ClvsMCE = plotSpeciesMCE(allBOTH.filter,'NR_Chloride_JIMENEZ','Cl','Chloride','Chl')
      SO4vsMCE = plotSpeciesMCE(allBOTH.filter,'Sulfate_JIMENEZ','Sulfate','Sulfate','SO4')
      NITvsMCE = plotSpeciesMCE(allBOTH.filter,'Nitrate_JIMENEZ','Nitrate','Nitrate','NO3')
      NH4vsMCE = plotSpeciesMCE(allBOTH.filter,'Ammonium_JIMENEZ','Ammonium','Ammonium','NH4')
      KvsMCE = plotSpeciesMCE(allBOTH.filter,'Potassium_JIMENEZ','K','K','K')
      NCATsMCE = plotSpeciesMCE(allBOTH.filter,'C6H5NO4_JIMENEZ','4-Nitrocatechol','4-Nitrocatechol','4-Nitrocatechol')
      allPM =ggarrange(OCvsMCE,LGvsMCE,ECvsMCE,ClvsMCE,NH4vsMCE,KvsMCE,
                common.legend = TRUE,nrow=2,ncol=3,
                labels = c("a)","b)","c)","d)","e)","f)"),
                hjust = c(-5,-5,-6,-6))
      ggsave('allPM.ps',allPM,width = 7*1.25*1.25, height=7*1.25)
      
      ind = which(allBOTH.filter$variable == 'OC_JIMENEZ' & allBOTH.filter$fuel != 'forest' )#& allBOTH.filter$fuel == 'corn')
      OC = allBOTH.filter$FinalEF[ind]
      ind = which(allBOTH.filter$variable == 'BC_SCHWARZ' & allBOTH.filter$fuel != 'forest' )#& allBOTH.filter$fuel == 'corn')
      EC = allBOTH.filter$FinalEF[ind]
      ggplot(allBOTH.filter[ind,])+geom_point(aes(x=MCE, y=EC/OC, col=fuel))+theme_classic()
        geom_point(data=xiaoxi,aes(x=xiaoxi$MCE, xiaoxi$BC/xiaoxi$OC), col='green')
      
        # ------- Figure  total VOCs ----------------
        ind = which(allBOTH.filter$variable == 'Short-lived NMVOC' & allBOTH.filter$fuel == 'corn')
      cor.test(allBOTH.filter$MCE[ind], allBOTH.filter$FinalEF[ind])
      ind = which(allBOTH.filter$variable == 'Long-lived NMVOC' & allBOTH.filter$fuel == 'corn')
      cor.test(allBOTH.filter$MCE[ind], allBOTH.filter$FinalEF[ind])
      # Figure ShortVOC
      ShortVOCvsMCE = plotSpeciesMCE(allBOTH.filter,'Short-lived NMVOC','Short-lived VOC','NMVOC','NMVOC')
      LongVOCvsMCE = plotSpeciesMCE(allBOTH.filter,'Long-lived NMVOC','Long-lived VOC','NMVOC','NMVOC')
      # ------- Figure individual VOCs ----------------
      #AcetaldehydevsMCE = plotSpeciesMCE(allBOTH.filter,'CH3CHO_NOAAPTR_ppbv_WARNEKE','Acetaldehyde','SO2','SO2')
      FormaldehydevsMCE    = plotSpeciesMCE(allBOTH.filter,'Formaldehyde_FRIED_HANISCO','Formaldehyde','Formaldehyde','Formaldehyde')
      FormaldehydevsMCE2    = plotSpeciesMCE(allBOTH.filter,'CH2O_ISAF_HANISCO','Formaldehyde','Formaldehyde','Formaldehyde')
      FormaldehydevsMCE3    = plotSpeciesMCE(allBOTH.filter,'CH2O_CAMS_pptv_FRIED','Formaldehyde','Formaldehyde','Formaldehyde')
      C2H4vsMCE    = plotSpeciesMCE(allBOTH.filter,'Ethene_GILMAN_BLAKE','C2H4','Ethene','Ethene')
      AcetvsMCE    = plotSpeciesMCE(allBOTH.filter,'Acetone_NOAAiWAS_GILMAN','Acetone','Acetone','Acetone')
      AcrvsMCE    = plotSpeciesMCE(allBOTH.filter,'Acrolein_NOAAPTR_ppbv_WARNEKE','Acrolein','Acrolein','Acrolein')
      EthanevsMCE    = plotSpeciesMCE(allBOTH.filter,'C2H6_CAMS_pptv_FRIED','Ethane','Ethane','Ethane')
      TOLUvsMCE    = plotSpeciesMCE(allBOTH.filter,'Toluene_NOAAPTR_ppbv_WARNEKE','Toluene','Toluene','Toluene')
      #Depolymerization in lignin (300–500C)produces guaiacols, (iso)eugenol, and syringol.
      SYRvsMCE    = plotSpeciesMCE(allBOTH.filter,'Syringol_NOAAPTR_ppbv','Syringol','Syringol','Syringol')
      GuaiacolvsMCE    = plotSpeciesMCE(allBOTH.filter,'Guaiacol_NOAAPTR_ppbv_WARNEKE','Guaiacol','Guaiacol','Guaiacol')
      #Furans and furfurals are dominantly formed from cellulose and hemicellulose (300–400oC).
      FuranvsMCE    = plotSpeciesMCE(allBOTH.filter,'Furan_ppt_BLAKE_','Furan','Furan','Furan')
      Furan2MevsMCE    = plotSpeciesMCE(allBOTH.filter,'2-Methylfuran_BLAKE_ppt','2-Methylfuran','2-Methylfuran','2-Methylfuran')
      Furan3MevsMCE    = plotSpeciesMCE(allBOTH.filter,'3-Methylfuran_BLAKE_ppt','3-Methylfuran','3-Methylfuran','3-Methylfuran')
     # FurfuralvsMCE    = plotSpeciesMCE(allBOTH.filter,'Furfural_ppt','Furfural','Furfural','Furfural')
      FurfuralvsMCE    = plotSpeciesMCE(allBOTH.filter,'Furfural_NOAAPTR_ppbv_WARNEKE','Furfural','Furfural','Furfural')
      #Higher temperatures allow reaction #of functional groups and covalent bonds in polymers and
      #monomers. The resulting fragmentation emits various VOCs:for example, hydroxyacetone, acetaldehyde, and acetic acid
      #from depolymerization of cellulose and/or hemicellulose
      AcetaldehydevsMCE    = plotSpeciesMCE(allBOTH.filter,'CH3CHO_NOAAPTR_ppbv_WARNEKE','Acetaldehyde','Acetaldehyde','Acetaldehyde')
      AceticvsMCE    = plotSpeciesMCE(allBOTH.filter,'GlycolaldehydeCH3COOH_NOAAPTR_ppbv_WARNEKE','GlycoaldehydeCH3COOH','GlycoaldehydeCH3COOH','GlycoaldehydeCH3COOH')
      HACvsMCE    = plotSpeciesMCE(allBOTH.filter,'C3H6O2_NOAAPTR_ppbv_WARNEKE','Hydroxyacetone','Hydroxyacetone','Hydroxyacetone')
      
      #Higher temperature pyrolysis breaks progressively stronger bonds in char (> 500C).
      #This aromatization process gives off aromatic compounds with short substituents (e.g., phenol), 
      #nonsubstituted aromatics (e.g.,benzene), and polycyclic aromatic hydrocarbons (PAHs) such as naphthalene).
      BenzenevsMCE    = plotSpeciesMCE(allBOTH.filter,'Benzene_NOAAPTR_ppbv_WARNEKE','Benzene','Benzene','Benzene')
      StyrenevsMCE    = plotSpeciesMCE(allBOTH.filter,'Styrene_NOAAPTR_ppbv_WARNEKE','Styrene','Styrene','Styrene')
      C8vsMCE    = plotSpeciesMCE(allBOTH.filter,'C8Aromatics_NOAAPTR_ppbv_WARNEKE','Xylenes','Xylenes','Xylenes')
      C9vsMCE    = plotSpeciesMCE(allBOTH.filter,'C9Aromatics_NOAAPTR_ppbv_WARNEKE','C9','C9','C9')
      NaphthalenevsMCE    = plotSpeciesMCE(allBOTH.filter,'Naphthalene_NOAAPTR_ppbv_WARNEKE','Naphthalene','Naphthalene','Naphthalene')
      PhenolvsMCE    = plotSpeciesMCE(allBOTH.filter,'PHENOL_WENNBERG','Phenol2','Phenol2','Phenol2')
      PhenolvsMCE2    = plotSpeciesMCE(allBOTH.filter,'Phenol_NOAAPTR_ppbv_WARNEKE','Phenol','Phenol','Phenol')
      CycHexvsMCE= plotSpeciesMCE(allBOTH.filter,'CycHexane_WAS_BLAKE','Cyclohexane','Cyclohexane','Cyclohexane')
    # Toluene outlier
      aroms = ggarrange(BenzenevsMCE,StyrenevsMCE, CycHexvsMCE, TOLUvsMCE+ylim(limits=c(0,1.1)), C8vsMCE, C9vsMCE, common.legend = TRUE)
      # Ethene_GILMAN_
     # PhenolvsMCE2    = plotSpeciesMCE(allBOTH.filter,'Phenol_NOAAPTR_ppbv_WARNEKE','Phenol','SO2','SO2')
      #?
      #MonovsMCE    = plotSpeciesMCE(allBOTH.filter,'beta-Pinene/Myrcene_BLAKE_ppt','beta-Pinene/Myrcene','SO2','SO2')
      MonovsMCE    = plotSpeciesMCE(allBOTH.filter,'Monoterpenes_NOAAPTR_ppbv_WARNEKE','Terpenes','Terpenes','Monoterpenes')
      CatvsMCE = plotSpeciesMCE(allBOTH.filter,'Catecholx5MeFurfural_NOAAPTR_ppbv_WARNEKE','Catechol/5MeFurfural','Catechol','Catechol')
    
      EOHvsMCE = plotSpeciesMCE(allBOTH.filter,'C2H5OH_NOAAPTR_ppbv_WARNEKE','Ethanol','Ethanol','Ethanol')
      CH3OHvsMCE = plotSpeciesMCE(allBOTH.filter,'CH3OH_NOAAPTR_ppbv_WARNEKE','Methanol','Methanol','Methanol')
      GLYOXALvsMCE = plotSpeciesMCE (allBOTH.filter,'CHOCHO_ACES_WOMACK','Glyoxal','Glyoxal','Glyoxal')
      ISOPvsMCE = plotSpeciesMCE (allBOTH.filter,'Isoprene_NOAAPTR_ppbv_WARNEKE','Isoprene','Isoprene','Isoprenehydroperoxyaldehydes')
      HCOOHvsMCE = plotSpeciesMCE (allBOTH.filter,'Formic acid_VERES_WARNEKE','Formic Acid','Formic acid','HCOOH')
     
      #ind = which(allBOTH.filter$kind == 'aromatic')
      #aromsplot = unique(allBOTH.filter$variable[ind])

      #for (i in 1:length(aromsplot)){
      #  plotSpeciesMCE(allBOTH.filter,aromsplot[i],aromsplot[i],aromsplot[i],aromsplot[i])
      #}
       ggarrange(ShortVOCvsMCE,LongVOCvsMCE,
                common.legend = TRUE,nrow=1,ncol=2,
                labels = c("a)","b)"),
                hjust = c(-5,-5))
      allVOCplot = ggarrange(FormaldehydevsMCE,AcetaldehydevsMCE,CatvsMCE, 
                GuaiacolvsMCE,FuranvsMCE,
                HACvsMCE, PhenolvsMCE,NaphthalenevsMCE, BenzenevsMCE,
                MonovsMCE,EOHvsMCE, CH3OHvsMCE,
                GLYOXALvsMCE, ISOPvsMCE, HCOOHvsMCE,
                common.legend = TRUE,nrow=5,ncol=3,
                labels = c("a)","b)","c)","d)","e)","f)","g)","h)","i)","j)","k)","l)","m)","n)","o)"),
                hjust = c(-8))
      ggsave("allVOCplot.pdf", allVOCplot, width = 7*1.25*2, height=7*1.25*2)
      # Figure NH3
      NOvsMCE = plotSpeciesMCE(allBOTH.filter,'Nitrogen oxide_ROLLINS_RYERSON','NO','NO','NO')
      NO2vsMCE = plotSpeciesMCE(allBOTH.filter,'Nitrogen dioxide_WOMACK_RYERSON','NO2','NO2','NO2')
      HNO2vsMCE = plotSpeciesMCE(allBOTH.filter,'Nitrous acid_WOMACK_VERES','HNO2','HNO2','HNO2')
      CH3CNvsMCE = plotSpeciesMCE(allBOTH.filter,'CH3CN_NOAAPTR_ppbv_WARNEKE','Acetonitrile','Acetonitrile','Acetonitrile')
      HCNvsMCE = plotSpeciesMCE(allBOTH.filter,'Hydrogen cyanide_VERES_WENNBERG_','HCN','CH3CN','CH3CN')
      NH3vsMCE = plotSpeciesMCE(allBOTH.filter,'NH3_WISTHALER','Ammonia','NH3','NH3')
      HNCOvsMCE = plotSpeciesMCE(allBOTH.filter,'Isocyanic acid_VERES_WARNEKE','HNCO','HNCO','HNCO')
      HCNvsMCE = plotSpeciesMCE(allBOTH.filter,'Hydrogen cyanide_VERES_WENNBERG_','HydrogenCyanide','Hydrogen cyanide','HCN')
      CH3NO3vsMCE = plotSpeciesMCE(allBOTH.filter,'Methyl nitrate_BLAKE_ppt','MethylNitrate','MethylNitrate','MethylNitrate')
      ggarrange(NH3vsMCE,NOvsMCE,NO2vsMCE, HNO2vsMCE,HCNvsMCE,CH3CNvsMCE, common.legend = TRUE,nrow=2,
                ncol=3,labels = c("a)","b)","c)","d)","e)","f)"),
                hjust = c(-5,-5,-6,-6,-5,-8))
      # Figure CH3Cl
      CH3ClvsMCE = plotSpeciesMCE(allBOTH.filter,'CH3Cl_WAS_BLAKE','Chloromethane','CH3Cl','Chloromethane')
      CH3IvsMCE = plotSpeciesMCE(allBOTH.filter,'Methyl Iodide_BLAKE_ppt','MethylIodide','CH3I','Methyl Iodide')
      CH3BrvsMCE = plotSpeciesMCE(allBOTH.filter,'Bromomethane_ppt_BLAKE','Bromomethane','CH3Br','Bromomethane')
      CH2Cl2vsMCE = plotSpeciesMCE(allBOTH.filter,'Dichloromethane_BLAKE_ppt','Dichloromethane','Dichloromethane','Dichloromethane')
      
      ind = which(was.all$fuel != 'forest' & was.all$fuel != 'house' & was.all$fuel != '?')
      ggplot(was.all[ind,])+geom_point(aes(x=CO_DACOM_DISKIN_BLAKE, y=CH2Cl2_WAS_BLAKE, col=fuel), size=4)+theme_classic()+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)
      
      ind = which(toga.all$fuel != 'forest' & toga.all$fuel != 'house' & toga.all$fuel != '?')
      ggplot(toga.all[ind,])+geom_point(aes(x=as.numeric(CO_DACOM_DISKIN_BECKY), y=CH2Cl2_ppt, col=fuel), size=4)+theme_classic()+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)
      
      #CH3ClvsCO = plotSpeciesCO(was.all,'CH3Cl_WAS_BLAKE','CO_DACOM_DISKIN_BLAKE','Chloromethane','Chloromethane','Chloromethane')
      # CH3COOClvsCO = plotSpeciesCO(allfires.1hz,'CH3COOCl_NOAACIMS_VERES','CO_DACOM_DISKIN','Chloroacetic acid','Chloromethane','Chloromethane')
      #C2H5ClvsCO = plotSpeciesCO(was.all,'C2H5Cl_WAS_BLAKE','CO_DACOM_DISKIN_BLAKE','Chloromethane','Chloromethane','Chloromethane')
      #C6H5ClvsCO = plotSpeciesCO(toga.all,'ClBenzene_ppt','CO_DACOM_DISKIN','Chlorobenzene','Chloromethane','Chloromethane')
      ###CH2Br2vsCO = plotSpeciesCO(was.all,'CH2Br2_WAS_BLAKE','CO_DACOM_DISKIN_BLAKE','Dibromomethane','Chloromethane','Chloromethane')
      #CH3BrvsCO = plotSpeciesCO(was.all,'CH3Br_WAS_BLAKE','CO_DACOM_DISKIN_BLAKE','Bromomethane','Chloromethane','Chloromethane')
      #CH3BrvsCO = plotSpeciesCO(toga.all,'CH3Br_ppt','CO_DACOM_DISKIN','Bromomethane','Chloromethane','Chloromethane')
      # Figure halogen
      allHal = ggarrange(CH3ClvsMCE,CH3BrvsMCE,CH3IvsMCE, common.legend = TRUE,ncol=3,labels = c("a)","b)","c)"),
                hjust = c(-5,-8,-8))
      ggsave(allHal, file='allHal.ps',width = 7*1.25*2, height=7*1.25*2/3)
      # all variables
      
      # Andreae comparison
      
      vars = unique(allBOTH.filter$variable)
      for (i in 1:length(vars)){
        vv = vars[i]
        # need at least 25% of the data?
        ind = which(allBOTH.filter$variable == vv & is.finite(allBOTH.filter$FinalEF) & allBOTH.filter$USEME == 1)
        ind3 = which(allBOTH.filter$variable == vv )
        
        if (length(ind)/length(ind3) > 0.25){
          ff=plotSpeciesMCE(allBOTH.filter, vv,vv,vv)
          print(ff)
          ggsave(filename=paste(vv, '.ps', sep=''),ff, path='/Users/ktravis1/OneDrive - NASA/FIREX/MCEplots/')

        }
      }
      
      ind = which(allBOTH.filter$formula == 'HCN' & allBOTH.filter$USEME == 1 & allBOTH.filter$fuel != 'forest')
      ind2 = which(allBOTH.filter$formula == 'HCN' & allBOTH.filter$USEME == 1 & allBOTH.filter$fuel == 'slash')
      tt = lmodel2(allBOTH.filter$FinalEF[ind2]~ allBOTH.filter$MCE[ind2])
      ggplot(allBOTH.filter[ind,]) + geom_point(aes(x=MCE,y=FinalEF, col=fuel), size=2)+ theme_classic()+
        geom_abline(mapping=aes(intercept=tt$coefficients[1], slope=tt$coefficients[2]))
      CH4vsMCE = CH4vsMCE + 
       # geom_point(data=xiaoxi.avg, aes(x=mean[1], y=mean[indC]), pch=0,col='green',stroke=2,size=3)+
        #geom_errorbar(data=xiaoxi.avg,aes(xmin=mean[1]-sd[1],xmax=mean[1]+sd[1], y=mean[indC]),col='green', width=1.3*2, position=position_dodge(0.05))+
        #geom_errorbar(data=xiaoxi.avg,aes(x=mean[1], ymin=mean[indC] - sd[indC], ymax=mean[indC]+ sd[indC]),col='green', width=0.0025, position=position_dodge(0.05))+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_errorbar(data=andreae,aes(xmin=as.numeric(andreae$average[1])-as.numeric(andreae$std.dev.[1]),xmax=as.numeric(andreae$average[1])+as.numeric(andreae$std.dev.[1]), y=as.numeric(andreae$average[indB])),col='purple', width=1.3*2, position=position_dodge(0.05))+
        geom_errorbar(data=andreae,aes(x=as.numeric(andreae$average[1]), ymin=as.numeric(andreae$average[indB])-as.numeric(andreae$std.dev.[indB]),ymax=as.numeric(andreae$average[indB])+as.numeric(andreae$std.dev.[indB])),col='purple', width=0.0025, position=position_dodge(0.05))+
      
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)+
        geom_errorbar(data=akagi,aes(xmin=akagi$PastureEF[1]-as.numeric(akagi$PastureSD[1]),
                                     xmax=akagi$PastureEF[1]+as.numeric(akagi$PastureSD[1]), 
                                     y=akagi$PastureEF[indA]),col='pink', width=1.3*2, position=position_dodge(0.05))+
        geom_errorbar(data=akagi,aes(x=akagi$PastureEF[1],
                                     ymin=akagi$PastureEF[indA] - as.numeric(akagi$PastureSD[indA]),
                                     ymax=akagi$PastureEF[indA] + as.numeric(akagi$PastureSD[indA])),
                      col='pink', width=0.0025, position=position_dodge(0.05))+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_errorbar(data=akagi,aes(xmin=akagi$CropEF[1]-as.numeric(akagi$CropSD[1]),
                                       xmax=akagi$CropEF[1]+as.numeric(akagi$CropSD[1]), 
                                       y=akagi$CropEF[indA]),col='pink', width=1.3*2, position=position_dodge(0.05))+
        geom_errorbar(data=akagi,aes(x=akagi$CropEF[1],
                                       ymin=akagi$CropEF[indA] - as.numeric(akagi$CropSD[indA]),
                                       ymax=akagi$CropEF[indA] + as.numeric(akagi$CropSD[indA])),
                        col='pink', width=0.0025, position=position_dodge(0.05))+
        
        geom_point(data=mccarty, aes(x=mccarty$mce, y=mccarty$ch4, col='grey', size=5, shape=fuel, stroke=3))+
        geom_errorbar(data=mccarty,aes(x=mccarty$mce,
                                     ymin=mccarty$ch4 - mccarty$ch4SD ,
                                     ymax=mccarty$ch4 + mccarty$ch4SD ),
                      col='grey', width=0.0025, position=position_dodge(0.05))
      # add in my averages
      indK = which(allBOTH.filter.ag.avg$Group.1 == 'CH4_DACOM_DISKIN')
        CH4vsMCE+geom_point(data=allBOTH.filter.ag.avg[indK,],aes(x=mce.5hz,y=FinalEF), size=3, col='purple')+
          geom_errorbar(data=allBOTH.filter.ag.avg[indK,],aes(x=allBOTH.filter.ag.avg$mce.5hz[indK],
                                       ymin=allBOTH.filter.ag.avg$FinalEF[indK] - allBOTH.filter.ag.avg$FinalEF.ag.sd[indK],
                                       ymax=allBOTH.filter.ag.avg$FinalEF[indK] + allBOTH.filter.ag.avg$FinalEF.ag.sd[indK]),
                        col='purple', width=0.0025, position=position_dodge(0.05))+
          geom_errorbar(data=allBOTH.filter.ag.avg[indK,],aes(xmin=allBOTH.filter.ag.avg$mce.5hz[indK]-allBOTH.filter.ag.avg$mce.ag.sd[indK],
                                                              xmax=allBOTH.filter.ag.avg$mce.5hz[indK]+allBOTH.filter.ag.avg$mce.ag.sd[indK],
                                                              y=allBOTH.filter.ag.avg$FinalEF[indK]) ,
                        col='purple',  width=1.3*2, position=position_dodge(0.05))+
          geom_point(data=allBOTH.filter.sc.avg[indK,],aes(x=mce.5hz,y=FinalEF), size=3, col='pink')+   
          geom_errorbar(data=allBOTH.filter.sc.avg[indK,],aes(x=allBOTH.filter.sc.avg$mce.5hz[indK],
                                                              ymin=allBOTH.filter.sc.avg$FinalEF[indK] - allBOTH.filter.sc.avg$FinalEF.sc.sd[indK],
                                                              ymax=allBOTH.filter.sc.avg$FinalEF[indK] + allBOTH.filter.sc.avg$FinalEF.sc.sd[indK]),
                        col='pink', width=0.0025, position=position_dodge(0.05))+
          geom_errorbar(data=allBOTH.filter.sc.avg[indK,],aes(xmin=allBOTH.filter.sc.avg$mce.5hz[indK]-allBOTH.filter.sc.avg$mce.sc.sd[indK],
                                                              xmax=allBOTH.filter.sc.avg$mce.5hz[indK]+allBOTH.filter.sc.avg$mce.sc.sd[indK],
                                                              y=allBOTH.filter.sc.avg$FinalEF[indK]) ,
                        col='pink',  width=1.3*2, position=position_dodge(0.05))+
          theme_classic()
      
      allplot = as.data.frame(cbind(mce=c(akagi$CropEF[1],akagi$PastureEF[1],as.numeric(andreae$average[1]), xiaoxi.avg$mean[1], mccarty$mce, allBOTH.filter.CO.avg$mce.5hz),
                      co = (c(akagi$CropEF[indA],akagi$PastureEF[indA],as.numeric(andreae$average[indB]), xiaoxi.avg$mean[indC], mccarty$co, allBOTH.filter.CO.avg$FinalEF)),
                      coSD = (c(akagi$CropSD[indA],akagi$PastureSD[indA],as.numeric(andreae$std.dev.[indB]), xiaoxi.avg$sd[indC], mccarty$coSD, allBOTH.filter.CO.avg$FinalEF.sd)),
                      fuel = c("agricultural residue","pasture maintenance","agricultural residue","rice",mccarty$fuel,allBOTH.filter.CO.avg$fuel),
                      study = c("Akagi","Akagi","Andreae","Liu","McCarty","McCarty","McCarty","McCarty","McCarty","McCarty","McCarty","This Study","This Study","This Study","This Study","This Study","This Study","This Study","This Study")))
      
      # Basic box plot
      ind = which(allBOTH.filter.ag$variable == 'CO_DACOM_DISKIN')
      ind2 = which(allBOTH.filter.sc$variable == 'CO_DACOM_DISKIN')
      p <- ggplot() + 
        geom_boxplot(data=allBOTH.filter.ag[ind,], aes(x=1, y=as.numeric(FinalEF)))+
        geom_boxplot(data=allBOTH.filter.sc[ind2,], aes(x=2, y=as.numeric(FinalEF)))+
        geom_boxplot(data=xiaoxi,aes(x=3,y=xiaoxi$CO))+ylab('CO, g/kg')+
        geom_boxplot(data=akagi[indB,],aes(x=4,y=CropEF))
      
      ggplot(allplot, aes(x=as.numeric(mce), y=as.numeric(co))) + 
        geom_bar(position="dodge", stat="identity") + theme_classic()+ 
       geom_text(aes(x =as.numeric(mce),y = 0,label =paste(fuel,study)), size=5,
                  vjust = -1,  hjust = 0,angle = 90, nudge_y = 0.01, nudge_x = 0.002)+
        ylab('CO EF, g/kg')
      ggplot(allplot) + geom_point(aes(x=as.numeric(mce), y=as.numeric(co),col=study, size=fuel), size=5) + theme_classic()
      
      ggsave(COvsMCE, file='COvsMCE.ps')
      sz=25
      
      # ------ Figure 2 - NO fire -------------
      
      ind =which(allBOTH.filter$variable == 'NO_LIF_ROLLINS' )
      ind2 =which(allBOTH.filter$variable == 'NO_RYERSON' )
      plot(allBOTH.filter$FinalEF[ind], allBOTH.filter$FinalEF[ind2], xlab='Rollins NO', ylab='Ryerson NO')
      abline(lm(allBOTH.filter$FinalEF[ind2]~allBOTH.filter$FinalEF[ind]))
      abline(0,1,lty=2)
      ind =which(allBOTH.filter$variable == 'NO2_ACES_WOMACK' )
      ind2 =which(allBOTH.filter$variable == 'NO2_RYERSON' )
      plot(allBOTH.filter$FinalEF[ind], allBOTH.filter$FinalEF[ind2], xlab='Womack NO2', ylab='Ryerson NO2')
      abline(lm(allBOTH.filter$FinalEF[ind2]~allBOTH.filter$FinalEF[ind]))
      abline(0,1,lty=2)
      
      
      abline(0,1, lty=2)
      Makeplots(allBOTH.filter,'NO_RYERSON','NO','NO','NO')
      ind = which(allBOTH.filter$variable == 'NO_LIF_ROLLINS' )
      ind2 = which(allBOTH.filter$variable == 'NO_LIF_ROLLINS' & allBOTH.filter$mce.5hz >= 0.09)
      par(mfrow=c(2,1),cex=1.5,oma=c(.1,.1,.1,.1))
      plot(allBOTH.filter$mce.5hz[ind], allBOTH.filter$FinalEF[ind], xlab='', ylab=expression(paste('NO, g/kg')), main='Corn')
      points(allBOTH.filter$mce.5hz[ind2], allBOTH.filter$FinalEF[ind2], pch=19, col='orange')
      abline(lm(allBOTH.filter$FinalEF[ind]~allBOTH.filter$mce.5hz[ind]))
      cor.test(allBOTH.filter$mce.5hz[ind], allBOTH.filter$FinalEF[ind])
      ind = which(allBOTH.filter$variable == 'NO2_ACES_WOMACK' )
      ind2 = which(allBOTH.filter$variable == 'NO2_ACES_WOMACK' & allBOTH.filter$mce.5hz >= 0.09)
      plot(allBOTH.filter$mce.5hz[ind], allBOTH.filter$FinalEF[ind], xlab='MCE', ylab=expression(paste('NO'[2],', g/kg')), main='')
      points(allBOTH.filter$mce.5hz[ind2], allBOTH.filter$FinalEF[ind2], pch=19, col='orange')
      abline(lm(allBOTH.filter$FinalEF[ind]~allBOTH.filter$mce.5hz[ind]))
      cor.test(allBOTH.filter$mce.5hz[ind], allBOTH.filter$FinalEF[ind])
      
      
      ind = which(allBOTH.filter$variable == 'NO_LIF_ROLLINS' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz >R2filter)
      tmp.avg = aggregate(allBOTH.filter[ind,], by=list(allBOTH.filter$fuel[ind]), FUN='mean', na.rm=TRUE)
      indBW = which(allBOTH.blackwater.filter$variable == 'NO_LIF_ROLLINS' & is.finite(allBOTH.blackwater.filter$FinalEF))
                    
      #akagi NO: crop, 2.06 @ 0.925 MCE; NO2
      indA = which(akagi$Species == 'NO')
      indC  = which(xiaoxi.avg$var == 'NO')
      NOfireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel), size=4)+
        ylab("NO EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = c(0.19, 0.7))+
        theme(text = element_text(size = sz)) +
        geom_point(data=xiaoxi, aes(x=MCE, y=NO), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)+
        geom_point(data=allBOTH.blackwater.filter[indBW,],aes(x=mce.5hz,y=FinalEF), col='black', size=5, stroke=3)#+
       # geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)
      
      
      indA = which(akagi$Species == 'NO')
      NOfireEFavg = ggplot(tmp.avg) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
       # geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, size=fuel))+
        geom_point(data=tmp.avg,aes(x=mce.5hz,y=FinalEF,  shape=Group.1,col=Group.1), size=8, stroke=1)+
        ylab("NO EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = c(0.19, 0.7))+
        theme(text = element_text(size = sz)) +
        geom_point(data=xiaoxi.avg, aes(x=mean[1], y=mean[indC]), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        # annotate(geom="text", x=0.83, y=0.8, label="Liu et al., 2016",size=6, col='green')+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)
      
      NOfireER = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalERtoCO*1E3, col=fuel,shape=fuel), size=4)+
        # geom_point(data=tmp.avg,aes(x=mce.5hz,y=FinalEF,  shape=Group.1),col='grey', size=8, stroke=1)+
        ylab("NO ER, ppt/ppb CO") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        labs(col="",shape="")
      
      # correlation with MCE
      cor.test(allBOTH.filter$FinalEF[ind], allBOTH.filter$mce.5hz[ind])
      # ------ Figure 2B - NO2 fire -------------
      ind = which(allBOTH.filter$variable == 'NO2_ACES_WOMACK' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      indA = which(akagi$Species == 'NO2')
      indBW = which(allBOTH.blackwater.filter$variable == 'NO2_ACES_WOMACK' & is.finite(allBOTH.blackwater.filter$FinalEF))
      
      NO2fireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel,shape=fuel),size=4)+
        ylab(expression(paste("NO"[2]," EF, g/kg"))) + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        geom_point(data=xiaoxi, aes(x=MCE, y=NO2), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
      geom_point(data=allBOTH.blackwater.filter[indBW,],aes(x=mce.5hz,y=FinalEF), col='black', size=5, stroke=3)#+
      
      NO2fireER = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalERtoCO*1E3, col=fuel, shape=fuel),size=4)+
        # geom_point(data=tmp.avg,aes(x=mce.5hz,y=FinalEF,  shape=Group.1),col='grey', size=8, stroke=1)+
        ylab(expression(paste("NO"[2]," ER, ppt/ppb CO"))) + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        labs(col="",shape="")
      allNO = ggarrange(NOfireEF, NOfireER, NO2fireEF, NO2fireER, common.legend = TRUE)
      ggsave(AllNO, file='AllNO.ps')
      
      # correlation with MCE
      cor.test(allBOTH.filter$FinalEF[ind], allBOTH.filter$mce.5hz[ind])
      # ------ Figure 2C - HNO2 fire -------------
      ind = which(allBOTH.filter$variable == 'HNO2_ACES_WOMACK' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      indB = which(andreae$Species == 'HONO')
      
      HNO2fireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel),size=4, stroke=3)+
        ylab("HONO EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+#ylim(c(0,1.5))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +  labs(col="",shape="")+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)
        
      HNO2fireER = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalERtoCO*1E3, col=fuel, shape=fuel),size=4)+
        # geom_point(data=tmp.avg,aes(x=mce.5hz,y=FinalEF,  shape=Group.1),col='grey', size=8, stroke=1)+
        ylab("HONO ER, ppt/ppb CO") + xlab("MCE")  + xlim(c(0.8, 1))+#ylim(c(0,6.5))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        labs(col="",shape="")
      ind = which(allBOTH.filter$variable == 'HNO2_NOAACIMS_VERES' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      HNO2fireEF2 = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel),size=4, stroke=3)+
        ylab("HONO EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+#ylim(c(0,1.5))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +  labs(col="",shape="") +
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)
      
      HNO2fireER2 = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalERtoCO*1E3, col=fuel, shape=fuel),size=4)+
        ylab("HONO ER, ppt/ppb CO") + xlab("MCE")  + xlim(c(0.8, 1))+#ylim(c(0,6.5))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+
        theme(text = element_text(size = sz)) +
        labs(col="",shape="")#+@theme(legend.text = element_text(color = cbp1))
      
      allHONO=ggarrange(HNO2fireEF,  HNO2fireER,HNO2fireEF2, HNO2fireER2,
                common.legend = TRUE,labels = c("a) Womack","b) Womack","c) Veres","d) Veres"),
                hjust =-2)
      ggsave(AllHONO, file='AllHONO.ps')
      # theme(text = element_text(size = 20)) +
      #scale_shape_manual(values=c(19,15,17,18,7,8,4))+
      #  annotate(geom="text", x=0.83, y=1.1, label="Liu et al., 2016",size=6, col='green')
      
      # correlation with MCE
      cor.test(allBOTH.filter$FinalEF[ind], allBOTH.filter$mce.5hz[ind])
      
      # Which fires am I not using
      allfires=unique(allBOTH$fire)
      usingfires=unique(allBOTH.filter$fire)
      ind = which(allfires %nin% usingfires)
      
      # ------ Figure 2D - CH3CN fire -------------
      ind = which(allBOTH.filter$variable == 'CH3CN_NOAAPTR_ppbv_WARNEKE' & 
                    is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.5hz > R2filter)
      indB = which(andreae$Species == 'Acetonitrile')
      indA = which(akagi$Species == 'Acetonitrile')
      indK = which(allBOTH.filter.ag.avg$Group.1 == 'CH3CN_NOAAPTR_ppbv_WARNEKE')
      CH3CNfireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel,shape=fuel), size=4, stroke=3)+
        ylab("CH3CN EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        geom_point(data=xiaoxi, aes(x=MCE, y=Acetonitrile), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      #+
      #  geom_point(data=allBOTH.filter.ag.avg[indK,], aes(x=mce.5hz, y=FinalEF), col='black',size=5,shape=3, stroke=3)+
      #  geom_point(data=allBOTH.filter.sc.avg[indK,], aes(x=mce.5hz, y=FinalEF), col='black',size=5,shape=4, stroke=3)
      
      CH3CNfireER = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalERtoCO*1E3, col=fuel,shape=fuel), size=4, stroke=3)+
        ylab("CH3CN ER, ppt/ppbCO") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        #geom_point(data=xiaoxi, aes(x=MCE, y=Acetonitrile), pch=0,col='green',stroke=2)+
        labs(col="",shape="")
      +
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)+
        geom_point(data=allBOTH.filter.ag.avg[indK,], aes(x=mce.5hz, y=FinalEF), col='black',size=5,shape=3, stroke=3)+
        geom_point(data=allBOTH.filter.sc.avg[indK,], aes(x=mce.5hz, y=FinalEF), col='black',size=5,shape=4, stroke=3)
      
      
      ind = which(allBOTH.filter$variable == 'HCN_NOAAPTR_ppbv_WARNEKE' & 
                    is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.5hz > R2filter)
      indB = which(andreae$Species == 'HCN')
      indA = which(akagi$Species == 'HydrogenCyanide')
      indK = which(allBOTH.filter.ag.avg$Group.1 == 'HCN_NOAAPTR_ppbv_WARNEKE')
      HCNfireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel,shape=fuel), size=4, stroke=3)+
        ylab("HCN EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        geom_point(data=xiaoxi, aes(x=MCE, y=HCN), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)+
        geom_point(data=allBOTH.filter.ag.avg[indK,], aes(x=mce.5hz, y=FinalEF), col='black',size=5,shape=3, stroke=3)+
        geom_point(data=allBOTH.filter.sc.avg[indK,], aes(x=mce.5hz, y=FinalEF), col='black',size=5,shape=4, stroke=3)
      
      
      # theme(text = element_text(size = 20)) +
      #scale_shape_manual(values=c(19,15,17,18,7,8,4))+
      #annotate(geom="text", x=0.83, y=1.1, label="Liu et al., 2016",size=6, col='green')
      HCNfireER = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalERtoCO*1E3, col=fuel, shape=fuel))+
        # geom_point(data=tmp.avg,aes(x=mce.5hz,y=FinalEF,  shape=Group.1),col='grey', size=8, stroke=1)+
        ylab("HCN ER, ppt/ppb CO") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        geom_point(data=xiaoxi, aes(x=MCE, y=Acetonitrile), pch=0,col='green',stroke=2)+
        labs(col="",shape="")
      
      AllCH3CN=ggarrange(CH3CNfireEF, CH3CNfireER,common.legend = TRUE)
      ggsave(AllCH3CN, file='AllCH3CN.ps')
      # theme(text = element_text(size = 20)) +
      #scale_shape_manual(values=c(19,15,17,18,7,8,4))+
      #  annotate(geom="text", x=0.83, y=1.1, label="Liu et al., 2016",size=6, col='green')
      
      # correlation with MCE
      cor.test(allBOTH.filter$FinalEF[ind], allBOTH.filter$mce.5hz[ind])
      # ------ Figure 2H - SO2 fire -------------
      ind = which(allBOTH.filter$variable == 'SO2_LIF_ROLLINS' )
      ind2 = which(allBOTH.filter$variable == 'SO2_LIF_ROLLINS' & allBOTH.filter$mce.5hz >= 0.93)
      plot(allBOTH.filter$mce.5hz[ind], allBOTH.filter$FinalEF[ind], xlab='MCE', ylab=expression(paste('SO'[2],', g/kg')), main='Corn')
      points(allBOTH.filter$mce.5hz[ind2], allBOTH.filter$FinalEF[ind2], pch=19, col='orange')
      abline(lm(allBOTH.filter$FinalEF[ind2]~allBOTH.filter$mce.5hz[ind2]))
      cor.test(allBOTH.filter$mce.5hz[ind2], allBOTH.filter$FinalEF[ind2])
      indBW = which(allBOTH.blackwater.filter$variable == 'SO2_LIF_ROLLINS' & is.finite(allBOTH.blackwater.filter$FinalEF))
      
      ind = which(allBOTH.filter$variable == 'SO2_LIF_ROLLINS' & allBOTH.filter$fuel == 'slash')
      cor.test(allBOTH.filter$mce.5hz[ind], allBOTH.filter$FinalEF[ind])
      
      
      ind = which(allBOTH.filter$variable == 'SO2_LIF_ROLLINS' & 
                    is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.5hz > R2filter)
      indB = which(andreae$Species == 'SO2')
      indA = which(akagi$Species == 'SO2')
      
      SO2fireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel,shape=fuel), size=4)+
        ylab(expression(paste("SO"[2]," EF, g/kg"))) + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ))+
        scale_shape_manual(values =fuelshapes)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        geom_point(data=xiaoxi, aes(x=MCE, y=SO2), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)#+
      #  geom_point(data=allBOTH.blackwater.filter[indBW,],aes(x=mce.5hz,y=FinalEF), col='black', size=5, stroke=3)#+
      SO2fireER = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalERtoCO*1E3, col=fuel, shape=fuel), size=4)+
        # geom_point(data=tmp.avg,aes(x=mce.5hz,y=FinalEF,  shape=Group.1),col='grey', size=8, stroke=1)+
        ylab(expression(paste("SO"[2]," ER, ppt/ppb CO")))+ xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+#theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        labs(col="",shape="")
      AllSO2=ggarrange(SO2fireEF, SO2fireER,common.legend = TRUE)
      ggsave(AllSO2, file='AllSO2.ps')
      # theme(text = element_text(size = 20)) +
      #scale_shape_manual(values=c(19,15,17,18,7,8,4))+
      #  annotate(geom="text", x=0.83, y=1.1, label="Liu et al., 2016",size=6, col='green')
      
      # correlation with MCE
      cor.test(allBOTH.filter$FinalEF[ind], allBOTH.filter$mce.5hz[ind])
      
      # ------ Figure 2K - Nitrate fire -------------
      ind = which(allBOTH.filter$variable == 'Nitrate_JIMENEZ' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      #indB = which(andreae$Species == 'BC or EC')
      indA = which(akagi$Species == 'Nitrate')
      NitratefireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel), size=4, stroke=3)+
        ylab("Nitrate EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        geom_point(data=xiaoxi, aes(x=MCE, y=NO3), pch=0,col='green',stroke=2)+
        labs(col="",shape="") +
       # geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      # ------ Figure 2J - BC fire -------------
      ind = which(allBOTH.filter$variable == 'BC_SCHWARZ' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      indB = which(andreae$Species == 'BC or EC')
      indA = which(akagi$Species == 'BlackCarbon')
      BCfireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel), size=4, stroke=2)+
        ylab("BC EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        geom_point(data=xiaoxi, aes(x=MCE, y=BC), pch=0,col='green',stroke=2)+
        labs(col="",shape="") +
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      BCfireER = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalERtoCO*1E3, col=fuel, shape=fuel), size=4, stroke=2)+
        # geom_point(data=tmp.avg,aes(x=mce.5hz,y=FinalEF,  shape=Group.1),col='grey', size=8, stroke=1)+
        ylab("BC ER, ppt/ppb CO") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+#theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        labs(col="",shape="")
      
      AllBC=ggarrange(BCfireEF, BCfireER,common.legend = TRUE)
      ggsave(AllBC, file='AllBC.ps')
      # ------ Figure 2L - NH3 fire -------------
      ind = which(allBOTH.filter$variable == 'NH3_WISTHALER' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      indB = which(andreae$Species == 'NH3')
      indA = which(akagi$Species == 'Ammonia')
      NH3fireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel), size=4)+
        ylab("NH3 EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        labs(col="",shape="") +
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      NH3fireER = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalERtoCO*1E3, col=fuel, shape=fuel), size=4)+
        # geom_point(data=tmp.avg,aes(x=mce.5hz,y=FinalEF,  shape=Group.1),col='grey', size=8, stroke=1)+
        ylab("NH3 ER, ppt/ppb CO") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+#theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        labs(col="",shape="")
      
      AllNH3=ggarrange(NH3fireEF, NH3fireER,common.legend = TRUE)
      ggsave(AllNH3, file='AllNH3.ps')
      
      # correlation with MCE
      cor.test(allBOTH.filter$FinalEF[ind], allBOTH.filter$mce.5hz[ind])
      
      # ------ Figure 2I - Cl2 fire -------------
      ind = which(allBOTH.filter$variable == 'Cl2_NOAACIMS_VERES' & 
                    is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      indB = which(andreae$Species == 'Cl2')
      indA = which(akagi$Species == 'Cl2')
      
      Cl2fireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel,shape=fuel))+
        ylab(expression(paste("Cl2 EF, g/kg"))) + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        #geom_point(data=xiaoxi, aes(x=MCE, y=Acetonitrile), pch=0,col='green',stroke=2)+
        labs(col="",shape="")#+
      #  geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
      #  geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
       # geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      Cl2fireER = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalERtoCO*1E3, col=fuel, shape=fuel))+
        # geom_point(data=tmp.avg,aes(x=mce.5hz,y=FinalEF,  shape=Group.1),col='grey', size=8, stroke=1)+
        ylab(expression(paste("Cl2 ER, ppt/ppb CO")))+ xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+#theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        labs(col="",shape="")
      ggarrange(Cl2fireEF, Cl2fireER,common.legend = TRUE)
      ggsave(AllCl2, file='AllCl2.ps')
      # theme(text = element_text(size = 20)) +
      #scale_shape_manual(values=c(19,15,17,18,7,8,4))+
      #  annotate(geom="text", x=0.83, y=1.1, label="Liu et al., 2016",size=6, col='green')
      
      # correlation with MCE
      cor.test(allBOTH.filter$FinalEF[ind], allBOTH.filter$mce.5hz[ind])
      
      ind = which(allBOTH.filter$variable == 'OC_JIMENEZ')
      tmpOC = allBOTH.filter$FinalERtoCO[ind]
      ind = which(allBOTH.filter$variable == 'BC_SCHWARZ')
      tmpBC = allBOTH.filter$FinalERtoCO[ind]
      mce = allBOTH.filter$mce.5hz[ind]
      fuel = allBOTH.filter$fuel[ind]
      dat = as.data.frame(cbind(tmpOC, tmpBC, mce,fuel))
      dat$BC_OC = as.numeric(dat$tmpBC)/as.numeric(dat$tmpOC)
       ggplot(dat) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=as.numeric(mce), y=BC_OC, col=fuel, shape=fuel), size=4,stroke=2)+
        ylab("BC/OC") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) + ylim(c(0,0.3))
         
        geom_point(data=xiaoxi, aes(x=MCE, y=OA/1.8), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      
      # ------ Figure 2E - OA fire -------------
      ind = which(allBOTH.filter$variable == 'OC_JIMENEZ' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      indB = which(andreae$Species == 'OC')
      indA = which(akagi$Species == 'OrganicCarbon')
       OCfireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel), size=4,stroke=2)+
        ylab("OC EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        geom_point(data=xiaoxi, aes(x=MCE, y=OA/1.8), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      OCfireER = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalERtoCO*1E3, col=fuel, shape=fuel), size=4,stroke=2)+
        # geom_point(data=tmp.avg,aes(x=mce.5hz,y=FinalEF,  shape=Group.1),col='grey', size=8, stroke=1)+
        ylab("OC ER, ppt/ppb CO") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+#theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        labs(col="",shape="")
      
      ggarrange(OCfireEF, OCfireER,BCfireEF, BCfireER,common.legend = TRUE)
      
      # correlation with MCE
      cor.test(allBOTH.filter$FinalEF[ind], allBOTH.filter$mce.5hz[ind])
      
      # ------ Figure 2I - Syringol fire -------------
      ind = which(allBOTH.filter$variable == 'Guaiacol_NOAAPTR_ppbv_WARNEKE' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      #indB = which(andreae$Species == 'Glycolaldehyde_acetic acid')
      #indA = which(akagi$Species == 'Glycolaldehyde')
      
      GuaiacolfireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel), size=4, stroke=2)+
        ylab("Guaiacol EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        # geom_point(data=xiaoxi, aes(x=MCE, y=Formaldehyde), pch=0,col='green',stroke=2)+
        labs(col="",shape="")#+
      #geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
      #geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
      #geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      # ------ Figure 2N - Hydroxyacetone fire -------------
      ind = which(allBOTH.filter$variable == 'C3H6O2_NOAAPTR_ppbv_WARNEKE' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      indB = which(andreae$Species == 'Acetol (hydroxyacetone) ')
      indB=indB[1]
      indA = which(akagi$Species == 'Acetol')
      
      C3H6O2fireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel), size=4, stroke=2)+
        ylab("C3H6O2 EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        ggtitle('(sum of methyl acetate ethyl formate and hydroxyacetone)')+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
       geom_point(data=xiaoxi, aes(x=MCE, y=Hydroxyacetone), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      
      # ------ Figure 2I - Syringol fire -------------
      ind = which(allBOTH.filter$variable == 'Syringol_NOAAPTR_ppbv' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      #indB = which(andreae$Species == 'Glycolaldehyde_acetic acid')
      #indA = which(akagi$Species == 'Glycolaldehyde')
      
      SyringolfireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel), size=4, stroke=2)+
        ylab("Syringol EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        # geom_point(data=xiaoxi, aes(x=MCE, y=Formaldehyde), pch=0,col='green',stroke=2)+
        labs(col="",shape="")#+
        #geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        #geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        #geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      # ------ Figure 2Z - Napthalene fire -------------
      ind = which(allBOTH.filter$variable == 'Naphthalene_NOAAPTR_ppbv_WARNEKE' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      indB = which(andreae$Species == 'PAHs')
      #indA = which(akagi$Species == 'Ethylene')
      
      PAHfireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel), size=4, stroke=2)+
        ylab("Naphthalene EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
      #  geom_point(data=xiaoxi, aes(x=MCE, y=Monoterpenes), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)#+
      #  geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
      #  geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      # ------ Figure 2Z - Monoterpenes fire -------------
      ind = which(allBOTH.filter$variable == 'Monoterpenes_NOAAPTR_ppbv_WARNEKE' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter )
      indB = which(andreae$Species == 'Terpenes')
      #indA = which(akagi$Species == 'Ethylene')
      
      TerpenesfireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel), size=4, stroke=2)+
        ylab("Terpenes EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        geom_point(data=xiaoxi, aes(x=MCE, y=Monoterpenes), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)#+
      #  geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
      #  geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      # ------ Figure 2Q - Isoprene fire -------------
      ind = which(allBOTH.filter$variable == 'Isoprene_NOAAPTR_ppbv_WARNEKE' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter )
      indB = which(andreae$Species == 'Isoprene')
      indA = which(akagi$Species == 'Isoprene')
      
      IsoprenefireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel), size=4, stroke=2)+
        ylab("Isoprene EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        geom_point(data=xiaoxi, aes(x=MCE, y=Isoprenepentadienescyclopentenefuran), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      
      # ------ Figure 2Z - Isoprene fire -------------
      ind = which(allBOTH.filter$variable == 'Isoprene_ppt' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      indB = which(andreae$Species == 'Isoprene')
      indA = which(akagi$Species == 'Isoprene')
      
      IsoprenefireEF2 = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel), size=4, stroke=2)+ylim(c(0,4))+
        ylab("Isoprene (TOGA) EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        # geom_point(data=xiaoxi, aes(x=MCE, y=Formaldehyde), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      # ------ Figure 2W - Ethene fire -------------
      ind = which(allBOTH.filter$variable == 'Ethene_WAS_BLAKE' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      indB = which(andreae$Species == 'C2H4')
      indA = which(akagi$Species == 'Ethylene')
      
      EthenefireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel), size=4, stroke=2)+
        ylab("Ethene EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        # geom_point(data=xiaoxi, aes(x=MCE, y=Formaldehyde), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      
      # ------ Figure 2G - GLYC fire -------------
      ind = which(allBOTH.filter$variable == 'GlycolaldehydeCH3COOH_NOAAPTR_ppbv_WARNEKE' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      indB = which(andreae$Species == 'Glycolaldehyde_acetic acid')
      indA = which(akagi$Species == 'Glycolaldehyde')
      
      GLYCfireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel), size=4, stroke=2)+
        ylab("Glycoaldehyde/Acetic Acid EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
       # geom_point(data=xiaoxi, aes(x=MCE, y=Formaldehyde), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      # ------ Figure 2F - Benz fire -------------
      ind = which(allBOTH.filter$variable == 'Toluene_NOAAPTR_ppbv_WARNEKE' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      indB = which(andreae$Species == 'Toluene')
      indA = which(akagi$Species == 'Toluene')
      
      TOLUfireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel), size=4, stroke=2)+
        ylab("Toluene EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        geom_point(data=xiaoxi, aes(x=MCE, y=Toluene), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      
      # ------ Figure 2F - Benz fire -------------
      ind = which(allBOTH.filter$variable == 'Benzene_NOAAPTR_ppbv_WARNEKE' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      indB = which(andreae$Species == 'Benzene')
      indA = which(akagi$Species == 'Benzene')
      
      BENZfireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel), size=4, stroke=2)+
        ylab("Benzene EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
         geom_point(data=xiaoxi, aes(x=MCE, y=Benzene), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      
      # ------ Figure 2F - CH2O fire -------------
      ind = which(allBOTH.filter$variable == 'CH2O_ISAF_HANISCO' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      indB = which(andreae$Species == 'Formaldehyde')
      indA = which(akagi$Species == 'Formaldehyde')
      
      CH2OfireEF = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel), size=4, stroke=2)+
        ylab("CH2O EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        geom_point(data=xiaoxi, aes(x=MCE, y=Formaldehyde), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      CH2OfireER = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalERtoCO*1E3, col=fuel, shape=fuel), size=4)+
        # geom_point(data=tmp.avg,aes(x=mce.5hz,y=FinalEF,  shape=Group.1),col='grey', size=8, stroke=1)+
        ylab("CH2O ER, ppt/ppb CO") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+#theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        labs(col="",shape="")
      ind = which(allBOTH.filter$variable == 'CH2O_CAMS_pptv_FRIED' & is.finite(allBOTH.filter$FinalEF) &
                    allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
      
      CH2OfireEF2 = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel,shape=fuel), size=4, stroke=2)+
        ylab("CH2O EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        geom_point(data=xiaoxi, aes(x=MCE, y=Formaldehyde), pch=0,col='green',stroke=2)+
        labs(col="",shape="")+
        geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
        geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
        geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
      
      CH2OfireER2 = ggplot(allBOTH.filter[ind,]) + 
        theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
        geom_point(aes(x=mce.5hz, y=FinalERtoCO*1E3, col=fuel, size=fuel))+
        # geom_point(data=tmp.avg,aes(x=mce.5hz,y=FinalEF,  shape=Group.1),col='grey', size=8, stroke=1)+
        ylab("CH2O ER, ppt/ppb CO") + xlab("MCE")  + xlim(c(0.8, 1))+
        scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
        scale_shape_manual(values =fuelshapes, limits=fuellimits)+
        theme(legend.background=element_blank())+#theme(legend.position = "none")+
        theme(text = element_text(size = sz)) +
        labs(col="",shape="")
      
      ggarrange(CH2OfireEF, CH2OfireEF2, common.legend = TRUE,
                labels = c("a) Hanisco ","b) Fried"),
                hjust =-2)
    }
    
   
    
    par(cex=1.2)
    ind = which(VOCsum$fuel == 'corn')
    ggplot(VOCsum[ind,])+geom_point(aes(x=mce.5hz, y=ERtoCO.5hz, col=NN),size=5)+
         xlab('MCE')+ ylab('VOC, ppt/ppb CO')+ggtitle('Corn')
    tl = lm(VOCsum$ERtoCO.5hz[ind]~VOCsum$mce.5hz[ind])
    abline(tl)
    # get pie chart
    ind = which(is.finite(VOCsum$FinalEF))
    forpie = which(test$uniqueid %in% unique(VOCsum$ff[ind]))
    testpie = test[forpie,]
    #testpie.avg = aggregate(testpie)
    
    # ------ Figure 2Q - CHOCHO fire -------------
    ind = which(allBOTH.filter$variable == 'CHOCHO_ACES_WOMACK' & is.finite(allBOTH.filter$FinalEF) &
                  allBOTH.filter$FinalR2 > R2filter & allBOTH.filter$R2toCO.1hz > R2filter)
    indB = which(andreae$Species == 'Glyoxal')
    #indA = which(akagi$Species == 'Glyoxal')
    zER = c(0.0014,0.0011,0.0012,0.0008,0.0012,0.0003,0.0014,0.0016,0.0020,0.0021,0.0009, 0.0015,0.0010, 0.0028,0.0014,0.0018, 0.0044)
    zMCE = c(0.93,0.917,  0.927,0.945,NaN,0.919,0.92,NaN,0.974,0.962,0.941,0.936,  0.93,NaN, 0.948,NaN, 0.969)
    zERCO2CO = zMCE/(1-zMCE)
    zEF = 0.41*1E3*mWCHOCHO/12 * zER/(1+zERCO2CO)
    zarzana =as.data.frame(cbind(zER,zMCE,zEF))
    
    CHOCHOfireEF = ggplot(allBOTH.filter[ind,]) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=MCE, y=FinalEF, col=fuel, shape=fuel), size=4, stroke=2)+
      ylab("CHOCHO EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
      scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
      scale_shape_manual(values =fuelshapes, limits=fuellimits)+
      theme(legend.background=element_blank())+theme(legend.position = "none")+
      theme(text = element_text(size = sz)) +
      labs(col="",shape="")+
      geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
      geom_point(data=zarzana, aes(x=zMCE,y=zEF))
    #      geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
#      geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)
    
    CHOCHOfireER = ggplot(allBOTH.filter[ind,]) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=MCE, y=FinalERtoCO*1E3, col=fuel, shape=fuel), size=4, stroke=2)+
      # geom_point(data=tmp.avg,aes(x=mce.5hz,y=FinalEF,  shape=Group.1),col='grey', size=8, stroke=1)+
      ylab("CHOCHO ER, ppt/ppb CO") + xlab("MCE")  + xlim(c(0.8, 1))+
      scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
      scale_shape_manual(values =fuelshapes, limits=fuellimits)+
      theme(legend.background=element_blank())+#theme(legend.position = "none")+
      theme(text = element_text(size = sz)) +
      labs(col="",shape="")+
      geom_point(data=zarzana, aes(x=zMCE,y=zER*1E3))
    
    ggarrange(CHOCHOfireEF, CHOCHOfireER, common.legend = TRUE,
              labels = c("a) EF ","b) ER"),
              hjust =-2)

    # ---------------------------------actually make it-------
    xiaoxi$VOC = xiaoxi$Formaldehyde + xiaoxi$Methanol + xiaoxi$Hydroxyacetone + xiaoxi$Acetaldehyde + xiaoxi$MVKMACRcrotonaldehyde+
      xiaoxi$Isoprenehydroperoxyaldehydes+xiaoxi$Isoprenepentadienescyclopentenefuran+xiaoxi$Benzene+xiaoxi$Monoterpenes+xiaoxi$Toluene
    indB = which(andreae$Species == 'Total NMOG, including unidentifiedd'); indB=indB[1]
    indA = which(akagi$Species == 'NMOC') # identified, identified + unidenentifed, use #2
    indA = indA[2]
    VOCfireEF = ggplot(VOCsum) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce.5hz, y=FinalEF, col=fuel, shape=fuel),size=5, stroke=2)+
      ylab("VOC EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
      scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
      scale_shape_manual(values =fuelshapes, limits=fuellimits)+
      theme(legend.background=element_blank())+theme(legend.position = "none")+
      theme(text = element_text(size = sz)) +
      geom_point(data=xiaoxi, aes(x=MCE, y=VOC), pch=0,col='green',stroke=2)+
      geom_point(data=andreae, aes(x=as.numeric(andreae$average[1]), y=as.numeric(andreae$average[indB])), col='purple',stroke=3, size=5, shape=3)+
      geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3, stroke=3)+
      geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4, stroke=3)+
      labs(col="",shape="")#+
    VOCfireER = ggplot(VOCsum) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce.5hz, y=ERtoCO.5hz*1E3, col=fuel, shape=fuel),size=5, stroke=2)+
      # geom_point(data=tmp.avg,aes(x=mce.5hz,y=FinalEF,  shape=Group.1),col='grey', size=8, stroke=1)+
      ylab("VOC ER, ppt/ppb CO") + xlab("MCE")  + xlim(c(0.8, 1))+
      scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
      scale_shape_manual(values =fuelshapes, limits=fuellimits)+
      theme(legend.background=element_blank())+#theme(legend.position = "none")+
      theme(text = element_text(size = sz)) +
      labs(col="",shape="")
    
    
    VOCfireEFshort = ggplot(VOCsum) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce.5hz, y=FinalEFshort, col=fuel, shape=fuel),size=5, stroke=2)+
      ylab("Short-lived VOC EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
      scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
      scale_shape_manual(values =fuelshapes, limits=fuellimits)+
      theme(legend.background=element_blank())+theme(legend.position = "none")+
      theme(text = element_text(size = sz)) +
      #  geom_point(data=xiaoxi, aes(x=MCE, y=Formaldehyde), pch=0,col='green',stroke=2)+
      labs(col="",shape="")#+
    VOCfireERshort = ggplot(VOCsum) + 
      theme_classic()+
      geom_point(aes(x=mce.5hz, y=ERtoCO.5hzshort*1E3, col=fuel),size=5)+
      ylab("VOC ER, ppt/ppb CO") + xlab("MCE")  + xlim(c(0.8, 1))+
      scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
      scale_shape_manual(values =fuelshapes, limits=fuellimits)+
      theme(legend.background=element_blank())+#theme(legend.position = "none")+
      theme(text = element_text(size = sz)) +
      labs(col="",shape="")
    
    VOCfireEFlong = ggplot(VOCsum) + 
      theme_classic()+#
      geom_point(aes(x=mce.5hz, y=FinalEFlong,  col=fuel, shape=fuel),size=5, stroke=2)+
      ylab("Long-lived VOC EF, g/kg") + xlab("MCE")  + xlim(c(0.8, 1))+
      scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
      scale_shape_manual(values =fuelshapes, limits=fuellimits)+
      theme(legend.background=element_blank())+theme(legend.position = "none")+
      theme(text = element_text(size = sz)) +
      labs(col="",shape="")#+
    VOCfireERlong = ggplot(VOCsum) + 
      theme_classic()+
      geom_point(aes(x=mce.5hz, y=ERtoCO.5hzlong*1E3, col=fuel),size=5)+
      ylab("VOC ER, ppt/ppb CO") + xlab("MCE")  + xlim(c(0.8, 1))+
      scale_color_manual(values = c(cbp1 ), limits=fuellimits)+
      scale_shape_manual(values =fuelshapes, limits=fuellimits)+
      theme(legend.background=element_blank())+#theme(legend.position = "none")+
      theme(text = element_text(size = sz)) +
      labs(col="",shape="")
    ggarrange(VOCfireEFshort, VOCfireERshort,VOCfireEFlong, VOCfireERlong,common.legend = TRUE, labels=c("Short-lived EF","Short-lived ER","Long-lived EF","Long-lived ER"))
    ggarrange(VOCfireEF, VOCfireER)
    # which fires have a relationship of VOC ER with MCE
    ind = which(VOCsum$fuel == 'corn')
    cor.test(VOCsum$mce.5hz[ind], VOCsum$ERtoCO.5hz[ind])
    ind = which(VOCsum$fuel == 'rice')
    cor.test(VOCsum$mce.5hz[ind], VOCsum$ERtoCO.5hz[ind])
    ind = which(VOCsum$fuel == 'soybean')
    cor.test(VOCsum$mce.5hz[ind], VOCsum$ERtoCO.5hz[ind])
             
    ind = which(VOCsum$fuel == 'slash')
    cor.test(VOCsum$mce.5hz[ind], VOCsum$ERtoCO.5hz[ind])
    ind = which(VOCsum$fuel == 'grass')
    cor.test(VOCsum$mce.5hz[ind], VOCsum$ERtoCO.5hz[ind])
    ind = which(VOCsum$fuel == 'pile')
    cor.test(VOCsum$mce.5hz[ind], VOCsum$ERtoCO.5hz[ind])
    
    # ----------- Average agricultural table ---------- -----
    
    
    # -------- Is corn distinct from others? ----------
    allBOTH.filter.avg$CornvsFuel1hz = NaN
    allBOTH.filter.avg$CornvsFuel5hz = NaN
    for (i in 1:length(allBOTH.filter.avg$Group.1)){
      # Get everything for corn
      ind = which(allBOTH.filter$fuel == 'corn' & allBOTH.filter$variable == allBOTH.filter.avg$Group.1[i] )
      # Get everything for other fuels
      ind2 = which(allBOTH.filter$fuel != 'corn' & allBOTH.filter$variable == allBOTH.filter.avg$Group.1[i])
      # Maake sure there is enough corn data
      tti = which(!is.nan(allBOTH.filter$EF1.1hz[ind]) & !is.na(allBOTH.filter$EF1.1hz[ind]))
      tti2 = which(!is.nan(allBOTH.filter$EF1.1hz[ind2]) & !is.na(allBOTH.filter$EF1.1hz[ind2]))
      if (length(tti) > 1 & length(tti2) > 1){
        tcornvsother.1hz = t.test(allBOTH.filter$EF1.1hz[ind], allBOTH.filter$EF1.1hz[ind2])
        allBOTH.filter.avg$CornvsFuel1hz[i] = tcornvsother.1hz$p.value
      }  
      tti = which(!is.nan(allBOTH.filter$EF1.5hz[ind]) & !is.na(allBOTH.filter$EF1.5hz[ind]))
      tti2 = which(!is.nan(allBOTH.filter$EF1.5hz[ind2]) & !is.na(allBOTH.filter$EF1.5hz[ind2]))
      if (length(tti) > 1 & length(tti2) > 1){
        tcornvsother.5hz = t.test(allBOTH.filter$EF1.5hz[ind], allBOTH.filter$EF1.5hz[ind2])
        allBOTH.filter.avg$CornvsFuel5hz[i] = tcornvsother.5hz$p.value
      }
    }
    
    ind = which(allBOTH.filter.avg$CornvsFuel5hz < 0.05) # 43/60
    ind = which(allBOTH.filter.avg$CornvsFuel1hz < 0.05) # 101/262
    
    
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
    
    # ---------- Average Table by Fuel --------
    # ----------- Different fuels
    allBOTH.filter.avg.fuel$FinalERtoCO.sd = allBOTH.filter.sd.fuel$FinalERtoCO
    allBOTH.filter.avg.fuel$FinalEF.sd = allBOTH.filter.sd.fuel$FinalEF
    
    ind = which(allBOTH.filter.avg.fuel$Group.1 == 'corn')
    tmpCORN = allBOTH.filter.avg.fuel[ind,] ; tmpCORN$Group.1 = tmpCORN$Group.2
    tmp2 = allBOTH.filter.sd.fuel[ind,]
    ind2 = which(allBOTH.filter.ag$fuel == 'corn')
    tmpCORN  = getplumesANDmce(tmpCORN, allBOTH.filter.ag[ind2,] )
    tmpCORN$kind = ''
    for (i in 1:length(tmpCORN$kind.1hz)){
      ind = which(tmpCORN$Group.1[i] == allBOTH.filter$variable)
      tmpCORN$kind[i] = allBOTH.filter$kind[ind[1]]
      tmpCORN$Category.5hz[i] = as.numeric(allBOTH.filter$Category.5hz[ind[1]])
      tmpCORN$Category.1hz[i] = as.numeric(allBOTH.filter$Category.1hz[ind[1]])
    }
    ind = which(allBOTH.filter.sd.fuel$variable == 'corn')
    AverageTableCORN = getAverageTable(tmpCORN, allBOTH.filter.sd.fuel[ind,])
    ind = which(tmpCORN$Group.1 == "CO_DACOM_DISKIN")
    tmpCORN$mce.5hz[ind]
    ind = which(tmp2$Group.2 == "CO_DACOM_DISKIN")
    tmp2$mce.5hz[ind]
    write.csv(AverageTableCORN,file="AverageTableCorn.csv")
    ind = which(allBOTH.filter.avg.fuel$Group.1 == 'rice')
    tmpRICE = allBOTH.filter.avg.fuel[ind,] ; tmpRICE$Group.1 = tmpRICE$Group.2
    tmp2 = allBOTH.filter.sd.fuel[ind,]
    ind2 = which(allBOTH.filter.ag$fuel == 'rice')
    tmpRICE  = getplumesANDmce(tmpRICE, allBOTH.filter.ag[ind2,] )
    AverageTableRICE = getAverageTable(tmpRICE, allBOTH.filter.sd.fuel[ind,])
    ind = which(tmpRICE$Group.1 == "CO_DACOM_DISKIN")
    tmpRICE$mce.5hz[ind]
    ind = which(tmp2$Group.2 == "CO_DACOM_DISKIN")
    tmp2$mce.5hz[ind]
    
    ind = which(allBOTH.filter.avg.fuel$Group.1 == 'soybean')
    tmpSOYBEAN = allBOTH.filter.avg.fuel[ind,] ; tmpSOYBEAN$Group.1 = tmpSOYBEAN$Group.2
    tmp2 = allBOTH.filter.sd.fuel[ind,]
    ind2 = which(allBOTH.filter.ag$fuel == 'soybean')
    tmpSOYBEAN  = getplumesANDmce(tmpSOYBEAN, allBOTH.filter.ag[ind2,] )
    AverageTableSOYBEAN = getAverageTable(tmpSOYBEAN, allBOTH.filter.sd.fuel[ind,])
    ind = which(tmpSOYBEAN$Group.1 == "CO_DACOM_DISKIN")
    tmpSOYBEAN$mce.5hz[ind]
    ind = which(tmp2$Group.2 == "CO_DACOM_DISKIN")
    tmp2$mce.5hz[ind]
    
    ind = which(allBOTH.filter.avg.fuel$Group.1 == 'winter wheat')
    tmpWHEAT = allBOTH.filter.avg.fuel[ind,] ; tmpWHEAT$Group.1 = tmpWHEAT$Group.2
    tmp2 = allBOTH.filter.sd.fuel[ind,]
    ind2 = which(allBOTH.filter.ag$fuel == 'winter wheat')
    tmpWHEAT  = getplumesANDmce(tmpWHEAT, allBOTH.filter.ag[ind2,] )
    AverageTableWHEAT = getAverageTable(tmpWHEAT, allBOTH.filter.sd.fuel[ind,])
    ind = which(tmpWHEAT$Group.1 == "CO_DACOM_DISKIN")
    tmpWHEAT$mce.5hz[ind]
    ind = which(tmp2$Group.2 == "CO_DACOM_DISKIN")
    tmp2$mce.5hz[ind]
    
    ind = which(allBOTH.filter.avg.fuel$Group.1 == 'grass')
    tmpGRASS = allBOTH.filter.avg.fuel[ind,] ; tmpGRASS$Group.1 = tmpGRASS$Group.2
    tmp2 = allBOTH.filter.sd.fuel[ind,]
    ind2 = which(allBOTH.filter.sc$fuel == 'grass')
    tmpGRASS  = getplumesANDmce(tmpGRASS, allBOTH.filter.sc[ind2,] )
    AverageTableGRASS = getAverageTable(tmpGRASS, allBOTH.filter.sd.fuel[ind,])
    ind = which(tmpGRASS$Group.1 == "CO_DACOM_DISKIN")
    tmpGRASS$mce.5hz[ind]
    ind = which(tmp2$Group.2 == "CO_DACOM_DISKIN")
    tmp2$mce.5hz[ind]
    
    ind = which(allBOTH.filter.avg.fuel$Group.1 == 'pile')
    tmpPILE = allBOTH.filter.avg.fuel[ind,] ; tmpPILE$Group.1 = tmpPILE$Group.2
    tmp2 = allBOTH.filter.sd.fuel[ind,]
    ind2 = which(allBOTH.filter.sc$fuel == 'pile')
    tmpPILE  = getplumesANDmce(tmpPILE, allBOTH.filter.sc[ind2,] )
    AverageTablePILE = getAverageTable(tmpPILE, allBOTH.filter.sd.fuel[ind,])
    ind = which(tmpPILE$Group.1 == "CO_DACOM_DISKIN")
    tmpPILE$mce.5hz[ind]
    ind = which(tmp2$Group.2 == "CO_DACOM_DISKIN")
    tmp2$mce.5hz[ind]
    
    ind = which(allBOTH.filter.avg.fuel$Group.1 == 'slash')
    tmpSLASH = allBOTH.filter.avg.fuel[ind,] ; tmpSLASH$Group.1 = tmpSLASH$Group.2
    tmp2 = allBOTH.filter.sd.fuel[ind,]
    ind2 = which(allBOTH.filter.sc$fuel == 'slash')
    tmpSLASH  = getplumesANDmce(tmpSLASH, allBOTH.filter.sc[ind2,] )
    AverageTableSLASH = getAverageTable(tmpSLASH, allBOTH.filter.sd.fuel[ind,])
    ind = which(tmpSLASH$Group.1 == "CO_DACOM_DISKIN")
    tmpSLASH$mce.5hz[ind]
    ind = which(tmp2$Group.2 == "CO_DACOM_DISKIN")
    tmp2$mce.5hz[ind]
    
    test = cbind(AverageTableCORN, AverageTableSOYBEAN, AverageTableRICE,AverageTableWHEAT,AverageTableAG,
                 AverageTableGRASS, AverageTablePILE, AverageTableSLASH, AverageTableSC)
    ind = which(is.finite(test[,2]))
    test = test[ind,]
    write.csv(test, file='AllEFs.csv')
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
    
    # ----- Willow
    ind = which(allBOTH.filter$fire == 'Willow')
    allBOTH.filter.Willow = allBOTH.filter[ind,]
    ind = which(allBOTH.filter.Willow$variable == 'CO_DACOM_DISKIN')
    ggplot(allBOTH.filter.Willow[ind,]) + geom_point(aes(x=mce.5hz, y=EF1.5hz, col=MAtoF.5hz))+
      scale_color_viridis()
    
    # monoterpene outliers
    ind = which(allBOTH.filter$variable == 'Monoterpenes_NOAAPTR_ppbv_WARNEKE')
    ggplot(allBOTH.filter[ind,]) + 
      theme_classic()+geom_label(aes(x=mce.5hz, y=EF1.5hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce.5hz, y=EF1.5hz, col=fuel),size=5)
    
    # HONO
    ind = which(allBOTH.filter$variable == 'HNO2_NOAACIMS_VERES' & is.finite(allBOTH.filter$ChosenEF.1hz))
    varHNO2 = as.data.frame(allBOTH.filter$ChosenEF.1hz[ind])
    colnames(varHNO2) = 'HNO2_NOAACIMS_VERES' 
    varHNO2$mce = allBOTH.filter$mce.1hz[ind]
    varHNO2$EF1.1hz = allBOTH.filter$EF1.1hz[ind]
    varHNO2$EF1CO.1hz = allBOTH.filter$EF1CO.1hz[ind]
    varHNO2$EF1int.1hz = allBOTH.filter$EF1int.1hz[ind]
    varHNO2$EF1COint.1hz = allBOTH.filter$EF1COint.1hz[ind]
    
    varHNO2$mce5 = allBOTH.filter$mce.5hz[ind]
    goodpasses2 = allBOTH.filter$uniqueid[ind]
    tmp = which(allBOTH.filter$uniqueid %in% goodpasses2)
    test = allBOTH.filter[tmp,]
    ind = which(test$variable == 'CO_DACOM_DISKIN')
    varHNO2$maxCO = test$maxval.5hz[ind]
    varHNO2$fuel = test$fuel[ind]
    
    # HCHO
    ind = which(allBOTH.filter$variable == "CH2O_ISAF_HANISCO" & is.finite(allBOTH.filter$FinalEF))
    varCH2O = as.data.frame(allBOTH.filter$FinalEF[ind])
    colnames(varCH2O) = 'CH2O_ISAF_HANISCO' 
    varCH2O$mce = allBOTH.filter$mce.5hz[ind]
    varCH2O$EF1.5hz = allBOTH.filter$EF1.5hz[ind]
    varCH2O$EF1CO.5hz = allBOTH.filter$EF1CO.5hz[ind]
    varCH2O$EF1int.5hz = allBOTH.filter$EF1int.5hz[ind]
    varCH2O$EF1COint.5hz = allBOTH.filter$EF1COint.5hz[ind]
    varCH2O$mce5 = allBOTH.filter$mce.5hz[ind]
    goodpasses2 = allBOTH.filter$uniqueid[ind]
    tmp = which(allBOTH.filter$uniqueid %in% goodpasses2)
    test = allBOTH.filter[tmp,]
    ind = which(test$variable == 'CO_DACOM_DISKIN')
    varCH2O$maxCO = test$maxval.5hz[ind]
    varCH2O$fuel = test$fuel[ind]
    indB = which(andreae$Species == 'Formaldehyde')
    indA = which(akagi$Species == 'Formaldehyde')
    CH2Ofire=ggplot(varCH2O) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce5, y=CH2O_ISAF_HANISCO, col=fuel,size=maxCO))+
      ylab("CH2O (ISAF) EF, g/kg") + xlab("MCE")+
      theme(legend.background=element_blank())+theme(legend.position = c(0.18, 0.7))+
      theme(text = element_text(size = sz))+
      scale_shape_manual(values=c(19,15,17,18,7,8,4))+
      scale_color_manual(values = c(cbp1 ),name="")+xlim(c(0.8,1))+  
      geom_point(data=xiaoxi, aes(x=MCE, y=Formaldehyde), pch=0,col='green')+
      scale_shape_manual(values=c(19,15,17,18,7,8,4))+
      scale_color_manual(values = c(cbp1 ))+xlim(c(0.8,1))+
      geom_point(data=andreae, aes(x=andreae$average[1], y=andreae$average[49]), col='purple', size=5, shape=3)+
      geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[42]), col='pink', size=5, shape=3)+
      geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[42]), col='pink', size=5, shape=4)
    
    # Toluene
    ind = which(allBOTH.filter$variable == 'Toluene_NOAAPTR_ppbv_WARNEKE' & is.finite(allBOTH.filter$ChosenEF.5hz))
    varTOLUENE = as.data.frame(allBOTH.filter$ChosenEF.5hz[ind])
    colnames(varTOLUENE) = 'Toluene_NOAAPTR_ppbv_WARNEKE' 
    varTOLUENE$mce = allBOTH.filter$mce.5hz[ind]
    varTOLUENE$EF1.5hz = allBOTH.filter$EF1.5hz[ind]
    varTOLUENE$EF1CO.5hz = allBOTH.filter$EF1CO.5hz[ind]
    varTOLUENE$EF1int.5hz = allBOTH.filter$EF1int.5hz[ind]
    varTOLUENE$EF1COint.5hz = allBOTH.filter$EF1COint.5hz[ind]
    varTOLUENE$mce5 = allBOTH.filter$mce.5hz[ind]
    goodpasses2 = allBOTH.filter$uniqueid[ind]
    tmp = which(allBOTH.filter$uniqueid %in% goodpasses2)
    test = allBOTH.filter[tmp,]
    ind = which(test$variable == 'CO_DACOM_DISKIN')
    varTOLUENE$maxCO = test$maxval.5hz[ind]
    varTOLUENE$fuel = test$fuel[ind]
    indB = which(andreae$Species == 'Toluene')
    indA = which(akagi$Species == 'Toluene')
    TOLUENEfire=ggplot(varTOLUENE) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce5, y=Toluene_NOAAPTR_ppbv_WARNEKE, col=fuel,size=maxCO))+
      ylab("Toluene EF, g/kg") + xlab("MCE")+
      theme(legend.background=element_blank())+theme(legend.position = c(0.18, 0.7))+
      theme(text = element_text(size = sz))+
      scale_shape_manual(values=c(19,15,17,18,7,8,4))+
      scale_color_manual(values = c(cbp1 ),name="")+xlim(c(0.8,1))+  
      geom_point(data=xiaoxi, aes(x=MCE, y=Toluene), pch=0,col='green')+
      scale_shape_manual(values=c(19,15,17,18,7,8,4))+
      scale_color_manual(values = c(cbp1 ))+xlim(c(0.8,1))+
      geom_point(data=andreae, aes(x=andreae$average[1], y=andreae$average[indB]), col='purple', size=5, shape=3)+
      geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[indA]), col='pink', size=5, shape=3)+
      geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[indA]), col='pink', size=5, shape=4)
    
    # Glyoxal
    ind = which(allBOTH.filter$names == 'Glyoxal' & is.finite(allBOTH.filter$FinalERtoCO))
    varCHOCHO = as.data.frame(allBOTH.filter$FinalERtoCO[ind])
    colnames(varCHOCHO) = 'CHOCHO_ACES_WOMACK' 
    varCHOCHO$mce = allBOTH.filter$MCE[ind]
    varCHOCHO$FinalEF = allBOTH.filter$FinalEF[ind]
    varCHOCHO$FinalERtoCO = allBOTH.filter$FinalERtoCO[ind]
    goodpasses2 = allBOTH.filter$uniqueid[ind]
    tmp = which(allBOTH.filter$uniqueid %in% goodpasses2)
    test = allBOTH.filter[tmp,]
    ind = which(test$variable == 'CO_DACOM_DISKIN')
    varCHOCHO$maxCO = test$maxval.5hz[ind]
    varCHOCHO$fuel = test$fuel[ind]
    
    CH3COCHOfire=ggplot(varCH3COCHO) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce5, y=EF1.1hz, col=fuel,size=maxCO))+
      # geom_point(aes(x=mce5, y=CH3COCHO_ACES_WOMACK, col='red',size=maxCO))+
      ylab("Methylglyoxal EF, g/kg") + xlab("MCE")+
      theme(legend.background=element_blank())+theme(legend.position = c(0.18, 0.7))+
      theme(text = element_text(size = sz))+
      scale_shape_manual(values=c(19,15,17,18,7,8,4))+
      scale_color_manual(values = c(cbp1 ),name="")+xlim(c(0.8,1))+
      geom_point(data=andreae, aes(x=andreae$average[1], y=andreae$average[53]), col='purple', size=5, shape=3)
    
    cor.test(varCH3COCHO$CH3COCHO_ACES_WOMACK, varCH3COCHO$mce5)
    ind = which(varCH3COCHO$fuel == 'corn')
    cor.test(varCH3COCHO$CH3COCHO_ACES_WOMACK[ind], varCH3COCHO$mce5[ind])
    
    
    # Methylglyoxal
    ind = which(allBOTH.filter$variable == 'CH3COCHO_ACES_WOMACK' & is.finite(allBOTH.filter$ChosenEF.1hz))
    varCH3COCHO = as.data.frame(allBOTH.filter$ChosenEF.1hz[ind])
    colnames(varCH3COCHO) = 'CH3COCHO_ACES_WOMACK' 
    varCH3COCHO$mce = allBOTH.filter$mce.1hz[ind]
    varCH3COCHO$EF1.1hz = allBOTH.filter$EF1.1hz[ind]
    varCH3COCHO$EF1CO.1hz = allBOTH.filter$EF1CO.1hz[ind]
    varCH3COCHO$mce5 = allBOTH.filter$mce.5hz[ind]
    goodpasses2 = allBOTH.filter$uniqueid[ind]
    tmp = which(allBOTH.filter$uniqueid %in% goodpasses2)
    test = allBOTH.filter[tmp,]
    ind = which(test$variable == 'CO_DACOM_DISKIN')
    varCH3COCHO$maxCO = test$maxval.5hz[ind]
    varCH3COCHO$fuel = test$fuel[ind]
    
    CH3COCHOfire=ggplot(varCH3COCHO) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce5, y=EF1.1hz, col=fuel,size=maxCO))+
      # geom_point(aes(x=mce5, y=CH3COCHO_ACES_WOMACK, col='red',size=maxCO))+
      ylab("Methylglyoxal EF, g/kg") + xlab("MCE")+
      theme(legend.background=element_blank())+theme(legend.position = c(0.18, 0.7))+
      theme(text = element_text(size = sz))+
      scale_shape_manual(values=c(19,15,17,18,7,8,4))+
      scale_color_manual(values = c(cbp1 ),name="")+xlim(c(0.8,1))+
      geom_point(data=andreae, aes(x=andreae$average[1], y=andreae$average[53]), col='purple', size=5, shape=3)
    
    cor.test(varCH3COCHO$CH3COCHO_ACES_WOMACK, varCH3COCHO$mce5)
    ind = which(varCH3COCHO$fuel == 'corn')
    cor.test(varCH3COCHO$CH3COCHO_ACES_WOMACK[ind], varCH3COCHO$mce5[ind])
    
    Fig2 = ggarrange(CH3OHfire,CH2Ofire,common.legend = TRUE)
    ggsave(Fig2, file='Fig2.ps',width = 7*1.25*2, height=7*1.25*2/3)
    # methylglyoxal outliers
    ind = which(allBOTH.filter$variable == 'CH3COCHO_ACES_WOMACK' & is.finite(allBOTH.filter$EF1.1hz))
    tmpEF = allBOTH.filter[ind,]
    goodpasses = tmpEF$uniqueid
    ggplot(allBOTH.filter[ind,]) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce.5hz, y=EF1.1hz, col=fuel),size=5)+#xlim()
      ylab("Methylglyoxal EF, g/kg") + xlab("MCE")+
      theme(legend.background=element_blank())+theme(legend.position = c(0.8, 0.93))+
      theme(text = element_text(size = sz)) 
    cor.test(allBOTH.filter$EF1.1hz[ind], allBOTH.filter$mce.5hz[ind])
    ind = which(allBOTH.filter$variable == 'CH3COCHO_ACES_WOMACK'& allBOTH.filter$fuel == 'corn')
    cor.test(allBOTH.filter$EF1.1hz[ind], allBOTH.filter$mce.5hz[ind])
    
    ind = which(allBOTH.filter$variable == 'OC_JIMENEZ')
    ggplot(allBOTH.filter[ind,]) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce.1hz, y=EF1.1hz, col=fuel),size=5)+
      ylab("OC EF, g/kg") + xlab("MCE")+
      theme(legend.background=element_blank())+theme(legend.position = c(0.2, 0.3))+
      theme(text = element_text(size = sz)) 
    cor.test(allBOTH.filter$EF1.1hz[ind], allBOTH.filter$mce.1hz[ind])
    ind = which(allBOTH.filter$variable == 'OC_JIMENEZ' & allBOTH.filter$fuel == 'corn')
    cor.test(allBOTH.filter$EF1.1hz[ind], allBOTH.filter$mce.1hz[ind])
    
    ind = which(allBOTH.filter$variable == 'CH3OH_TOGA_APEL' & allBOTH.filter$EF1.1hz <= 230)# &
    #              allBOTH.filter$fuel == 'corn')
    ggplot(allBOTH.filter[ind,]) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce.5hz, y=EF1.1hz, col=fuel),size=5)+
      ylab("Acetaldehyde EF, g/kg") + xlab("MCE")+
      theme(legend.background=element_blank())+theme(legend.position = c(0.2, 0.3))+
      theme(text = element_text(size = sz)) 
    cor.test(allBOTH.filter$EF1.1hz[ind], allBOTH.filter$mce.5hz[ind])
    
    ind = which(allBOTH.filter$variable == 'BC_SCHWARZ')
    ggplot(allBOTH.filter[ind,]) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce.1hz, y=EF1.1hz, col=fuel),size=5)+
      ylab("BC EF, g/kg") + xlab("MCE")+
      theme(legend.background=element_blank())+theme(legend.position = c(0.2, 0.3))+
      theme(text = element_text(size = sz)) 
    cor.test(allBOTH.filter$EF1.1hz[ind], allBOTH.filter$mce.1hz[ind])
    ind = which(allBOTH.filter$variable == 'BC_SCHWARZ' & allBOTH.filter$fuel == 'corn')
    cor.test(allBOTH.filter$EF1.1hz[ind], allBOTH.filter$mce.1hz[ind])
    
    ind = which(allBOTH.filter$variable == 'All5HzVOC')
    ggplot(allBOTH.filter[ind,]) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce.5hz, y=EF1.5hz, col=fuel),size=5)+
      ylab("VOC, g/kg") + xlab("MCE")+
      theme(legend.background=element_blank())+theme(legend.position = c(0.2, 0.3))+
      theme(text = element_text(size = sz)) 
    cor.test(allBOTH.filter$EF1.1hz[ind], allBOTH.filter$mce.1hz[ind])
    ind = which(allBOTH.filter$variable == 'BC_SCHWARZ' & allBOTH.filter$fuel == 'corn')
    cor.test(allBOTH.filter$EF1.1hz[ind], allBOTH.filter$mce.1hz[ind])
    
    # ------------- Methanol with CO size -------
    ind = which(allBOTH.filter$variable == 'CH3OH_NOAAPTR_ppbv_WARNEKE')# & allBOTH.filter$ChosenEF.R2.5hz < R2filter)# allBOTH.filter$EF1.5hz < 7)
    toinvestigate = as.data.frame(cbind(fire=allBOTH.filter$fire[ind], 
                                        pass=allBOTH.filter$pass[ind], R2=allBOTH.filter$R2toX.5hz[ind],
                                        R2CO=allBOTH.filter$R2toCO.5hz[ind],
                                        fuel=allBOTH.filter$fuel[ind],
                                        mce = allBOTH.filter$mce.5hz[ind]))
    varCH3OH = as.data.frame(allBOTH.filter$ChosenEF.5hz[ind] )
    colnames(varCH3OH) = 'CH3OH_NOAAPTR_ppbv_WARNEKE'
    varCH3OH$EF1.5hz = allBOTH.filter$EF1.5hz[ind]
    varCH3OH$mce5 = allBOTH.filter$mce.5hz[ind]
    varCH3OH$mce = allBOTH.filter$mce.1hz[ind]
    goodpasses2 = allBOTH.filter$uniqueid[ind]
    tmp = which(allBOTH.filter$uniqueid %in% goodpasses2)
    test = allBOTH.filter[tmp,]
    ind = which(test$variable == 'CO_DACOM_DISKIN')
    varCH3OH$maxCO = test$maxval.5hz[ind]
    varCH3OH$fuel = test$fuel[ind]
    indA = which(akagi$Species == 'Methanol')
    CH3OHfire = ggplot(varCH3OH) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce5, y=CH3OH_NOAAPTR_ppbv_WARNEKE, col=fuel,size=maxCO))+
      ylab("Methanol EF, g/kg") + xlab("MCE")+
      theme(legend.background=element_blank())+theme(legend.position = c(0.18, 0.7))+
      theme(text = element_text(size = sz)) +
      geom_point(data=xiaoxi, aes(x=MCE, y=Methanol), pch=0,col='green')+
      scale_shape_manual(values=c(19,15,17,18,7,8,4))+
      scale_color_manual(values = c(cbp1 ))+xlim(c(0.8,1))+
      geom_point(data=andreae, aes(x=andreae$average[1], y=andreae$average[42]), col='purple', size=5, shape=3)+
      geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[40]), col='pink', size=5, shape=3)+
      geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[40]), col='pink', size=5, shape=4)
    
    cor.test(varCH3OH$CH3OH_NOAAPTR_ppbv_WARNEKE, varCH3OH$mce5)
    ind = which(varCH3OH$fuel == 'corn')
    cor.test(varCH3OH$CH3OH_NOAAPTR_ppbv_WARNEKE[ind], varCH3OH$mce5[ind])
    
    # investigate 1hz MCE outliers
    diff=(allBOTH.filter$mce.1hz - allBOTH.filter$mce.5hz)/allBOTH.filter$mce.5hz
    ind = which(abs(diff) > 0.05)
    iqq = unique(allBOTH.filter$uniqueid[ind])
    ind = which(allBOTH.filter$uniqueid %in% iqq)
    toinvestigate = allBOTH.filter[ind,]
    ind = which(toinvestigate$variable == 'CO_DACOM_DISKIN')
    toinvestigate = toinvestigate[ind,]
    
    
    # ------ OC fire ------
    OCfire1 = ggplot(allBOTH.filter[ind,]) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce.1hz, y=EF1.1hz*2, col=fuel),size=5)+
      ylab("OA EF, g/kg") + xlab("MCE") + xlim(c(0.8, 1))+
      theme(legend.background=element_blank())+theme(legend.position = c(0.18, 0.7))+
      theme(text = element_text(size = sz))  + geom_point(data=xiaoxi, aes(x=mce, y=oa), pch=0,col='green')
    ind = which(allBOTH.filter$variable == 'OC_JIMENEZ' & is.finite(allBOTH.filter$EF1.1hz))
    
    # ------------- Furfural with CO size -------
    ind = which(allBOTH.filter$variable == 'Furfural_NOAAPTR_ppbv_WARNEKE')# & allBOTH.filter$ChosenEF.R2.5hz < R2filter)# allBOTH.filter$EF1.5hz < 7)
    toinvestigate = as.data.frame(cbind(fire=allBOTH.filter$fire[ind], 
                                        pass=allBOTH.filter$pass[ind], R2=allBOTH.filter$R2toX.5hz[ind],
                                        R2CO=allBOTH.filter$R2toCO.5hz[ind],
                                        fuel=allBOTH.filter$fuel[ind],
                                        mce = allBOTH.filter$mce.5hz[ind]))
    varFURFURAL = as.data.frame(allBOTH.filter$ChosenEF.5hz[ind] )
    colnames(varFURFURAL) = 'FURFURAL_NOAAPTR_ppbv_WARNEKE'
    varFURFURAL$EF1.5hz = allBOTH.filter$EF1.5hz[ind]
    varFURFURAL$mce5 = allBOTH.filter$mce.5hz[ind]
    varFURFURAL$mce = allBOTH.filter$mce.1hz[ind]
    goodpasses2 = allBOTH.filter$uniqueid[ind]
    tmp = which(allBOTH.filter$uniqueid %in% goodpasses2)
    test = allBOTH.filter[tmp,]
    ind = which(test$variable == 'CO_DACOM_DISKIN')
    varFURFURAL$maxCO = test$maxval.5hz[ind]
    varFURFURAL$fuel = test$fuel[ind]
    indA = which(akagi$Species == 'Methanol')
    FURFURALfire = ggplot(varFURFURAL) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce5, y=FURFURAL_NOAAPTR_ppbv_WARNEKE, col=fuel,size=maxCO))+
      ylab("Methanol EF, g/kg") + xlab("MCE")+
      theme(legend.background=element_blank())+theme(legend.position = c(0.18, 0.7))+
      theme(text = element_text(size = sz)) +
      geom_point(data=xiaoxi, aes(x=MCE, y=Methanol), pch=0,col='green')+
      scale_shape_manual(values=c(19,15,17,18,7,8,4))+
      scale_color_manual(values = c(cbp1 ))+xlim(c(0.8,1))+
      geom_point(data=andreae, aes(x=andreae$average[1], y=andreae$average[42]), col='purple', size=5, shape=3)+
      geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[40]), col='pink', size=5, shape=3)+
      geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[40]), col='pink', size=5, shape=4)
    
    cor.test(varFURFURAL$FURFURAL_NOAAPTR_ppbv_WARNEKE, varFURFURAL$mce5)
    ind = which(varFURFURAL$fuel == 'corn')
    cor.test(varFURFURAL$FURFURAL_NOAAPTR_ppbv_WARNEKE[ind], varFURFURAL$mce5[ind])
    
    # ----------------- NO with MCE -------
    ind = which(allBOTH.filter$variable == 'NO_LIF_ROLLINS')
    varNO = as.data.frame(allBOTH.filter$ChosenEF.5hz[ind])
    colnames(varNO) = 'NO_LIF_ROLLINS' 
    varNO$EF1.5hz = allBOTH.filter$EF1.5hz[ind]
    varNO$EF1.1hz = allBOTH.filter$EF1.1hz[ind]
    varNO$EF1CO.5hz = allBOTH.filter$EF1CO.5hz[ind]
    
    varNO$mce5 = allBOTH.filter$mce.5hz[ind]
    varNO$mce = allBOTH.filter$mce.1hz[ind]
    goodpasses2 = allBOTH.filter$uniqueid[ind]
    tmp = which(allBOTH.filter$uniqueid %in% goodpasses2)
    test = allBOTH.filter[tmp,]
    ind = which(test$variable == 'CO_DACOM_DISKIN')
    varNO$maxCO = test$maxval.5hz[ind]
    varNO$fuel = test$fuel[ind]
    ind = which(test$variable == 'NO2_ACES_WOMACK')
    varNO$NO2_ACES_WOMACK = test$EF1.1hz[ind]
    ind = which(test$variable == 'NO_RYERSON')
    varNO$NO_RYERSON = test$EF1.1hz[ind]
    ind = which(test$variable == 'NO_RYERSON')
    varNO$NO_RYERSON = test$EF1.1hz[ind]
    ind = which(test$variable == 'NO2_RYERSON')
    varNO$NO2_RYERSON = test$EF1.1hz[ind]
    varNO$NOx = varNO$EF1.5hz + varNO$NO2_ACES_WOMACK
    
    NOfire = ggplot(varNO) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce5, y=NOx, col=fuel,size=maxCO))+
      ylab("NOx EF, g/kg") + xlab("MCE")+
      theme(legend.background=element_blank())+theme(legend.position = c(0.18, 0.7))+
      theme(text = element_text(size = sz)) +
      geom_point(data=xiaoxi, aes(x=MCE, y=NO+NO2), pch=0,col='green')+
      scale_shape_manual(values=c(19,15,17,18,7,8,4))+
      scale_color_manual(values = c(cbp1 ))+xlim(c(0.8,1))+
      geom_point(data=andreae, aes(x=andreae$average[1], y=andreae$average[97]), col='purple', size=5, shape=3)+
      geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[100]), col='pink', size=5, shape=3)+
      geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[100]), col='pink', size=5, shape=4)
    
    cor.test(varNO$NO_LIF_ROLLINS, varNO$mce5)
    ind = which(varNO$fuel == 'corn' & varNO$mce5 > 0.95)
    cor.test(varNO$NO_LIF_ROLLINS[ind], varNO$mce5[ind])
    # ------------- CH4 with CO size -------
    ind = which(allBOTH.filter$variable == 'CH4_DACOM_DISKIN')
    
    varCH4 = as.data.frame(allBOTH.filter$ChosenEF.5hz[ind] ) 
    colnames(varCH4) = 'CH4_DACOM_DISKIN' 
    varCH4$EF1.5hz = allBOTH.filter$EF1.5hz[ind]
    varCH4$EF1CO.5hz = allBOTH.filter$EF1CO.5hz[ind]
    varCH4$EF1intfill.5hz = allBOTH.filter$EF1intfill.5hz[ind]
    varCH4$EF1COintfill.5hz = allBOTH.filter$EF1COintfill.5hz[ind]
    varCH4$mce5 = allBOTH.filter$mce.5hz[ind]
    varCH4$mce = allBOTH.filter$mce.5hz[ind]
    goodpasses2 = allBOTH.filter$uniqueid[ind]
    tmp = which(allBOTH.filter$uniqueid %in% goodpasses2)
    test = allBOTH.filter[tmp,]
    ind = which(test$variable == 'CO_DACOM_DISKIN')
    varCH4$maxCO = test$maxval.5hz[ind]
    varCH4$fuel = test$fuel[ind]
    indA = which(akagi$Species == 'Methane')
    CH4fire = ggplot(varCH4) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce5, y=EF1CO.5hz, col=fuel,size=maxCO))+
      ylab("CH4 EF, g/kg") + xlab("MCE")+
      theme(legend.background=element_blank())+theme(legend.position = c(0.18, 0.7))+
      theme(text = element_text(size = sz)) +
      #geom_point(data=xiaoxi, aes(x=MCE, y=OA/2), pch=0,col='green')+
      scale_shape_manual(values=c(19,15,17,18,7,8,4))+
      scale_color_manual(values = c(cbp1 ))+xlim(c(0.8,1))  +
      geom_point(data=andreae, aes(x=andreae$average[1], y=andreae$average[4]), col='purple', size=5, shape=3)+
      geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[4]), col='pink', size=5, shape=3)+
      geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[4]), col='pink', size=5, shape=4)
    
    cor.test(varCH4$CH4_JIMENEZ, varCH4$mce5)
    ind = which(varCH4$fuel == 'corn')
    cor.test(varCH4$CH4_JIMENEZ[ind], varCH4$mce5[ind])
    
    # ------------- OC with CO size -------
    ind = which(allBOTH.filter$variable == 'OC_JIMENEZ' & is.finite(allBOTH.filter$EF1.1hz))
    
    varOC = as.data.frame(allBOTH.filter$ChosenEF.1hz[ind] ) # roughly the OA/OC
    colnames(varOC) = 'OC_JIMENEZ' 
    varOC$EF1.1hz = allBOTH.filter$EF1.1hz[ind]
    varOC$EF1CO.1hz = allBOTH.filter$EF1CO.1hz[ind]
    varOC$EF1intfill.1hz = allBOTH.filter$EF1intfill.1hz[ind]
    varOC$EF1COintfill.1hz = allBOTH.filter$EF1COintfill.1hz[ind]
    varOC$mce5 = allBOTH.filter$mce.5hz[ind]
    varOC$mce = allBOTH.filter$mce.1hz[ind]
    goodpasses2 = allBOTH.filter$uniqueid[ind]
    tmp = which(allBOTH.filter$uniqueid %in% goodpasses2)
    test = allBOTH.filter[tmp,]
    ind = which(test$variable == 'CO_DACOM_DISKIN')
    varOC$maxCO = test$maxval.5hz[ind]
    varOC$fuel = test$fuel[ind]
    indA = which(akagi$Species == 'Organic Carbon')
    OCfire = ggplot(varOC) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce5, y=OC_JIMENEZ, col=fuel,size=maxCO))+
      ylab("OC EF, g/kg") + xlab("MCE")+
      theme(legend.background=element_blank())+theme(legend.position = c(0.18, 0.7))+
      theme(text = element_text(size = sz)) +
      geom_point(data=xiaoxi, aes(x=MCE, y=OA/2), pch=0,col='green')+
      scale_shape_manual(values=c(19,15,17,18,7,8,4))+
      scale_color_manual(values = c(cbp1 ))+xlim(c(0.8,1))  +
      geom_point(data=andreae, aes(x=andreae$average[1], y=andreae$average[115]), col='purple', size=5, shape=3)+
      geom_point(data=akagi, aes(x=akagi$CropEF[1], y=akagi$CropEF[110]), col='pink', size=5, shape=3)+
      geom_point(data=akagi, aes(x=akagi$PastureEF[1], y=akagi$PastureEF[110]), col='pink', size=5, shape=4)
    
    cor.test(varOC$OC_JIMENEZ, varOC$mce5)
    ind = which(varOC$fuel == 'corn')
    cor.test(varOC$OC_JIMENEZ[ind], varOC$mce5[ind])
    
    Fig3 = ggarrange(OCfire,NOfire,common.legend = TRUE)
    ggsave(Fig3, file='Fig3.ps',width = 7*1.25*2, height=7*1.25*2/3)
    
    ind = which(allBOTH.filter$variable == 'All5HzVOC'& is.finite(allBOTH.filter$EF1.5hz))
    ggplot(allBOTH.filter[ind,]) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce.5hz, y=EF1.5hz, col=fuel),size=5)+
      ylab("VOC EF, g/kg") + xlab("MCE") +
      theme(legend.background=element_blank())+theme(legend.position = c(0.18, 0.7))+
      theme(text = element_text(size = sz)) 
    
    ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' )
    ggplot(allBOTH.filter[ind,]) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      #geom_point(aes(x=mce.5hz, y=EF1.5hz, col=fuel),size=5)+
      geom_point(aes(x=mce_int.1hz, y=EF1int.1hz, col=fuel),size=5)+
      
      ylab("CO EF, g/kg") + xlab("MCE") +
      theme(legend.background=element_blank())+theme(legend.position = c(0.18, 0.7))+
      theme(text = element_text(size = sz)) +
      geom_point(data=xiaoxi, aes(x=mce, y=co), col='green')
    
    ind = which(allBOTH.filter$variable == 'CH4_DACOM_DISKIN' )
    ggplot(allBOTH.filter[ind,]) + 
      theme_classic()+#geom_label(aes(x=mce.1hz, y=EF1.1hz, label=fire),nudge_x = .005,nudge_y = .05)+
      geom_point(aes(x=mce.5hz, y=EF1.5hz, col=fuel),size=5)+
      ylab("CH4 EF, g/kg") + xlab("MCE") +
      theme(legend.background=element_blank())+theme(legend.position = c(0.18, 0.27))+
      theme(text = element_text(size = sz)) 
    ind = which(allBOTH.filter$variable == 'CH4_DACOM_DISKIN' & allBOTH.filter$fuel == 'corn')
    cor.test(allBOTH.filter$EF1.5hz[ind], allBOTH.filter$mce.5hz[ind])
              
    ind = which(Ricearoni.1hz.EF$variable == 'CH2O_ppt')
    
              
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
      
      # ---------- OC 
      ind = which(allBOTH.filter.allfuels$variable == 'OC_JIMENEZ' & allBOTH.filter.allfuels$fire != 'BlackwaterRiver' )
      averagebyfuel = aggregate(allBOTH.filter.allfuels[ind,], by=list(allBOTH.filter.allfuels$fuel[ind]), FUN='mean', na.rm=TRUE)
      averagebyfuelsd = aggregate(allBOTH.filter.allfuels[ind,], by=list(allBOTH.filter.allfuels$fuel[ind]), FUN='sd', na.rm=TRUE)
      ind = which(allBOTH.blackwater.filter$variable == 'OC_JIMENEZ' )
      averagebyfuelB = aggregate(allBOTH.blackwater.filter[ind,], by=list(allBOTH.blackwater.filter$fuel[ind]), FUN='mean', na.rm=TRUE)
      averagebyfuelBsd = aggregate(allBOTH.blackwater.filter[ind,], by=list(allBOTH.blackwater.filter$fuel[ind]), FUN='sd', na.rm=TRUE)
      
      iA = which(akagi$Species == 'OrganicCarbon')
      fuel = c(averagebyfuel$Group.1, averagebyfuelB$Group.1,'Akagi: savannah','Akagi: pasture maintenance','Akagi: crop residue')
      EF = c(averagebyfuel$FinalEF, averagebyfuelB$FinalEF,akagi$Savannah[iA],akagi$PastureEF[iA],akagi$CropEF[iA])
      mce = c(averagebyfuel$mce.5hz, averagebyfuelB$mce.5hz, 0.94,0.88,0.91)
      EFsd = c(averagebyfuelsd$FinalEF, averagebyfuelBsd$FinalEF,akagi$`Savannah SD`[iA],akagi$PastureSD[iA],akagi$CropSD[iA])
      mcesd = c(averagebyfuelsd$mce.5hz, averagebyfuelBsd$mce.5hz,NaN,NaN,NaN)
      averagebyfuel2 = as.data.frame(cbind(fuel,EF,mce, EFsd, mcesd))
      averagebyfuel2$EF = as.numeric(averagebyfuel2$EF);averagebyfuel2$EFsd = as.numeric(averagebyfuel2$EFsd)
      averagebyfuel2$mce = as.numeric(averagebyfuel2$mce);averagebyfuel2$mcesd = as.numeric(averagebyfuel2$mcesd)
      ind = which(is.finite(averagebyfuel2$EF))
      ggplot(averagebyfuel2[ind,]) + geom_point(aes(x=mce, y=EF, col=fuel), size=5)+
        theme_classic()+xlab("MCE")+ylab("EF OC, g/kg")+
        geom_errorbar(data=averagebyfuel2[ind,],aes(xmin=mce-mcesd,xmax=mce+mcesd, y=EF,col=fuel),  position=position_dodge(0.05))+
        geom_errorbar(data=averagebyfuel2[ind,],aes(x=mce,ymin=EF-EFsd, ymax=EF+EFsd,col=fuel), position=position_dodge(0.05))+
        theme(legend.background=element_blank())+ theme(text = element_text(size = 20)) +
        scale_colour_hue()
      
      # ---------- C2H6 
      ind = which(allBOTH.filter.allfuels$variable == 'C2H6_CAMS_pptv_FRIED' & allBOTH.filter.allfuels$fire != 'BlackwaterRiver' )
      averagebyfuel = aggregate(allBOTH.filter.allfuels[ind,], by=list(allBOTH.filter.allfuels$fuel[ind]), FUN='mean', na.rm=TRUE)
      averagebyfuelsd = aggregate(allBOTH.filter.allfuels[ind,], by=list(allBOTH.filter.allfuels$fuel[ind]), FUN='sd', na.rm=TRUE)
      ind = which(allBOTH.blackwater.filter$variable == 'C2H6_CAMS_pptv_FRIED' )
      averagebyfuelB = aggregate(allBOTH.blackwater.filter[ind,], by=list(allBOTH.blackwater.filter$fuel[ind]), FUN='mean', na.rm=TRUE)
      averagebyfuelBsd = aggregate(allBOTH.blackwater.filter[ind,], by=list(allBOTH.blackwater.filter$fuel[ind]), FUN='sd', na.rm=TRUE)
      
      iA = which(akagi$Species == 'Ethane')
      fuel = c(averagebyfuel$Group.1, averagebyfuelB$Group.1,'Akagi: savannah','Akagi: pasture maintenance','Akagi: crop residue')
      EF = c(averagebyfuel$FinalEF, averagebyfuelB$FinalEF,akagi$Savannah[iA],akagi$PastureEF[iA],akagi$CropEF[iA])
      mce = c(averagebyfuel$mce.5hz, averagebyfuelB$mce.5hz, 0.94,0.88,0.91)
      EFsd = c(averagebyfuelsd$FinalEF, averagebyfuelBsd$FinalEF,akagi$`Savannah SD`[iA],akagi$PastureSD[iA],akagi$CropSD[iA])
      mcesd = c(averagebyfuelsd$mce.5hz, averagebyfuelBsd$mce.5hz,NaN,NaN,NaN)
      averagebyfuel2 = as.data.frame(cbind(fuel,EF,mce, EFsd, mcesd))
      averagebyfuel2$EF = as.numeric(averagebyfuel2$EF);averagebyfuel2$EFsd = as.numeric(averagebyfuel2$EFsd)
      averagebyfuel2$mce = as.numeric(averagebyfuel2$mce);averagebyfuel2$mcesd = as.numeric(averagebyfuel2$mcesd)
      ind = which(is.finite(averagebyfuel2$EF))
      ggplot(averagebyfuel2[ind,]) + geom_point(aes(x=mce, y=EF, col=fuel), size=5)+
        theme_classic()+xlab("MCE")+ylab("EF C2H6, g/kg")+
        geom_errorbar(data=averagebyfuel2[ind,],aes(xmin=mce-mcesd,xmax=mce+mcesd, y=EF,col=fuel),  position=position_dodge(0.05))+
        geom_errorbar(data=averagebyfuel2[ind,],aes(x=mce,ymin=EF-EFsd, ymax=EF+EFsd,col=fuel), position=position_dodge(0.05))+
        theme(legend.background=element_blank())+ theme(text = element_text(size = 20)) +
        scale_colour_hue()
      
      
      # ---------- CO
      ind = which(allBOTH.filter.allfuels$variable == 'CO_DACOM_DISKIN' & allBOTH.filter.allfuels$fire != 'BlackwaterRiver' )
      averagebyfuel = aggregate(allBOTH.filter.allfuels[ind,], by=list(allBOTH.filter.allfuels$fuel[ind]), FUN='mean', na.rm=TRUE)
      averagebyfuelsd = aggregate(allBOTH.filter.allfuels[ind,], by=list(allBOTH.filter.allfuels$fuel[ind]), FUN='sd', na.rm=TRUE)
      ind = which(allBOTH.blackwater.filter$variable == 'CO_DACOM_DISKIN' )
      averagebyfuelB = aggregate(allBOTH.blackwater.filter[ind,], by=list(allBOTH.blackwater.filter$fuel[ind]), FUN='mean', na.rm=TRUE)
      averagebyfuelBsd = aggregate(allBOTH.blackwater.filter[ind,], by=list(allBOTH.blackwater.filter$fuel[ind]), FUN='sd', na.rm=TRUE)
      
      iA = which(akagi$Species == 'Carbon Monoxide')
      fuel = c(averagebyfuel$Group.1, averagebyfuelB$Group.1,'Akagi: savannah','Akagi: pasture maintenance','Akagi: crop residue','Akagi: temperature')
      EF = c(averagebyfuel$FinalEF, averagebyfuelB$FinalEF,akagi$Savannah[iA],akagi$PastureEF[iA],akagi$CropEF[iA],akagi$TemperateEF[iA])
      mce = c(averagebyfuel$mce.5hz, averagebyfuelB$mce.5hz, 0.94,0.88,0.91,0.95)
      EFsd = c(averagebyfuelsd$FinalEF, averagebyfuelBsd$FinalEF,akagi$`Savannah SD`[iA],akagi$PastureSD[iA],akagi$CropSD[iA],akagi$`Temperate SD`[iA])
      mcesd = c(averagebyfuelsd$mce.5hz, averagebyfuelBsd$mce.5hz,NaN,NaN,NaN,NaN)
      averagebyfuel2 = as.data.frame(cbind(fuel,EF,mce, EFsd, mcesd))
      averagebyfuel2$EF = as.numeric(averagebyfuel2$EF);averagebyfuel2$EFsd = as.numeric(averagebyfuel2$EFsd)
      averagebyfuel2$mce = as.numeric(averagebyfuel2$mce);averagebyfuel2$mcesd = as.numeric(averagebyfuel2$mcesd)
      ind = which(is.finite(averagebyfuel2$EF))
      ggplot(averagebyfuel2[ind,]) + geom_point(aes(x=mce, y=EF, col=fuel), size=5)+
        theme_classic()+xlab("MCE")+ylab("EF CO, g/kg")+
        geom_errorbar(data=averagebyfuel2[ind,],aes(xmin=mce-mcesd,xmax=mce+mcesd, y=EF,col=fuel),  position=position_dodge(0.05))+
        geom_errorbar(data=averagebyfuel2[ind,],aes(x=mce,ymin=EF-EFsd, ymax=EF+EFsd,col=fuel), position=position_dodge(0.05))+
        theme(legend.background=element_blank())+ theme(text = element_text(size = 20)) +
        scale_colour_hue()
      
      # ---------- Benzene 
      ind = which(allBOTH.filter.allfuels$variable == 'Benzene_NOAAPTR_ppbv_WARNEKE' & allBOTH.filter.allfuels$fire != 'BlackwaterRiver' )
      averagebyfuel = aggregate(allBOTH.filter.allfuels[ind,], by=list(allBOTH.filter.allfuels$fuel[ind]), FUN='mean', na.rm=TRUE)
      averagebyfuelsd = aggregate(allBOTH.filter.allfuels[ind,], by=list(allBOTH.filter.allfuels$fuel[ind]), FUN='sd', na.rm=TRUE)
      ind = which(allBOTH.blackwater.filter$variable == 'Benzene_NOAAPTR_ppbv_WARNEKE' )
      averagebyfuelB = aggregate(allBOTH.blackwater.filter[ind,], by=list(allBOTH.blackwater.filter$fuel[ind]), FUN='mean', na.rm=TRUE)
      averagebyfuelBsd = aggregate(allBOTH.blackwater.filter[ind,], by=list(allBOTH.blackwater.filter$fuel[ind]), FUN='sd', na.rm=TRUE)
      
      iA = which(akagi$Species == 'Benzene')
      fuel = c(averagebyfuel$Group.1, averagebyfuelB$Group.1,'Akagi: savannah','Akagi: pasture maintenance','Akagi: crop residue')
      EF = c(averagebyfuel$FinalEF, averagebyfuelB$FinalEF,akagi$Savannah[iA],akagi$PastureEF[iA],akagi$CropEF[iA])
      mce = c(averagebyfuel$mce.5hz, averagebyfuelB$mce.5hz, 0.94,0.88,0.91)
      EFsd = c(averagebyfuelsd$FinalEF, averagebyfuelBsd$FinalEF,akagi$`Savannah SD`[iA],akagi$PastureSD[iA],akagi$CropSD[iA])
      mcesd = c(averagebyfuelsd$mce.5hz, averagebyfuelBsd$mce.5hz,NaN,NaN,NaN)
      averagebyfuel2 = as.data.frame(cbind(fuel,EF,mce, EFsd, mcesd))
      averagebyfuel2$EF = as.numeric(averagebyfuel2$EF);averagebyfuel2$EFsd = as.numeric(averagebyfuel2$EFsd)
      averagebyfuel2$mce = as.numeric(averagebyfuel2$mce);averagebyfuel2$mcesd = as.numeric(averagebyfuel2$mcesd)
      ind = which(is.finite(averagebyfuel2$EF))
      ggplot(averagebyfuel2[ind,]) + geom_point(aes(x=mce, y=EF, col=fuel), size=5)+
        theme_classic()+xlab("MCE")+ylab("EF Benzene, g/kg")+
        geom_errorbar(data=averagebyfuel2[ind,],aes(xmin=mce-mcesd,xmax=mce+mcesd, y=EF,col=fuel),  position=position_dodge(0.05))+
        geom_errorbar(data=averagebyfuel2[ind,],aes(x=mce,ymin=EF-EFsd, ymax=EF+EFsd,col=fuel), position=position_dodge(0.05))+
        theme(legend.background=element_blank())+ theme(text = element_text(size = 20)) +
        scale_colour_hue()
      
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

# TOGA/WAS meacetate
plot( toga.all$C3H6O2_NOAAPTR_TM,toga.all$MeAcetate_ppt, pch=19, xlab='Methyl acetate, ppt', ylab='C3H6O2 (NOAA PTRMS), ppt')
xx=toga.all$C3H6O2_NOAAPTR_TM
yy=as.numeric(toga.all$MeAcetate_ppt)
tt = lm(yy~xx+0)
abline(tt, col='black')
text(4000,6000,paste("TOGA slope = ", round(tt$coefficients[1], digits = 2)))

points(toga.all$C3H6O2_NOAAPTR_TM, toga.all$MeAcetate_WAS_TM, pch=19, col='red')
yy=as.numeric(toga.all$MeAcetate_WAS_TM)
tt2 = lm(yy~xx+0)
abline(tt2, col='red')
text(4000,5700,paste("WAS slope = ", round(tt2$coefficients[1], digits = 2)), col='red')

# TOGA pyrrole
plot( toga.all$C4H5N_NOAAPTR_TM, toga.all$Pyrrole_ppt, pch=19, xlab='C4H5N (NOAA PTRMS), ppt', ylab='Pyrrole, ppt')
xx=as.numeric(toga.all$C4H5N_NOAAPTR_TM)
yy=as.numeric(toga.all$Pyrrole_ppt)
ind = which(xx > 400 & yy < 100) # don't include outliers in slope
xx[ind]= NaN; yy[ind] = NaN
tt = lm(yy~xx+0)
abline(tt, col='black')
text(300,600,paste("TOGA slope = ", round(tt$coefficients[1], digits = 2)))

# Monoterpenes
tmp = as.data.frame(cbind(toga.all$aPinene_ppt,toga.all$bPineneMyrcene_ppt , toga.all$Camphene_ppt ,toga.all$Tricyclene_ppt))
ind = which(is.finite(toga.all$bPineneMyrcene_ppt) & is.finite(toga.all$Camphene_ppt)) # don't include outliers in slope
plot( toga.all$Monoterpenes_NOAAPTR_TM[ind],rowSums(tmp[ind,], na.rm=TRUE), pch=19, xlab='Monoterpenes, ppt', ylab='aPinene+bPinene+Myrcene+Camphene+Tricyclene, ppt')
xx=as.numeric(toga.all$Monoterpenes_NOAAPTR_TM)
yy=as.numeric(rowSums(tmp, na.rm=TRUE))

xx=xx[ind]; yy=yy[ind]
tt = lm(yy~xx+0)
abline(tt, col='black')
text(400,600,paste("TOGA slope = ", round(tt$coefficients[1], digits = 2)))

# C9 aromatics
tmp = as.data.frame(cbind(as.numeric(toga.all$x124rimeBenzene_WAS_TM), as.numeric(toga.all$x135rimeBenzene_WAS_TM), as.numeric(toga.all$iPropBenzene_WAS_TM),
                          as.numeric(toga.all$nPropBenzene_WAS_TM), as.numeric(toga.all$x2EthToluene_WAS_TM), as.numeric(toga.all$x3EthToluene_WAS_TM), 
                          as.numeric(toga.all$x4EthToluene_WAS_TM)))
ind = which(is.finite(as.numeric(toga.all$x124rimeBenzene_WAS_TM)) &is.finite(as.numeric(toga.all$x3EthToluene_WAS_TM)) &
        is.finite( as.numeric(toga.all$x4EthToluene_WAS_TM)) & is.finite(as.numeric(toga.all$nPropBenzene_WAS_TM)) & is.finite(as.numeric(toga.all$x2EthToluene_WAS_TM)))
plot( as.numeric(toga.all$C9Aromatics_NOAAPTR_TM[ind]),rowSums(tmp[ind,], na.rm=TRUE), pch=19, xlab='C9 aromatics, ppt', ylab='WAS C9 aromatics, ppt', xlim=c(0,800), ylim=c(0,800))
xx=as.numeric(toga.all$C9Aromatics_NOAAPTR_TM[ind])
yy=as.numeric(rowSums(tmp[ind,], na.rm=TRUE))
tt = lm(yy~xx+0)
abline(tt, col='black')
text(200,250,paste("WAS slope = ", round(tt$coefficients[1], digits = 2)))

# Isoprene ---
plot(toga.all$Isoprene_NOAAPTR_TM, toga.all$Isoprene_ppt, ylab='Isoprene, ppt', xlab='NOAA PTRMS Isoprene, ppt', pch=19, xlim=c(0,10E3), ylim=c(0,10E3))
tt = lm(toga.all$Isoprene_ppt ~toga.all$Isoprene_NOAAPTR_TM)
abline(tt, col='grey')
abline(0,1, lty=1)
points(toga.all$Isoprene_NOAAPTR_TM, toga.all$Isoprene_NOAAiWAS_TM, col='red', pch=19)
points(toga.all$Isoprene_NOAAPTR_TM, toga.all$Isoprene_WAS_TM, col='blue', pch=19)

