# Run Everything
print('Calculate Emission Factors') #---------- Calculate EFs ----------
source("CalculateEmissionFactorsv2.R")
print('Emissions Analysis') #---------- Analyze EFs ----------
source("FIREXEmissionsAnalysis.R")
print('Table 1') #---------- Table 1 ----------
source("Table1.R") # Tble of fire and plume counts
print("EF Table") #---------- Table 4 ----------
source("TestTableJustEFv2.R")
print("Calculate MCE dependence")#---------- MCE dependence ----------
source("MCEFig.R")
print('Slope differences') # ---------- Slope differences Where FinalEF_MCE92 gets calculated  --------- 
source("SlopeDifferences.R")
#print('Supplement Xioaxi')
#source("Supplement_Xiaoxi.R")
print('Unique species')
source("UniqueSpeciesv2.R")

print(' Run plots')
source("RunPlotsv2.R")
#print(' Andreae Supplement')
#source("Supplement_Andreae.R")

ind = which(allBOTH$variable == 'CO_DACOM_DISKIN'  & is.finite(allBOTH$EF1.5hz))
ind = which(allBOTH$variable == 'CO_DACOM_DISKIN' & allBOTH$fuel == 'corn' & 
              as.numeric(allBOTH$R2toX.5hz) >= R2filterCO & as.numeric(allBOTH$maxval.5hz) > COcutoff) 

ind = which(allBOTH.filter$variable =='Acrolein_NOAAPTR_ppbv_WARNEKE' & allBOTH.filter$fuel == 'corn'  )#is.finite(allBOTH$EF1.5hz))
tmpACROLEIN = allBOTH.filter[ind,]
ind = which(allBOTH.filter$variable =='CO_DACOM_DISKIN' & allBOTH.filter$fuel == 'corn'  )#is.finite(allBOTH$EF1.5hz))
tmpCO = allBOTH.filter[ind,]

# 99 EFs
# For corn fires, 112 that were analyzed
# 1 fire had a bad correlation : Kingpin pass 2

allBOTH.filter.mean = aggregate(allBOTH.filter, by=list(allBOTH.filter$variable), FUN='mean', na.rm=TRUE)

for (i in 1:length(allBOTH.filter.mean$kind)){
  ind = which(allBOTH.filter$variable == allBOTH.filter.mean$Group.1[i])
  allBOTH.filter.mean$kind[i] = allBOTH.filter$kind[ind[1]]
  allBOTH.filter.mean$formula[i] = allBOTH.filter$formula[ind[1]]
  allBOTH.filter.mean$names[i] = allBOTH.filter$names[ind[1]]
  allBOTH.filter.mean$PI[i] = allBOTH.filter$PI[ind[1]]
  # allBOTH.filter.mean.fuelmajor$fuel[i] = allBOTH.filter$fuel[ind[1]]
  # print(c(allBOTH.filter$variable[ind[1]], allBOTH.filter.mean.fuelmajor$Group.2[i]))
}

propertiesSpecies = as.data.frame(cbind(allBOTH.filter.mean$Group.1,
                                        allBOTH.filter.mean$names,allBOTH.filter.mean$PI,
                                        allBOTH.filter.mean$mWs,allBOTH.filter.mean$OHrate.1hz, 
                                        allBOTH.filter.mean$OHrate.5hz, allBOTH.filter.mean$lifetime_jval,
                                        allBOTH.filter.mean$lifetime, allBOTH.filter.mean$LifetimeCat))
# ------ Fix instrument names -----------

ind = which(propertiesSpecies$V3 == 'BLAKE')
propertiesSpecies$V3[ind] = 'WAS'
propertiesSpecies$V3[3] = 'NIR spect.'
ind = which(propertiesSpecies$V3 == 'DISKIN')
propertiesSpecies$V3[ind] = 'DACOM'
ind = which(propertiesSpecies$V3 == 'JIMENEZ')
propertiesSpecies$V3[ind] = 'AMS'
ind = which(propertiesSpecies$V3 == 'ppt')
propertiesSpecies$V3[ind] = 'TOGA'
ind = which(propertiesSpecies$V3 == 'ppt+BLAKE+GILMAN')
propertiesSpecies$V3[ind] = 'TOGA, iWAS, WAS'
ind = which(propertiesSpecies$V3 == 'GILMAN+ppt')
propertiesSpecies$V3[ind] = 'TOGA, iWAS'
ind = which(propertiesSpecies$V3 == 'GILMAN+BLAKE')
propertiesSpecies$V3[ind] = 'iWAS, WAS'
ind = which(propertiesSpecies$V3 == 'BLAKE+GILMAN')
propertiesSpecies$V3[ind] = 'iWAS, WAS'
ind = which(propertiesSpecies$V3 == 'APEL+BLAKE')
propertiesSpecies$V3[ind] = 'TOGA, WAS'
ind = which(propertiesSpecies$V3 == 'BLAKE+ppt')
propertiesSpecies$V3[ind] = 'TOGA, WAS'
ind = which(propertiesSpecies$V3 == 'ppt+BLAKE')
propertiesSpecies$V3[ind] = 'TOGA, WAS'
ind = which(propertiesSpecies$V3 == 'ppt+BLAKE+WARNEKE')
propertiesSpecies$V3[ind] = 'NOAA PTR-ToF-MS, TOGA, WAS'
ind = which(propertiesSpecies$V3 == 'VERES')
propertiesSpecies$V3[ind] = 'NOAA CIMS'
ind = which(propertiesSpecies$V3 == 'WISTHALER')
propertiesSpecies$V3[ind] = 'OSLO PTRMS'
ind = which(propertiesSpecies$V3 == 'FRIED')
propertiesSpecies$V3[ind] = 'CAMS'
ind = which(propertiesSpecies$V3 == 'GILMAN')
propertiesSpecies$V3[ind] = 'iWAS'
ind = which(propertiesSpecies$V3 == 'FRIED+HANISCO')
propertiesSpecies$V3[ind] = 'CAMS, ISAF'
ind = which(propertiesSpecies$V3 == 'GILMAN+FRIED')
propertiesSpecies$V3[ind] = 'iWAS, CAMS'
ind = which(propertiesSpecies$V3 == 'VERES+WARKEKE')
propertiesSpecies$V3[ind] = 'NOAA CIMS, NOAA PTR-ToF-MS'
ind = which(propertiesSpecies$V3 == 'WARNEKE' | propertiesSpecies$V3 == 'ppbv')
propertiesSpecies$V3[ind] = 'NOAA PTR-ToF-MS'
ind = which(propertiesSpecies$V3 == 'WOMACK')
propertiesSpecies$V3[ind] = 'ACES'
ind = which(propertiesSpecies$V3 == 'WOMACK+RYERSON')
propertiesSpecies$V3[ind] = 'ACES, NOAA NOyO3'
ind = which(propertiesSpecies$V3 == 'WOMACK+VERES')
propertiesSpecies$V3[ind] = 'ACES, NOAA CIMS'
ind = which(propertiesSpecies$V3 == 'Warneke w/Blake+GILMAN+APEL')
propertiesSpecies$V3[ind] = 'NOAA PTR-ToF-MS, TOGA+WAS+iWAS spec.'
ind = which(propertiesSpecies$V3 == 'Warneke w/APEL')
propertiesSpecies$V3[ind] = 'NOAA PTR-ToF-MS, TOGA spec.'
ind = which(propertiesSpecies$V3 == 'WENNBERG')
propertiesSpecies$V3[ind] = 'CIT-CIMS'
ind = which(propertiesSpecies$V3 == 'WARNEKE+ppt+GILMAN')
propertiesSpecies$V3[ind] = 'NOAA PTR-ToF-MS, TOGA, iWAS'
ind = which(propertiesSpecies$V3 == 'WARNEKE+BLAKE+ppt')
propertiesSpecies$V3[ind] = 'NOAA PTR-ToF-MS, TOGA, WAS'
ind = which(propertiesSpecies$V3 == 'VERES+WENNBERG+WARNEKE+ppt')
propertiesSpecies$V3[ind] = 'NOAA CIMS, CIT-CIMS, NOAA PTR-ToF-MS, TOGA'
ind = which(propertiesSpecies$V3 == 'VERES+WARNEKE')
propertiesSpecies$V3[ind] = 'NOAA CIMS, NOAA PTR-ToF-MS'
ind = which(propertiesSpecies$V3 == 'ROLLINS+RYERSON')
propertiesSpecies$V3[ind] = 'NOAA LIF, NOAA NOyO3'
ind = which(propertiesSpecies$V3 == 'ROLLINS')
propertiesSpecies$V3[ind] = 'NOAA LIF'
ind = which(propertiesSpecies$V3 == 'RYERSON')
propertiesSpecies$V3[ind] = 'NOAA NOyO3'
ind = which(propertiesSpecies$V3 == 'SCHWARZ')
propertiesSpecies$V3[ind] = 'NOAA SP2'
ind = which(propertiesSpecies$V3 == 'APEL')
propertiesSpecies$V3[ind] = 'TOGA'
ind = which(propertiesSpecies$V3 == 'BLAKE+GILMAN+FRIED')
propertiesSpecies$V3[ind] = 'CAMS, WAS, iWAS'
ind = which(propertiesSpecies$V3 == 'VERES+WENNBERG+ppt')
propertiesSpecies$V3[ind] = 'NOAA CIMS, CIT-CIMS, TOGA'
ind = which(propertiesSpecies$V3 == 'WARNEKE+ppt+GILMAN+BLAKE')
propertiesSpecies$V3[ind] = 'NOAA PTR-ToF-MS, TOGA, WAS, iWAS'



ind = which(allBOTH.filter.mean$USEME != 0)
write.csv(propertiesSpecies[ind,],file='propertiesSpecies.csv')


# 1hz vs 5Hz
ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' & is.finite(allBOTH.filter$FinalEF))
fuelshapes = c(19,0,2,18,7,8,4,6,1,12,13)
fuellimits = c("corn","soybean","rice","winter wheat","grass","pile","slash","shrub","Blackwater")
cbp1 = c('#377eb8','#e41a1c','#4daf4a','#a65628','#984ea3','#ff7f00',"#000000")

ind = which(allBOTH.filter$variable == 'CH4_DACOM_DISKIN'  & allBOTH.filter$fuel != 'coniferous/decidous' & allBOTH.filter$fuel != 'forest'& 
              is.finite(allBOTH.filter$FinalEF) & is.finite(allBOTH.filter$EF1COintfill.5hz))
ggplot(allBOTH.filter[ind,])+geom_point(aes(x=EF1CO.1hz, y=EF1intfill.1hz, col=fuel2, shape=fuelORIG), size=3)+ 
  theme_classic()   + scale_color_manual(values = c(cbp1 ))+
  scale_shape_manual(values =fuelshapes)+ylab('Slope method') + xlab('Discrete method')+
  geom_abline(slope=1,intercept=0)
c( mean(allBOTH.filter$FinalEF[ind]),mean(allBOTH.filter$EF1COintfill.5hz[ind], na.rm=TRUE))

# ----- MA/F Aging -------
ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' & allBOTH.filter$MAtoF.5hz < 0.2 & 
              is.finite(allBOTH.filter$MAtoF.5hz) & allBOTH.filter$fuel2 == 'agriculture')
ag = allBOTH.filter[ind,]
ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' & allBOTH.filter$MAtoF.5hz < 0.2 &
              is.finite(allBOTH.filter$MAtoF.5hz) &
              allBOTH.filter$fuel2 == 'prescribed' & allBOTH.filter$fire != 'BlackwaterRiver')
presc = allBOTH.filter[ind,]
t.test(ag$MAtoF.5hz, presc$MAtoF.5hz)

fuelshapes = c(19,0,2,18,7,8,4,6,1,12,13)
fuellimits = c("corn","soybean","rice","winter wheat","grass","pile","slash","shrub","Blackwater")
cbp1 = c('#377eb8','#e41a1c','#4daf4a','#a65628','#984ea3','#ff7f00',"#000000")

ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN'  & allBOTH.filter$fuel != 'coniferous/decidous' & allBOTH.filter$fuel != 'forest')
ggplot(allBOTH.filter[ind,])+geom_point(aes(x=MAtoF.1hz, y=age.1hz/60/60, col=fuel2, shape=fuelORIG), size=3)+ 
  ylim(c(0,4)) + theme_classic()   + scale_color_manual(values = c(cbp1 ))+
  scale_shape_manual(values =fuelshapes)+ylab('Smoke age (h)') + xlab('MA/F')

mean(ag$OAtoOC.5hz)
# ------- OA frac of PM1 -------
ind = which(allBOTH.filter$names == 'Organic Carbon')
OC = allBOTH.filter[ind,]
OC$FinalEF = OC$FinalEF*OC$OAtoOC.5hz
ind = which(allBOTH.filter$names == 'PM1')
PM1 = allBOTH.filter[ind,]
PM1$FinalEFfrac = OC$FinalEF/PM1$FinalEF
ind = which(allBOTH.filter$names == 'Levoglucosan')
Lg = allBOTH.filter[ind,]

ind = which(PM1$fuel2 != 'coniferous/decidous' & is.finite(PM1$FinalEFfrac))
ggplot(PM1[ind,])+geom_point(aes(x=MCE, y=FinalEFfrac, col=fuel2, shape=fuelORIG), size=3)+ 
  ylim(c(0,4)) + theme_classic()   + scale_color_manual(values = c(cbp1 ))+
  scale_shape_manual(values =fuelshapes)+ylab('OA/PM1') + xlab('MCE')+ylim(c(0.4,1.2))
# ------- BC binning -------
ind = which(allBOTH.filter$names == 'Black Carbon')
BC = allBOTH.filter[ind,]
require(creditmodel)
bins=10
bin = cut_equal(BC$MCE, g = bins, sp_values = NULL, cut_bin = "equal_depth")
bin = bin[which(is.finite(bin))]

alt_group = ( floor(BC$FinalEF *2) )
dothis=1
if (dothis == 1){
  for (i in 1:(length(bin))){
    if (i == 1){
      ind = which(BC$MCE <= bin[i])
      alt_group[ind] = i
      counts = length(ind)
    } else{
      ind = which(BC$MCE > bin[i-1] & BC$MCE <= bin[i])
      counts = c(counts, length(ind))
      alt_group[ind] = i
    }
    #print(c(i,min(aircraft.m.olympic$ALTP[ind]), max(aircraft.m.olympic$ALTP[ind])))
  }
  # last bin
  ind = which(BC$MCE > bin[length(bin)])
  counts = c(counts, length(ind))
  alt_group[ind] = i+1
}  
BC$alt_group = alt_group
BC.mean = aggregate(BC,by=list(BC$fuel2,BC$alt_group), 'mean', na.rm=TRUE)
ind = which(BC.mean$Group.1 == 'agriculture')
plot(BC.mean$MCE[ind], BC.mean$FinalEF[ind], pch=19)
cor.test(BC.mean$MCE[ind], BC.mean$FinalEF[ind])

# ----- Monoterpenes ------
ind = which(allBOTH.filter$fuel2 == 'agriculture' & allBOTH.filter$names == 'Monoterpenes')
ag =allBOTH.filter[ind,]
ind = which(allBOTH.filter$fuel2 == 'prescribed' & allBOTH.filter$names == 'Monoterpenes' & allBOTH.filter$fire != 'BlackwaterRiver')
presc =allBOTH.filter[ind,]
t.test(ag$FinalEF, presc$FinalEF)

ind = which(allBOTH.filter$fuel2 == 'agriculture' & allBOTH.filter$names == 'Furan' & allBOTH.filter$USEME == 1)
ag =allBOTH.filter[ind,]
ind = which(allBOTH.filter$fuel2 == 'prescribed' & allBOTH.filter$names == 'Furan' & allBOTH.filter$USEME == 1)
presc =allBOTH.filter[ind,]
t.test(ag$FinalEF, presc$FinalEF)
 # ------- RGF --------
ind = which(allBOTH.filter$fuel2 == 'prescribed' & allBOTH.filter$names == 'Glyoxal' )
ind2 = which(allBOTH.filter$fuel2 == 'prescribed' & allBOTH.filter$names == 'Formaldehyde' & allBOTH.filter$USEME ==1 )
RGF_presc = allBOTH.filter$FinalERtoCO[ind]/allBOTH.filter$FinalERtoCO[ind2]
ind = which(allBOTH.filter$fuel2 == 'agriculture' & allBOTH.filter$names == 'Glyoxal')
ind2 = which(allBOTH.filter$fuel2 == 'agriculture' & allBOTH.filter$names == 'Formaldehyde' & allBOTH.filter$USEME ==1 )
RGF_ag = allBOTH.filter$FinalERtoCO[ind]/allBOTH.filter$FinalERtoCO[ind2]

