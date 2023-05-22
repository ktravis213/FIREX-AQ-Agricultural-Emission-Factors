# Run Everything
require(pals)

print('Calculate Emission Factors') #---------- Calculate EFs ----------
source("CalculateEmissionFactorsv2.R")
print('Emissions Analysis') #---------- Analyze EFs ----------
source("FIREXEmissionsAnalysis.R")
print('Table 1') #---------- Table 1 ----------
source("Table1.R") # Tble of fire and plume counts
print("EF Table") #---------- Table 4 ----------
source("TestTableJustEFv2add.R")
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

# ---- CH4 1hz vs 5Hz -----
ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' & is.finite(allBOTH.filter$FinalEF))
fuelshapes = c(19,0,2,18,7,8,4,6,1,12,13)
fuellimits = c("corn","soybean","rice","winter wheat","grass","pile","slash","shrub","Blackwater")
cbp1 = c('#377eb8','#e41a1c','#4daf4a','#a65628','#984ea3','#ff7f00',"#000000")

ind = which(allBOTH.filter$variable == 'CH4_DACOM_DISKIN'  & 
              allBOTH.filter$fuel != 'coniferous/decidous' & allBOTH.filter$fuel != 'forest'& 
              is.finite(allBOTH.filter$FinalEF) & is.finite(allBOTH.filter$EF1COintfill.1hz))
ggplot(allBOTH.filter[ind,])+geom_point(aes(x=EF1CO.1hz, y=EF1intfill.1hz, col=fuel2, shape=fuelORIG), size=3)+ 
  theme_classic()   + scale_color_manual(values = c(cbp1 ))+
  scale_shape_manual(values =fuelshapes)+ylab('Slope method') + xlab('Discrete method')+
  geom_abline(slope=1,intercept=0)
c( mean(allBOTH.filter$FinalEF[ind]),mean(allBOTH.filter$EF1COintfill.1hz[ind], na.rm=TRUE))

# --- MCE corn vs. others ----
ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' & is.finite(allBOTH.filter$FinalEF) & allBOTH.filter$fuelORIG == 'corn')
corn = allBOTH.filter$MCE[ind]
ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' & is.finite(allBOTH.filter$FinalEF) & allBOTH.filter$fuelORIG == 'rice')
rice = allBOTH.filter$MCE[ind]
ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' & is.finite(allBOTH.filter$FinalEF) & allBOTH.filter$fuelORIG == 'soybean')
soy = allBOTH.filter$MCE[ind]
ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' & is.finite(allBOTH.filter$FinalEF) & allBOTH.filter$fuelORIG == 'pile')
pile = allBOTH.filter$MCE[ind]
ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' & is.finite(allBOTH.filter$FinalEF) & allBOTH.filter$fuelORIG == 'slash')
slash = allBOTH.filter$MCE[ind]

# ----- MA/F Aging -------
ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' & allBOTH.filter$MAtoF.5hz < 0.2 & is.finite(allBOTH.filter$FinalEF) &
              is.finite(allBOTH.filter$MAtoF.5hz) & allBOTH.filter$fuel2 == 'agriculture')
ag = allBOTH.filter[ind,]
ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' & allBOTH.filter$MAtoF.5hz < 0.2 & is.finite(allBOTH.filter$FinalEF) &
              is.finite(allBOTH.filter$MAtoF.5hz) &
              allBOTH.filter$fuel2 == 'prescribed' & allBOTH.filter$fire != 'BlackwaterRiver')
presc = allBOTH.filter[ind,]
t.test(ag$MAtoF.5hz, presc$MAtoF.5hz)
t.test(ag$OAtoOC.5hz, presc$OAtoOC.5hz)

fuelshapes = c(19,0,2,18,7,8,4,6,1,12,13)
fuellimits = c("corn","soybean","rice","winter wheat","grass","pile","slash","shrub","Blackwater")
cbp1 = c('#377eb8','#e41a1c','#4daf4a','#a65628','#984ea3','#ff7f00',"#000000")

ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN'  & allBOTH.filter$fuel != 'coniferous/decidous' & allBOTH.filter$fuel != 'forest')
ggplot(allBOTH.filter[ind,])+geom_point(aes(x=MAtoF.1hz, y=age.1hz/60/60, col=fuel2, shape=fuelORIG), size=3)+ 
  ylim(c(0,4)) + theme_classic()   + scale_color_manual(values = c(cbp1 ))+
  scale_shape_manual(values =fuelshapes)+ylab('Smoke age (h)') + xlab('MA/F')

mean(ag$OAtoOC.5hz)
# ------- OA frac of PM1 -------
ind = which(allBOTH.filter$names == 'Organic carbon')
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
ind = which(allBOTH.filter$names == 'Black carbon')
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

# ---- speciation ------
ind = which(allBOTH.filter$fuel2 == 'agriculture' | allBOTH.filter$fuel2 == 'prescribed' | allBOTH.filter$fuel2 == 'grass' | allBOTH.filter$fuel2 == 'Blackwater')
tmp = allBOTH.filter[ind,]
ind = which(tmp$variable == 'MVKMAC_NOAAPTR_ppbv_WARNEKE')
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'Methyl vinyl ketone' & tmp$USEME == 2)
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'Methacrolein' & tmp$USEME == 2)
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == '2-Butenals' & tmp$USEME == 2)
mean(tmp$FinalEF[ind], na.rm=TRUE)

ind = which(tmp$variable == 'Monoterpenes_NOAAPTR_ppbv_WARNEKE')
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'Tricyclene' & tmp$USEME == 1)
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'Camphene' & tmp$USEME == 1)
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'Myrcene' & tmp$USEME == 1)
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == '⍺-Pinene' & tmp$USEME == 1)
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'β-Pinene' & tmp$USEME == 1)
mean(tmp$FinalEF[ind], na.rm=TRUE)

ind = which(tmp$names == 'Acetone/Propanal' & tmp$USEME == 1)
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'Acetone' )
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'Propanal')
mean(tmp$FinalEF[ind], na.rm=TRUE)

ind = which(tmp$variable == 'C4Carbonyls_NOAAPTR_ppbv_WARNEKE')
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'Methyl ethyl ketone'& tmp$USEME == 1 & tmp$PI != "BLAKE" & tmp$PI != "GILMAN")
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'Isobutanal' & tmp$USEME == 1 )
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'Butanal' & tmp$USEME == 1 )
mean(tmp$FinalEF[ind], na.rm=TRUE)

ind = which(tmp$variable == 'C8Aromatics_NOAAPTR_ppbv_WARNEKE' & tmp$USEME == 1)
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'Ethylbenzene' )
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'o-Xylene')
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'm,p-Xylene')
mean(tmp$FinalEF[ind], na.rm=TRUE)

ind = which(tmp$names == 'Pyrrole/Butenenitrile')
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'Pyrrole')
mean(tmp$FinalEF[ind], na.rm=TRUE)

ind = which(tmp$names == 'Methyl acetate/Ethyl formate/Hydroxyacetone')
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'Methyl acetate')
mean(tmp$FinalEF[ind], na.rm=TRUE)

ind = which(tmp$variable == 'C9Aromatics_NOAAPTR_ppbv_WARNEKE')
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == '1,2,4-Trimethylbenzene')
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == '1,3,5-Trimethylbenzene')
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == '2-Ethyltoluene')
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == '3-Ethyltoluene')
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == '4-Ethyltoluene')
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'i-Propylbenzene')
mean(tmp$FinalEF[ind], na.rm=TRUE)
ind = which(tmp$names == 'n-Propylbenzene')
mean(tmp$FinalEF[ind], na.rm=TRUE)

# ---- NOy/VOC
ind = which(allBOTH.filter$names== 'NOy' & 
              allBOTH.filter$fuel2 != 'coniferous/decidous' &
              allBOTH.filter$fuel2 != "forest" )
noy = allBOTH.filter$FinalERtoCO[ind]

ind = which(allBOTH.filter$names== 'Short-lived VOC' & 
              allBOTH.filter$fuel2 != 'coniferous/decidous' &
              allBOTH.filter$fuel2 != "forest" )
short = allBOTH.filter$FinalERtoCO[ind]
shortEF = allBOTH.filter$FinalEF[ind]

ind = which(allBOTH.filter$names== 'Long-lived VOC' & 
              allBOTH.filter$fuel2 != 'coniferous/decidous' &
              allBOTH.filter$fuel2 != "forest" )
long  = allBOTH.filter$FinalERtoCO[ind]
longEF  = allBOTH.filter$FinalEF[ind]

xx = allBOTH.filter$MCE[ind]
yy= noy/(short + long)
model <- lm(log(yy)~ log(xx))
log(y) = -0.741 + 16.546*log(xx)
y = 0.476637*xx^16.546 
plot(allBOTH.filter$MCE[ind],noy/(short+long), xlab='MCE', ylab='NOy/NMVOCs')
y = 1.72*xx^31 
points(xx,y, col='red')

plot(allBOTH.filter$MCE[ind], short+long)
