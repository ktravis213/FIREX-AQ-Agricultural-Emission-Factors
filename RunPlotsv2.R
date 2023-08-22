# Run Plots
# Figure 1 is from Emily Gargulinski
# Figure 3 is in CalculateEmissionFactorsv2.R - Chips
source("plotSpeciesMCE.R"); source("xioaxi.R"); require(CropScapeR); require(pals)
options(scipen=0, digits=7)
print('Fig 2')
# ---------- Figure 2: Map and Histogram --------------------
fuelshapes = c(19,15,2,18,7,8,4,2,16,12,13)
fuellimits = c("corn","soybean","rice","winter wheat","grass","pile","slash","shrub","BlackwaterRiver")

cbp1a <- c( "#E69F00", "#56B4E9", "#009E73",
            "#F0E442", "#0072B2", "#D55E00", 
            "#CC79A7","#000000","#999999","#CC0000")
# colorblind safe option
cbp1b = c('#377eb8','#e41a1c','#4daf4a','#f781bf','#a65628','#984ea3','#ff7f00','#ffff33')
cbp1c=c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9',
        '#74add1','#4575b4','#313695')
cbp1d <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
           "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbp1 = cbp1a

domap=0
if (domap == 1){
  # Example 3. Retrieve data for a rectangle box defined by four corner points in 2018.
  #data <- GetCDLImage(aoi = c(130783,2203171,153923,2217961), year1 = 2019, year2 = 2019, type = 'b')
  #head(data, 5)
  
  map.us <- get_stamenmap(c(-98.5,30,-82,41), zoom=5)
  ind = which(allBOTH.filter$fire == 'BlackwaterRiver'); allBOTH.filter$fuel[ind] = 'BlackwaterRiver'
  ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' & allBOTH.filter$fuel != 'coniferous/decidous' & allBOTH.filter$fuel != 'forest')
  all5hz.map = allBOTH.filter[ind,]
  mapfire = ggmap(map.us) +
    geom_point(data = all5hz.map, aes(x = lon.5hz,y = lat.5hz,colour =fuelORIG, shape=fuelORIG), size=6,stroke=2,alpha=0.65)+ 
    ggtitle("Ag Fires, 8/21-9/03 ")+  scale_color_manual(values = c(cbp1 ), 
                                                         limits=fuellimits)+
    scale_shape_manual(values =fuelshapes, limits=fuellimits)+labs(col="",shape="")
}
#ggmap(map.us)+geom_point(data=allBOTH.filter[ind,],aes(x=lon.5hz, y=lat.5hz,col=MCE),size=5)+
#  theme_classic()+scale_color_viridis_b()
# 
# ind = which(all5hz.map$fuel == 'corn' | all5hz.map$fuel == 'rice' | 
#               all5hz.map$fuel == 'soybean' | all5hz.map$fuel == 'winter wheat')   
# ind = which(all5hz.map$fuel == 'grass' | all5hz.map$fuel == 'shrub')   
# ind = which(all5hz.map$fuel == 'slash' | all5hz.map$fuel == 'pile')   
# 
# ggmap(map.us) +
#   geom_point(data = all5hz.map[ind,], aes(x = lon,y = lat,colour =fuel, size=catCO))

# ---- MCE  histogram
ind = which(allBOTH.filter$variable == 'CO_DACOM_DISKIN' & 
              allBOTH.filter$fuel2 != 'coniferous/decidous' &
              allBOTH.filter$fuel2 != "forest" &
             is.finite(allBOTH.filter$FinalEF))
tmpMCE = allBOTH.filter[ind,]
MCEhist = ggplot(tmpMCE, aes(x=MCE, col=fuelORIG, fill=fuelORIG)) + geom_histogram()+theme_classic()+
  scale_color_manual(values = c(cbp1 ),limits=fuellimits) + scale_fill_manual(values = c(cbp1 ),
                                                                              limits=fuellimits) +
  theme(legend.position = "top")+labs(fill="")
# ----- Add GG ------
ggMCE = 0.900; ggSD = 0.023
ggdat = as.data.frame(cbind(ggMCE,ggSD))
MCEhist + geom_point(data=ggdat,mapping = aes(x=ggMCE,y=0))
ggsave(filename = 'MCEhist.pdf',MCEhist,width = 6,height = 6,units = 'in')

if (domap == 1){ggsave(filename = 'Figure1_FireMap.pdf',mapfire,width = 6,height = 6,units = 'in')}

# 
print('Fig 3')
# ---- Figure 3: Plot CH4 EF vs. MCE for separated fuel types -----
CH4vsMCE = plotSpeciesMCE(allBOTH.filter,'CH4_DACOM_DISKIN','Methane','CH4','CH4')
# Figure CH4
nP = 1
aspR = 1.75
ww = 10
ggsave('CH4vsMCE.pdf',CH4vsMCE,width = ww, height=ww/aspR)

print('Fig 4')
# --------- Figure 4: VOC contribution pie chart -------
ind = which(allBOTH.filter$names == ' Furan and fragments')
allBOTH.filter$USEME[ind] =0
ind = which(allBOTH.filter$names == 'Phenol' & allBOTH.filter$PI == 'WARNEKE')
allBOTH.filter$USEME[ind] =0
ind = which(allBOTH.filter$names == 'Phenol' & allBOTH.filter$PI == 'WARNEKE')
allBOTH.filter$USEME[ind] =0
ind = which(allBOTH.filter$fuelORIG == 'corn' | allBOTH.filter$fuelORIG == 'rice' | allBOTH.filter$fuelORIG == 'soybean' |
              allBOTH.filter$fuelORIG == 'winter wheat' | allBOTH.filter$fuelORIG == 'pile' | allBOTH.filter$fuelORIG == 'slash' |
              allBOTH.filter$fuelORIG == 'grass' | allBOTH.filter$fuelORIG == 'shrub' | 
              allBOTH.filter$fuelORIG == 'Blackwater')
allBOTH.filter.ag = allBOTH.filter[ind,]
ind = which(allBOTH.filter.ag$Category == 1 & allBOTH.filter.ag$USEME == 1 )
allBOTH.filter.ag = allBOTH.filter.ag[ind,]
allBOTH.filter.ag.avg = aggregate(allBOTH.filter.ag, by=list(allBOTH.filter.ag$names, allBOTH.filter.ag$PI), FUN='mean', na.rm=TRUE)

tmp = as.data.frame(cbind(names=allBOTH.filter.ag.avg$Group.1,PI=allBOTH.filter.ag.avg$Group.2, EF=allBOTH.filter.ag.avg$FinalEF, lifetime=allBOTH.filter.ag.avg$LifetimeCat,
                          MW=allBOTH.filter.ag.avg$mWs, kOH=allBOTH.filter.ag.avg$OHrate.1hz, nC=allBOTH.filter.ag.avg$nCs))
# Calculate reactivity
tmp$dummyppt = as.numeric(tmp$EF)*6.022E23*1E12/1E15/as.numeric(tmp$MW)/2.69E19# 1km3, 273K, 1atm, 1kg fuel
tmp$dummyohr = tmp$dummyppt*1E-12*as.numeric(tmp$kOH)*2.69E19

write.csv(tmp, file='VOCcontributions.csv')

print('Fig PM2.5')
# --------- Figure X: PM2.5 -------
ind = which(allBOTH.filter$names == ' Furan and fragments')
allBOTH.filter$USEME[ind] =0
ind = which(allBOTH.filter$fuel == 'corn' | allBOTH.filter$fuel == 'rice' | allBOTH.filter$fuel == 'soybean' |
              allBOTH.filter$fuel == 'winter wheat' | allBOTH.filter$fuel == 'pile' | allBOTH.filter$fuel == 'slash' |
              allBOTH.filter$fuel == 'grass' | allBOTH.filter$fuel == 'shrub' | 
              allBOTH.filter$fuel == 'BlackwaterRiver')
allBOTH.filter.ag = allBOTH.filter[ind,]
ind = which(allBOTH.filter.ag$Category == 1 & allBOTH.filter.ag$USEME == 1 )
allBOTH.filter.ag = allBOTH.filter.ag[ind,]
allBOTH.filter.ag.avg = aggregate(allBOTH.filter.ag, by=list(allBOTH.filter.ag$names), FUN='median', na.rm=TRUE)

tmp = as.data.frame(cbind(names=allBOTH.filter.ag.avg$Group.1, EF=allBOTH.filter.ag.avg$FinalEF, allBOTH.filter.ag.avg$lifetime))

write.csv(tmp, file='VOCcontributions.csv')

print('Fig 5')
# ------- Figure  5: total VOCs ----------------
ind = which(allBOTH.filter$fuel2 == 'agriculture' & allBOTH.filter$names == 'Short-lived VOC')
short.ag = allBOTH.filter[ind,]
cor.test(allBOTH.filter$MCE[ind], allBOTH.filter$FinalEF[ind])
ind = which(allBOTH.filter$fuel2 == 'agriculture' & allBOTH.filter$names == 'Long-lived VOC')
long.ag = allBOTH.filter[ind,]

cor.test(allBOTH.filter$MCE[ind], allBOTH.filter$FinalEF[ind])
ind = which(allBOTH.filter$fuel2 == 'agriculture' & allBOTH.filter$names == 'Methane')
ch4.ag = allBOTH.filter[ind,]
ch4.ag$FinalEF = ch4.ag$FinalEF + short.ag$FinalEF + long.ag$FinalEF

ind = which(allBOTH.filter$fuel2 == 'prescribed' & allBOTH.filter$names == 'Short-lived VOC')
cor.test(allBOTH.filter$MCE[ind], allBOTH.filter$FinalEF[ind])
ind = which(allBOTH.filter$fuel2 == 'prescribed' & allBOTH.filter$names == 'Long-lived VOC')
cor.test(allBOTH.filter$MCE[ind], allBOTH.filter$FinalEF[ind])

# Figure ShortVOC
ShortVOCvsMCE = plotSpeciesMCE(allBOTH.filter,'Short-lived VOC','Short-lived VOC','VOC','VOC')
LongVOCvsMCE = plotSpeciesMCE(allBOTH.filter,'Long-lived VOC','Long-lived VOC','VOC','VOC')

#  NOx/VOC?
ind = which(allBOTH.filter$variable == 'NOx (as NO)')
NOx = allBOTH.filter[ind,]
ind = which(allBOTH.filter$variable == 'VOC')
VOC = allBOTH.filter[ind,]
NOx$FinalEF = NOx$FinalEF/VOC$FinalEF
NOx$variable = 'NOx_VOC'
allBOTH.filter = rbind(allBOTH.filter,NOx)

NOxVOCvsMCE = plotSpeciesMCE(allBOTH.filter,'NOx_VOC','NOx','NOx','NOx')
shortlong= ggarrange(ShortVOCvsMCE,LongVOCvsMCE, NOxVOCvsMCE,
                     common.legend = TRUE,nrow=1,ncol=3,
                     labels = c("a)","b)","c)"),
                     hjust = c(-5,-5,-5))

ggsave("Shortlong.pdf", shortlong, width = ww*3, height=ww/aspR)


print('Fig 6')
# ------ Figure 6: Selection of VOCs -------
FormaldehydevsMCE    = plotSpeciesMCE(allBOTH.filter,'Formaldehyde_FRIED_HANISCO','Formaldehyde','Formaldehyde','Formaldehyde')
FormaldehydevsMCE_HANISCO   = plotSpeciesMCE(allBOTH.filter,'CH2O_ISAF_HANISCO','Formaldehyde','Formaldehyde','Formaldehyde')
FormaldehydevsMCE_FRIED   = plotSpeciesMCE(allBOTH.filter,'CH2O_CAMS_pptv_FRIED','Formaldehyde','Formaldehyde','Formaldehyde')
ind = which(allBOTH.filter$formula == 'CH2O' & allBOTH.filter$PI == 'HANISCO' & allBOTH.filter$fuel2 != 'forest' & allBOTH.filter$fuel2 != "coniferous/decidous")
tmpHANISCO = allBOTH.filter[ind,]
ind = which(allBOTH.filter$formula == 'CH2O' & allBOTH.filter$PI == 'FRIED' & allBOTH.filter$fuel2 != 'forest' & allBOTH.filter$fuel2 != "coniferous/decidous")
tmpFRIED = allBOTH.filter[ind,]
t.test(tmpHANISCO$FinalEF, tmpFRIED$FinalEF)

GLYCvsMCE    = plotSpeciesMCE(allBOTH.filter,'GlycolaldehydeCH3COOH_NOAAPTR_ppbv_WARNEKE','GlycoaldehydeCH3COOH','GlycoaldehydeCH3COOH','GlycoaldehydeCH3COOH')

ggarrange(FormaldehydevsMCE, FormaldehydevsMCE_FRIED, FormaldehydevsMCE_HANISCO,ncol = 3,labels=c("Average","Fried","Hanisco"))
AcetaldehydevsMCE    = plotSpeciesMCE(allBOTH.filter,'CH3CHO_NOAAPTR_ppbv_WARNEKE','Acetaldehyde','Acetaldehyde','Acetaldehyde')
GuaiacolvsMCE    = plotSpeciesMCE(allBOTH.filter,'Guaiacol_NOAAPTR_ppbv_WARNEKE','Guaiacol','Guaiacol','Guaiacol')
SYRvsMCE    = plotSpeciesMCE(allBOTH.filter,'Syringol_NOAAPTR_ppbv','Syringol','Syringol','Syringol')
FuranvsMCE    = plotSpeciesMCE(allBOTH.filter,'Furan_ppt','Furan','Furan','Furan')
MeFuranvsMCE    = plotSpeciesMCE(allBOTH.filter,'x2MeFuranx3MeFuran_NOAAPTR_ppbv_WARNEKE','23Methylfuran','23Methylfuran','23Methylfuran')
CatecholvsMCE    = plotSpeciesMCE(allBOTH.filter,'Catecholx5MeFurfural_NOAAPTR_ppbv_WARNEKE','Catechol','Catechol','Catechol')
HACvsMCE    = plotSpeciesMCE(allBOTH.filter,'C3H6O2_NOAAPTR_ppbv_WARNEKE','Hydroxyacetone','Hydroxyacetone','Hydroxyacetone')
PhenolvsMCE    = plotSpeciesMCE(allBOTH.filter,'PHENOL_WENNBERG','Phenol2','Phenol2','Phenol2') # EF is SO much lower than Akagi etc
PropeneHPvsMCE    = plotSpeciesMCE(allBOTH.filter,'PROPENE.HP_WENNBERG','Propene hydroxyperoxide','PropeneHP','PropeneHP')

BenzenevsMCE    = plotSpeciesMCE(allBOTH.filter,'Benzene_NOAAPTR_ppbv_WARNEKE','Benzene','Benzene','Benzene')
NaphthalenevsMCE    = plotSpeciesMCE(allBOTH.filter,'Naphthalene_NOAAPTR_ppbv_WARNEKE','Naphthalene','Naphthalene','Naphthalene')
C2H2vsMCE    = plotSpeciesMCE(allBOTH.filter,'Ethyne_GILMAN_BLAKE','Acetylene','Ethyne','Ethyne')
MonovsMCE    = plotSpeciesMCE(allBOTH.filter,'Monoterpenes_NOAAPTR_ppbv_WARNEKE','Terpenes','Terpenes','Monoterpenes')
CH3OHvsMCE = plotSpeciesMCE(allBOTH.filter,'CH3OH_NOAAPTR_ppbv_WARNEKE','Methanol','Methanol','Methanol')
GLYOXALvsMCE = plotSpeciesMCE(allBOTH.filter,'CHOCHO_ACES_WOMACK','Glyoxal','Glyoxal','Glyoxal')
ISOPvsMCE = plotSpeciesMCE(allBOTH.filter,'Isoprene_BLAKE_ppt_','Isoprene','Isoprene','Isoprenehydroperoxyaldehydes')
HCOOHvsMCE = plotSpeciesMCE(allBOTH.filter,'Formic acid_VERES_WARNEKE','Formic Acid','Formic acid','HCOOH')

#C2H4vsMCE    = plotSpeciesMCE(allBOTH.filter,'Ethene_GILMAN_BLAKE','C2H4','Ethene','Ethene')
#AcetvsMCE    = plotSpeciesMCE(allBOTH.filter,'Acetone_NOAAiWAS_GILMAN','Acetone','Acetone','Acetone')
#AcrvsMCE    = plotSpeciesMCE(allBOTH.filter,'Acrolein_NOAAPTR_ppbv_WARNEKE','Acrolein','Acrolein','Acrolein')
#EthanevsMCE    = plotSpeciesMCE(allBOTH.filter,'C2H6_CAMS_pptv_FRIED','Ethane','Ethane','Ethane')
#TOLUvsMCE    = plotSpeciesMCE(allBOTH.filter,'Toluene_NOAAPTR_ppbv_WARNEKE','Toluene','Toluene','Toluene')
#Depolymerization in lignin (300–500C)produces guaiacols, (iso)eugenol, and syringol.
#Furans and furfurals are dominantly formed from cellulose and hemicellulose (300–400oC).
#Furan2MevsMCE    = plotSpeciesMCE(allBOTH.filter,'2-Methylfuran_BLAKE_ppt','2-Methylfuran','2-Methylfuran','2-Methylfuran')
#Furan3MevsMCE    = plotSpeciesMCE(allBOTH.filter,'3-Methylfuran_BLAKE_ppt','3-Methylfuran','3-Methylfuran','3-Methylfuran')
# FurfuralvsMCE    = plotSpeciesMCE(allBOTH.filter,'Furfural_ppt','Furfural','Furfural','Furfural')
#FurfuralvsMCE    = plotSpeciesMCE(allBOTH.filter,'Furfural_NOAAPTR_ppbv_WARNEKE','Furfural','Furfural','Furfural')
#Higher temperatures allow reaction #of functional groups and covalent bonds in polymers and
#monomers. The resulting fragmentation emits various VOCs:for example, hydroxyacetone, acetaldehyde, and acetic acid
#from depolymerization of cellulose and/or hemicellulose
#AceticvsMCE    = plotSpeciesMCE(allBOTH.filter,'GlycolaldehydeCH3COOH_NOAAPTR_ppbv_WARNEKE','GlycoaldehydeCH3COOH','GlycoaldehydeCH3COOH','GlycoaldehydeCH3COOH')

#Higher temperature pyrolysis breaks progressively stronger bonds in char (> 500C).
#This aromatization process gives off aromatic compounds with short substituents (e.g., phenol), 
#nonsubstituted aromatics (e.g.,benzene), and polycyclic aromatic hydrocarbons (PAHs) such as naphthalene).
# PhenolvsMCE2    = plotSpeciesMCE(allBOTH.filter,'Phenol_NOAAPTR_ppbv_WARNEKE','Phenol','SO2','SO2')
#?
#MonovsMCE    = plotSpeciesMCE(allBOTH.filter,'beta-Pinene/Myrcene_BLAKE_ppt','beta-Pinene/Myrcene','SO2','SO2')
#CatvsMCE = plotSpeciesMCE(allBOTH.filter,'Catecholx5MeFurfural_NOAAPTR_ppbv_WARNEKE','Catechol/5MeFurfural','Catechol','Catechol')
#EOHvsMCE = plotSpeciesMCE(allBOTH.filter,'C2H5OH_NOAAPTR_ppbv_WARNEKE','Ethanol','Ethanol','Ethanol')


allVOCplot = ggarrange(FormaldehydevsMCE,AcetaldehydevsMCE,MonovsMCE, 
                       GuaiacolvsMCE,CatecholvsMCE,PhenolvsMCE,
                       BenzenevsMCE,NaphthalenevsMCE, GLYOXALvsMCE,
                       common.legend = TRUE,nrow=3,ncol=3,
                       labels = c("a)","b)","c)","d)","e)","f)","g)","h)","i)"),
                       hjust = c(-8))
                       
#SYRvsMCE,  HACvsMCE, 
#C2H2vsMCE,  CH3OHvsMCE,
#ISOPvsMCE, HCOOHvsMCE,

ggsave("allVOCplot.pdf", allVOCplot, width = ww*3, height=ww*3/aspR)

# ------ AROMATICs for supplement -----
StyrenevsMCE    = plotSpeciesMCE(allBOTH.filter,'Styrene_NOAAPTR_ppbv_WARNEKE','Styrene','Styrene','Styrene')
TOLUvsMCE    = plotSpeciesMCE(allBOTH.filter,'Toluene_NOAAPTR_ppbv_WARNEKE','Toluene','Toluene','Toluene')
C8vsMCE    = plotSpeciesMCE(allBOTH.filter,'C8Aromatics_NOAAPTR_ppbv_WARNEKE','Xylenes','Xylenes','Xylenes')
C9vsMCE    = plotSpeciesMCE(allBOTH.filter,'C9Aromatics_NOAAPTR_ppbv_WARNEKE','C9','C9','C9')
#PhenolvsMCE2    = plotSpeciesMCE(allBOTH.filter,'Phenol_NOAAPTR_ppbv_WARNEKE','Phenol','Phenol','Phenol')
CycHexvsMCE= plotSpeciesMCE(allBOTH.filter,'CycHexane_WAS_BLAKE','Cyclohexane','Cyclohexane','Cyclohexane')
# Toluene outlier
aroms = ggarrange(StyrenevsMCE,
                  TOLUvsMCE, C8vsMCE, C9vsMCE, common.legend = TRUE)

print('Fig 7')
# ------ Figure 7: Nitrogen species -----
CresolvsMCE = plotSpeciesMCE(allBOTH.filter,'CRESOL_WENNBERG','Cresol','NO4','NO4')
NCvsMCE = plotSpeciesMCE(allBOTH.filter,'NITROCRESOL_WENNBERG','Nitrocresol','NO4','NO4')
NPvsMCE = plotSpeciesMCE(allBOTH.filter,'NITROPHENOL_WENNBERG','Nitrophenol','NO4','NO4')
NCTvsMCE = plotSpeciesMCE(allBOTH.filter,'NITROCATECHOL_WENNBERG','Nitrocatechol','NO4','NO4')
NMCTvsMCE = plotSpeciesMCE(allBOTH.filter,'NITROMETHYLCATECHOL_WENNBERG','Nitromethylcatechol','NO4','NO4')

NOvsMCE = plotSpeciesMCE(allBOTH.filter,'Nitrogen oxide_ROLLINS_RYERSON','NO','NO','NO')
NO2vsMCE = plotSpeciesMCE(allBOTH.filter,'Nitrogen dioxide_WOMACK_RYERSON','NO2','NO2','NO2')
NOxvsMCE = plotSpeciesMCE(allBOTH.filter,'NOx (as NO)','NOx (as NO)','NOx (as NO)','NOxasNO')
HNO2vsMCE = plotSpeciesMCE(allBOTH.filter,'Nitrous acid_WOMACK_VERES','HNO2','HNO2','HNO2')
HNO2vsMCE2 = plotSpeciesMCE(allBOTH.filter,'HNO2_ACES_WOMACK','HNO2','HNO2','HNO2')
HNO2vsMCE3 = plotSpeciesMCE(allBOTH.filter,'HNO2_NOAACIMS_VERES','HNO2','HNO2','HNO2')
CH3CNvsMCE = plotSpeciesMCE(allBOTH.filter,'CH3CN_NOAAPTR_ppbv_WARNEKE','Acetonitrile','Acetonitrile','Acetonitrile')
HCNvsMCE = plotSpeciesMCE(allBOTH.filter,'Hydrogen cyanide_VERES_WENNBERG__','Hydrogen cyanide','Hydrogen cyanide','HCN')
HCNvsMCE1 = plotSpeciesMCE(allBOTH.filter,'HCN_NOAAPTR_ppbv_WARNEKE','Hydrogen cyanide','HCN','HCN')
HCNvsMCE2 = plotSpeciesMCE(allBOTH.filter,'HCN_WENNBERG','Hydrogen cyanide','Hydrogen cyanide','HCN')
HCNvsMCE3 = plotSpeciesMCE(allBOTH.filter,'HCN_NOAACIMS_VERES','Hydrogen cyanide','HCN','HCN')
HNO3vsMCE = plotSpeciesMCE(allBOTH.filter,'HNO3_WENNBERG','Nitric acid','HNO3','HNO3')
PANvsMCE = plotSpeciesMCE(allBOTH.filter,'PAN_HUEY','PAN','PAN','PAN')

#NH3vsMCE = plotSpeciesMCE(allBOTH.filter,'NH3_WISTHALER','Ammonia','NH3','NH3')
HNCOvsMCE = plotSpeciesMCE(allBOTH.filter,'Isocyanic acid_VERES_WARNEKE','HNCO','HNCO','HNCO')
CH3NO3vsMCE = plotSpeciesMCE(allBOTH.filter,'Methyl nitrate_BLAKE_ppt','MethylNitrate','MethylNitrate','MethylNitrate')
CH3NO2vsMCE = plotSpeciesMCE(allBOTH.filter, 'CH3NO2_NOAAPTR_ppbv_WARNEKE','Nitromethane','Nitromethane','CH3NO2')
allNOx=ggarrange(NOxvsMCE, HNO2vsMCE,HCNvsMCE, HNCOvsMCE,  common.legend = TRUE,nrow=2,
                 ncol=2,labels = c("a)","b)","c)","d)"), hjust = c(-5,-5,-6,-6))
ggsave("allNOxplot.pdf", allNOx, width = ww*2, height=ww*2/aspR)

# ------ Figure 8: Halogenated species -------
# Figure CH3Cl
CH3ClvsMCE = plotSpeciesMCE(allBOTH.filter,'CH3Cl_WAS_BLAKE','Chloromethane','CH3Cl','Chloromethane')
CH3IvsMCE = plotSpeciesMCE(allBOTH.filter,'Methyl Iodide_BLAKE_ppt','MethylIodide','CH3I','Methyl Iodide')
CH3BrvsMCE = plotSpeciesMCE(allBOTH.filter,'Methyl bromide_ppt_BLAKE','Bromomethane','CH3Br','Bromomethane')
CH2Cl2vsMCE = plotSpeciesMCE(allBOTH.filter,'Dichloromethane_BLAKE_ppt','Dichloromethane','Dichloromethane','Dichloromethane')
allHal = ggarrange(CH3ClvsMCE,CH3BrvsMCE,CH3IvsMCE, common.legend = TRUE,ncol=3,labels = c("a)","b)","c)"),
                   hjust = c(-5,-8,-8))
ggsave('allHal.pdf',allHal, width = ww*3, height=ww/aspR)

print('Fig 9')
# ------ Figure 9: Aerosols ------

PM1vsMCE = plotSpeciesMCE(allBOTH.filter,'PM1','PM1','PM1','PM1')
ECvsMCE = plotSpeciesMCE(allBOTH.filter,'BC_SCHWARZ','BC','Black Carbon','BC')
#ind = which(allBOTH.filter$variable == 'CNgt3nm')
#allBOTH.filter$FinalEF[ind] = allBOTH.filter$FinalEF[ind]*1E15
#CNvsMCE = plotSpeciesMCE(allBOTH.filter,'CNgt3nm','CN','CN','CN')

xiaoxi$OC = xiaoxi$OA/2
OCvsMCE = plotSpeciesMCE(allBOTH.filter,'OC_JIMENEZ','OC','Organic Carbon','OC')
LGvsMCE = plotSpeciesMCE(allBOTH.filter,'C6H10O5_JIMENEZ','Levoglucosan','Levoglucosan','Levoglucosan')

ClvsMCE = plotSpeciesMCE(allBOTH.filter,'NR_Chloride_JIMENEZ','Cl','Chloride','Chl')
SO4vsMCE = plotSpeciesMCE(allBOTH.filter,'Sulfate_JIMENEZ','Sulfate','Sulfate','SO4')
NITvsMCE = plotSpeciesMCE(allBOTH.filter,'Nitrate_JIMENEZ','Nitrate','Nitrate','NO3')
NH4vsMCE = plotSpeciesMCE(allBOTH.filter,'Ammonium_JIMENEZ','Ammonium','Ammonium','NH4')
KvsMCE = plotSpeciesMCE(allBOTH.filter,'Potassium_JIMENEZ','K','K','K')
NCATvsMCE = plotSpeciesMCE(allBOTH.filter,'C6H5NO4_JIMENEZ','4-Nitrocatechol','4-Nitrocatechol','4-Nitrocatechol')
allPM =ggarrange(PM1vsMCE, OCvsMCE,LGvsMCE,ECvsMCE,NITvsMCE, ClvsMCE,NH4vsMCE,KvsMCE,
                 common.legend = TRUE,nrow=2,ncol=4,
                 labels = c("a)","b)","c)","d)","e)","f)","g)","h)"),
                 hjust = c(-5,-5,-6,-6,-5,-5))
ggsave('allPM.pdf',allPM,width = ww*4, height=ww*2/aspR)

print('Fig 10')
# ------- Figure  10: SO2 + H2O2-------------
# Other
H2O2vsMCE = plotSpeciesMCE(allBOTH.filter,'H2O2_WENNBERG','H2O2','H2O2','H2O2')
SO2vsMCE = plotSpeciesMCE(allBOTH.filter,'SO2_LIF_ROLLINS','Sulfur dioxide','SO2','SO2')
CH3SHvsMCE = plotSpeciesMCE(allBOTH.filter,'CH3SH_ppt','CH3SH','CH3SH','CH3SH')
OCSvsMCE = plotSpeciesMCE(allBOTH.filter,'OCS_WAS_BLAKE','OCS','OCS','OCS')

otherFIG = ggarrange(SO2vsMCE,CH3SHvsMCE,H2O2vsMCE,
                     common.legend = TRUE,nrow=1,ncol=3,
                     labels = c("a)","b)","c)"),   hjust = c(-5,-5,-5))
ggsave('OTHERFIG.pdf',otherFIG,width = ww*3, height=ww*1/aspR)
ind = which(allBOTH.filter$formula == 'H2O2')
h2o2=allBOTH.filter[ind,]
h2o2$fuelORIG <- factor(h2o2$fuelORIG, levels = c("winter wheat", "soybean","rice","corn","grass","slash","pile","shrub","Blackwater"))
ind = which(is.finite(h2o2$MAtoF.5hz) & h2o2$MAtoF.5hz < 0.2 & h2o2$fuel2 != 'forest' & h2o2$fuel2 != 'coniferous/decidous')
fuelshapes = c(19,0,2,18,7,8,4,6,1,12,13)
h2o2age =ggplot(h2o2[ind,])+geom_point(aes(x=MAtoF.5hz,y=FinalERtoCO, shape=fuelORIG, col=fuel2), size=3)+theme_classic()+
  scale_color_manual(values = c(cbp1 ))+
  scale_shape_manual(values =fuelshapes)+
  geom_point(aes(x=0,y=0.0013), col='purple', size=3, shape=2)
ggsave('h2o2age.pdf',h2o2age,width = ww*1, height=ww*1/aspR)


# --- Figure 11: Correlation with MCE analysis ----- 
run12=1
if (run12 == 1){
  # USE varsALL from SlopeDifferences.R
  ind = which(varsALL$Category != 5 & varsALL$USEME > 0)
  namesTODO = varsALL[ind,]
  ind = which(namesTODO$Category == 1 & namesTODO$LifetimeCat == 1)
  namesTODO$Category[ind] = 1.5 # Short-lived?
  ## 0 is inorganic, , 1 is VOC, 2 is nitrogen containing, 3 is halogen, 4 aerosol
  
  namesTODO = namesTODO[order(namesTODO$Rval_ag),]
  ind = which(namesTODO$names != 'PM1' & namesTODO$names != 'Carbon Monoxide' &
                namesTODO$names != 'Methane' & namesTODO$kind != 'CO2' & namesTODO$names != 'NOx (as NO)')
  namesTODO = namesTODO[ind,]
  
  library(pals)
  namesTODO$R2_ag = namesTODO$Rval_ag^2
  ind = which(namesTODO$names == 'sum of monoterpenes')
  namesTODO$names[ind] = "Monoterpenes"
  ind = which(namesTODO$names == '2,5-Dimethylfuran/2-Ethylfuran/Other unidentified organic compounds')
  namesTODO$names[ind] = '2,5-Dimethylfuran + 2-Ethylfuran + unknown'
  ind = which(namesTODO$names == '2,3-Butanedione/2-Oxobutanal/1,4-Butanedial')
  namesTODO$names[ind] = '2,3-Butanedione + 2-Oxobutanal + 1,4-Butanedial'
  ind = which(namesTODO$names == '2-Methylphenol (=o-cresol)/Anisol')
  namesTODO$names[ind] = '2-Methylphenol (=o-cresol) + Anisol'
  ind = which(namesTODO$names == 'Acetic acid/Glycolaldehyde')
  namesTODO$names[ind] = 'Acetic acid + Glycolaldehyde'
  
  ind = which( namesTODO$names != " Furan and fragments" & namesTODO$names != 'Long-lived VOC' & namesTODO$names != 'VOC' & namesTODO$names != 'Short-lived VOC' &
                 namesTODO$names != 'Nitrate' & namesTODO$vars != 'Phenol_NOAAPTR_ppbv_WARNEKE' & 
                 namesTODO$names!= 'CCN_034_stdPT' & namesTODO$names != 'CN > 20nm' & namesTODO$names != 'CN > 6nm'& namesTODO$names != 'Ngt100nm_LAS_stdPT') #is.finite(namesTODO$corMCE_ag)  &
  tmp = namesTODO[ind,]
  ind = which(round(tmp$R2_ag,2) >= 0.5 | tmp$Rval_ag > 0)#.3)
  tmp = tmp[ind,]
  tmp$NUM = seq(1,length(tmp$kind))
  require(ggpattern)
  tmp$kind2 = 'NotNerd'
  ind = which(tmp$kind == 'oVOC')
  tmp$kind2[ind] = 'Nerd'
  # Change aerosols to a different pattern
  ind = which(tmp$kind == 'aerosol')
  tmp$kind2[ind] = 'DefNerd'
  ind = which(tmp$names == '4-nitrocatechol'); tmp$kind[ind] = 2
  ind = which(tmp$names == 'Chloride'); tmp$Category[ind] = 3
  ind = which(tmp$names == 'Sulfate'); tmp$Category[ind] = 6
  
  ind = which(tmp$names == 'HNO3_NO3'); tmp$Category[ind] = 2; tmp$names[ind] = 'Nitrate'
  ind = which(tmp$names == 'Organic Carbon'); tmp$Category[ind] = 1; tmp$kind2[ind] = 'Nerd'
  
  agEFMCE = ggplot(tmp , aes(y=round(Rval_ag, digits = 2), x=NUM, fill=factor(Category), pattern=factor(kind2))) + #ylim(0,-1)+
    scale_x_discrete(breaks=c("-1","0","1"),labels=c("-1","0","1"))+
    geom_bar_pattern(position="dodge", stat="identity") + theme_classic()+ 
    geom_text(aes(x =NUM,y = 0,label = names), size=5,
              vjust = 0,  hjust = 0,angle = 0, nudge_y = 0.02, nudge_x = -0.15)+
    ylab('R EF MCE') + xlab(c(""))+ylab("Correlation with MCE") +
    scale_pattern_manual(values = c(Nerd = "stripe", NotNerd = "none",DefNerd="wave")) +
    scale_fill_manual(values=as.vector(polychrome(26)))+ theme(legend.position = "top")+labs(fill="")
  
  agEFMCE = agEFMCE + coord_flip()
  
  ww=10
  ggsave(agEFMCE, file='AgEFMCE.pdf', width = ww, height=ww*1.25)
  #ggsave(prescEFMCE, file='PrescEFMCE.pdf', width = ww, height=ww/aspR)
  #ggsave(grassEFMCE, file='GrassEFMCE.pdf', width = ww, height=ww/aspR)
}

print('Fig 13')
# ------- Figure 13: ? Xiaoxi comparison ------
#source("TestTableJustEF_AgvsPresc.R")
#load('setupdataEF.RData') # from TestTableJustEFv2.R to get the averages accounting for fuel type fraction
outputdata = as.data.frame(cbind(new.data.frame$Category, new.data.frame$LifetimeCat,new.data.frame$mWs,
                                  new.data.frame$names,new.data.frame$formula,new.data.frame$PI,
                                  new.data.frame$FinalEF_mean, new.data.frame$FinalEF_sd, new.data.frame$COUNT_EFFINAL,
                                  new.data.frame.rice$FinalEF_mean, new.data.frame.rice$FinalEF_sd,         new.data.frame.rice$COUNT_EFFINAL,
                                  new.data.frame.soybean$FinalEF_mean, new.data.frame.soybean$FinalEF_sd, new.data.frame.soybean$COUNT_EFFINAL,  
                                  new.data.frame.ag$FinalEF_mean, new.data.frame.ag$FinalEF_sd, new.data.frame.ag$COUNT_EFFINAL,  
                                  
                                  
                                  new.data.frame.slash$FinalEF_mean,   new.data.frame.slash$FinalEF_sd,     new.data.frame.slash$COUNT_EFFINAL,
                                  new.data.frame.pile$FinalEF_mean,   new.data.frame.pile$FinalEF_sd,         new.data.frame.pile$COUNT_EFFINAL,
                                  new.data.frame.presc$FinalEF_mean,   new.data.frame.presc$FinalEF_sd,         new.data.frame.presc$COUNT_EFFINAL,
                                  
                                  
                                  new.data.frame.grass$FinalEF_mean,   new.data.frame.grass$FinalEF_sd,     new.data.frame.grass$COUNT_EFFINAL))

colnames(outputdata) = c("Category", 'LifetimeCat','mWs',  'names','formula','PI',
                                 'FinalEF_mean_corn',    'FinalEF_sd',      'COUNT_EFFINAL',
                                 'FinalEF_mean_rice',    'FinalEF_sd_rice', 'COUNT_EFFINAL_rice',
                                 'FinalEF_mean_soy',     'FinalEF_sd_soy',   'COUNT_EFFINAL_soy',  
                                 'FinalEF_mean_ag',      'FinalEF_sd_ag',    'COUNT_EFFINAL_ag',  
                                 'FinalEF_mean_slash',   'FinalEF_sd_slash','COUNT_EFFINAL_slash',
                                 'FinalEF_mean_pile',    'FinalEF_sd_pile',    'COUNT_EFFINAL_pile',
                                 'FinalEF_mean_presc',   'FinalEF_sd_presc',    'COUNT_EFFINAL_presc',
                                 'FinalEF_mean_grass',   'FinalEF_sd_grass',    'COUNT_EFFINAL_grass')

outputdata$kind = ''
for (i in 1:length(outputdata$Category)){
  ind = which(allBOTH.filter$names == outputdata$names[i])
  outputdata$kind[i] = allBOTH.filter$kind[ind[1]]
}
nP = 1
aspR = 1.75
ww = 10
# ---- EFs from Xiaoxi Liu (2016 ACP SEAC4RS 'Rice', straw)
source("xioaxi.R")
source("getERsv2.R")
require(robustbase)
xiaoxi2$EF = as.numeric(xiaoxi2$EF)
xiaoxi2.med = aggregate(xiaoxi2,by=list(xiaoxi2$name), FUN='mean', na.rm=TRUE)
xiaoxi2.sd = aggregate(xiaoxi2,by=list(xiaoxi2$name), FUN='sd', na.rm=TRUE)
xiaoxi2.med$EF_sd = xiaoxi2.sd$EF
xiaoxi2.med$name =xiaoxi2.med$Group.1
#xiaoxi$fuel = c('Corn','Corn','Corn', 'Rice','Rice', 'Rice','Rice','Rice','Rice','Rice','Rice','Soybean', 'Popcorn','DblCrpWinWht_Soy',   'WoodyWetlands')

ind = which(outputdata$names == 'Black carbon')
outputdata$names[ind] = 'BC'
ind = which(outputdata$names == 'Organic aerosol')
outputdata$names[ind] = "OA"; outputdata$formula[ind] = 'OA' #kludge
ind = which(outputdata$names == 'Organic carbon')
outputdata$names[ind] = "OC"
# ---
ind = which(outputdata$names == 'Chloride')
outputdata$formula[ind] = 'Chl'
ind = which(outputdata$names == 'sum of monoterpenes')
outputdata$formula[ind] = 'Monoterpenes'
ind = which(outputdata$names == 'MVK/MACR')
outputdata$formula[ind] = 'MVKMACRcrotonaldehyde'

ind = which(outputdata$names == 'Isoprene')
outputdata$formula[ind] = 'Isoprenepentadienescyclopentenefuran'  # unclear if I should add pentadienes, cyclopoentene, furan to this
ind = which(outputdata$names == "Methyl acetate/Ethyl formate/Hydroxyacetone")
outputdata$names[ind] = 'Hydroxyacetone'
outputdata$XiaoxiEF = NaN
outputdata$Xiaoxisd = NaN
outputdata$XiaoxiName = NaN
outputdata$XiaoxiN = NaN
cc = colnames(xiaoxi)
iK= which(cc == 'NOxasNO')
cc[iK] = 'NOx (as NO)'
for (i in 1:length(outputdata$XiaoxiEF)){
  ind = which(outputdata$formula[i] == xiaoxi2.med$name |
                outputdata$names[i] == xiaoxi2.med$name)
  if (length(ind) > 0){ outputdata$XiaoxiEF[i] = xiaoxi2.med$EF[ind]}
  if (length(ind) > 0){ outputdata$XiaoxiName[i] = xiaoxi2.med$name[ind]}
  if (length(ind) > 0){ outputdata$Xiaoxisd[i] = xiaoxi2.med$EF_sd[ind]}
  if (length(ind) > 0){ 
    iq = which(cc == xiaoxi2.med$name[ind] )
    iqn = which(is.finite(xiaoxi[,iq]))
    outputdata$XiaoxiN[i] = length(iqn)
  }
}

ind = which(is.finite(outputdata$XiaoxiEF))
outputdata=outputdata[ind,]

print('Fig Akagi')
# --------- Get Akagi/Andreae emission factors ------------
f2 = '/Users/ktravis1/OneDrive - NASA/FIREX/FinalAnalysisForGithub/InputFiles/OtherStudies/Andreae-BB-EMFactors-14Apr2021_justtable1.csv'
andreae = read.csv(f2)
require(readxl)
akagi=readxl::read_xlsx('/Users/ktravis1/OneDrive - NASA/FIREX/FinalAnalysisForGithub/InputFiles/OtherStudies/Akagi_acp-11-4039-2011-supplement/Tables 1-5_4.27.11.xlsx')
outputdata$AkagiName = outputdata$names
outputdata$AndreaeName = outputdata$Names
outputdata$AkagiEF = NaN; outputdata$AkagiSD=NaN
outputdata$AndreaeEF = NaN; outputdata$AndreaeSD=NaN

#ind = which(outputdata$names == 'OA')
#outputdata$names[ind] = 'OC'
fix=c()
for (i in 1:length(outputdata$names)){
  ind = which(outputdata$names[i] == akagi$Species)
  indB = which(outputdata$names[i] == andreae$Katie | outputdata$formula[i] == andreae$Katie | outputdata$formula[i] == andreae$Species)
  if (length(ind) == 0 ){fix=c(fix,outputdata$names[i])}
  
  if (length(ind) >2){print(">2")}
  if (length(ind) == 1){
    outputdata$AkagiName[i] = akagi$Species[ind]
    outputdata$AkagiEF[i] = akagi$CropEF[ind]
    outputdata$AkagiSD[i] = akagi$CropSD[ind]
    
    print(c(outputdata$names[i], akagi$AkagiName[ind]))
  }
  if (length(indB) >= 1 ){
   # outputdata$AkagiName[i] = akagi$Species[ind]
    outputdata$AndreaeEF[i] =andreae$average[indB[1]]
    outputdata$AndreaeSD[i] = andreae$std.dev.[indB[1]]
  }
}
#write.csv(outputdata, file='Xiaoxi_Matchup.csv')
outputdata$AndreaeEF = as.numeric(outputdata$AndreaeEF)
ind = which(outputdata$names != 'Nitrogen oxide' & outputdata$names != 'Nitrogen dioxide')
outputdata = outputdata[ind,]
outputdata = outputdata[order(outputdata$XiaoxiEF,decreasing = TRUE),]
outputdata$variableNUM = seq(1,length(outputdata$names))
outputdata$AkagiEF=as.numeric(outputdata$AkagiEF)
outputdata$AkagiSD=as.numeric(outputdata$AkagiSD)
outputdata$AndreaeEF=as.numeric(outputdata$AndreaeEF)
outputdata$AndreaeSD=as.numeric(outputdata$AndreaeSD)
ind = which(outputdata$names != 'Nitrate') # just using HNO3_NO3
outputdata = outputdata[ind,]
ind = which(outputdata$names == 'HNO3_NO3')
outputdata$names[ind] = 'Nitrate'
ind = which(outputdata$names == 'sum of monoterpenes')
outputdata$names[ind] = 'Monoterpenes'

ind1 = which(outputdata$kind == 'aerosol' | outputdata$kind == 'NOy' | outputdata$kind == 'Nitrogen-containing' | outputdata$kind == 'sulfur' | outputdata$kind == 'nitrogen')
ind2 = which(outputdata$kind == 'oVOC' | outputdata$kind == 'CO' | outputdata$kind == 'CO2' | 
               outputdata$kind == 'aromatic' | outputdata$kind == 'alkene' | outputdata$kind == 'CH2O')
tmp1 = outputdata[ind1,]; tmp2 = outputdata[ind2,]
tmp1$variableNUM = seq(1,length(tmp1$variableNUM));tmp2$variableNUM = seq(1,length(tmp2$variableNUM))
tmp1$AndreaeEF1 = tmp1$AndreaeEF - tmp1$AndreaeSD
tmp1$AkagiEF1 = tmp1$AkagiEF - tmp1$AkagiSD

ind = which(tmp1$AndreaeEF1 < 0)
if (length(ind)> 0){ tmp1$AndreaeEF1[ind] = 0.001}
ind = which(tmp1$AkagiEF1 < 0)
if (length(ind)> 0){ tmp1$AkagiEF1[ind] = 0.001}

tmp2$AndreaeEF1 = tmp2$AndreaeEF - tmp2$AndreaeSD
ind = which(tmp2$AndreaeEF1 < 0)
if (length(ind)> 0){ tmp2$AndreaeEF1[ind] = 0.001}

EFcomparisonInOrg=ggplot(tmp1) + 
  geom_point(aes(x=as.numeric(FinalEF_mean_ag),y=variableNUM), size=4,col='blue')+theme_classic()+ grids(linetype = "dashed")+
  geom_errorbar(aes(xmin=as.numeric(FinalEF_mean_ag)-as.numeric(FinalEF_sd_ag), xmax=as.numeric(FinalEF_mean_ag)+as.numeric(FinalEF_sd_ag),y=variableNUM), width=.2,size=2,
                position=position_dodge(.9), col='blue') +
  
  geom_point(aes(x=XiaoxiEF,y=variableNUM+0.1), size=4,col='black', show.legend=FALSE)+
  geom_errorbar(aes(xmin=XiaoxiEF-Xiaoxisd, xmax=XiaoxiEF+Xiaoxisd,y=variableNUM+0.1), width=.2,size=2,
                position=position_dodge(.9), col='black') +
  scale_x_continuous(trans='log10', limits=c(1E-03,1E+2))+
  
  geom_point(aes(x=AkagiEF,y=variableNUM+0.2), size=4,col='orange', show.legend=FALSE)+
  geom_errorbar(aes(xmin=(AkagiEF1), xmax=(AkagiEF+AkagiSD),y=variableNUM+0.2), width=.2,size=2,
                position=position_dodge(.9),col='orange') +
  
  geom_point(aes(x=AndreaeEF,y=variableNUM+0.4), size=4,col='purple', show.legend=FALSE)+
  geom_errorbar(aes(xmin=(AndreaeEF1), xmax=(AndreaeEF+AndreaeSD),y=variableNUM+0.4), width=.2,size=2,
                position=position_dodge(.9),col='purple') +
  #scale_color_manual(values = c("purple", "blue", "green","yellow","red"))+
  #guides(fill=guide_legend(title="nCs"))+
  theme(text = element_text(size=20))+
  geom_text(aes(x=0.01,y =variableNUM,label = names), size=5,
            vjust = 0,  hjust = 1,angle = 0, nudge_y = -.01, nudge_x = 0.2)+xlab(expression(paste('EF, g kg'^-1))) 
EFcomparisonOrg=ggplot(tmp2) + 
  geom_point(aes(x=as.numeric(FinalEF_mean_ag),y=variableNUM), size=4,col='blue')+theme_classic()+ grids(linetype = "dashed")+
  geom_errorbar(aes(xmin=as.numeric(FinalEF_mean_ag)-as.numeric(FinalEF_sd_ag), xmax=as.numeric(FinalEF_mean_ag)+as.numeric(FinalEF_sd_ag),y=variableNUM), width=.2,size=2,
                position=position_dodge(.9), col='blue') +
  
  geom_point(aes(x=XiaoxiEF,y=variableNUM+0.1), size=4,col='black', show.legend=FALSE)+
  geom_errorbar(aes(xmin=XiaoxiEF-Xiaoxisd, xmax=XiaoxiEF+Xiaoxisd,y=variableNUM+0.1), width=.2,size=2,
                position=position_dodge(.9), col='black') +
  scale_x_continuous(trans='log10', limits=c(1E-04,1E+4))+
  
  geom_point(aes(x=AkagiEF,y=variableNUM+0.2), size=4,col='orange', show.legend=FALSE)+
  geom_errorbar(aes(xmin=(AkagiEF-AkagiSD), xmax=(AkagiEF+AkagiSD),y=variableNUM+0.2), width=.2,size=2,
                position=position_dodge(.9),col='orange') +
  
  geom_point(aes(x=AndreaeEF,y=variableNUM+0.4), size=4,col='purple', show.legend=FALSE)+
  geom_errorbar(aes(xmin=(AndreaeEF1), xmax=(AndreaeEF+AndreaeSD),y=variableNUM+0.4), width=.2,size=2,
                position=position_dodge(.9),col='purple') +
  #scale_color_manual(values = c("purple", "blue", "green","yellow","red"))+
  #guides(fill=guide_legend(title="nCs"))+
  theme(text = element_text(size=20))+
  geom_text(aes(x=0.01,y =variableNUM,label = names), size=5,
            vjust = 0,  hjust = 1,angle = 0, nudge_y = -.01, nudge_x = 0.2)+xlab(expression(paste('EF, g kg'^-1))) 

# --- Normalize all data to the mean Liu et al value -----
xiaoxiNormalizeInorg = tmp1
xiaoxiNormalizeInorg = subset(xiaoxiNormalizeInorg, select = -c( FinalEF_mean_corn,FinalEF_sd ,COUNT_EFFINAL, FinalEF_mean_rice,   FinalEF_sd_rice,
                                                                 COUNT_EFFINAL_rice,  FinalEF_mean_soy,     FinalEF_sd_soy,
                                                                 COUNT_EFFINAL_soy,FinalEF_mean_slash ,  FinalEF_sd_slash,
                                                                 COUNT_EFFINAL_slash, FinalEF_mean_pile,   FinalEF_sd_pile,COUNT_EFFINAL_pile,
                                                                 FinalEF_mean_grass,FinalEF_sd_grass, COUNT_EFFINAL_grass) )
# just trying to show rangesm - forget it

xiaoxiNormalizeInorg$FinalEF_sd_ag = as.numeric(tmp1$FinalEF_sd_ag)/(as.numeric(tmp1$FinalEF_mean_ag))
xiaoxiNormalizeInorg$AkagiSD       = (tmp1$AkagiSD/tmp1$AkagiEF)
xiaoxiNormalizeInorg$AndreaeSD     = (tmp1$AndreaeSD/tmp1$AndreaeEF)
xiaoxiNormalizeInorg$Xiaoxisd      = as.numeric(tmp1$Xiaoxisd)/as.numeric(tmp1$XiaoxiEF)

xiaoxiNormalizeOrg = tmp2
xiaoxiNormalizeOrg$FinalEF_sd_ag = as.numeric(tmp2$FinalEF_mean_ag)*as.numeric(tmp2$FinalEF_sd_ag)/
  (as.numeric(tmp2$FinalEF_mean_ag))
xiaoxiNormalizeOrg$AkagiSD       = xiaoxiNormalizeOrg$AkagiEF*(tmp2$AkagiSD/tmp2$AkagiEF)
xiaoxiNormalizeOrg$AndreaeSD     = xiaoxiNormalizeOrg$AndreaeEF*(tmp2$AndreaeSD/tmp2$AndreaeEF)
xiaoxiNormalizeOrg$Xiaoxisd      = as.numeric(tmp2$XiaoxiEF) * (as.numeric(tmp2$Xiaoxisd)/
                                                                  as.numeric(tmp2$XiaoxiEF))

EFcomparisonNormInOrg=ggplot(xiaoxiNormalizeInorg)  + #scale_x_continuous(limits=c(0,6))+
  geom_point(aes(x=rep(1,length(ind1)),y=variableNUM+0.2 ),size=4,col='blue', show.legend=FALSE)+
  geom_errorbar(aes(xmin=as.numeric(1-FinalEF_sd_ag), xmax=as.numeric(1+FinalEF_sd_ag),y=variableNUM+0.2), col='blue',width=.2,size=2,
                position=position_dodge(.9)) +
  #scale_x_continuous(trans='log10')+
  
  geom_point(aes(x=rep(1,length(ind1)),y=variableNUM), size=4,col='black')+theme_classic()+ grids(linetype = "dashed")+
  geom_errorbar(aes(xmin=as.numeric(1-Xiaoxisd), xmax=as.numeric(1+Xiaoxisd),y=variableNUM), width=.2,size=2,
                position=position_dodge(.9), col='black') +

  theme(text = element_text(size=20))+
  geom_text(aes(x=-.1,y =variableNUM,label = names), size=5,
            vjust = 0,  hjust = 1,angle = 0, nudge_y = -.01, nudge_x = 0.2)+xlab(expression(paste('Normalized EF'))) +
#  annotate("text",x=4,y=3,label=paste("Andreae, 2019"), col='green')+
#  annotate("text",x=4,y=4,label=paste("Akagi et al., 2011"), col='orange')+
  annotate("text",x=4,y=5,label=paste("Xiaoxi et al., 2016"), col='black')+
  annotate("text",x=4,y=6,label=paste("FIREX-AQ (this study)"), col='blue')

EFcomparisonNormOrg=ggplot(xiaoxiNormalizeOrg) + scale_x_continuous(limits=c(0,6))+
  geom_point(aes(x=rep(1,length(ind2)),y=variableNUM+0.2 ),size=4,col='blue', show.legend=FALSE)+
  geom_errorbar(aes(xmin=as.numeric(V7), xmax=as.numeric(V8),y=variableNUM+0.2), col='blue',width=.2,size=2,
                position=position_dodge(.9)) +
  geom_point(aes(x=rep(1,length(ind2)),y=variableNUM),size=4, col='black')+theme_classic()+ grids(linetype = "dashed")+
  geom_errorbar(aes(xmin=as.numeric(Xiaoxi25), xmax=as.numeric(Xiaoxi75),y=variableNUM), width=.2,size=2,
                position=position_dodge(.9), col='black') +
  theme(text = element_text(size=20))+
  geom_text(aes(x=-.1,y =variableNUM,label = names), size=5,
            vjust = 0,  hjust = 1,angle = 0, nudge_y = -.01, nudge_x = 0.2)+xlab(expression(paste('Normalized EF'))) +
  annotate("text",x=4,y=5,label=paste("Xiaoxi et al., 2016"), col='black')+
  annotate("text",x=4,y=6,label=paste("FIREX-AQ (this study)"), col='blue')
#ggsave('EFComparisonNorm.pdf',EFcomparisonNorm,width = ww*1, height=ww*1/aspR)
tmp=ggarrange(EFcomparisonInOrg,EFcomparisonOrg)#,EFcomparisonNormInOrg, EFcomparisonNormOrg)
ggsave('EFComparisonBoth4.pdf',tmp,width = ww*2, height=ww*1/aspR)
ind = which(is.finite(outputdata$XiaoxiEF))
write.csv(outputdata[ind,],'Fig13.csv')
# ---- Figure S1 C9 aromatics --------
print("Fig. S1 C9 arom")
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

# ---- Figure S2 Monoterpenes ------
print("Fig. S2 Monoterpenes ")
tmp = as.data.frame(cbind(toga.all$aPinene_ppt,toga.all$bPineneMyrcene_ppt , toga.all$Camphene_ppt ,toga.all$Tricyclene_ppt))
ind = which(is.finite(toga.all$bPineneMyrcene_ppt) & is.finite(toga.all$Camphene_ppt)) # don't include outliers in slope
plot( toga.all$Monoterpenes_NOAAPTR_TM[ind],rowSums(tmp[ind,], na.rm=TRUE), pch=19, xlab='Monoterpenes, ppt', ylab='aPinene+bPinene+Myrcene+Camphene+Tricyclene, ppt')
xx=as.numeric(toga.all$Monoterpenes_NOAAPTR_TM)
yy=as.numeric(rowSums(tmp, na.rm=TRUE))

xx=xx[ind]; yy=yy[ind]
tt = lm(yy~xx+0)
abline(tt, col='black')
text(500,800,paste("TOGA slope = ", round(tt$coefficients[1], digits = 2)))
# ----- Figure S3 TOGA/WAS meacetate -----
print("Fig. S3 MeAcetate ")
plot( toga.all$C3H6O2_NOAAPTR_TM,toga.all$MeAcetate_ppt, pch=19,ylab='Methyl acetate, ppt', xlab='C3H6O2 (NOAA PTRMS), ppt')
xx=toga.all$C3H6O2_NOAAPTR_TM
yy=as.numeric(toga.all$MeAcetate_ppt)
tt = lm(yy~xx+0)
abline(tt, col='black')
text(4000,6000,paste("TOGA slope = ", round(tt$coefficients[1], digits = 2)))

#points(toga.all$C3H6O2_NOAAPTR_TM, toga.all$MeAcetate_WAS_TM, pch=19, col='red')
#yy=as.numeric(toga.all$MeAcetate_WAS_TM)
#tt2 = lm(yy~xx+0)
#abline(tt2, col='red')
#text(4000,5600,paste("WAS slope = ", round(tt2$coefficients[1], digits = 2)), col='red')

# ---- Figure S4 TOGA pyrrole -------
print('Fig S4, pyrrole')
plot( toga.all$C4H5N_NOAAPTR_TM, toga.all$Pyrrole_ppt, pch=19, xlab='C4H5N (NOAA PTRMS), ppt', ylab='Pyrrole, ppt')
xx=as.numeric(toga.all$C4H5N_NOAAPTR_TM)
yy=as.numeric(toga.all$Pyrrole_ppt)
ind = which(xx > 400 & yy < 100) # don't include outliers in slope
xx[ind]= NaN; yy[ind] = NaN
tt = lm(yy~xx+0)
abline(tt, col='black')
text(400,800,paste("TOGA slope = ", round(tt$coefficients[1], digits = 2)))

# Isoprene ---
print('Isoprene')
plot(toga.all$Isoprene_NOAAPTR_TM, toga.all$Isoprene_ppt, ylab='Isoprene, ppt', xlab='NOAA PTRMS Isoprene, ppt', pch=19, xlim=c(0,10E3), ylim=c(0,10E3))
tt = lm(toga.all$Isoprene_ppt ~toga.all$Isoprene_NOAAPTR_TM)
abline(tt, col='grey')
abline(0,1, lty=1)
points(toga.all$Isoprene_NOAAPTR_TM, toga.all$Isoprene_NOAAiWAS_TM, col='red', pch=19)
points(toga.all$Isoprene_NOAAPTR_TM, toga.all$Isoprene_WAS_TM, col='blue', pch=19)

# furan, butanal, isobutanal, methylfurans, 
plot(toga.all$Furan_NOAAPTR_TM, toga.all$Furan_ppt, xlab='NOAA PTR-ToF-MS Furan, ppt', ylab='Furan, ppt',xlim=c(0,max(toga.all$Furan_ppt, na.rm=TRUE)), ylim=c(0,max(toga.all$Furan_ppt, na.rm=TRUE)))
points(toga.all$Furan_NOAAPTR_TM, toga.all$Furan_NOAAiWAS_TM, col='red', pch=19)
points(toga.all$Furan_NOAAPTR_TM, toga.all$Furan_WAS_TM, col='blue', pch=19)
#abline(0,1, lty=1)
tt =lm(as.numeric(toga.all$Furan_ppt) ~ toga.all$Furan_NOAAPTR_TM)
abline(tt)
text(3e3,2e4,paste("TOGA slope = ", round(tt$coefficients[2], digits = 2)))
tt =lm(as.numeric(toga.all$Furan_NOAAiWAS_TM) ~ toga.all$Furan_NOAAPTR_TM)
abline(tt, col='red')
text(3e3,1.8e4,paste("iWAS slope = ", round(tt$coefficients[2], digits = 2)), col='red')
tt =lm(as.numeric(toga.all$Furan_WAS_TM) ~ toga.all$Furan_NOAAPTR_TM)
abline(tt, col='blue')
text(3e3,1.6e4,paste("WAS slope = ", round(tt$coefficients[2], digits = 2)), col='blue')

# ] methylfurans, 
plot(toga.all$x2MeFuranx3MeFuran_NOAAPTR_TM, toga.all$x2MeFuran_ppt+toga.all$x3MeFuran_ppt, xlab='NOAA PTR-ToF-MS Methylfuran, ppt',
     ylab='Methylfuran, ppt',xlim=c(0,max(toga.all$x2MeFuranx3MeFuran_NOAAPTR_TM, na.rm=TRUE)), ylim=c(0,max(toga.all$x2MeFuranx3MeFuran_NOAAPTR_TM, na.rm=TRUE)))
points(toga.all$x2MeFuranx3MeFuran_NOAAPTR_TM, as.numeric(toga.all$x2MeFuran_WAS_TM)+as.numeric(toga.all$x3MeFuran_WAS_TM), col='red', pch=19)
#abline(0,1, lty=1)
tt =lm(as.numeric(toga.all$Furan_ppt) ~ toga.all$Furan_NOAAPTR_TM)
abline(tt)
text(3e3,2e4,paste("TOGA slope = ", round(tt$coefficients[2], digits = 2)))
tt =lm(as.numeric(toga.all$Furan_NOAAiWAS_TM) ~ toga.all$Furan_NOAAPTR_TM)
abline(tt, col='red')
abline(tt, col='blue')
text(3e3,1.6e4,paste("WAS slope = ", round(tt$coefficients[2], digits = 2)), col='blue')

# furan, butanal, isobutanal, methylfurans, 
plot(toga.all$Butanal_ppt, as.numeric(toga.all$Butanal_WAS_TM), xlab='TOGA Butanal, ppt', ylab='WAS Butanal, ppt',xlim=c(0,max(toga.all$Butanal_ppt, na.rm=TRUE)),
     ylim=c(0,max(toga.all$Butanal_ppt, na.rm=TRUE)))
#abline(0,1, lty=1)
tt =lm(as.numeric(toga.all$Butanal_WAS_TM) ~ toga.all$Butanal_ppt)
abline(tt)
text(300,1000,paste("WAS slope = ", round(tt$coefficients[2], digits = 2)))

# furan, butanal, isobutanal, methylfurans, 
plot(toga.all$iButanal_ppt, as.numeric(toga.all$iButanal_WAS_TM), xlab='TOGA iButanal, ppt', ylab='WAS iButanal, ppt',xlim=c(0,max(toga.all$iButanal_ppt, na.rm=TRUE)),
     ylim=c(0,max(toga.all$iButanal_ppt, na.rm=TRUE)))
#abline(0,1, lty=1)
tt =lm(as.numeric(toga.all$iButanal_WAS_TM) ~ toga.all$iButanal_ppt)
abline(tt)
text(300,1000,paste("WAS slope = ", round(tt$coefficients[2], digits = 2)))

# ---- Fig. S7 -----
allBOTH.filter$fuel = allBOTH.filter$fuelORIG
ind = which(allBOTH.filter$variable == 'BC_SCHWARZ')
tmpBC = allBOTH.filter[ind,]
ind = which(allBOTH.filter$variable == 'OC_JIMENEZ')
tmpOC = allBOTH.filter[ind,]
tmpBC$FinalERtoCO = tmpBC$FinalERtoCO/tmpOC$FinalERtoCO
tmpBC$variable = 'BCtoOC'
ECvsMCE_ER = plotSpeciesCOMCE(tmpBC,'BCtoOC','BC/OC Ratio','Black Carbon2','BC2')
ggsave('ECvsMCE_ER.pdf',ECvsMCE_ER,width = ww*1, height=ww*1/aspR)

# ---
xx=toga.all$C4Carbonyls_NOAAPTR_TM/1E3
yy=(as.numeric(toga.all$MEK_NOAAiWAS_TM )+  as.numeric(toga.all$iButanal_ppt) +  as.numeric(toga.all$Butanal_ppt))/1E3
plot(toga.all$C4Carbonyls_NOAAPTR_TM/1E3, (as.numeric(toga.all$MEK_NOAAiWAS_TM )+  as.numeric(toga.all$iButanal_ppt) +  as.numeric(toga.all$Butanal_ppt))/1E3, xlim=c(0,11), ylim=c(0,11))
abline(0,1)
tt = lm(yy~xx)
abline(tt, lty=2)


# Ag vs. prescribed CH4 related to aging?

ind = which(allBOTH.filter$variable == 'CH4_DACOM_DISKIN' & allBOTH.filter$fuel2 == 'agriculture' & is.finite(allBOTH.filter$MAtoF.5hz))
ag = allBOTH.filter[ind,]

ind = which(allBOTH.filter$variable == 'CH4_DACOM_DISKIN' & allBOTH.filter$fuel2 == 'prescribed' & is.finite(allBOTH.filter$MAtoF.5hz))
presc = allBOTH.filter[ind,]
t.test(ag$MAtoF.5hz, presc$MAtoF.5hz, na.rm=TRUE)

