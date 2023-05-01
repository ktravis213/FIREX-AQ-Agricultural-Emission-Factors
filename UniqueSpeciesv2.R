#source("SlopeDifferences.R")
require(tidyverse); require(ggbreak)

# Which species have a large difference between agricultural and prescribed
ind = which(allBOTH.filter$fuel2 == 'agriculture' & allBOTH.filter$USEME != -1)
allBOTH.filter.ag = allBOTH.filter[ind,]

# Make average pile/slash
ind = which(allBOTH.filter$fuel2 == 'prescribed'& allBOTH.filter$USEME != -1 & allBOTH.filter$fuelORIG != 'Blackwater')
allBOTH.filter.presc = allBOTH.filter[ind,]

# Get rid of things we don't want
ind = which(allBOTH.filter$PI !='NOx/NOy' & allBOTH.filter$PI !='PAN/NOy' & allBOTH.filter$PI != 'MAtoF' & 
              allBOTH.filter$USEME != -1 & is.finite(allBOTH.filter$USEME))
varsDIFFS = as.data.frame(unique(allBOTH.filter$variable[ind]))
colnames(varsDIFFS) = 'vars'
# p values
varsDIFFS$AgvsPresc = NaN; varsDIFFS$AgCorrvsPrescCorr = NaN; varsDIFFS$AgCorrvsPresc = NaN; varsDIFFS$AgvsPrescCorr = NaN 
# averages
varsDIFFS$AvgAg = NaN;varsDIFFS$AvgPresc = NaN; varsDIFFS$AvgAgCorr = NaN;varsDIFFS$AvgPrescCorr = NaN
varsDIFFS$Category = NaN; varsDIFFS$OHrate = NaN; varsDIFFS$MWs = NaN; varsDIFFS$kind = ''
minD = 10
for (i in 1:length(varsDIFFS$vars)){
  print(varsDIFFS$vars[i])
  
  ind = which(allBOTH.filter.ag$variable == varsDIFFS$vars[i] & is.finite(allBOTH.filter.ag$FinalEF))
  ind2 = which(allBOTH.filter.presc$variable == varsDIFFS$vars[i]& is.finite(allBOTH.filter.presc$FinalEF))
  # corrected to MCE=0.92
  ind4 = which(allBOTH.filter.ag$variable == varsDIFFS$vars[i]& is.finite(allBOTH.filter.ag$FinalEF_MCE92))
  ind5 = which(allBOTH.filter.presc$variable == varsDIFFS$vars[i]& is.finite(allBOTH.filter.presc$FinalEF_MCE92))
  
  varsDIFFS$AvgAg[i]        = mean(allBOTH.filter.ag$FinalEF[ind])
  varsDIFFS$AvgPresc[i]     = mean(allBOTH.filter.presc$FinalEF[ind2])
  varsDIFFS$AvgAgCorr[i]    = mean(allBOTH.filter.ag$FinalEF_MCE92[ind])
  varsDIFFS$AvgPrescCorr[i] = mean(allBOTH.filter.presc$FinalEF_MCE92[ind2])
  
  varsDIFFS$Category[i] = allBOTH.filter.ag$Category[ind[1]]
  varsDIFFS$MWs[i] = allBOTH.filter.ag$mWs[ind[1]]
  varsDIFFS$kind[i] = allBOTH.filter.ag$kind[ind[1]]
  varsDIFFS$OHrate[i] = max(c(allBOTH.filter.ag$OHrate.1hz[ind[1]],allBOTH.filter.ag$OHrate.5hz[ind[1]]), na.rm=TRUE)
  # Is ag unique from presc
  if (length(ind) > minD & length(ind2) > minD & length(ind2)/length(ind) > 0.1){
    tmp= t.test(allBOTH.filter.ag$FinalEF[ind], allBOTH.filter.presc$FinalEF[ind2]) # Is ag unique from presc
   if (tmp$p.value < 0.05){ varsDIFFS$AgvsPresc[i] = tmp$p.value }
  }
  # Is corrected ag unique from presc
  if (length(ind4) > minD & length(ind2) > minD & length(ind2)/length(ind4) > 0.1){
    tmp= t.test(allBOTH.filter.ag$FinalEF_MCE92[ind4], allBOTH.filter.presc$FinalEF[ind2]) # Corrected ag vs. presc
    if (tmp$p.value < 0.05){varsDIFFS$AgCorrvsPresc[i] = tmp$p.value}
  }
  # Is ag unique from corr presc
  if (length(ind) > minD & length(ind5) > minD & length(ind5)/length(ind) > 0.1){
    tmp= t.test(allBOTH.filter.ag$FinalEF[ind], allBOTH.filter.presc$FinalEF_MCE92[ind5]) # Corrected presc vs. ag
    if (tmp$p.value < 0.05){varsDIFFS$AgvsPrescCorr[i] = tmp$p.value}
  }
  # Is corrected ag unique from corrected presc
  if (length(ind4) > minD & length(ind5) > minD & length(ind5)/length(ind4) > 0.1){
    tmp= t.test(allBOTH.filter.ag$FinalEF_MCE92[ind4], allBOTH.filter.presc$FinalEF_MCE92[ind5]) # Corrected presc vs. corrected ag
    if (tmp$p.value < 0.05){varsDIFFS$AgCorrvsPrescCorr[i] = tmp$p.value}
  }
}
varsDIFFS$names = ''; varsDIFFS$LifetimeCat = NaN
varsDIFFS$PI = ''; varsDIFFS$USEME=NaN
for (i in 1:length(varsDIFFS$vars)){
  ind = which(varsDIFFS$vars[i] == allBOTH.filter$variable)
  varsDIFFS$PI[i] = allBOTH.filter$PI[ind[1]]
  varsDIFFS$names[i] = allBOTH.filter$names[ind[1]]
  varsDIFFS$USEME[i] = allBOTH.filter$USEME[ind[1]]
  varsDIFFS$LifetimeCat[i] = allBOTH.filter$LifetimeCat[ind[1]]
}

# ----- Start with corrected ag and corrected presc -------
varsDIFFS$FinalAg = varsDIFFS$AvgAgCorr
varsDIFFS$FinalPresc = varsDIFFS$AvgPrescCorr
varsDIFFS$Finalpval = varsDIFFS$AgCorrvsPrescCorr
# ----- Both uncorrected? --------
ind = which(!is.finite(varsDIFFS$FinalAg) & !is.finite(varsDIFFS$FinalPresc))
varsDIFFS$Finalpval[ind] = varsDIFFS$AgvsPresc[ind]
# ----- Just Presc uncorrected? --------
ind = which(is.finite(varsDIFFS$FinalAg) & !is.finite(varsDIFFS$FinalPresc))
varsDIFFS$Finalpval[ind] = varsDIFFS$AgCorrvsPresc[ind]
# ----- Just Ag uncorrected? --------
ind = which(!is.finite(varsDIFFS$FinalAg) & is.finite(varsDIFFS$FinalPresc))
varsDIFFS$Finalpval[ind] = varsDIFFS$AgvsPrescCorr[ind]

# use uncorrected Presc if no corrected
varsDIFFS$FinalAg[ind] = varsDIFFS$AvgAg[ind]
varsDIFFS$FinalPresc[ind] = varsDIFFS$AvgPresc[ind]
ind = which(!is.finite(varsDIFFS$FinalPresc))
varsDIFFS$FinalPresc[ind] = varsDIFFS$AvgPresc[ind]
ind = which(!is.finite(varsDIFFS$FinalAg))
varsDIFFS$FinalAg[ind] = varsDIFFS$AvgAg[ind]


varsDIFFS$Diff = abs((varsDIFFS$FinalAg-varsDIFFS$FinalPresc)*100/varsDIFFS$FinalPresc)
varsDIFFS$Ratio = varsDIFFS$FinalAg/varsDIFFS$FinalPresc
# ***** # ***** Ok, need to make consistency now between Table 4, where averages are weighted, and here. So replace ratios with Table 4 ratios but only for species not corrected for MCE.# ***** # ***** 
outputdata = as.data.frame(cbind(new.data.frame$Category, new.data.frame$LifetimeCat,new.data.frame$mWs,
                                 new.data.frame$names,new.data.frame$formula,new.data.frame$PI,
                                 new.data.frame.ag$FinalEF_mean, new.data.frame.ag$FinalEF_sd, new.data.frame.ag$COUNT_EFFINAL,  
                                 new.data.frame.presc$FinalEF_mean,   new.data.frame.presc$FinalEF_sd,         new.data.frame.presc$COUNT_EFFINAL))

colnames(outputdata) = c("Category", 'LifetimeCat','mWs',  'names','formula','PI', 
                         'FinalEF_mean_ag',      'FinalEF_sd_ag',    'COUNT_EFFINAL_ag',
                         'FinalEF_mean_presc',   'FinalEF_sd_presc',    'COUNT_EFFINAL_presc')

ind = which(varsDIFFS$names == 'Methyl Bromide')
varsDIFFS$names[ind] = 'Methyl bromide'
varsDIFFS$TableRatio = NaN; varsDIFFS$TableDiff = NaN
for (i in 1:length(varsDIFFS$Category)){
  # So only if we DIDNT correct for MCE replace with Table ratio
  if (is.nan(varsDIFFS$AvgAgCorr[i]) & is.nan(varsDIFFS$AvgPrescCorr[i])){
    ind = which(outputdata$names == varsDIFFS$names[i])
    if (length(ind) > 1  & varsDIFFS$names[i] == 'Phenol'){ind = ind[2]} # Phenol, only want 2nd one which is Wennberg
    if (length(ind) > 1  & varsDIFFS$names[i] == 'Methyl acetate'){ind = ind[1]} # Methyl acetate, only want 1st one which is Toga
    if (length(ind) > 0){varsDIFFS$TableRatio[i] = as.numeric(outputdata$FinalEF_mean_ag[ind])/as.numeric(outputdata$FinalEF_mean_presc[ind])}
    if (length(ind) > 0){varsDIFFS$TableDiff[i] = abs((as.numeric(outputdata$FinalEF_mean_ag[ind])-as.numeric(outputdata$FinalEF_mean_presc[ind]))*100/
                                                        as.numeric(outputdata$FinalEF_mean_presc[ind]))}
  }
}
# Use regular ratio if TableRatio is NaN
ind = which(is.nan(varsDIFFS$TableRatio))
varsDIFFS$TableRatio[ind] = varsDIFFS$Ratio[ind]
ind = which(is.nan(varsDIFFS$TableDiff))
varsDIFFS$TableDiff[ind] = varsDIFFS$Diff[ind]
# ***** # ***** # ***** # ***** # ***** # ***** # ***** # ***** # ***** # ***** # ***** # ***** # ***** # ***** # ***** # ***** # ***** # ***** # ***** # ***** 


indD = which(abs(varsDIFFS$TableDiff) >= 50)
# kludge for NaN pval
ind = which(!is.finite(varsDIFFS$Finalpval))
varsDIFFS$Finalpval[ind] = 1
indD = which(round(varsDIFFS$Finalpval,2) < 0.05 & abs(varsDIFFS$TableDiff) >= 50 )
test = varsDIFFS[indD,]
varsDIFFS$USEFORPLOT = 0
varsDIFFS$USEFORPLOT[indD] = 1
# species where p-values are not NaN
ind = which(varsDIFFS$USEME != 0)
write.csv(varsDIFFS[ind,], 'differentEFs.csv')
# ----- Fix names
test$kindORIG = test$kind
ind = which(test$names == 'sum of monoterpenes')
test$names[ind] = 'Monoterpenes'
ind = which(test$names == 'Methyl acetate/Ethyl formate/Hydroxyacetone')
test$names[ind] = 'Hydroxyacetone + Methyl acetate + Ethyl formate'
ind = which(test$names == '2,3-Butanedione/2-Oxobutanal/1,4-Butanedial')
test$names[ind] = '2,3-Butanedione + 2-Oxobutanal + 1,4-Butanedial'
# ------ Consolidate kinds
ind = which(test$Category == 1 & test$LifetimeCat == 1)
test$kind[ind] = 'Short-lived NMVOC'
ind = which(test$Category == 1 & test$LifetimeCat == 2)
test$kind[ind] = 'Long-lived NMVOC'
ind = which(test$Category == 2)
test$kind[ind] = 'Nitrogen-containing'
ind = which(test$names == 'Sulfur dioxide' | test$kind == 'sulfur')
test$kind[ind] = 'Sulfur-containing'
ind = which(test$kind == 'halocarbon' | test$kind == 'halogen' | test$kind == 'CH3Br_ppt')
test$kind[ind] = 'Halogen-containing species'
ind = which(test$kind == 'aerosol' | test$names == 'Methane')
test$kind[ind] = 'Other'
#ind = which(test$kind == 'OVOC' | test$kind == 'oVOC')
#test$kind2[ind] = 'NotNerd'
# Re-order
test = test[order(test$TableRatio,decreasing=TRUE),]
# ---- New order
test$kind2 = "NotNerd" # Not aerosol # DefNerd = aerosol
indA = which(test$names == 'Chloride'); test$kind2[indA] = "DefNerd"; test$kind[indA] = "Halogen-containing species"
indA = which(test$names == 'Potassium'); test$kind2[indA] = "DefNerd"; test$kind[indA] = "Other"
indA = which(test$names == 'Ammonium'); test$kind2[indA] = "DefNerd"; test$kind[indA] = "Other"
indA = which(test$names == 'PM1'); test$kind2[indA] = "DefNerd"; test$kind[indA] = "Other"
indA = which(test$names == 'Organic carbon'); test$kind2[indA] = "DefNerd"; test$kind[indA] = "Other"
ind = which(test$kindORIG == 'oVOC'); test$kind2[ind] = "MyNerd"

ind = which(test$names != '⍺-Pinene' & test$names != 'Camphene' &  test$names != 'beta-Pinene' &
              test$names != 'beta-Pinene/Myrcene' & test$names != 'Peroxyacryloyl nitrate' & 
              test$names != 'Peroxyacetyl nitrate' &
              test$names != 'CCN_034' & test$USEME > 0 & test$names != 'NOx (as NO)')
test=test[ind,]
test$NUM = seq(1,length(test$vars))
test<- test[order(-test$TableRatio),]

uPLOT2=ggplot(test, aes(x=NUM,y=TableRatio, label=names)) + geom_point(aes(shape=kind, col=kind), size=5) + geom_hline(aes(yintercept=1), lty=2)+
  theme_classic() +ylab('Ratio') +
  ylab('Ratio of Agricultural EF to Prescribed EF') +
  theme(axis.text.x = element_text(angle =45, vjust =1, hjust=1))+
  theme(text = element_text(size = 20)) + theme(legend.position = c(.9, 0.8)) +
 # ylim(c(0,10))+
  geom_text(aes(x =NUM,y = 0,label = names), size=5,
            vjust = 0,  hjust = 1,angle = 45, nudge_y = 0.01, nudge_x = 0.2)+
#  scale_pattern_manual(values = c(Nerd = "stripe", NotNerd = "none")) +
  scale_color_manual(values=as.vector(polychrome(26)))+ theme(legend.position = "top")+labs(fill="")

# ---- Split y-axis
  #scale_y_cut(breaks=c(8,20), )+
  #scale_y_cut(breaks=c(6,11))+
#  scale_y_break(c(11,13), scale=.05)+
#  scale_y_break(c(13,22), scale=.05)+
  #scale_y_break(c(1, 6), scale=1)+ 
  #scale_y_break(c(0, 1), scale=.1) +
  #scale_y_cut(breaks=c(.3,1,2,8,20),which=c(1,1,1,1,1),scales=c(.1,1,1,1,1),space=.2)+
# ---------- Figure 13
ggsave(uPLOT2, file='RatioPLOT.pdf',width = 7*1.25*1.75, height=7*1.25)

