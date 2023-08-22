allBOTH.filter.mean.fire = aggregate(allBOTH.filter, by=list(allBOTH.filter$fire,allBOTH.filter$fuel,allBOTH.filter$variable),   FUN='mean', na.rm=TRUE)
allBOTH.filter.sd.fire = aggregate(allBOTH.filter, by=list(allBOTH.filter$fire,allBOTH.filter$fuel,allBOTH.filter$variable), FUN='sd', na.rm=TRUE)
allBOTH.filter.mean.fire$FinalEF_sd = allBOTH.filter.sd.fire$FinalEF
allBOTH.filter.mean.fire$FinalERtoCO_sd = allBOTH.filter.sd.fire$FinalERtoCO

ind = which(allBOTH.filter.mean.fire$USEME > 0 & allBOTH.filter.mean.fire$Category != 5 & is.finite(allBOTH.filter.mean.fire$Category))
tmp = allBOTH.filter.mean.fire[ind,]
cc = colnames(tmp)
ind = which(cc == 'Group.1' | cc == 'Group.2' | cc == 'Group.3' | cc == 'FinalEF' | cc == 'FinalEF_sd' | cc == 'MCE' |
              cc == 'FinalERtoCO' | cc == 'FinalERtoCO_sd')
write.csv(tmp[,ind], 'ForGLENN.csv')
