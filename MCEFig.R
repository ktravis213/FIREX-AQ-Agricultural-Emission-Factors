# MCE dependence binned by 0.01
vars = unique(allBOTH.filter$variable)
MCE = seq(84,98,by=1)
MCE = MCE/100
fuels = c("agriculture","prescribed","grassland")
EFs = rep(NaN, length(MCE))
require(creditmodel)
bins=40
bin = cut_equal(allBOTH.filter$MCE, g = bins, sp_values = NULL, cut_bin = "equal_depth")
bin = bin[which(is.finite(bin))]
alt_group = ( floor(allBOTH.filter$FinalEF *2) )
for (i in 1:(length(bin))){
  if (i == 1){
    ind = which(allBOTH.filter$MCE <= bin[i])
    alt_group[ind] = i
    counts = length(ind)
   } else{
     ind = which(allBOTH.filter$MCE > bin[i-1] &allBOTH.filter$MCE <= bin[i])
     counts = c(counts, length(ind))
    alt_group[ind] = i
   }
}
# last bin
ind = which(allBOTH.filter$MCE > bin[length(bin)])
counts = c(counts, length(ind))
alt_group[ind] = i+1
allBOTH.filter$alt_group = alt_group

# Don't include blackwater
ind = which(allBOTH.filter$fuelORIG != 'Blackwater')
allBOTH.filter.bin = aggregate(allBOTH.filter[ind,], by=list(allBOTH.filter$variable[ind], allBOTH.filter$fuel2[ind], allBOTH.filter$alt_group[ind]), FUN='mean', na.rm=TRUE)
for (i in 1:length(allBOTH.filter.bin$Group.1)){
  ind = which(allBOTH.filter$variable==allBOTH.filter.bin$Group.1[i] )
  allBOTH.filter.bin$names[i] = allBOTH.filter$names[ind[1]]
  
}
varsALL = as.data.frame(cbind(vars))
varsALL$Rval_ag = NaN; varsALL$Slope_ag = NaN; varsALL$Slopestd_ag = NaN; varsALL$Intercept_ag = NaN
varsALL$Rval_presc = NaN; varsALL$Slope_presc = NaN; varsALL$Slopestd_presc = NaN; varsALL$Intercept_presc = NaN
varsALL$Rval_grass = NaN; varsALL$Slope_grass = NaN; varsALL$Slopestd_grass = NaN; varsALL$Intercept_grass = NaN

for (i in 1:length(varsALL$vars)){
  # ---------------------- Agriculture --------------------
  ind = which(allBOTH.filter.bin$Group.1 == varsALL$vars[i] & allBOTH.filter.bin$Group.2 == 'agriculture' & is.finite(allBOTH.filter.bin$FinalEF))
  if (length(ind) > 4){
    tmp = lm(allBOTH.filter.bin$FinalEF[ind]~allBOTH.filter.bin$MCE[ind])
    sL = summary(tmp)
    sQ = cor.test(allBOTH.filter.bin$FinalEF[ind],allBOTH.filter.bin$MCE[ind])
    if (sQ$p.value < 0.05){
      varsALL$Rval_ag[i]      = sQ$estimate
      varsALL$Slope_ag[i]     = tmp$coefficients[2]
      varsALL$Slopestd_ag[i]  = sL$coefficients[2,2]
      varsALL$Intercept_ag[i] = tmp$coefficients[1]
    }
  }
  # ---------------------- Prescribed --------------------
  ind = which(allBOTH.filter.bin$Group.1 == varsALL$vars[i] & allBOTH.filter.bin$Group.2 == 'prescribed' & is.finite(allBOTH.filter.bin$FinalEF))
  if (length(ind) > 4){
    tmp = lm(allBOTH.filter.bin$FinalEF[ind]~allBOTH.filter.bin$MCE[ind])
    sL = summary(tmp)
    sQ = cor.test(allBOTH.filter.bin$FinalEF[ind],allBOTH.filter.bin$MCE[ind])
    if (sQ$p.value < 0.05){
      varsALL$Rval_presc[i]      = sQ$estimate
      varsALL$Slope_presc[i]     = tmp$coefficients[2]
      varsALL$Slopestd_presc[i]  = sL$coefficients[2,2]
      varsALL$Intercept_presc[i] = tmp$coefficients[1]
    }
  }
  # ---------------------- Grassland --------------------
  ind = which(allBOTH.filter.bin$Group.1 == varsALL$vars[i] & allBOTH.filter.bin$Group.2 == 'grass' & is.finite(allBOTH.filter.bin$FinalEF))
  if (length(ind) > 4){
    tmp = lm(allBOTH.filter.bin$FinalEF[ind]~allBOTH.filter.bin$MCE[ind])
    sL = summary(tmp)
    sQ = cor.test(allBOTH.filter.bin$FinalEF[ind],allBOTH.filter.bin$MCE[ind])
    if (sQ$p.value < 0.05){
      varsALL$Rval_grass[i]      = sQ$estimate
      varsALL$Slope_grass[i]     = tmp$coefficients[2]
      varsALL$Slopestd_grass[i]  = sL$coefficients[2,2]
      varsALL$Intercept_grass[i] = tmp$coefficients[1]
    }
  }
}
varsALL$names = ''; varsALL$kind = ''; varsALL$USEME = NaN; varsALL$Category = NaN; varsALL$LifetimeCat = NaN
for (i in 1:length(varsALL$vars)){
  ind = which(allBOTH.filter$variable == varsALL$vars[i])
  varsALL$USEME[i] = max(allBOTH.filter$USEME[ind],na.rm=TRUE)
  varsALL$Category[i] = allBOTH.filter$Category[ind[1]]
  varsALL$LifetimeCat[i] = allBOTH.filter$LifetimeCat[ind[1]]
  varsALL$names[i] = allBOTH.filter$names[ind[1]]
  varsALL$kind[i] = allBOTH.filter$kind[ind[1]]
}

# ------- Plots ------
doplots=0
for (i in 1:length(varsALL$vars)){
  
  print(varsALL$vars[i])
  ind = which(allBOTH.filter$variable == vars[i] & allBOTH.filter$fuelORIG != 'forest' & 
                allBOTH.filter$fuelORIG != 'coniferous/decidous' & allBOTH.filter$fuelORIG != 'Blackwater' )
  fuelshapes = c(19,15,1)
  fuellimits = c("agriculture","grass","prescribed")
  
  cbp1 <- c( "#E69F00", "#56B4E9", "#009E73")
  allBOTH.filter.CH4 = allBOTH.filter[ind,]
  #ind = which(allBOTH.filter.CH4$fuel == 'corn' | allBOTH.filter.CH4$fuel == 'rice' | allBOTH.filter.CH4$fuel == 'soybean' |
  #              allBOTH.filter.CH4$fuel == 'winter wheat')
  #allBOTH.filter.CH4$fuel[ind] = 'agricultural'
  
  #ind = which(allBOTH.filter.CH4$fuel == 'pile' | allBOTH.filter.CH4$fuel == 'slash' | allBOTH.filter.CH4$fuel == 'shrub' )
  #allBOTH.filter.CH4$fuel[ind] = 'prescribed'
  allBOTH.filter.CH4$fuel2 <- factor(allBOTH.filter.CH4$fuel2,
                                     levels = c("agriculture","grass","prescribed","Blackwater"))
  
  xx = c(0.84,0.98)
  CH4vsMCE = ggplot(allBOTH.filter.CH4)+
    geom_point(aes(x=MCE,y=FinalEF, col=fuel2, shape=fuel2), stroke=2,size=3)+
    theme_classic()+
    xlab("MCE")+ylab(paste(vars[i],", g/kg",sep=''))+labs(col="",shape="")+
    theme(legend.background=element_blank())+
    # theme(legend.position = c(0.2, .6upSD*max(allBOTH.filter.CH4.avg$FinalEF, na.rm=TRUE)))+
    theme(text = element_text(size = 20)) +
    # theme(legend.direction="horizontal")+
    #scale_shape_manual(values=c(19,15,17,18,7,8,4))+
    scale_color_manual(values = c(cbp1 ))+
    scale_shape_manual(values =fuelshapes)+
    xlim(xx)
  CH4vsMCE = CH4vsMCE + 
    #        geom_text(x = 0.98, y=.4,label='Akagi et al., 2011 (Pasture maintenance)',hjust = 0, size = 6, col='light pink') +
    coord_cartesian(xlim = xx, # This focuses the x-axis on the range of interest
                    clip = 'off') +   # This keeps the labels from disappearing
    theme(plot.margin = unit(c(1,3,1,1), "lines")) +# This widens the right margin
    # annotate(geom="shape",shape="cross",  color="pink", fill="pink")+
    
    #  geom_text(x = 0.98, y=.35, label='Akagi et al., 2011 (Crop residue)',hjust = 0, size = 6, col='pink') +
    coord_cartesian(xlim = xx, # This focuses the x-axis on the range of interest
                    clip = 'off') +   # This keeps the labels from disappearing
    theme(plot.margin = unit(c(1,3,1,1), "lines")) # This widens the right margin
  
  scaleFUN <- function(x) sprintf("%.2f", x)
  
  #make all y-axis same decimal points
  CH4vsMCE=CH4vsMCE + scale_y_continuous(labels=scaleFUN)
  # Get fuel specific slopes
  # agr
  CH4vsMCE = CH4vsMCE + geom_abline(slope=varsALL$Slope_ag[i], intercept = varsALL$Intercept_ag[i], col= cbp1[1])
  # grass
  CH4vsMCE = CH4vsMCE + geom_abline(slope=varsALL$Slope_grass[i], intercept = varsALL$Intercept_grass[i], col= cbp1[2])
  # presc    
  CH4vsMCE = CH4vsMCE + geom_abline(slope=varsALL$Slope_presc[i], intercept = varsALL$Intercept_presc[i], col= cbp1[3])
  
  tmp = vars[i]
  tmp = str_replace(tmp,"/","")
  tmp = str_replace(tmp,"-","")
  if (doplots == 1){ggsave(CH4vsMCE,filename = paste0('FiguresMCE/',tmp,".pdf"))}
}

# for (i in 1:length(varsALL$vars)){
#   nnM = varsALL$names[i]
#   ind = which(allBOTH.filter$names == nnM & allBOTH.filter$USEME == 1 &
#                 allBOTH.filter$fuel2 == 'agriculture' & is.finite(allBOTH.filter$FinalEF))
#   if (length(ind) > 0){
#     plot(allBOTH.filter$MCE[ind], allBOTH.filter$FinalEF[ind], pch=19, ylab=paste0('EF ',nnM, ' g/kg'), xlab='MCE')
#     cor.test(allBOTH.filter$MCE[ind], allBOTH.filter$FinalEF[ind])
#     ind = which(allBOTH.filter.bin$names == nnM & allBOTH.filter.bin$USEME == 1 & 
#                   allBOTH.filter.bin$Group.2 == 'agriculture')
#     points(allBOTH.filter.bin$MCE[ind], allBOTH.filter.bin$FinalEF[ind], col='red', lwd=3, type='o')
#     cor.test(allBOTH.filter.bin$MCE[ind], allBOTH.filter.bin$FinalEF[ind])
#   }
# }
