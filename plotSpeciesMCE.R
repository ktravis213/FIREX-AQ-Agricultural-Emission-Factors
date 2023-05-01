# ---- Figure1: Plot CH4 EF vs. MCE for separated fuel types -----
plotSpeciesMCE = function(allBOTH.filterIN,variableN,akagiN,andreaeN,xiaoxiN){
  doavgshape = 0
  
  xiaoxi$fuel = c('Corn','Corn','Corn', 'Rice','Rice',
                  'Rice','Rice','Rice','Rice',
                  'Rice','Rice','Soybean',
                  'Popcorn','DblCrpWinWht_Soy',
                  'WoodyWetlands')
  ind = which(xiaoxi$fuel == 'Corn' | xiaoxi$fuel == 'Popcorn'); xiaoxi$fuel[ind] = 'corn'
  ind = which(xiaoxi$fuel == 'Rice'); xiaoxi$fuel[ind] = 'rice'
  ind = which(xiaoxi$fuel == 'Soybean'); xiaoxi$fuel[ind] = 'soybean'
  ind = which(xiaoxi$fuel == 'DblCrpWinWht_Soy'); xiaoxi$fuel[ind] = 'soybean'
  ind = which(xiaoxi$fuel == 'WoodyWetlands'); xiaoxi$fuel[ind] = 'slash'
  indA = which(akagi$Species == akagiN)
  indB = which(andreae$Katie == andreaeN)
  cc = colnames(xiaoxi)
  indC  = which(cc == xiaoxiN)
  print(c(indA,indB,indC))
  ind = which(allBOTH.filterIN$fire == 'BlackwaterRiver')
  allBOTH.filterIN$fuelORIG[ind] = 'BlackwaterRiver'
  ind = which(allBOTH.filterIN$variable == variableN & allBOTH.filterIN$fire != 'Vivian' &
                allBOTH.filterIN$fire != 'Copper Breaks' &
                allBOTH.filterIN$fire != 'Invictus' &
                allBOTH.filterIN$fuelORIG != 'coniferous/decidous')
  fuelshapes = c(19,0,2,18,7,8,4,6,1,12,13)
  fuellimits = c("corn","soybean","rice","winter wheat","grass","pile","slash","shrub","BlackwaterRiver")
  
  cbp1a <- c( "#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", 
              "#CC79A7","#000000","#999999","#CC0000")
  # colorblind safe option
  #cbp1b = c('#377eb8','#e41a1c','#4daf4a','#f781bf','#a65628','#984ea3','#ff7f00','#ffff33')
  cbp1b = c('#377eb8','#e41a1c','#4daf4a','#a65628','#984ea3','#ff7f00',"#000000")
  cbp1c=c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9',
          '#74add1','#4575b4','#313695')
  cbp1d <- c( "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000","#E69F00", "#56B4E9", "#009E73")
  
  cbp1 = cbp1b
  allBOTH.filter.CH4 = allBOTH.filterIN[ind,]
  allBOTH.filter.CH4.avg = aggregate(allBOTH.filter.CH4, by=list(allBOTH.filter.CH4$fuelORIG), FUN='median', na.rm=TRUE)
  #allBOTH.filter.CH4.sd = aggregate(allBOTH.filter.CH4, by=list(allBOTH.filter.CH4$fuelORIG), FUN='sd', na.rm=TRUE)
  allBOTH.filter.CH4.25 = aggregate(allBOTH.filter.CH4$FinalEF, by=list(allBOTH.filter.CH4$fuelORIG), 'quantile',probs=c(q1),na.rm=TRUE)
  allBOTH.filter.CH4.75 = aggregate(allBOTH.filter.CH4$FinalEF, by=list(allBOTH.filter.CH4$fuelORIG), 'quantile',probs=c(q2),na.rm=TRUE)
  allBOTH.filter.CH4.25.mce = aggregate(allBOTH.filter.CH4$MCE, by=list(allBOTH.filter.CH4$fuelORIG), 'quantile',probs=c(q1),na.rm=TRUE)
  allBOTH.filter.CH4.75.mce = aggregate(allBOTH.filter.CH4$MCE, by=list(allBOTH.filter.CH4$fuelORIG), 'quantile',probs=c(q2),na.rm=TRUE)
  allBOTH.filter.CH4.avg = getplumesANDmcebyfuel1VAR(allBOTH.filter.CH4.avg, allBOTH.filter.CH4 )
  
  allBOTH.filter.CH4.avg$FinalEF25 = allBOTH.filter.CH4.25$x
  allBOTH.filter.CH4.avg$FinalEF75 = allBOTH.filter.CH4.75$x
  allBOTH.filter.CH4.avg$MCE.25 = allBOTH.filter.CH4.25.mce$x
  allBOTH.filter.CH4.avg$MCE.75 = allBOTH.filter.CH4.75.mce$x
 
  allBOTH.filter.CH4$fuelORIG <- factor(allBOTH.filter.CH4$fuelORIG,
                                    levels = c("winter wheat", "soybean","rice","corn","grass","slash","pile","shrub","BlackwaterRiver"))
  allBOTH.filter.CH4.avg$fuelORIG = allBOTH.filter.CH4.avg$Group.1
  allBOTH.filter.CH4.avg$fuelORIG <- factor(allBOTH.filter.CH4.avg$fuelORIG,
                                        levels = c("winter wheat", "soybean","rice","corn","grass","slash","pile","shrub","BlackwaterRiver"))
  
  # Slopes
  ind = which(allBOTH.filter.CH4$fuelORIG == 'corn'& is.finite(allBOTH.filter.CH4$FinalEF))
  if (length(ind) > 2){
    ttC = cor.test(allBOTH.filter.CH4$FinalEF[ind], allBOTH.filter.CH4$MCE[ind])
    ttCORN = lm(allBOTH.filter.CH4$FinalEF[ind] ~ allBOTH.filter.CH4$MCE[ind])
    ttC2 = summary(ttCORN)
    if (ttC$p.value > 0.05) {ttCORN$coefficients[1:2] = NaN}
  }
  
  ind = which(allBOTH.filter.CH4$fuelORIG == 'soybean'& is.finite(allBOTH.filter.CH4$FinalEF))
  if (length(ind) > 2){
    ttC = cor.test(allBOTH.filter.CH4$FinalEF[ind], allBOTH.filter.CH4$MCE[ind])
    ttSOY = lm(allBOTH.filter.CH4$FinalEF[ind] ~ allBOTH.filter.CH4$MCE[ind])
    ttSY2 = summary(ttSOY)
    if (ttC$p.value > 0.05) {ttSOY$coefficients[1:2] = NaN}
  } else{
    ttSOY = data.frame('coefficients'=c(NaN,NaN))
    tmp = lm(c(1,1,2,5,9)~c(1,1,1,1,5))
    ttSY2 =summary(tmp)
    ttSY2$coefficients[,] = NaN
    
  }
  
  ind = which(allBOTH.filter.CH4$fuelORIG == 'rice'& is.finite(allBOTH.filter.CH4$FinalEF))
  if (length(ind)> 2){
     ttC = cor.test(allBOTH.filter.CH4$FinalEF[ind], allBOTH.filter.CH4$MCE[ind])
    ttRICE= lm(allBOTH.filter.CH4$FinalEF[ind] ~ allBOTH.filter.CH4$MCE[ind])
    ttR2 = summary(ttRICE)
    if (ttC$p.value > 0.05) {ttRICE$coefficients[1:2] = NaN}
  }else{
    ttRICE= data.frame('coefficients'=c(NaN,NaN))
    tmp = lm(c(1,1,2,5,9)~c(1,1,1,1,5))
    ttR2 =summary(tmp)
    ttR2$coefficients[,] = NaN
    
  }
  
  ind = which(allBOTH.filter.CH4$fuelORIG == 'pile'& is.finite(allBOTH.filter.CH4$FinalEF))
  if (length(ind) > 2){
    ttC = cor.test(allBOTH.filter.CH4$FinalEF[ind], allBOTH.filter.CH4$MCE[ind])
    ttPILE= lm(allBOTH.filter.CH4$FinalEF[ind] ~ allBOTH.filter.CH4$MCE[ind])
    ttP2 = summary(ttPILE)
    if (ttC$p.value > 0.05) {ttPILE$coefficients[1:2] = NaN}
  } else{
    ttPILE= data.frame('coefficients'=c(NaN,NaN))
    tmp = lm(c(1,1,2,5,9)~c(1,1,1,1,5))
    ttP2 =summary(tmp)
    ttP2$coefficients[,] = NaN
    
  }
  
  ind = which(allBOTH.filter.CH4$fuelORIG == 'slash'& is.finite(allBOTH.filter.CH4$FinalEF))
  if (length(ind) > 2){
    ttC = cor.test(allBOTH.filter.CH4$FinalEF[ind], allBOTH.filter.CH4$MCE[ind])
    ttSLASH= lm(allBOTH.filter.CH4$FinalEF[ind] ~ allBOTH.filter.CH4$MCE[ind])
    ttS2 = summary(ttSLASH)
    if (ttC$p.value > 0.05) {ttSLASH$coefficients[1:2] = NaN}
  } else{
    ttSLASH= data.frame('coefficients'=c(NaN,NaN))
    tmp = lm(c(1,1,2,5,9)~c(1,1,1,1,5))
    ttS2 =summary(tmp)
    ttS2$coefficients[,] = NaN
  }
  ind = which(allBOTH.filter.CH4$fuelORIG == 'grass' & is.finite(allBOTH.filter.CH4$FinalEF))
  if (length(ind) > 2){
    ttC = cor.test(allBOTH.filter.CH4$FinalEF[ind], allBOTH.filter.CH4$MCE[ind])
    ttGRASS= lm(allBOTH.filter.CH4$FinalEF[ind] ~ allBOTH.filter.CH4$MCE[ind])
    ttG2 = summary(ttGRASS)
    if (ttC$p.value > 0.05) {ttGRASS$coefficients[1:2] = NaN}
  }  else{
    ttGRASS= data.frame('coefficients'=c(NaN,NaN))
    tmp = lm(c(1,1,2,5,9)~c(1,1,1,1,5))
    ttG2 =summary(tmp)
    ttG2$coefficients[,] = NaN
  }
  
  slopes = c(ttCORN$coefficients[2], ttSOY$coefficients[2], ttRICE$coefficients[2], ttPILE$coefficients[2], ttSLASH$coefficients[2], ttGRASS$coefficients[2])
  slopesE = c(ttC2$coefficients[2,2], ttSY2$coefficients[2,2], ttR2$coefficients[2,2], ttP2$coefficients[2,2], ttS2$coefficients[2,2], ttG2$coefficients[2,2])
  ints = c(ttCORN$coefficients[1], ttSOY$coefficients[1], ttRICE$coefficients[1], ttPILE$coefficients[1], ttSLASH$coefficients[1], ttGRASS$coefficients[1])
  fuels = c('corn','soybean','rice','pile','slash','grass')
  alldat = as.data.frame(cbind(as.numeric(slopes),as.numeric(slopesE),as.numeric(ints),fuels))
  colnames(alldat) = c("slopes","slopeError","ints","fuels")
  print(alldat)
  ierr = which(allBOTH.filter.CH4.avg$COUNT_EFFINAL > 1)
  require(ggplot2)
  szz=4 ;fz=30
  xx = c(0.84,0.98)
  yyK=max(allBOTH.filter.CH4$FinalEF, na.rm=TRUE)*0.2
  yLL=(bquote(paste(.(akagiN)," EF, g kg"^-1,sep='')))
  if (akagiN == 'NO2'){ yLL=expression(paste('NO'[2]," EF, g kg"^-1,sep=''))}
  if (akagiN == 'SO2'){ yLL=expression(paste('SO'[2]," EF, g kg"^-1,sep=''))}
  if (akagiN == 'H2O2'){ yLL=expression(paste('H'[2],'O'[2], " EF, g kg"^-1,sep=''))}
  
  #if (akagiN == 'NOx (as NO)'){ yLL=expression(paste('NO'[x],"(as NO) EF, g kg"^-1,sep=''))}
  
  CH4vsMCE = ggplot(allBOTH.filter.CH4)+
    geom_point(aes(x=MCE,y=FinalEF, col=fuel2, shape=fuelORIG), stroke=2,size=szz)+
    theme_classic()+
    xlab("MCE")+labs(col="",shape="")+
    ylab(yLL) +
    theme(legend.background=element_blank())+
   # theme(legend.position = c(0.2, .65*max(allBOTH.filter.CH4.avg$FinalEF, na.rm=TRUE)))+
    theme(text = element_text(size = fz)) +
   # theme(legend.direction="horizontal")+
    #scale_shape_manual(values=c(19,15,17,18,7,8,4))+
    scale_color_manual(values = c(cbp1 ))+
    scale_shape_manual(values =fuelshapes)+
    xlim(xx)
  if (doavgshape == 1){
    CH4vsMCE = CH4vsMCE + geom_errorbar(data=allBOTH.filter.CH4.avg[ierr,],aes(xmin=MCE.25,xmax=MCE.75, y=FinalEF, col=fuelORIG),  position=position_dodge(0.05))+
    geom_errorbar(data=allBOTH.filter.CH4.avg[ierr,],aes(x=MCE, ymin=FinalEF25, ymax=FinalEF75,col=fuelORIG), width=0.0025, position=position_dodge(0.05))+
    geom_point(data=allBOTH.filter.CH4.avg,aes(x=MCE,y=FinalEF, shape=fuelORIG), col='white', size=12, stroke = 1)+
    geom_point(data=allBOTH.filter.CH4.avg,aes(x=MCE,y=FinalEF, col=fuelORIG,shape=fuelORIG), size=10, stroke = 2)
  }
  
  if (length(indC) > 0){
      CH4vsMCE = CH4vsMCE + 
#          geom_text(x = 1, y=.5,label='Liu et al., 2017',hjust = 0, size = 6, col=cbp1[7]) +
          coord_cartesian(xlim = xx, # This focuses the x-axis on the range of interest
                          clip = 'off') +   # This keeps the labels from disappearing
          theme(plot.margin = unit(c(1,3,1,1), "lines")) +# This widens the right margin
        geom_point(data=xiaoxi, aes(x=MCE, y=xiaoxi[,indC]),pch=0,col=cbp1[7],stroke=2,size=szz)
          
        #geom_errorbar(data=xiaoxi.avg,aes(xmin=mean[1]-sd[1],xmax=mean[1]+sd[1], y=mean[indC]),col=cbp1[7], width=1.3*2, position=position_dodge(0.05))+
        #geom_errorbar(data=xiaoxi.avg,aes(x=mean[1], ymin=mean[indC] - sd[indC], ymax=mean[indC]+ sd[indC]),col=cbp1[7], width=0.0025, position=position_dodge(0.05))
  }
  cropShape = 23
  savannahShape = 12
  temperateShape = 6
  if (length(indB) > 0){
    CH4vsMCE = CH4vsMCE + 
      #geom_text(x = 1, y=.45,label='Andreae, 2019',hjust = 0, size = 6, col='purple') +
      coord_cartesian(xlim = xx, # This focuses the x-axis on the range of interest
                      clip = 'off') +   # This keeps the labels from disappearing
      theme(plot.margin = unit(c(1,3,1,1), "lines")) +# This widens the right margin
      geom_point(data=andreae, aes(x=as.numeric(average[1]), y=as.numeric(average[indB])), col='#bc5fc9',stroke=3, size=5, shape=cropShape)+
      geom_point(data=andreae, aes(x=as.numeric(Savannah[1]), y=as.numeric(Savannah[indB])), col='#bc5fc9',stroke=3, size=5, shape=savannahShape)+
      geom_point(data=andreae, aes(x=as.numeric(Temperate[1]), y=as.numeric(Temperate[indB])), col='#bc5fc9',stroke=3, size=5, shape=temperateShape)+
      annotate(geom="text", x=.88, y=yyK, label="Andreae (Crop)",fontface="italic",size=5, col='#bc5fc9')+
      annotate(geom="text", x=.88, y=yyK*0.95, label="Andreae (Savannah)",fontface="italic",size=5, col='#bc5fc9')+
      annotate(geom="text", x=.88, y=yyK*0.9, label="Andreae (Temperate)",fontface="italic",size=5, col='#bc5fc9')
    #+
    #      geom_errorbar(data=andreae,aes(xmin=as.numeric(average[1])-as.numeric(andreae$std.dev.[1]),xmax=as.numeric(average[1])+as.numeric(andreae$std.dev.[1]), y=as.numeric(average[indB])),col='purple', width=1.3*2, position=position_dodge(0.05))+
#      geom_errorbar(data=andreae,aes(x=as.numeric(average[1]), ymin=as.numeric(average[indB])-as.numeric(andreae$std.dev.[indB]),ymax=as.numeric(average[indB])+as.numeric(andreae$std.dev.[indB])),col='purple', width=0.0025, position=position_dodge(0.05))
  }
  if (length(indA) > 0){
    CH4vsMCE = CH4vsMCE + 
#        geom_text(x = 0.98, y=.4,label='Akagi et al., 2011 (Pasture maintenance)',hjust = 0, size = 6, col='light cbp1c[10]') +
        coord_cartesian(xlim = xx, # This focuses the x-axis on the range of interest
                        clip = 'off') +   # This keeps the labels from disappearing
        theme(plot.margin = unit(c(1,3,1,1), "lines")) +# This widens the right margin
       # annotate(geom="shape",shape="cross",  color="cbp1c[10]", fill="cbp1c[10]")+
    
      #  geom_text(x = 0.98, y=.35, label='Akagi et al., 2011 (Crop residue)',hjust = 0, size = 6, col='cbp1c[10]') +
        coord_cartesian(xlim = xx, # This focuses the x-axis on the range of interest
                        clip = 'off') +   # This keeps the labels from disappearing
        theme(plot.margin = unit(c(1,3,1,1), "lines")) +# This widens the right margin
     # annotate(geom="shape",shape="plus", xmin=c(2,4), xmax=c(3,5), ymin=c(20,10) , ymax=c(30,20), alpha=0.2, color="cbp1c[10]", fill="cbp1c[10]")+
      
      geom_point(data=akagi, aes(x=as.numeric(Savannah[1]), y=as.numeric(Savannah[indA])), col=cbp1[6], size=5, shape=savannahShape, stroke=3)+
    #  geom_errorbar(data=akagi,aes(xmin=Savannah[1]-as.numeric(SavannahSD[1]),
    #                               xmax=Savannah[1]+as.numeric(SavannahSD[1]), 
    #                               y=Savannah[indA]),col='cbp1c[10]', width=1.3*2, position=position_dodge(0.05))+
   
      geom_point(data=akagi, aes(x=CropEF[1], y=CropEF[indA]), col= cbp1[6], size=5, shape=cropShape, stroke=3)+
      geom_point(data=akagi, aes(x=TemperateEF[1], y=TemperateEF[indA]), col= cbp1[6], size=5, shape=temperateShape, stroke=3)+
      annotate(geom="text", x=.88, y=yyK*0.9, label="Akagi (Crop)",fontface="italic",size=5, col=cbp1[6])+
      annotate(geom="text", x=.88, y=yyK*0.85, label="Akagi (Savannah)",fontface="italic",size=5, col=cbp1[6])+
    annotate(geom="text", x=.88, y=yyK*0.8, label="Akagi (Temperate)",fontface="italic",size=5, col=cbp1[6])
    
    # geom_errorbar(data=akagi,aes(xmin=CropEF[1]-as.numeric(CropSD[1]),
    #                               xmax=CropEF[1]+as.numeric(CropSD[1]), 
    #                               y=CropEF[indA]),col='cbp1c[10]', width=1.3*2, position=position_dodge(0.05))+
    #  geom_errorbar(data=akagi,aes(x=CropEF[1],
    #                               ymin=CropEF[indA] - as.numeric(CropSD[indA]),
     #                              ymax=CropEF[indA] + as.numeric(CropSD[indA])),
     #               col='cbp1c[10]', width=0.0025, position=position_dodge(0.05))
  }
    
  #ind = which( allBOTH.filter.CH4$MCE > 0.92)
  ind = which(  allBOTH.filter.CH4$fuelORIG == 'corn')# & allBOTH.filter.CH4$MCE > 0.92)
  print(cor.test(allBOTH.filter.CH4$FinalEF[ind], allBOTH.filter.CH4$MCE[ind]))
  
  #Our transformation function
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  #make all y-axis same decimal points
  CH4vsMCE=CH4vsMCE + scale_y_continuous(labels=scaleFUN) + theme(legend.position="none")
  doslopes=0
  if (doslopes == 1){
    CH4vsMCE=CH4vsMCE + geom_abline(data=alldat, aes(intercept=as.numeric(ints),slope=as.numeric(slopes), col=fuels))
    
  }
  return(CH4vsMCE)
}
# ----  ER vs. MCE -----
plotSpeciesCOMCE = function(allBOTH.filterIN,variableN,akagiN,andreaeN,xiaoxiN){
  doavgshape = 0
  
  xiaoxi$fuel = c('Corn','Corn','Corn', 'Rice','Rice',
                  'Rice','Rice','Rice','Rice',
                  'Rice','Rice','Soybean',
                  'Popcorn','DblCrpWinWht_Soy',
                  'WoodyWetlands')
  ind = which(xiaoxi$fuel == 'Corn' | xiaoxi$fuel == 'Popcorn'); xiaoxi$fuel[ind] = 'corn'
  ind = which(xiaoxi$fuel == 'Rice'); xiaoxi$fuel[ind] = 'rice'
  ind = which(xiaoxi$fuel == 'Soybean'); xiaoxi$fuel[ind] = 'soybean'
  ind = which(xiaoxi$fuel == 'DblCrpWinWht_Soy'); xiaoxi$fuel[ind] = 'soybean'
  ind = which(xiaoxi$fuel == 'WoodyWetlands'); xiaoxi$fuel[ind] = 'slash'
  indA = which(akagi$Species == akagiN)
  indB = which(andreae$Katie == andreaeN)
  cc = colnames(xiaoxi)
  indC  = which(cc == xiaoxiN)
  print(c(indA,indB,indC))
  ind = which(allBOTH.filterIN$fire == 'BlackwaterRiver')
  allBOTH.filterIN$fuelORIG[ind] = 'BlackwaterRiver'
  ind = which(allBOTH.filterIN$variable == variableN & allBOTH.filterIN$fire != 'Vivian' &
                allBOTH.filterIN$fire != 'Copper Breaks' &
                allBOTH.filterIN$fire != 'Invictus' &
                allBOTH.filterIN$fuel != 'coniferous/decidous')
  fuelshapes = c(19,0,2,18,7,8,4,6,1,12,13)
  fuellimits = c("corn","soybean","rice","winter wheat","grass","pile","slash","shrub","BlackwaterRiver")
  
  cbp1a <- c( "#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", 
              "#CC79A7","#000000","#999999","#CC0000")
  # colorblind safe option
  #cbp1b = c('#377eb8','#e41a1c','#4daf4a','#f781bf','#a65628','#984ea3','#ff7f00','#ffff33')
  cbp1b = c('#377eb8','#e41a1c','#4daf4a','#a65628','#984ea3','#ff7f00',"#000000")
  cbp1c=c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9',
          '#74add1','#4575b4','#313695')
  cbp1d <- c( "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000","#E69F00", "#56B4E9", "#009E73")
  
  cbp1 = cbp1b
  allBOTH.filter.CH4 = allBOTH.filterIN[ind,]
  allBOTH.filter.CH4.avg = aggregate(allBOTH.filter.CH4, by=list(allBOTH.filter.CH4$fuel), FUN='median', na.rm=TRUE)
  #allBOTH.filter.CH4.sd = aggregate(allBOTH.filter.CH4, by=list(allBOTH.filter.CH4$fuel), FUN='sd', na.rm=TRUE)
  allBOTH.filter.CH4.25 = aggregate(allBOTH.filter.CH4$FinalERtoCO, by=list(allBOTH.filter.CH4$fuel), 'quantile',probs=c(q1),na.rm=TRUE)
  allBOTH.filter.CH4.75 = aggregate(allBOTH.filter.CH4$FinalERtoCO, by=list(allBOTH.filter.CH4$fuel), 'quantile',probs=c(q2),na.rm=TRUE)
  allBOTH.filter.CH4.25.mce = aggregate(allBOTH.filter.CH4$MCE, by=list(allBOTH.filter.CH4$fuel), 'quantile',probs=c(q1),na.rm=TRUE)
  allBOTH.filter.CH4.75.mce = aggregate(allBOTH.filter.CH4$MCE, by=list(allBOTH.filter.CH4$fuel), 'quantile',probs=c(q2),na.rm=TRUE)
  allBOTH.filter.CH4.avg = getplumesANDmcebyfuel1VAR(allBOTH.filter.CH4.avg, allBOTH.filter.CH4 )
  
  allBOTH.filter.CH4.avg$FinalERtoCO25 = allBOTH.filter.CH4.25$x
  allBOTH.filter.CH4.avg$FinalERtoCO75 = allBOTH.filter.CH4.75$x
  allBOTH.filter.CH4.avg$MCE.25 = allBOTH.filter.CH4.25.mce$x
  allBOTH.filter.CH4.avg$MCE.75 = allBOTH.filter.CH4.75.mce$x
  
  allBOTH.filter.CH4$fuel <- factor(allBOTH.filter.CH4$fuel,
                                    levels = c("winter wheat", "soybean","rice","corn","grass","slash","pile","shrub","BlackwaterRiver"))
  allBOTH.filter.CH4.avg$fuel = allBOTH.filter.CH4.avg$Group.1
  allBOTH.filter.CH4.avg$fuel <- factor(allBOTH.filter.CH4.avg$fuel,
                                        levels = c("winter wheat", "soybean","rice","corn","grass","slash","pile","shrub","BlackwaterRiver"))
  
  # Slopes
  ind = which(allBOTH.filter.CH4$fuel == 'corn'& is.finite(allBOTH.filter.CH4$FinalERtoCO))
  if (length(ind) > 2){
    ttC = cor.test(allBOTH.filter.CH4$FinalERtoCO[ind], allBOTH.filter.CH4$MCE[ind])
    ttCORN = lm(allBOTH.filter.CH4$FinalERtoCO[ind] ~ allBOTH.filter.CH4$MCE[ind])
    ttC2 = summary(ttCORN)
    if (ttC$p.value > 0.05) {ttCORN$coefficients[1:2] = NaN}
  }
  
  ind = which(allBOTH.filter.CH4$fuel == 'soybean'& is.finite(allBOTH.filter.CH4$FinalERtoCO))
  if (length(ind) > 2){
    ttC = cor.test(allBOTH.filter.CH4$FinalERtoCO[ind], allBOTH.filter.CH4$MCE[ind])
    ttSOY = lm(allBOTH.filter.CH4$FinalERtoCO[ind] ~ allBOTH.filter.CH4$MCE[ind])
    ttSY2 = summary(ttSOY)
    if (ttC$p.value > 0.05) {ttSOY$coefficients[1:2] = NaN}
  } else{
    ttSOY = data.frame('coefficients'=c(NaN,NaN))
    tmp = lm(c(1,1,2,5,9)~c(1,1,1,1,5))
    ttSY2 =summary(tmp)
    ttSY2$coefficients[,] = NaN
    
  }
  
  ind = which(allBOTH.filter.CH4$fuel == 'rice'& is.finite(allBOTH.filter.CH4$FinalERtoCO))
  if (length(ind)> 2){
    ttC = cor.test(allBOTH.filter.CH4$FinalERtoCO[ind], allBOTH.filter.CH4$MCE[ind])
    ttRICE= lm(allBOTH.filter.CH4$FinalERtoCO[ind] ~ allBOTH.filter.CH4$MCE[ind])
    ttR2 = summary(ttRICE)
    if (ttC$p.value > 0.05) {ttRICE$coefficients[1:2] = NaN}
  }else{
    ttRICE= data.frame('coefficients'=c(NaN,NaN))
    tmp = lm(c(1,1,2,5,9)~c(1,1,1,1,5))
    ttR2 =summary(tmp)
    ttR2$coefficients[,] = NaN
    
  }
  
  ind = which(allBOTH.filter.CH4$fuel == 'pile'& is.finite(allBOTH.filter.CH4$FinalERtoCO))
  if (length(ind) > 2){
    ttC = cor.test(allBOTH.filter.CH4$FinalERtoCO[ind], allBOTH.filter.CH4$MCE[ind])
    ttPILE= lm(allBOTH.filter.CH4$FinalERtoCO[ind] ~ allBOTH.filter.CH4$MCE[ind])
    ttP2 = summary(ttPILE)
    if (ttC$p.value > 0.05) {ttPILE$coefficients[1:2] = NaN}
  } else{
    ttPILE= data.frame('coefficients'=c(NaN,NaN))
    tmp = lm(c(1,1,2,5,9)~c(1,1,1,1,5))
    ttP2 =summary(tmp)
    ttP2$coefficients[,] = NaN
    
  }
  
  ind = which(allBOTH.filter.CH4$fuel == 'slash'& is.finite(allBOTH.filter.CH4$FinalERtoCO))
  if (length(ind) > 2){
    ttC = cor.test(allBOTH.filter.CH4$FinalERtoCO[ind], allBOTH.filter.CH4$MCE[ind])
    ttSLASH= lm(allBOTH.filter.CH4$FinalERtoCO[ind] ~ allBOTH.filter.CH4$MCE[ind])
    ttS2 = summary(ttSLASH)
    if (ttC$p.value > 0.05) {ttSLASH$coefficients[1:2] = NaN}
  }
  ind = which(allBOTH.filter.CH4$fuel == 'grass' & is.finite(allBOTH.filter.CH4$FinalERtoCO))
  if (length(ind) > 2){
    ttC = cor.test(allBOTH.filter.CH4$FinalERtoCO[ind], allBOTH.filter.CH4$MCE[ind])
    ttGRASS= lm(allBOTH.filter.CH4$FinalERtoCO[ind] ~ allBOTH.filter.CH4$MCE[ind])
    ttG2 = summary(ttGRASS)
    if (ttC$p.value > 0.05) {ttGRASS$coefficients[1:2] = NaN}
  }  else{
    ttGRASS= data.frame('coefficients'=c(NaN,NaN))
    tmp = lm(c(1,1,2,5,9)~c(1,1,1,1,5))
    ttG2 =summary(tmp)
    ttG2$coefficients[,] = NaN
  }
  
  slopes = c(ttCORN$coefficients[2], ttSOY$coefficients[2], ttRICE$coefficients[2], ttPILE$coefficients[2], ttSLASH$coefficients[2], ttGRASS$coefficients[2])
  slopesE = c(ttC2$coefficients[2,2], ttSY2$coefficients[2,2], ttR2$coefficients[2,2], ttP2$coefficients[2,2], ttS2$coefficients[2,2], ttG2$coefficients[2,2])
  ints = c(ttCORN$coefficients[1], ttSOY$coefficients[1], ttRICE$coefficients[1], ttPILE$coefficients[1], ttSLASH$coefficients[1], ttGRASS$coefficients[1])
  fuels = c('corn','soybean','rice','pile','slash','grass')
  alldat = as.data.frame(cbind(as.numeric(slopes),as.numeric(slopesE),as.numeric(ints),fuels))
  colnames(alldat) = c("slopes","slopeError","ints","fuels")
  print(alldat)
  ierr = which(allBOTH.filter.CH4.avg$COUNT_EFFINAL > 1)
  require(ggplot2)
  szz=4 ;fz=30
  xx = c(0.84,0.98)
  CH4vsMCE = ggplot(allBOTH.filter.CH4)+
    geom_point(aes(x=MCE,y=FinalERtoCO, col=fuel2, shape=fuel), stroke=2,size=szz)+
    theme_classic()+
    xlab("MCE")+labs(col="",shape="")+
    ylab(bquote(paste(.(akagiN)," EF, g kg"^-1,sep='')))+
    theme(legend.background=element_blank())+
    # theme(legend.position = c(0.2, .65*max(allBOTH.filter.CH4.avg$FinalERtoCO, na.rm=TRUE)))+
    theme(text = element_text(size = fz)) +
    # theme(legend.direction="horizontal")+
    #scale_shape_manual(values=c(19,15,17,18,7,8,4))+
    scale_color_manual(values = c(cbp1 ))+
    scale_shape_manual(values =fuelshapes)+
    xlim(xx)
  if (doavgshape == 1){
    CH4vsMCE = CH4vsMCE + geom_errorbar(data=allBOTH.filter.CH4.avg[ierr,],aes(xmin=MCE.25,xmax=MCE.75, y=FinalERtoCO, col=fuel),  position=position_dodge(0.05))+
      geom_errorbar(data=allBOTH.filter.CH4.avg[ierr,],aes(x=MCE, ymin=FinalERtoCO25, ymax=FinalERtoCO75,col=fuel), width=0.0025, position=position_dodge(0.05))+
      geom_point(data=allBOTH.filter.CH4.avg,aes(x=MCE,y=FinalERtoCO, shape=fuel), col='white', size=12, stroke = 1)+
      geom_point(data=allBOTH.filter.CH4.avg,aes(x=MCE,y=FinalERtoCO, col=fuel,shape=fuel), size=10, stroke = 2)
  }
  
  if (length(indC) > 0){
    CH4vsMCE = CH4vsMCE + 
      #          geom_text(x = 1, y=.5,label='Liu et al., 2017',hjust = 0, size = 6, col=cbp1[7]) +
      coord_cartesian(xlim = xx, # This focuses the x-axis on the range of interest
                      clip = 'off') +   # This keeps the labels from disappearing
      theme(plot.margin = unit(c(1,3,1,1), "lines")) +# This widens the right margin
      geom_point(data=xiaoxi, aes(x=MCE, y=xiaoxi[,indC]),pch=0,col=cbp1[7],stroke=2,size=szz)
    
    #geom_errorbar(data=xiaoxi.avg,aes(xmin=mean[1]-sd[1],xmax=mean[1]+sd[1], y=mean[indC]),col=cbp1[7], width=1.3*2, position=position_dodge(0.05))+
    #geom_errorbar(data=xiaoxi.avg,aes(x=mean[1], ymin=mean[indC] - sd[indC], ymax=mean[indC]+ sd[indC]),col=cbp1[7], width=0.0025, position=position_dodge(0.05))
  }
  if (length(indB) > 0){
    CH4vsMCE = CH4vsMCE + 
      #geom_text(x = 1, y=.45,label='Andreae, 2019',hjust = 0, size = 6, col='purple') +
      coord_cartesian(xlim = xx, # This focuses the x-axis on the range of interest
                      clip = 'off') +   # This keeps the labels from disappearing
      theme(plot.margin = unit(c(1,3,1,1), "lines")) +# This widens the right margin
      geom_point(data=andreae, aes(x=as.numeric(average[1]), y=as.numeric(average[indB])), col='#bc5fc9',stroke=3, size=5, shape=3)#+
    #      geom_errorbar(data=andreae,aes(xmin=as.numeric(average[1])-as.numeric(andreae$std.dev.[1]),xmax=as.numeric(average[1])+as.numeric(andreae$std.dev.[1]), y=as.numeric(average[indB])),col='purple', width=1.3*2, position=position_dodge(0.05))+
    #      geom_errorbar(data=andreae,aes(x=as.numeric(average[1]), ymin=as.numeric(average[indB])-as.numeric(andreae$std.dev.[indB]),ymax=as.numeric(average[indB])+as.numeric(andreae$std.dev.[indB])),col='purple', width=0.0025, position=position_dodge(0.05))
  }
  if (length(indA) > 0){
    CH4vsMCE = CH4vsMCE + 
      #        geom_text(x = 0.98, y=.4,label='Akagi et al., 2011 (Pasture maintenance)',hjust = 0, size = 6, col='light cbp1c[10]') +
      coord_cartesian(xlim = xx, # This focuses the x-axis on the range of interest
                      clip = 'off') +   # This keeps the labels from disappearing
      theme(plot.margin = unit(c(1,3,1,1), "lines")) +# This widens the right margin
      # annotate(geom="shape",shape="cross",  color="cbp1c[10]", fill="cbp1c[10]")+
      
      #  geom_text(x = 0.98, y=.35, label='Akagi et al., 2011 (Crop residue)',hjust = 0, size = 6, col='cbp1c[10]') +
      coord_cartesian(xlim = xx, # This focuses the x-axis on the range of interest
                      clip = 'off') +   # This keeps the labels from disappearing
      theme(plot.margin = unit(c(1,3,1,1), "lines")) +# This widens the right margin
      # annotate(geom="shape",shape="plus", xmin=c(2,4), xmax=c(3,5), ymin=c(20,10) , ymax=c(30,20), alpha=0.2, color="cbp1c[10]", fill="cbp1c[10]")+
      
      #   geom_point(data=akagi, aes(x=PastureEF[1], y=PastureEF[indA]), col=cbp1[6], size=5, shape=4, stroke=3)+
      # geom_errorbar(data=akagi,aes(xmin=PastureEF[1]-as.numeric(PastureSD[1]),
      #                               xmax=PastureEF[1]+as.numeric(PastureSD[1]), 
      #                               y=PastureEF[indA]),col='cbp1c[10]', width=1.3*2, position=position_dodge(0.05))+
      # geom_errorbar(data=akagi,aes(x=PastureEF[1],
      #                               ymin=PastureEF[indA] - as.numeric(PastureSD[indA]),
      #                               ymax=PastureEF[indA] + as.numeric(PastureSD[indA])),
      #                col='cbp1c[10]', width=0.0025, position=position_dodge(0.05))+
      geom_point(data=akagi, aes(x=CropEF[1], y=CropEF[indA]), col= cbp1[6], size=5, shape=3, stroke=3)
    # geom_errorbar(data=akagi,aes(xmin=CropEF[1]-as.numeric(CropSD[1]),
    #                               xmax=CropEF[1]+as.numeric(CropSD[1]), 
    #                               y=CropEF[indA]),col='cbp1c[10]', width=1.3*2, position=position_dodge(0.05))+
    #  geom_errorbar(data=akagi,aes(x=CropEF[1],
    #                               ymin=CropEF[indA] - as.numeric(CropSD[indA]),
    #                              ymax=CropEF[indA] + as.numeric(CropSD[indA])),
    #               col='cbp1c[10]', width=0.0025, position=position_dodge(0.05))
  }
  
  #ind = which( allBOTH.filter.CH4$MCE > 0.92)
  ind = which(  allBOTH.filter.CH4$fuel == 'corn')# & allBOTH.filter.CH4$MCE > 0.92)
  print(cor.test(allBOTH.filter.CH4$FinalERtoCO[ind], allBOTH.filter.CH4$MCE[ind]))
  
  #Our transformation function
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  #make all y-axis same decimal points
  CH4vsMCE=CH4vsMCE + scale_y_continuous(labels=scaleFUN) + theme(legend.position="none")
  doslopes=0
  if (doslopes == 1){
    CH4vsMCE=CH4vsMCE + geom_abline(data=alldat, aes(intercept=as.numeric(ints),slope=as.numeric(slopes), col=fuels))
    
  }
  return(CH4vsMCE)
}



