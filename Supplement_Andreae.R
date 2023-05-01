 # --------- Supplement ----------
f2 = '/Users/ktravis1/OneDrive - NASA/ForGITHUB/InputFiles/OtherStudies/Andreae-BB-EMFactors-14Apr2021_justtable1.csv'
andreae = read.csv(f2)
require(readxl)
akagi=readxl::read_xlsx('/Users/ktravis1/OneDrive - NASA/ForGITHUB/InputFiles/OtherStudies/Akagi_acp-11-4039-2011-supplement/Tables 1-5_4.27.11.xlsx')


# ------- Need to apply Table 4 --------
outputdata = as.data.frame(cbind(new.data.frame$Category, new.data.frame$LifetimeCat,new.data.frame$mWs,
                                 new.data.frame$names,new.data.frame$formula,new.data.frame$PI,
                                 new.data.frame.ag$FinalEF_mean, new.data.frame.ag$FinalEF_sd, new.data.frame.ag$COUNT_EFFINAL,  
                                 new.data.frame.presc$FinalEF_mean,   new.data.frame.presc$FinalEF_sd,         new.data.frame.presc$COUNT_EFFINAL,
                                 new.data.frame.grass$FinalEF_mean,   new.data.frame.grass$FinalEF_sd,     new.data.frame.grass$COUNT_EFFINAL,
                                 new.data.frame.blackwater$FinalEF_mean, new.data.frame.blackwater$FinalEF_sd, new.data.frame.blackwater$COUNT_EFFINAL))

colnames(outputdata) = c("Category", 'LifetimeCat','mWs',  'names','formula','PI',
                        'FinalEF_mean_ag',      'FinalEF_sd_ag',    'COUNT_EFFINAL_ag',  
                         'FinalEF_mean_presc',   'FinalEF_sd_presc',    'COUNT_EFFINAL_presc',
                         'FinalEF_mean_grass',   'FinalEF_sd_grass',    'COUNT_EFFINAL_grass',
                        'FinalEF_mean_brsf',   'FinalEF_sd_brsf',    'COUNT_EFFINAL_brsf')

outputdata$kind = ''
for (i in 1:length(outputdata$Category)){
  ind = which(allBOTH.filter$names == outputdata$names[i])
  outputdata$kind[i] = allBOTH.filter$kind[ind[1]]
}

ind = which(outputdata$names == 'OA')
outputdata$names[ind] = 'OC'
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

# Sum lumped species
ind = which(outputdata$names == 'o-Xylene')
ind2 = which(outputdata$names == 'm,p-Xylene')
tmp = outputdata[ind,]
tmp$FinalEF_mean_ag = as.numeric(outputdata$FinalEF_mean_ag[ind])       +  as.numeric(outputdata$FinalEF_mean_ag[ind2])
tmp$FinalEF_mean_presc = as.numeric(outputdata$FinalEF_mean_presc[ind]) +  as.numeric(outputdata$FinalEF_mean_presc[ind2])
tmp$FinalEF_mean_grass = as.numeric(outputdata$FinalEF_mean_grass[ind]) +  as.numeric(outputdata$FinalEF_mean_grass[ind2])
tmp$FinalEF_mean_brsf = as.numeric(outputdata$FinalEF_mean_brsf[ind])   +  as.numeric(outputdata$FinalEF_mean_brsf[ind2])
tmp$names  = 'Xylenes'
outputdata = rbind(outputdata,tmp)

ind = which(outputdata$names== '1,3-Butadiene' )
ind2 = which(outputdata$names == '1,2-Butadiene')
tmp = outputdata[ind,]
tmp$FinalEF_mean_ag = as.numeric(outputdata$FinalEF_mean_ag[ind])       +  as.numeric(outputdata$FinalEF_mean_ag[ind2])
tmp$FinalEF_mean_presc = as.numeric(outputdata$FinalEF_mean_presc[ind]) +  as.numeric(outputdata$FinalEF_mean_presc[ind2])
tmp$FinalEF_mean_grass = as.numeric(outputdata$FinalEF_mean_grass[ind]) +  as.numeric(outputdata$FinalEF_mean_grass[ind2])
tmp$FinalEF_mean_brsf = as.numeric(outputdata$FinalEF_mean_brsf[ind])   +  as.numeric(outputdata$FinalEF_mean_brsf[ind2])
tmp$names = 'Butadiene'
outputdata = rbind(outputdata,tmp)

ind = which(outputdata$names== 'trans-2-Pentene' )
ind2 = which(outputdata$names == 'cis-2-Pentene')
tmp = outputdata[ind,]
tmp$FinalEF_mean_ag = as.numeric(outputdata$FinalEF_mean_ag[ind])       +  as.numeric(outputdata$FinalEF_mean_ag[ind2])
tmp$FinalEF_mean_presc = as.numeric(outputdata$FinalEF_mean_presc[ind]) +  as.numeric(outputdata$FinalEF_mean_presc[ind2])
tmp$FinalEF_mean_grass = as.numeric(outputdata$FinalEF_mean_grass[ind]) +  as.numeric(outputdata$FinalEF_mean_grass[ind2])
tmp$FinalEF_mean_brsf = as.numeric(outputdata$FinalEF_mean_brsf[ind])   +  as.numeric(outputdata$FinalEF_mean_brsf[ind2])
tmp$names = '2-pentenes'
outputdata = rbind(outputdata,tmp)

ind = which(outputdata$names== '2-Methyl-1-butene' )
ind2 = which(outputdata$names == '3-Methyl-1-butene')
ind3 = which(outputdata$names == '2-Methyl-2-butene' )
tmp = outputdata[ind,]
tmp$FinalEF_mean_ag = as.numeric(outputdata$FinalEF_mean_ag[ind])       +  as.numeric(outputdata$FinalEF_mean_ag[ind2])+  as.numeric(outputdata$FinalEF_mean_ag[ind3])
tmp$FinalEF_mean_presc = as.numeric(outputdata$FinalEF_mean_presc[ind]) +  as.numeric(outputdata$FinalEF_mean_presc[ind2])+  as.numeric(outputdata$FinalEF_mean_presc[ind3])
tmp$FinalEF_mean_grass = as.numeric(outputdata$FinalEF_mean_grass[ind]) +  as.numeric(outputdata$FinalEF_mean_grass[ind2]) +  as.numeric(outputdata$FinalEF_mean_grass[ind3])
tmp$FinalEF_mean_brsf = as.numeric(outputdata$FinalEF_mean_brsf[ind])   +  as.numeric(outputdata$FinalEF_mean_brsf[ind2]) +  as.numeric(outputdata$FinalEF_mean_brsf[ind3])
tmp$names = 'Methyl-butenes'
outputdata = rbind(outputdata,tmp)

# --------- Get Andreae emission factors ------------
outputdata$AndreaeEF = NaN
outputdata$AndreaeEFsd = NaN
outputdata$AndreaeName = NaN
outputdata$AndreaeNN = NaN
outputdata$AndreaeEF_temp = NaN
outputdata$AndreaeEFsd_temp = NaN
outputdata$AndreaeNN_temp = NaN
outputdata$AndreaeEF_sav = NaN
outputdata$AndreaeEFsd_sav = NaN
outputdata$AndreaeNN_sav = NaN
fix=c()
for (i in 1:length(outputdata$AndreaeEF)){
  ind2 = which( outputdata$names[i] == andreae$Katie)
  if (length(ind2) == 0){
    fix=c(fix,outputdata$names[i])
    } #else{  print(c(outputdata$names[i], andreae$Katie[ind]))}
  if (length(ind2) >1){
#    outputdata$AndreaeEF[i] = andreae$average[ind2]
#    outputdata$AndreaeEFsd[i] = andreae$std.dev.[ind2]
#    outputdata$AndreaeName[i] = andreae$Species[ind2]
#    outputdata$AndreaeNN[i] = andreae$N[ind2]
    print(">2")
    }
  
  if (length(ind2) ==1){
    outputdata$AndreaeEF[i] = andreae$average[ind2]
    outputdata$AndreaeEFsd[i] = andreae$std.dev.[ind2]
    outputdata$AndreaeName[i] = andreae$Species[ind2]
    outputdata$AndreaeNN[i] = andreae$N[ind2]
    outputdata$AndreaeEF_sav[i] = andreae$Savannah[ind2]
    outputdata$AndreaeEFsd_sav[i] = andreae$SavannahSD[ind2]
    outputdata$AndreaeNN_sav[i] = andreae$SavannahN[ind2]
    
    outputdata$AndreaeEF_temp[i] = andreae$Temperate[ind2]
    outputdata$AndreaeEFsd_temp[i] = andreae$TemperateSD[ind2]
    outputdata$AndreaeNN_temp[i] = andreae$TemperateN[ind2]
  }
}

# ---------- Make a comparison plot to Andreae, 2019 - ------
# ---- recreate Gigi's plot
doplot = 0
if (doplot == 1){
  library(viridis)
  # ----- plot ----
  
  #outputdata= subset(outputdata, select = -c( Group.1))
  ind = which(outputdata$variable != 'NMVOC')
  outputdata= outputdata[order(outputdata$FinalEF,decreasing = FALSE),]
  outputdata$variableNUM = seq(1,length(outputdata$variable))
  outputdata$AndreaeEF = as.numeric(outputdata$AndreaeEF)
  outputdata$AndreaeEFsd = as.numeric(outputdata$AndreaeEFsd)
  # ---- sort from highest to lowest EF
  outputdata =outputdata[order(outputdata$FinalEF,decreasing = TRUE),]
  require(ggbreak)
  ind = which(is.finite(outputdata$FinalEF) & outputdata$variable != 'CO2_ppb'
              & outputdata$variable != 'CO_DACOM_DISKIN' & 
                outputdata$variable != 'CH4_DACOM_DISKIN'  &
                outputdata$variable != 'NMVOC'  &
                outputdata$USEME == 1 &
               is.finite(outputdata$AndreaeEF))
  tmp = outputdata[ind,]
  tmp$variableNUM = seq(1,length(tmp$variable))
  tmp$AndreaeNN = as.numeric(tmp$AndreaeNN)
  ind = which(tmp$FinalEF > .1)
  andreaePLOT1 = ggplot(tmp[ind,], aes(fill=AndreaeNN,y=FinalEF, x=variableNUM)) + 
    geom_bar(position="dodge", stat="identity") + theme_classic()+ # ylim(-3,7)+
    geom_errorbar(aes(x=variableNUM,ymin=FinalEF_25, ymax=FinalEF_75), width=.2,
                  position=position_dodge(.9), col='grey') +
    geom_point(aes(x=variableNUM, y=AndreaeEF), show.legend=FALSE)+
    geom_errorbar(aes(x=variableNUM,ymin=AndreaeEF-AndreaeEFsd, ymax=AndreaeEF+AndreaeEFsd), width=.2,
                  position=position_dodge(.9)) +
    scale_fill_viridis()+
    #scale_fill_manual(values = c("purple", "blue", "green","yellow","red"))+
    #scale_fill_manual(values=as.vector(polychrome(26)))+ theme(legend.position = "top")+labs(fill="")+
    guides(fill=guide_legend(title="#EFs"))+
    theme(text = element_text(size=20))  +
    geom_text(aes(x =variableNUM,y = 0,label = names), size=5,
              vjust = 0,  hjust = 1,angle = 90, nudge_y = -.01, nudge_x = 0.2)+ylab('EF, g/kg') 
  
  ind = which(tmp$FinalEF <100)
  andreaePLOT2 = ggplot(tmp[ind,], aes(fill=AndreaeNN,y=FinalEF, x=variableNUM)) + 
    geom_bar(position="dodge", stat="identity") + theme_classic()+ # ylim(-3,7)+
    geom_errorbar(aes(x=variableNUM,ymin=FinalEF_25, ymax=FinalEF_75), width=.2,
                  position=position_dodge(.9), col='grey') +
    geom_point(aes(x=variableNUM, y=AndreaeEF), show.legend=FALSE)+
    geom_errorbar(aes(x=variableNUM,ymin=AndreaeEF-AndreaeEFsd, ymax=AndreaeEF+AndreaeEFsd), width=.2,
                  position=position_dodge(.9)) +
    scale_fill_viridis()+
    #scale_fill_manual(values = c("purple", "blue", "green","yellow","red"))+
    #scale_fill_manual(values=as.vector(polychrome(26)))+ theme(legend.position = "top")+labs(fill="")+
    guides(fill=guide_legend(title="#EFs"))+
    theme(text = element_text(size=20))  +
    geom_text(aes(x =variableNUM,y = 0,label = names), size=5,
              vjust = 0,  hjust = 1,angle = 90, nudge_y = -.01, nudge_x = 0.2)+ylab('EF, g/kg') +
    scale_y_cut(breaks=0.5,1)
  ##ggsave(plot = andreaePLOT2, file='AndreaePLOT.pdf',width = 7*1.25*1.75*2, height=7*1.25*2)
}  
# --------------- Make agricultural output -----

ind = which(outputdata$PI == 'BLAKE')
outputdata$PI[ind] = 'WAS'

ind = which(outputdata$PI == 'DISKIN')
outputdata$PI[ind] = 'DACOM'
ind = which(outputdata$PI == 'DISKIN' & outputdata$formula == 'CH4')
outputdata$PI[ind] = 'NIR spect.'
ind = which(outputdata$PI == 'JIMENEZ')
outputdata$PI[ind] = 'AMS'
ind = which(outputdata$PI == 'ppt')
outputdata$PI[ind] = 'TOGA'
ind = which(outputdata$PI == 'ppt+BLAKE+GILMAN')
outputdata$PI[ind] = 'TOGA, iWAS, WAS'
ind = which(outputdata$PI == 'GILMAN+ppt')
outputdata$PI[ind] = 'TOGA, iWAS'
ind = which(outputdata$PI == 'GILMAN+BLAKE')
outputdata$PI[ind] = 'iWAS, WAS'
ind = which(outputdata$PI == 'BLAKE+GILMAN')
outputdata$PI[ind] = 'iWAS, WAS'
ind = which(outputdata$PI == 'APEL+BLAKE')
outputdata$PI[ind] = 'TOGA, WAS'
ind = which(outputdata$PI == 'BLAKE+ppt')
outputdata$PI[ind] = 'TOGA, WAS'
ind = which(outputdata$PI == 'ppt+BLAKE')
outputdata$PI[ind] = 'TOGA, WAS'
ind = which(outputdata$PI == 'ppt+BLAKE+WARNEKE')
outputdata$PI[ind] = 'NOAA PTRMS, TOGA, WAS'
ind = which(outputdata$PI == 'VERES')
outputdata$PI[ind] = 'NOAA CIMS'
ind = which(outputdata$PI == 'WISTHALER')
outputdata$PI[ind] = 'OSLO PTRMS'
ind = which(outputdata$PI == 'FRIED')
outputdata$PI[ind] = 'CAMS'
ind = which(outputdata$PI == 'GILMAN')
outputdata$PI[ind] = 'iWAS'
ind = which(outputdata$PI == 'FRIED+HANISCO')
outputdata$PI[ind] = 'CAMS, ISAF'
ind = which(outputdata$PI == 'GILMAN+FRIED')
outputdata$PI[ind] = 'iWAS, CAMS'
ind = which(outputdata$PI == 'VERES+WARKEKE')
outputdata$PI[ind] = 'NOAA CIMS, NOAA PTR-TOF-MS'
ind = which(outputdata$PI == 'WARNEKE' | outputdata$PI == 'ppbv')
outputdata$PI[ind] = 'NOAA PTR-ToF-MS'
ind = which(outputdata$PI == 'WOMACK')
outputdata$PI[ind] = 'ACES'
ind = which(outputdata$PI == 'WOMACK+RYERSON')
outputdata$PI[ind] = 'ACES, NOAA NOyO3'
ind = which(outputdata$PI == 'WOMACK+VERES')
outputdata$PI[ind] = 'ACES, NOAA CIMS'
ind = which(outputdata$PI == 'Warneke w/Blake+GILMAN+APEL')
outputdata$PI[ind] = 'NOAA PTR-TOF-MS, TOGA+WAS+iWAS spec.'
ind = which(outputdata$PI == 'Warneke w/APEL')
outputdata$PI[ind] = 'NOAA PTR-TOF-MS, TOGA spec.'
ind = which(outputdata$PI == 'WENNBERG')
outputdata$PI[ind] = 'CIT-CIMS'
ind = which(outputdata$PI == 'WARNEKE+ppt+GILMAN')
outputdata$PI[ind] = 'NOAA PTR-TOF-MS, TOGA, iWAS'
ind = which(outputdata$PI == 'WARNEKE+BLAKE+ppt')
outputdata$PI[ind] = 'NOAA PTR-TOF-MS, TOGA, WAS'
ind = which(outputdata$PI == 'VERES+WENNBERG+WARNEKE+ppt')
outputdata$PI[ind] = 'NOAA CIMS, CIT-CIMS, NOAA PTR-TOF-MS, TOGA'
ind = which(outputdata$PI == 'VERES+WARNEKE')
outputdata$PI[ind] = 'NOAA CIMS, NOAA PTR-TOF-MS'
ind = which(outputdata$PI == 'ROLLINS+RYERSON')
outputdata$PI[ind] = 'NOAA LIF, NOAA NOyO3'
ind = which(outputdata$PI == 'ROLLINS')
outputdata$PI[ind] = 'NOAA LIF'
ind = which(outputdata$PI == 'RYERSON')
outputdata$PI[ind] = 'NOAA NOyO3'
ind = which(outputdata$PI == 'SCHWARZ')
outputdata$PI[ind] = 'NOAA SP2'
ind = which(outputdata$PI == 'APEL')
outputdata$PI[ind] = 'TOGA'

  
write.csv(outputdata, file='Andreae_Matchup.csv')

