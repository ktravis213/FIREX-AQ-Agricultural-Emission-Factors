# Function to speciate PTRMS data based on iWAS, WAS, WAS etc
speciateSpecies = function(speciateData){
  
  # kludge for different aggregation (Group.1 or Group.2)
  cc = colnames(speciateData)
  ind = which(cc == 'Group.2')
  if (length(ind) == 0){
    speciateData$Group.2 = speciateData$Group.1
  }
  ind = which(speciateData$fuel != 'coniferous/decidous' & speciateData$fuel != 'forest')
  if (length(ind) > 0){speciateData = speciateData[ind,]}
  doxylene =1; doacetone = 1; domvk =1; domek =1
  if (doxylene == 1){
    # --------- Sum of m-xylene p-xylene o-xylene and ethyl benzene --------------
    # ugg need to get rid of "coniferous/decidous
  
    ind1 = which(speciateData$formula== 'C8H10' & speciateData$names == 'Ethylbenzene' & 
                   speciateData$PI == 'ppt+BLAKE+GILMAN')
    ind2 = which(speciateData$formula== 'C8H10' & speciateData$names == 'o-Xylene' & 
                   speciateData$PI == 'ppt+BLAKE+GILMAN')
    ind3 = which(speciateData$formula== 'C8H10' & speciateData$names == 'm,p-Xylene' & 
                   speciateData$PI == 'ppt+BLAKE+GILMAN')
    ind = which(speciateData$formula== 'C8H10')
    speciateData$USEME[ind] = 0
    
    ind = which(speciateData$formula== 'C8H10' & speciateData$PI == 'WARNEKE')
    fracebz =  speciateData$FinalEF[ind1]/(speciateData$FinalEF[ind1] + speciateData$FinalEF[ind2]+ + speciateData$FinalEF[ind3])
    fracox =  speciateData$FinalEF[ind2]/(speciateData$FinalEF[ind1] + speciateData$FinalEF[ind2]+ + speciateData$FinalEF[ind3])
    fracmpx =  speciateData$FinalEF[ind3]/(speciateData$FinalEF[ind1] + speciateData$FinalEF[ind2]+ + speciateData$FinalEF[ind3])
    # check for negs
    iw = which(!is.finite(fracebz))
    if (length(iw) != 0){fracebz[iw] = mean(fracebz, na.rm=TRUE)}
    iw = which(!is.finite(fracox))
    if (length(iw) != 0){fracox[iw] = mean(fracox, na.rm=TRUE)}
    iw = which(!is.finite(fracmpx))
    if (length(iw) != 0){fracmpx[iw] = mean(fracmpx, na.rm=TRUE)}
    
    #  ----- Ethylbenzene -----
    newline = speciateData[ind,]
    newline$PI = 'Warneke w/Blake+GILMAN+APEL'
    #newline$Group.2 = 'Warneke w/Blake+GILMAN+APEL'
    newline$names = 'Ethylbenzene'; newline$variable = 'Ethylbenzene'
    newline$FinalEF = newline$FinalEF*fracebz;  newline$FinalEF_mean = newline$FinalEF_mean*fracebz
    newline$FinalERtoCO = newline$FinalERtoCO*fracebz ; newline$FinalERtoCO_mean = newline$FinalERtoCO_mean*fracebz
    newline$FinalEF_sd = newline$FinalEF_sd*fracebz
   # newline$FinalEF_75 = newline$FinalEF_75*fracebz
    newline$USEME = 1
    speciateData = rbind(speciateData,newline)
    
    # ------- o-Xylene -----
    newline = speciateData[ind,]
    newline$PI = 'Warneke w/Blake+GILMAN+APEL'
    #newline$Group.2 = 'Warneke w/Blake+GILMAN+APEL'
    newline$names = 'o-Xylene'; newline$variable = 'o-Xylene'
    newline$FinalEF = newline$FinalEF*fracox;  newline$FinalEF_mean = newline$FinalEF_mean*fracox
    newline$FinalERtoCO = newline$FinalERtoCO*fracox; newline$FinalERtoCO_mean = newline$FinalERtoCO_mean*fracox
    newline$FinalEF_sd = newline$FinalEF_sd*fracox
    #newline$FinalEF_75 = newline$FinalEF_75*fracox
    newline$USEME = 1
    speciateData = rbind(speciateData,newline)
    
    # ------- m,p-Xylene -----
    newline = speciateData[ind,]
    newline$PI = 'Warneke w/Blake+GILMAN+APEL'
    #newline$Group.2 = 'Warneke w/Blake+GILMAN+APEL'
    newline$names = 'm,p-Xylene'; newline$variable = 'm,p-Xylene'
    newline$FinalEF = newline$FinalEF*fracmpx; newline$FinalEF_mean = newline$FinalEF_mean*fracmpx
    newline$FinalERtoCO = newline$FinalERtoCO*fracmpx; newline$FinalERtoCO_mean = newline$FinalERtoCO_mean*fracmpx
    newline$FinalEF_sd = newline$FinalEF_sd*fracmpx
    #newline$FinalEF_75 = newline$FinalEF_75*fracmpx
    newline$USEME = 1
    speciateData = rbind(speciateData,newline)
  }
  if (doacetone == 1){
    # -------- AcetonePropanal --------------
    ind = which(speciateData$Group.2== 'AcetonePropanal_NOAAPTR_ppbv_WARNEKE')
    speciateData$USEME[ind] = 0
    ind2 = which(speciateData$Group.2 == 'Propanal_ppt')
    speciateData$USEME[ind2] = 0
    ind3 = which(speciateData$Group.2== 'Acetone_ppt')
    speciateData$USEME[ind3] = 0
    
    fracacetone =  speciateData$FinalEF[ind3]/(speciateData$FinalEF[ind3] + speciateData$FinalEF[ind2])
    # check for negs
    iw = which(!is.finite(fracacetone))
    fracacetone[iw] = mean(fracacetone, na.rm=TRUE)
    fracpropanal = 1-fracacetone
    
    # ---- Acetone -----
    newline = speciateData[ind,]
    newline$PI = 'Warneke w/APEL'
    #newline$Group.2 = 'Warneke w/APEL'
    newline$names = 'Acetone'; newline$variable = 'Acetone'
     
    newline$FinalEF = newline$FinalEF*fracacetone;  newline$FinalEF_mean = newline$FinalEF_mean*fracacetone
    newline$FinalERtoCO = newline$FinalERtoCO*fracacetone; newline$FinalERtoCO_mean = newline$FinalERtoCO_mean*fracacetone
    newline$FinalEF_sd = newline$FinalEF_sd*fracacetone
   ## newline$FinalEF_75 = newline$FinalEF_75*fracacetone
    newline$USEME = 1
    newline$lifetime = speciateData$lifetime[ind3]
    newline$LifetimeCat = speciateData$LifetimeCat[ind3]
    speciateData = rbind(speciateData,newline)
    
    # ---- Propanal ------
    newline = speciateData[ind,]
    newline$PI = 'Warneke w/APEL'
    newline$names = 'Propanal'; newline$variable = 'Propanal'
    newline$FinalEF = newline$FinalEF*fracpropanal; newline$FinalEF_mean = newline$FinalEF_mean*fracpropanal
    newline$FinalERtoCO = newline$FinalERtoCO*fracpropanal; newline$FinalERtoCO_mean = newline$FinalERtoCO_mean*fracpropanal
    newline$FinalEF_sd = newline$FinalEF_sd*fracpropanal
  #  newline$FinalEF_75 = newline$FinalEF_75*fracpropanal
    newline$lifetime = speciateData$lifetime[ind2]
    newline$LifetimeCat = speciateData$LifetimeCat[ind2]
    
    newline$USEME = 1
    speciateData = rbind(speciateData,newline)
  }
  if (domvk == 1){
    # -------- MVK/MACR/2Butenals-------
    ind = which(speciateData$Group.2== 'MVKMAC_NOAAPTR_ppbv_WARNEKE')
    speciateData$USEME[ind] = 0
    ind2 = which(speciateData$Group.2 == 'MVK_ppt')
    speciateData$USEME[ind2] = 0
    ind3 = which(speciateData$Group.2== 'MAC_ppt')
    speciateData$USEME[ind3] = 0
    ind4 = which(speciateData$Group.2== 'x2Butenals_ppt')
    speciateData$USEME[ind4] = 0
    
    fracmvk=  speciateData$FinalEF[ind2]/
      (speciateData$FinalEF[ind2] + speciateData$FinalEF[ind3] + speciateData$FinalEF[ind4])
    fracmacr =  speciateData$FinalEF[ind3]/
      (speciateData$FinalEF[ind2] + speciateData$FinalEF[ind3] + speciateData$FinalEF[ind4])
    iw = which(!is.finite(fracmvk))
    fracmvk[iw] = mean(fracmvk, na.rm=TRUE)
    iw = which(!is.finite(fracmacr))
    fracmacr[iw] = mean(fracmacr, na.rm=TRUE)
    fracbutenal= 1-fracmvk - fracmacr
    
    # ---- MVK ----
    newline = speciateData[ind,]
    newline$PI = 'Warneke w/APEL'
    #newline$Group.2 = 'Warneke w/APEL'
    newline$names = speciateData$names[ind2[1]]
    newline$lifetime = speciateData$lifetime[ind2[1]]
    newline$LifetimeCat = speciateData$LifetimeCat[ind2[1]]
    
    newline$FinalEF = newline$FinalEF*fracmvk; newline$FinalEF_mean = newline$FinalEF_mean*fracmvk
    newline$FinalERtoCO = newline$FinalERtoCO*fracmvk; newline$FinalERtoCO_mean = newline$FinalERtoCO_mean*fracmvk
    newline$FinalEF_sd = newline$FinalEF_sd*fracmvk
   # newline$FinalEF_75 = newline$FinalEF_75*fracmvk
    newline$USEME = 1
    newline$lifetime = speciateData$lifetime[ind2[1]]
    newline$LifetimeCat = speciateData$LifetimeCat[ind2[1]]
    
    speciateData = rbind(speciateData,newline)
    
    # ---- MACR ------
    newline = speciateData[ind,]
    newline$PI = 'Warneke w/APEL'
    newline$names = speciateData$names[ind3[1]]
    newline$FinalEF = newline$FinalEF*fracmacr; newline$FinalEF_mean = newline$FinalEF_mean*fracmacr
    newline$FinalERtoCO = newline$FinalERtoCO*fracmacr; newline$FinalERtoCO_mean = newline$FinalERtoCO_mean*fracmacr
    newline$FinalEF_sd = newline$FinalEF_sd*fracmacr
  #  newline$FinalEF_75 = newline$FinalEF_75*fracmacr
    #newline$Group.2 = 'MAC_ppt'
    newline$USEME = 1
    newline$lifetime = speciateData$lifetime[ind3[1]]
    newline$LifetimeCat = speciateData$LifetimeCat[ind3[1]]
    
    speciateData = rbind(speciateData,newline)
    
    # ---- Butenal ------
    newline = speciateData[ind,]
    newline$PI = 'Warneke w/APEL'
    newline$names = speciateData$names[ind4[1]]
    newline$FinalEF = newline$FinalEF*fracbutenal; newline$FinalEF_mean = newline$FinalEF_mean*fracbutenal
    newline$FinalERtoCO = newline$FinalERtoCO*fracbutenal; newline$FinalERtoCO_mean = newline$FinalERtoCO_mean*fracbutenal
    newline$FinalEF_sd = newline$FinalEF_sd*fracbutenal
   # newline$FinalEF_75 = newline$FinalEF_75*fracbutenal
    newline$USEME = 1
    #newline$Group.2 = 'x2Butenals_ppt'
    newline$lifetime = speciateData$lifetime[ind4[1]]
    newline$LifetimeCat = speciateData$LifetimeCat[ind4[1]]
    speciateData = rbind(speciateData,newline)
  }
  if (domek == 1){
    # -------- MEK + ibutanal + butanal -------
    ind = which(speciateData$Group.2== 'C4Carbonyls_NOAAPTR_ppbv_WARNEKE' )
    speciateData$USEME[ind] = 0
    ind2 = which(speciateData$Group.2 == 'Butanal_ppt')
    speciateData$USEME[ind2] = 0
    ind3 = which(speciateData$Group.2== 'MEK_ppt')
    speciateData$USEME[ind3] = 0
    ind4 = which(speciateData$Group.2== 'iButanal_ppt')
    speciateData$USEME[ind4] = 0
    
    fracbutanal=  speciateData$FinalEF[ind2]/
      (speciateData$FinalEF[ind2] + speciateData$FinalEF[ind3] + speciateData$FinalEF[ind4])
    fracibutanal =  speciateData$FinalEF[ind4]/
      (speciateData$FinalEF[ind2] + speciateData$FinalEF[ind3] + speciateData$FinalEF[ind4])
    iw = which(!is.finite(fracbutanal))
    fracbutanal[iw] = mean(fracbutanal, na.rm=TRUE)
    iw = which(!is.finite(fracibutanal))
    fracibutanal[iw] = mean(fracibutanal, na.rm=TRUE)
    fracMEK= 1-fracibutanal - fracbutanal
    
    # ---------   Butanal -----
    newline = speciateData[ind,]
    newline$PI = 'Warneke w/APEL'
    newline$names = speciateData$names[ind2[1]]
    newline$lifetime = speciateData$lifetime[ind2[1]]
    newline$LifetimeCat = speciateData$LifetimeCat[ind2[1]]
    
    newline$FinalEF = newline$FinalEF*fracbutanal;  newline$FinalEF_mean = newline$FinalEF_mean*fracbutanal
    newline$FinalERtoCO = newline$FinalERtoCO*fracbutanal; newline$FinalERtoCO_mean = newline$FinalERtoCO_mean*fracbutanal
    newline$FinalEF_sd = newline$FinalEF_sd*fracbutanal
   # newline$FinalEF_75 = newline$FinalEF_75*fracbutanal
    newline$USEME = 1
    newline$lifetime = speciateData$lifetime[ind2[1]]
    newline$LifetimeCat = speciateData$LifetimeCat[ind2[1]]
    
    speciateData = rbind(speciateData,newline)
    
    # ------- MEK -------
    newline = speciateData[ind,]
    newline$PI = 'Warneke w/APEL'
    newline$names = speciateData$names[ind3[1]]
    newline$FinalEF = newline$FinalEF*fracMEK; newline$FinalEF_mean = newline$FinalEF_mean*fracMEK
    newline$FinalERtoCO = newline$FinalERtoCO*fracMEK; newline$FinalERtoCO_mean = newline$FinalERtoCO_mean*fracMEK
    newline$FinalEF_sd = newline$FinalEF_sd*fracMEK
  #  newline$FinalEF_75 = newline$FinalEF_75*fracMEK
    #newline$Group.2 = 'MAC_ppt'
    newline$USEME = 1
    newline$lifetime = speciateData$lifetime[ind3[1]]
    newline$LifetimeCat = speciateData$LifetimeCat[ind3[1]]
    
    speciateData = rbind(speciateData,newline)
    
    # ------------- iButanal -----------
    newline = speciateData[ind,]
    newline$PI = 'Warneke w/APEL'
    newline$names = speciateData$names[ind4[1]]
    newline$FinalEF = newline$FinalEF*fracibutanal;  newline$FinalEF_mean = newline$FinalEF_mean*fracibutanal
    newline$FinalERtoCO = newline$FinalERtoCO*fracibutanal; newline$FinalERtoCO_mean = newline$FinalERtoCO_mean*fracibutanal
    newline$FinalEF_sd = newline$FinalEF_sd*fracibutanal
   # newline$FinalEF_75 = newline$FinalEF_75*fracibutanal
    newline$USEME = 1
    newline$lifetime = speciateData$lifetime[ind4[1]]
    newline$LifetimeCat = speciateData$LifetimeCat[ind4[1]]
    speciateData = rbind(speciateData,newline)
   }
  return(speciateData)
  # -------- Methyl acetate/Ethyl formate/Hydroxyacetone -------?
}