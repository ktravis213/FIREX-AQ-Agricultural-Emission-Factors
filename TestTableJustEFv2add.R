# Do ER or do EF
library(htmlTable); library(plyr);library(magrittr)

source("Table1.R")
source("speciateSpecies.R")
require(flextable)
q1=0.25; q2=0.75
`%notin%` <- Negate(`%in%`)
doER = 0
doEF = 1
dowheatshrub =1
options(scipen = 1, digits=3)

cornfrac = table1$fires[1]/sum(table1$fires[1:4])
ricefrac = table1$fires[2]/sum(table1$fires[1:4])
soyfrac = table1$fires[3]/sum(table1$fires[1:4])
wheatfrac = table1$fires[4]/sum(table1$fires[1:4])

slashfrac = table1$fires[6]/sum(table1$fires[6:8])
pilefrac = table1$fires[7]/sum(table1$fires[6:8])
shrubfrac = table1$fires[8]/sum(table1$fires[6:8])

# 
allBOTH.filter.median = aggregate(allBOTH.filter, by=list(allBOTH.filter$variable), FUN='median', na.rm=TRUE)
propertiesSpecies = as.data.frame(cbind(allBOTH.filter.median$Group.1,
                                        allBOTH.filter.median$names,allBOTH.filter.median$PI,
                                        allBOTH.filter.median$mWs,allBOTH.filter.median$OHrate.1hz, 
                                        allBOTH.filter.median$OHrate.5hz, allBOTH.filter.median$lifetime_jval,allBOTH.filter.median$lifetime))
ind = which(allBOTH.filter.median$USEME != 0)
write.csv(propertiesSpecies[ind,],file='propertiesSpecies.csv')
library(htmlTable); library(plyr)
library(magrittr)
# --------------

ind = which(allBOTH.filter$fire == 'BlackwaterRiver') 
allBOTH.filter$fuel[ind] = 'Blackwater'
# Make larger fuel categories
allBOTH.filter$fuelORIG = allBOTH.filter$fuel
allBOTH.filter$fuel2 = allBOTH.filter$fuel
ind = which(allBOTH.filter$fuel2 == 'corn' | allBOTH.filter$fuel2 == 'soybean' | allBOTH.filter$fuel2 == 'rice' |
              allBOTH.filter$fuel2 == 'winter wheat')
allBOTH.filter$fuel2[ind] = 'agriculture'
ind = which(allBOTH.filter$fuel2 == 'pile' | allBOTH.filter$fuel2 == 'slash' | allBOTH.filter$fuel2 == 'shrub')
allBOTH.filter$fuel2[ind] = 'prescribed'
#ind = which(allBOTH.filter$fire == 'Blackwater')
#allBOTH.filter$fuel2[ind] = 'Blackwater'

#ind = which(allBOTH.filter.median.fuel$fire == 'BlackwaterRiver') 
#allBOTH.filter.median.fuel$fuel[ind] = 'Blackwater'
allBOTH.filter.median.fuel = aggregate(allBOTH.filter, by=list(allBOTH.filter$fuel,allBOTH.filter$variable), FUN='median', na.rm=TRUE)
allBOTH.filter.mean.fuel = aggregate(allBOTH.filter, by=list(allBOTH.filter$fuel,allBOTH.filter$variable),   FUN='mean', na.rm=TRUE)
allBOTH.filter.sd.fuel = aggregate(allBOTH.filter, by=list(allBOTH.filter$fuel,allBOTH.filter$variable), FUN='sd', na.rm=TRUE)

allBOTH.filter.median.fuel$FinalEF_mean = allBOTH.filter.mean.fuel$FinalEF
allBOTH.filter.median.fuel$FinalEF_sd = allBOTH.filter.sd.fuel$FinalEF
allBOTH.filter.median.fuel$FinalERtoCO_mean = allBOTH.filter.mean.fuel$FinalERtoCO*1E3
allBOTH.filter.median.fuel$FinalERtoCO_sd = allBOTH.filter.sd.fuel$FinalERtoCO*1E3
allBOTH.filter.25.fuel = aggregate(allBOTH.filter$FinalEF, by=list(allBOTH.filter$fuel,allBOTH.filter$variable), 'quantile',probs=c(q1),na.rm=TRUE)
allBOTH.filter.75.fuel = aggregate(allBOTH.filter$FinalEF, by=list(allBOTH.filter$fuel,allBOTH.filter$variable), 'quantile',probs=c(q2),na.rm=TRUE)
allBOTH.filter.25.fuelER = aggregate(allBOTH.filter$FinalERtoCO, by=list(allBOTH.filter$fuel,allBOTH.filter$variable), 'quantile',probs=c(q1),na.rm=TRUE)
allBOTH.filter.75.fuelER = aggregate(allBOTH.filter$FinalERtoCO, by=list(allBOTH.filter$fuel,allBOTH.filter$variable), 'quantile',probs=c(q2),na.rm=TRUE)

allBOTH.filter.median.fuel$FinalEF_25 = allBOTH.filter.25.fuel$x
allBOTH.filter.median.fuel$FinalEF_75 = allBOTH.filter.75.fuel$x
allBOTH.filter.median.fuel$FinalERtoCO_25 = allBOTH.filter.25.fuelER$x*1E3
allBOTH.filter.median.fuel$FinalERtoCO_75 = allBOTH.filter.75.fuelER$x*1E3

# need to recover fire, kind, formula, and names
allBOTH.filter.median.fuel$fuel = allBOTH.filter.median.fuel$Group.1
allBOTH.filter.median.fuel$variable = allBOTH.filter.median.fuel$Group.2
for (i in 1:length(allBOTH.filter.median.fuel$kind)){
  ind = which(allBOTH.filter$variable == allBOTH.filter.median.fuel$Group.2[i])
  allBOTH.filter.median.fuel$kind[i] = allBOTH.filter$kind[ind[1]]
  allBOTH.filter.median.fuel$formula[i] = allBOTH.filter$formula[ind[1]]
  allBOTH.filter.median.fuel$names[i] = allBOTH.filter$names[ind[1]]
  allBOTH.filter.median.fuel$PI[i] = allBOTH.filter$PI[ind[1]]
  # print(c(allBOTH.filter$variable[ind[1]], allBOTH.filter.median.fuel$Group.2[i]))
}
allBOTH.filter.median.fuel= getplumesANDmcebyfuel(allBOTH.filter.median.fuel, allBOTH.filter )
#allBOTH.filter.median.fuel = speciateSpecies(allBOTH.filter.median.fuel)

tmpTABLE = speciateSpecies(allBOTH.filter.median.fuel)
# Here add kludge for headers
newLine = tmpTABLE

# ok want to put C9 aromatics and monoterpenes together regardless of lifetime
ind = which(tmpTABLE$formula == 'C9H12')
tmpTABLE$LifetimeCat[ind] = 1
ind = which(tmpTABLE$formula == 'C10H16')
tmpTABLE$LifetimeCat[ind] = 1
ind =which(tmpTABLE$LifetimeCat == 2 & tmpTABLE$Category == 1)
tmpTABLE$Category[ind] = 1.5
# get rid of stuff without a category, category = 5 (oxidation products) and with USEME = 1
# decided to report H2O2
ind = which(tmpTABLE$formula == 'H2O2')
tmpTABLE$Category[ind] = 7

ind = which(tmpTABLE$formula == 'NOy')
tmpTABLE$Category[ind] = 2
ind = which(tmpTABLE$PI == 'MOORE')
tmpTABLE$Category[ind] = 4.5

ind = which(is.finite(tmpTABLE$Category) & tmpTABLE$Category != 5 & tmpTABLE$USEME > 0)
tmpTABLE=tmpTABLE[ind,]

ind = which(tmpTABLE$variable == 'Short-lived NMVOC' | tmpTABLE$variable == 'Long-lived NMVOC')
tmpTABLE$Category[ind] = 1.5
# ----- get fuel specific species -------
ind = which(tmpTABLE$Group.1 == 'corn'  )
new.data.frame = tmpTABLE[ind,]
ind = which(tmpTABLE$Group.1 == 'rice' )
new.data.frame.rice = tmpTABLE[ind,]
ind = which(tmpTABLE$Group.1 == 'soybean' )
new.data.frame.soybean = tmpTABLE[ind,]
ind = which(tmpTABLE$Group.1 == 'winter wheat' )
new.data.frame.wheat = tmpTABLE[ind,]
ind = which(tmpTABLE$Group.1 == 'grass' )
new.data.frame.grass = tmpTABLE[ind,]
ind = which(tmpTABLE$Group.1 == 'slash' )
new.data.frame.slash = tmpTABLE[ind,]
ind = which(tmpTABLE$Group.1 == 'pile' )
new.data.frame.pile = tmpTABLE[ind,]
ind = which(tmpTABLE$Group.1 == 'shrub' )
new.data.frame.shrub = tmpTABLE[ind,]
ind = which(tmpTABLE$Group.1 == 'Blackwater' )
new.data.frame.blackwater = tmpTABLE[ind,]

# ------- Make average ag --------
# but account for zeros by replacing with corn which has complete EFs
dothis=1
new.data.frame.ag=new.data.frame
new.data.frame.presc=new.data.frame
if (dothis == 1){
  new.data.frame.ag = new.data.frame
  new.data.frame.wheat2 = new.data.frame.wheat
  ind = which(!is.finite(new.data.frame.wheat2$FinalEF_mean))
  new.data.frame.wheat2$FinalEF_mean[ind] = new.data.frame$FinalEF_mean[ind]
  new.data.frame.wheat2$FinalERtoCO_mean[ind] = new.data.frame$FinalERtoCO_mean[ind]
  ind = which(!is.finite(new.data.frame.wheat2$FinalEF_mean) | !is.finite(new.data.frame.wheat$FinalEF_sd))
  new.data.frame.wheat2$FinalEF_sd[ind] = new.data.frame$FinalEF_sd[ind]
  new.data.frame.wheat2$FinalERtoCO_sd[ind] = new.data.frame$FinalERtoCO_sd[ind]
  
  new.data.frame.soybean2 = new.data.frame.soybean
  ind = which(!is.finite(new.data.frame.soybean2$FinalEF_mean))
  new.data.frame.soybean2$FinalEF_mean[ind] = new.data.frame$FinalEF_mean[ind]
  new.data.frame.soybean2$FinalERtoCO_mean[ind] = new.data.frame$FinalERtoCO_mean[ind]
  ind = which(!is.finite(new.data.frame.soybean2$FinalEF_mean) | !is.finite(new.data.frame.soybean2$FinalEF_sd))
  new.data.frame.soybean2$FinalEF_sd[ind] = new.data.frame$FinalEF_sd[ind]
  new.data.frame.soybean2$FinalERtoCO_sd[ind] = new.data.frame$FinalERtoCO_sd[ind]
  
  new.data.frame.rice2 = new.data.frame.rice
  ind = which(!is.finite(new.data.frame.rice2$FinalEF_mean))
  new.data.frame.rice2$FinalEF_mean[ind] = new.data.frame$FinalEF_mean[ind]
  new.data.frame.rice2$FinalERtoCO_mean[ind] = new.data.frame$FinalEF_mean[ind]
  ind = which(!is.finite(new.data.frame.rice2$FinalEF_mean) |!is.finite(new.data.frame.rice2$FinalEF_sd) )
  new.data.frame.rice2$FinalEF_sd[ind] = new.data.frame$FinalEF_sd[ind]
  new.data.frame.rice2$FinalERtoCO_sd[ind] = new.data.frame$FinalERtoCO_sd[ind]
  
  new.data.frame.ag$FinalEF_mean = new.data.frame$FinalEF_mean*cornfrac  + new.data.frame.rice2$FinalEF_mean*ricefrac + new.data.frame.soybean2$FinalEF_mean*soyfrac + new.data.frame.wheat2$FinalEF_mean*wheatfrac
  new.data.frame.ag$FinalERtoCO_mean = new.data.frame$FinalERtoCO_mean*cornfrac  + new.data.frame.rice2$FinalERtoCO_mean*ricefrac  + new.data.frame.soybean2$FinalERtoCO_mean*soyfrac + new.data.frame.wheat2$FinalERtoCO_mean*wheatfrac
  new.data.frame.ag$COUNT_EFFINAL    = new.data.frame$COUNT_EFFINAL+new.data.frame.rice2$COUNT_EFFINAL + new.data.frame.soybean2$COUNT_EFFINAL + new.data.frame.wheat2$COUNT_EFFINAL
  new.data.frame.ag$FinalEF_sd   = sqrt((((cornfrac*new.data.frame$FinalEF_sd)^2 +
                                            (ricefrac*new.data.frame.rice2$FinalEF_sd)^2  +
                                                (soyfrac*new.data.frame.soybean2$FinalEF_sd)^2 + 
                                            (wheatfrac*new.data.frame.wheat2$FinalEF_sd)^2)))
  new.data.frame.ag$FinalERtoCO_sd   = sqrt((((cornfrac*new.data.frame$FinalERtoCO_sd)^2 + (ricefrac*new.data.frame.rice2$FinalERtoCO_sd)^2  +
                                               (soyfrac*new.data.frame.soybean2$FinalERtoCO_sd)^2 + (wheatfrac*new.data.frame.wheat2$FinalERtoCO_sd)^2)))
  # ---- Make average presc ------------
  new.data.frame.presc = new.data.frame.pile
  # fill missing pile
  new.data.frame.pile2 = new.data.frame.pile
  ind = which(!is.finite(new.data.frame.pile$FinalEF_mean) & is.finite(new.data.frame.slash$FinalEF_mean))
  new.data.frame.pile2$FinalEF_mean[ind] = new.data.frame.slash$FinalEF_mean[ind]
  new.data.frame.pile2$FinalERtoCO_mean[ind] = new.data.frame.slash$FinalERtoCO_mean[ind]
  ind = which(!is.finite(new.data.frame.pile$FinalEF_sd) & is.finite(new.data.frame.slash$FinalEF_sd))
  new.data.frame.pile2$FinalEF_sd[ind] = new.data.frame.slash$FinalEF_sd[ind]
  new.data.frame.pile2$FinalERtoCO_sd[ind] = new.data.frame.slash$FinalERtoCO_sd[ind]
  # fill missing slash
  new.data.frame.slash2 = new.data.frame.slash
  ind = which(is.finite(new.data.frame.pile2$FinalEF_mean) & !is.finite(new.data.frame.slash$FinalEF_mean))
  new.data.frame.slash2$FinalEF_mean[ind] = new.data.frame.pile$FinalEF_mean[ind]
  new.data.frame.slash2$FinalERtoCO_mean[ind] = new.data.frame.pile$FinalERtoCO_mean[ind]
  ind = which(is.finite(new.data.frame.pile2$FinalEF_sd) & !is.finite(new.data.frame.slash$FinalEF_sd))
  new.data.frame.slash2$FinalEF_sd[ind] = new.data.frame.pile$FinalEF_sd[ind]
  new.data.frame.slash2$FinalERtoCO_sd[ind] = new.data.frame.pile$FinalERtoCO_sd[ind]
  # fill missing shrub
  new.data.frame.shrub2 = new.data.frame.shrub
       # with average of pile + slash
  tmp = new.data.frame.pile2
  tmp$FinalEF_mean = (pilefrac*new.data.frame.pile2$FinalEF_mean +slashfrac*new.data.frame.slash2$FinalEF_mean)
  tmp$FinalERtoCO_mean = (pilefrac*new.data.frame.pile2$FinalERtoCO_mean +slashfrac*new.data.frame.slash2$FinalERtoCO_mean)
  tmp$FinalEF_sd= sqrt(((pilefrac*new.data.frame.pile2$FinalEF_sd)^2 +(slashfrac*new.data.frame.slash2$FinalEF_sd)^2))
  tmp$FinalERtoCO_sd= sqrt((pilefrac*new.data.frame.pile2$FinalERtoCO_sd)^2 + (slashfrac*new.data.frame.slash2$FinalERtoCO_sd)^2)
  
  ind = which(is.finite(tmp$FinalEF_mean) & !is.finite(new.data.frame.shrub2$FinalEF_mean))
  new.data.frame.shrub2$FinalEF_mean[ind] = (tmp$FinalEF_mean[ind])
  new.data.frame.shrub2$FinalERtoCO_mean[ind] = (tmp$FinalERtoCO_mean[ind])
  ind = which(is.finite(tmp$FinalEF_sd) & !is.finite(new.data.frame.shrub2$FinalEF_sd))
  new.data.frame.shrub2$FinalEF_sd[ind] = tmp$FinalEF_sd[ind]
  new.data.frame.shrub2$FinalERtoCO_sd[ind] = tmp$FinalERtoCO_sd[ind]

  new.data.frame.presc$FinalEF_mean = new.data.frame.pile2$FinalEF_mean*pilefrac+new.data.frame.slash2$FinalEF_mean*slashfrac + new.data.frame.shrub2$FinalEF_mean*shrubfrac
  new.data.frame.presc$FinalEF_sd = sqrt((((pilefrac*new.data.frame.pile2$FinalEF_sd)^2 +(slashfrac*new.data.frame.slash2$FinalEF_sd)^2  +
                                                 (shrubfrac*new.data.frame.shrub2$FinalEF_sd)^2 )))
  new.data.frame.presc$FinalERtoCO_mean = new.data.frame.pile2$FinalERtoCO_mean*pilefrac+new.data.frame.slash2$FinalERtoCO_mean*slashfrac + new.data.frame.shrub2$FinalERtoCO_mean*shrubfrac
  new.data.frame.presc$FinalERtoCO_sd = sqrt((((pilefrac*new.data.frame.pile2$FinalERtoCO_sd)^2 +
                                                 (slashfrac*new.data.frame.slash2$FinalERtoCO_sd)^2)  +
                                                 (shrubfrac*new.data.frame.shrub2$FinalERtoCO_sd)^2 ))
  new.data.frame.presc$COUNT_EFFINAL = new.data.frame.pile2$COUNT_EFFINAL+new.data.frame.slash2$COUNT_EFFINAL + new.data.frame.shrub2$COUNT_EFFINAL
}
#  --- fix PIs ------- 
ind= which(new.data.frame$PI == "ppt")
new.data.frame$PI[ind] = "APEL"
ind= which(new.data.frame$PI == "ppb")
new.data.frame$PI[ind] = "DISKIN"
setupdataEF = as.data.frame(cbind(new.data.frame$Category, new.data.frame$LifetimeCat,new.data.frame$mWs,
                                  new.data.frame$names,new.data.frame$formula,new.data.frame$PI,
                                  new.data.frame$FinalEF_mean, new.data.frame$FinalEF_sd, new.data.frame$COUNT_EFFINAL,
                                  new.data.frame.rice$FinalEF_mean, new.data.frame.rice$FinalEF_sd,         new.data.frame.rice$COUNT_EFFINAL,
                                  new.data.frame.soybean$FinalEF_mean, new.data.frame.soybean$FinalEF_sd, new.data.frame.soybean$COUNT_EFFINAL,  
                                  new.data.frame.ag$FinalEF_mean, new.data.frame.ag$FinalEF_sd, new.data.frame.ag$COUNT_EFFINAL,  
                                  
                                  
                                  new.data.frame.slash$FinalEF_mean,   new.data.frame.slash$FinalEF_sd,     new.data.frame.slash$COUNT_EFFINAL,
                                  new.data.frame.pile$FinalEF_mean,   new.data.frame.pile$FinalEF_sd,         new.data.frame.pile$COUNT_EFFINAL,
                                  new.data.frame.presc$FinalEF_mean,   new.data.frame.presc$FinalEF_sd,         new.data.frame.presc$COUNT_EFFINAL,
                                  
                                  
                                  new.data.frame.grass$FinalEF_mean,   new.data.frame.grass$FinalEF_sd,     new.data.frame.grass$COUNT_EFFINAL,
                                  new.data.frame.wheat$FinalEF_mean,   new.data.frame.wheat$FinalEF_sd,     new.data.frame.wheat$COUNT_EFFINAL,
                                  new.data.frame.shrub$FinalEF_mean,   new.data.frame.shrub$FinalEF_sd,     new.data.frame.shrub$COUNT_EFFINAL,
                                  new.data.frame.blackwater$FinalEF_mean,   new.data.frame.blackwater$FinalEF_sd,     new.data.frame.blackwater$COUNT_EFFINAL))


# Save for comparison to Akagi and Andreae
if (doEF == 1){save(setupdataEF,file='setupdataEF.RData')}
# Ok, don't want to order nitrogen compounds by lifetimecat
ind = which(setupdataEF$V1 == 2 | setupdataEF$V1 == 3)
setupdataEF$V2[ind] = 1
setupdataEF =setupdataEF[ with(setupdataEF, order(V1,V2, as.numeric(V3) )),]#mWs,names, )),]
#setupdataEF =setupdataEF[ with(setupdataEF, order(Category,LifetimeCat, -FinalEF )),]#mWs,names, )),]
setupdataEF= subset(setupdataEF, select = -c( V1, V2,V3))
# ------ ERs -----

setupdataER = as.data.frame(cbind(new.data.frame$Category, new.data.frame$LifetimeCat,new.data.frame$mWs,
                                  new.data.frame$names,new.data.frame$formula,new.data.frame$PI,
                                  new.data.frame$FinalERtoCO_mean, new.data.frame$FinalERtoCO_sd, new.data.frame$COUNT_EFFINAL,
                                  new.data.frame.rice$FinalERtoCO_mean, new.data.frame.rice$FinalERtoCO_sd,         new.data.frame.rice$COUNT_EFFINAL,
                                  new.data.frame.soybean$FinalERtoCO_mean, new.data.frame.soybean$FinalERtoCO_sd, new.data.frame.soybean$COUNT_EFFINAL,  
                                  new.data.frame.ag$FinalERtoCO_mean, new.data.frame.ag$FinalERtoCO_sd, new.data.frame.ag$COUNT_EFFINAL,  
                                  
                                  
                                  new.data.frame.slash$FinalERtoCO_mean,   new.data.frame.slash$FinalERtoCO_sd,     new.data.frame.slash$COUNT_EFFINAL,
                                  new.data.frame.pile$FinalERtoCO_mean,   new.data.frame.pile$FinalERtoCO_sd,         new.data.frame.pile$COUNT_EFFINAL,
                                  new.data.frame.presc$FinalERtoCO_mean,   new.data.frame.presc$FinalERtoCO_sd,         new.data.frame.presc$COUNT_EFFINAL,
                                  
                                  
                                  new.data.frame.grass$FinalERtoCO_mean,   new.data.frame.grass$FinalERtoCO_sd,     new.data.frame.grass$COUNT_EFFINAL,
                                  
                                  new.data.frame.wheat$FinalERtoCO_mean,   new.data.frame.wheat$FinalERtoCO_sd,     new.data.frame.wheat$COUNT_EFFINAL,
                                  new.data.frame.shrub$FinalERtoCO_mean,   new.data.frame.shrub$FinalERtoCO_sd,     new.data.frame.shrub$COUNT_EFFINAL,
                                  new.data.frame.blackwater$FinalERtoCO_mean,   new.data.frame.blackwater$FinalERtoCO_sd,     new.data.frame.blackwater$COUNT_EFFINAL))

# Ok, don't want to order nitrogen compounds by lifetimecat
ind = which(setupdataER$V1 == 2 | setupdataER$V1 == 3)
setupdataER$V2[ind] = 1
setupdataER =setupdataER[ with(setupdataER, order(V1,V2, as.numeric(V3) )),]#mWs,names, )),]
#setupdataER =setupdataER[ with(setupdataER, order(Category,LifetimeCat, -FinalER )),]#mWs,names, )),]
setupdataER= subset(setupdataER, select = -c( V1, V2,V3))

# ------------
#rownames(setupdata) = new.data.frame$Group.2
#colnames(setupdata)=  c("EF, g/kg","25%","75%","ERtoCO, ppt/ppb","n", "EF, g/kg","25%","75%","ERtoCO, ppt/ppb","n")
save(setupdataEF, file='setupdataEF.RData')
save(setupdataER, file='setupdataER.RData')

require(stringr)
if (doEF == 1){ setupdata = setupdataEF}
if (doER == 1){ setupdata = setupdataER}

setupdata2 = setupdata
dd=dim(setupdata2)
for (i in 1:length(setupdata2$V6)){
  for (j in 4:dd[2]){
    tt = as.numeric(setupdata2[i,j])
    if(is.finite(tt)){
      if (tt == 0){
        setupdata2[i,j] = "NA"
      }
    }
    if (j != 6 & j != 9 & j != 12 & j != 15 & j != 18 & j != 21 & j != 24 & j != 27 & j != 30 & j != 33 & j != 36){
      tt = as.numeric(setupdata2[i,j])
      if (is.na(tt)){
        setupdata2[i,j] = "NA"
        #  } else if {
        #    setupdata2[i,j]=sprintf('%1.3g',tt)
        # }
      } else if (tt >= 1E13){
        setupdata2[i,j] = sprintf('%1.2e',tt)
      } else if (tt >= 50 & tt <1E13){
        setupdata2[i,j] = sprintf('%1.0f',tt)
      } else if (tt >=1 & tt <= 50){
        setupdata2[i,j] = sprintf('%1.2f',tt)
      }else if (tt >= 0.001 & tt < 1){
        setupdata2[i,j] = sprintf('%0.3f',tt)
      }else {#if (tt <= 0.001){
        setupdata2[i,j] = sprintf('%0.2e',tt)
      } 
    }
    #print(c(i,j,tt, setupdata[i,j]))
  }
}
# make sure formats match for sd
for (i in 1:length(setupdata2$V7)){
  
  if (setupdata2$V8[i] != "NA" & setupdata2$V7[i] != "NA"){
    if (as.numeric(setupdata2$V8[i]) < 0.001 & as.numeric(setupdata2$V7[i]) >= 0.001){
      setupdata2$V8[i] = sprintf('%0.4f',as.numeric(setupdata2$V8[i]))
    }
  }
  if (setupdata2$V11[i] != "NA" & setupdata2$V10[i] != "NA"){
    if (as.numeric(setupdata2$V11[i]) < 0.001 & as.numeric(setupdata2$V10[i]) >= 0.001){
      setupdata2$V11[i] = sprintf('%0.4f',as.numeric(setupdata2$V11[i]))
    }
  }
  if (setupdata2$V14[i] != "NA" & setupdata2$V13[i] != "NA"){
    if (as.numeric(setupdata2$V14[i]) < 0.001 & as.numeric(setupdata2$V13[i]) >= 0.001){
      setupdata2$V14[i] = sprintf('%0.4f',as.numeric(setupdata2$V14[i]))
    }
  }
  if (setupdata2$V17[i] != "NA" & setupdata2$V16[i] != "NA"){
    if (as.numeric(setupdata2$V17[i]) < 0.001 & as.numeric(setupdata2$V16[i]) >= 0.001){
      setupdata2$V17[i] = sprintf('%0.4f',as.numeric(setupdata2$V17[i]))
    }
  }
  if (setupdata2$V20[i] != "NA" & setupdata2$V19[i] != "NA"){
    if (as.numeric(setupdata2$V20[i]) < 0.001 & as.numeric(setupdata2$V19[i]) >= 0.001){
      setupdata2$V20[i] = sprintf('%0.4f',as.numeric(setupdata2$V20[i]))
    }
  }
  if (setupdata2$V23[i] != "NA" & setupdata2$V22[i] != "NA"){
    if (as.numeric(setupdata2$V23[i]) < 0.001 & as.numeric(setupdata2$V22[i]) >= 0.001){
      setupdata2$V23[i] = sprintf('%0.4f',as.numeric(setupdata2$V23[i]))
    }
  }
  if (setupdata2$V25[i] != "NA" & setupdata2$V26[i] != "NA"){
    if (as.numeric(setupdata2$V26[i]) < 0.001 & as.numeric(setupdata2$V25[i]) >= 0.001){
      setupdata2$V26[i] = sprintf('%0.4f',as.numeric(setupdata2$V26[i]))
    } 
  }
  if (setupdata2$V28[i] != "NA" & setupdata2$V29[i] != "NA"){
    if (as.numeric(setupdata2$V29[i]) < 0.001 & as.numeric(setupdata2$V28[i]) >= 0.001){
      setupdata2$V29[i] = sprintf('%0.4f',as.numeric(setupdata2$V29[i]))
    }
  }
  if (setupdata2$V31[i] != "NA" & setupdata2$V32[i] != "NA"){
    if (as.numeric(setupdata2$V32[i]) < 0.001 & as.numeric(setupdata2$V31[i]) >= 0.001){
      setupdata2$V32[i] = sprintf('%0.4f',as.numeric(setupdata2$V32[i]))
    }
  }
  if (setupdata2$V34[i] != "NA" & setupdata2$V35[i] != "NA"){
    if (as.numeric(setupdata2$V35[i]) < 0.001 & as.numeric(setupdata2$V34[i]) >= 0.001){
      setupdata2$V35[i] = sprintf('%0.4f',as.numeric(setupdata2$V35[i]))
    }
  }
}
# Clean up format
dd=dim(setupdata2)
for (i in 1:length(setupdata2$V6)){
  for (j in 4:dd[2]){
    setupdata2[i,j]= str_replace_all(setupdata2[i,j],"0.0e+00","")
    setupdata2[i,j]= str_replace_all(setupdata2[i,j],"e-0","e-")
  }
}

# fixcolnames
tct = seq(3,dd[2]+1)
colnames(setupdata2) = paste0('V',tct)
ind = which(setupdata2$V6 != "NA" &  is.finite(as.numeric(setupdata2$V7) ))
setupdata2$V6[ind] = paste(as.numeric(setupdata2$V6[ind])," (",as.numeric(setupdata2$V7[ind]),")", sep='')
setupdata2= subset(setupdata2, select = -c(V7) )
# -------- V8 is count
ind = which(setupdata2$V9 != "NA" &  is.finite(as.numeric(setupdata2$V10) ))
setupdata2$V9[ind] = paste(setupdata2$V9[ind]," (",setupdata2$V10[ind],")", sep='')
setupdata2= subset(setupdata2, select = -c(V10) )
# -------- V11 is count
ind = which(setupdata2$V12 != "NA" &  is.finite(as.numeric(setupdata2$V13) ))
setupdata2$V12[ind] = paste(setupdata2$V12[ind]," (",setupdata2$V13[ind],")", sep='')
setupdata2= subset(setupdata2, select = -c(V13) )
# -------- V14 is count
ind = which(setupdata2$V15 != "NA" & is.finite(as.numeric(setupdata2$V16) ))
setupdata2$V15[ind] = paste(setupdata2$V15[ind]," (",setupdata2$V16[ind],")", sep='')
setupdata2= subset(setupdata2, select = -c(V16) )
# -------- V17 is count
ind = which(setupdata2$V18 != "NA" & is.finite(as.numeric(setupdata2$V19) ))
setupdata2$V18[ind] = paste(setupdata2$V18[ind]," (",setupdata2$V19[ind],")", sep='')
setupdata2= subset(setupdata2, select = -c(V19) )
# -------- V20 is count
ind = which(setupdata2$V21 != "NA" &  is.finite(as.numeric(setupdata2$V22) ))
setupdata2$V21[ind] = paste(setupdata2$V21[ind]," (",setupdata2$V22[ind], ")",sep='')
setupdata2= subset(setupdata2, select = -c(V22) )
# -------- V23 is count
ind = which(setupdata2$V24 != "NA" &  is.finite(as.numeric(setupdata2$V25) ))
setupdata2$V24[ind] = paste(setupdata2$V24[ind]," (",setupdata2$V25[ind], ")",sep='')
setupdata2= subset(setupdata2, select = -c(V25) )
# -------- V26 is count
ind = which(setupdata2$V27 != "NA" &  is.finite(as.numeric(setupdata2$V28) ))
setupdata2$V27[ind] = paste(setupdata2$V27[ind]," (",setupdata2$V28[ind], ")",sep='')
setupdata2= subset(setupdata2, select = -c(V28) )
# -------- V29 is count
ind = which(setupdata2$V30 != "NA" &  is.finite(as.numeric(setupdata2$V31) ))
setupdata2$V30[ind] = paste(setupdata2$V30[ind]," (",setupdata2$V31[ind], ")",sep='')
setupdata2= subset(setupdata2, select = -c(V31) )
# -------- V33 is count
ind = which(setupdata2$V33 != "NA" &  is.finite(as.numeric(setupdata2$V34) ))
setupdata2$V33[ind] = paste(setupdata2$V33[ind]," (",setupdata2$V34[ind], ")",sep='')
setupdata2= subset(setupdata2, select = -c(V34) )
# -------- V36 is count
ind = which(setupdata2$V36 != "NA" &  is.finite(as.numeric(setupdata2$V37) ))
setupdata2$V36[ind] = paste(setupdata2$V36[ind]," (",setupdata2$V37[ind], ")",sep='')
setupdata2= subset(setupdata2, select = -c(V37) )

setHtmlTableTheme(theme = "Google docs")

tt=unique(new.data.frame$Category)
test = aggregate(new.data.frame$Category, by=list(new.data.frame$Category), FUN='sum', na.rm=TRUE)
test$val = test$x/test$Group.1; test$val[1] = 3
output <- 
  matrix(nrow=length(new.data.frame$Group.2),ncol = 20, byrow = TRUE)
print('Fix instrument names ')
# ------ Fix instrument names -----------
ind = which(setupdata2$V5 == 'BLAKE')
setupdata2$V5[ind] = 'WAS'
setupdata2$V5[3] = 'NIR spect.'
ind = which(setupdata2$V5 == 'DISKIN')
setupdata2$V5[ind] = 'DACOM'
ind = which(setupdata2$V5 == 'JIMENEZ')
setupdata2$V5[ind] = 'AMS'
ind = which(setupdata2$V5 == 'ppt')
setupdata2$V5[ind] = 'TOGA'
ind = which(setupdata2$V5 == 'ppt+BLAKE+GILMAN')
setupdata2$V5[ind] = 'TOGA, iWAS, WAS'
ind = which(setupdata2$V5 == 'BLAKE+ppt+GILMAN')
setupdata2$V5[ind] = 'TOGA, iWAS, WAS'
ind = which(setupdata2$V5 == 'GILMAN+ppt')
setupdata2$V5[ind] = 'TOGA, iWAS'
ind = which(setupdata2$V5 == 'ppt+GILMAN')
setupdata2$V5[ind] = 'TOGA, iWAS'
ind = which(setupdata2$V5 == 'GILMAN+BLAKE')
setupdata2$V5[ind] = 'iWAS, WAS'
ind = which(setupdata2$V5 == 'BLAKE+GILMAN')
setupdata2$V5[ind] = 'iWAS, WAS'
ind = which(setupdata2$V5 == 'APEL+BLAKE')
setupdata2$V5[ind] = 'TOGA, WAS'
ind = which(setupdata2$V5 == 'BLAKE+ppt')
setupdata2$V5[ind] = 'TOGA, WAS'
ind = which(setupdata2$V5 == 'ppt+BLAKE')
setupdata2$V5[ind] = 'TOGA, WAS'
ind = which(setupdata2$V5 == 'ppt+BLAKE+WARNEKE')
setupdata2$V5[ind] = 'PTRMS*, TOGA, WAS'
ind = which(setupdata2$V5 == 'VERES')
setupdata2$V5[ind] = 'NOAA CIMS'
ind = which(setupdata2$V5 == 'WISTHALER')
setupdata2$V5[ind] = 'OSLO PTRMS'
ind = which(setupdata2$V5 == 'FRIED')
setupdata2$V5[ind] = 'CAMS'
ind = which(setupdata2$V5 == 'GILMAN')
setupdata2$V5[ind] = 'iWAS'
ind = which(setupdata2$V5 == 'FRIED+HANISCO')
setupdata2$V5[ind] = 'CAMS, ISAF'
ind = which(setupdata2$V5 == 'GILMAN+FRIED')
setupdata2$V5[ind] = 'iWAS, CAMS'
ind = which(setupdata2$V5 == 'VERES+WARKEKE')
setupdata2$V5[ind] = 'NOAA CIMS, PTRMS*'
ind = which(setupdata2$V5 == 'WARNEKE' | setupdata2$V5 == 'ppbv')
setupdata2$V5[ind] = 'PTRMS*'
ind = which(setupdata2$V5 == 'WOMACK')
setupdata2$V5[ind] = 'ACES'
ind = which(setupdata2$V5 == 'WOMACK+RYERSON')
setupdata2$V5[ind] = 'ACES, NOAA NOyO3'
ind = which(setupdata2$V5 == 'WOMACK+VERES')
setupdata2$V5[ind] = 'ACES, NOAA CIMS'
ind = which(setupdata2$V5 == 'Warneke w/Blake+GILMAN+APEL')
setupdata2$V5[ind] = 'PTRMS*, TOGA+WAS+iWAS spec.'
ind = which(setupdata2$V5 == 'Warneke w/APEL')
setupdata2$V5[ind] = 'PTRMS*, TOGA spec.'
ind = which(setupdata2$V5 == 'WENNBERG')
setupdata2$V5[ind] = 'CIT-CIMS'
ind = which(setupdata2$V5 == 'WARNEKE+ppt+GILMAN')
setupdata2$V5[ind] = 'PTRMS*, TOGA, iWAS'
ind = which(setupdata2$V5 == 'WARNEKE+BLAKE+ppt')
setupdata2$V5[ind] = 'PTRMS*, TOGA, WAS'
ind = which(setupdata2$V5 == 'VERES+WENNBERG+WARNEKE+ppt')
setupdata2$V5[ind] = 'NOAA CIMS, CIT-CIMS, PTRMS*, TOGA'
ind = which(setupdata2$V5 == 'VERES+WENNBERG+ppt+WARNEKE')
setupdata2$V5[ind] = 'NOAA CIMS, CIT-CIMS, PTRMS*, TOGA'
ind = which(setupdata2$V5 == 'VERES+WARNEKE')
setupdata2$V5[ind] = 'NOAA CIMS, PTRMS*'
ind = which(setupdata2$V5 == 'ROLLINS+RYERSON')
setupdata2$V5[ind] = 'NOAA LIF, NOAA NOyO3'
ind = which(setupdata2$V5 == 'ROLLINS')
setupdata2$V5[ind] = 'NOAA LIF'
ind = which(setupdata2$V5 == 'RYERSON')
setupdata2$V5[ind] = 'NOAA NOyO3'
ind = which(setupdata2$V5 == 'SCHWARZ')
setupdata2$V5[ind] = 'NOAA SP2'
ind = which(setupdata2$V5 == 'APEL')
setupdata2$V5[ind] = 'TOGA'
ind = which(setupdata2$V5 == 'BLAKE+GILMAN+FRIED')
setupdata2$V5[ind] = 'CAMS, WAS, iWAS'
ind = which(setupdata2$V5 == 'VERES+WENNBERG+ppt')
setupdata2$V5[ind] = 'NOAA CIMS, CIT-CIMS, TOGA'
ind = which(setupdata2$V5 == 'WARNEKE+ppt+GILMAN+BLAKE')
setupdata2$V5[ind] = 'PTRMS*, TOGA, WAS, iWAS'
headerR =  c("Names","Formula","Instrument","EF1, g/kg","n1", "EF2, g/kg","n2","EF3, g/kg","n3","EF4, g/kg","n4","EF5, g/kg","n5","EF6, g/kg","n6","EF7, g/kg","n7","EF8, g/kg", "n8","EF9, g/kg", "n9","EF10, g/kg", "n10","EF11, g/kg", "n11")
setupdata=setupdata2

colnames(setupdata) = headerR

# --------------------------------------- Main table ---------------
ind = which(setupdata$Names == 'Carbon Dioxide')
newrow = setupdata[ind,]; newrow$`EF1, g/kg`=''; newrow$`EF2, g/kg`='';newrow$`EF3, g/kg`=''; newrow$`EF4, g/kg`=''
newrow$`EF5, g/kg`='';newrow$`EF6, g/kg` = ''; newrow$`EF7, g/kg`='';newrow$`EF8, g/kg` = ''; newrow$Formula=''; newrow$Instrument='';newrow$Names = 'Short-lived NMVOC'
newrow$n1=''; newrow$n2=''; newrow$n3 = ''; newrow$n4=''; newrow$n5=''; newrow$n6 = '' ; newrow$n7=''; newrow$n8=''
setupdata = rbind(setupdata[1:ind,],newrow,setupdata[-(1:ind),])
ind = which(setupdata$Names == 'Ethyne'); ind = ind-1
newrow = setupdata[ind,]; newrow$`EF1, g/kg`=''; newrow$`EF2, g/kg`='';newrow$`EF3, g/kg`=''; newrow$`EF4, g/kg`=''
newrow$`EF5, g/kg`='';newrow$`EF6, g/kg` = ''; newrow$`EF7, g/kg`='';newrow$`EF8, g/kg` = '';  newrow$Formula=''; newrow$Instrument='';newrow$Names = 'Long-lived NMVOC'
newrow$n1=''; newrow$n2=''; newrow$n3 = ''; newrow$n4=''; newrow$n5=''; newrow$n6 = '' ; newrow$n7=''; newrow$n8=''
setupdata = rbind(setupdata[1:ind,],newrow,setupdata[-(1:ind),])
ind = which(setupdata$Names == 'Methanethiol'); ind = ind-1
newrow = setupdata[ind,]; newrow$`EF1, g/kg`=''; newrow$`EF2, g/kg`='';newrow$`EF3, g/kg`=''; newrow$`EF4, g/kg`=''
newrow$`EF5, g/kg`='';newrow$`EF6, g/kg` = ''; newrow$`EF7, g/kg`='';newrow$`EF8, g/kg` = '';  newrow$Formula=''; newrow$Instrument='';newrow$Names = 'Sulfur-containing Species'
newrow$n1=''; newrow$n2=''; newrow$n3 = ''; newrow$n4=''; newrow$n5=''; newrow$n6 = '' ; newrow$n7=''; newrow$n8=''
setupdata = rbind(setupdata[1:ind,],newrow,setupdata[-(1:ind),])
ind = which(setupdata$Names == 'Black carbon'); ind = ind-1
newrow = setupdata[ind,]; newrow$`EF1, g/kg`=''; newrow$`EF2, g/kg`='';newrow$`EF3, g/kg`=''; newrow$`EF4, g/kg`=''
newrow$`EF5, g/kg`='';newrow$`EF6, g/kg` = ''; newrow$`EF7, g/kg`='';newrow$`EF8, g/kg` = '';  newrow$Formula=''; newrow$Instrument='';newrow$Names = 'Aerosols'
newrow$n1=''; newrow$n2=''; newrow$n3 = ''; newrow$n4=''; newrow$n5=''; newrow$n6 = '' ; newrow$n7=''; newrow$n8=''
setupdata = rbind(setupdata[1:ind,],newrow,setupdata[-(1:ind),])
ind = which(setupdata$Names == 'Methyl chloride'); ind = ind-1
newrow = setupdata[ind,]; newrow$`EF1, g/kg`=''; newrow$`EF2, g/kg`='';newrow$`EF3, g/kg`=''; newrow$`EF4, g/kg`=''
newrow$`EF5, g/kg`='';newrow$`EF6, g/kg` = ''; newrow$`EF7, g/kg`='';newrow$`EF8, g/kg` = '';  newrow$Formula=''; newrow$Instrument='';newrow$Names = 'Halogenated Species'
newrow$n1=''; newrow$n2=''; newrow$n3 = ''; newrow$n4=''; newrow$n5=''; newrow$n6 = '' ; newrow$n7=''; newrow$n8=''
setupdata = rbind(setupdata[1:ind,],newrow,setupdata[-(1:ind),])
ind = which(setupdata$Names == 'Hydrogen cyanide'); ind = ind-1
newrow = setupdata[ind,]; newrow$`EF1, g/kg`=''; newrow$`EF2, g/kg`='';newrow$`EF3, g/kg`=''; newrow$`EF4, g/kg`=''
newrow$`EF5, g/kg`='';newrow$`EF6, g/kg` = ''; newrow$`EF7, g/kg`='';newrow$`EF8, g/kg` = '';  newrow$Formula=''; newrow$Instrument='';newrow$Names = 'Nitrogen-containing Species'
newrow$n1=''; newrow$n2=''; newrow$n3 = ''; newrow$n4=''; newrow$n5=''; newrow$n6 = '' ; newrow$n7=''; newrow$n8=''
setupdata = rbind(setupdata[1:ind,],newrow,setupdata[-(1:ind),])
# --------  Main table -----
ft <- flextable(setupdata)
ft <- theme_vanilla(ft)
ft <- add_header_row(ft,
                     colwidths = c(1,1,1,2,2,2,2,2,2,2,2,2,2,2),
                     values =  c("A","B","C","Corn", "Rice","Soybean","Ag","Slash","Pile","Presc","Grassland","Winter wheat","Shrubland","Blackwater"))
ft <- add_header_row(ft,
                     colwidths = c(1,1,1,8,8,6),
                     values = c('','','','Agricultural Residue','Land Clearing','Other'))

ft = fontsize(ft,part="body", size=6)
ft = fontsize(ft,part="header", size=6)
ft <- line_spacing(ft, space = 1, part = "all")
ft = align(ft,part="header",align="center")
ww=.7; ww2= 0.3
sc = 1

ft = width(ft,j=c(2,3,4,6,8,10,12,14,16,18),width = ww*sc)
ft = width(ft,j=c(5,7,9,11,13,15,17,19),width = ww2*sc)
ft = width(ft,j=1,width = 1)

#ft <-  autofit(ft)
#save_as_docx(ft,path='/Users/ktravis1/OneDrive - NASA/FIREX/test.docx')
#ft = fit_to_width(ft,max_width = 2, max_iter = 5)
sect_properties <- prop_section(
  page_size = page_size(orient = "landscape",
                        width = 8.3, height = 11.7),
  type = "continuous",
  page_margins = page_mar()
)

if (doEF == 1){ save_as_docx(`Emission Factors` = ft,  path ='testEFmeanadd.docx', pr_section = sect_properties)}
if (doER == 1){ save_as_docx(`Emission Ratios` = ft,  path ='testERmeanadd.docx', pr_section = sect_properties)}
