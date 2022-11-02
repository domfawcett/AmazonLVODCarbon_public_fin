#Author: D. Fawcett
#
#Description:
#This code calculates correlations between the L-VOD AGC and the modelled AGC change including different processes (deforestation, degradation, SF growth, IF change)
#
#Outputs:
#-Correlation statistics and errors for <90% old growth forest fraction grid cells
#-statistics reported in SI Table S6
#
#Date: 16/08/2022

library(mblm)
library(viridis)
library(rgdal)
library(Kendall)
library(trend)
library(scales)
library(raster)
library(ncdf4)
library(gridExtra)
library(rasterVis)
library(maptools)
library(tls)
library(epiR)
library(ggplot2)
library(rgdal) 
library(weights)
library(deming)
library(scales)

setwd("D:/RECCAP")


#read rasters of modelled and L-VOD derived AGC

deforestedBiomassmapbiomasstack <- stack('./CodeOutputs/deforestedbiomassstackV17.tif')
degradedBiomassstack <- stack('./CodeOutputs/degradedbiomassstackV17.tif')
edgeBiomassstack <-  stack('./CodeOutputs/edgebiomasschangestackV17.tif')
SFBiomassstack <- stack('./CodeOutputs/SFBiomassgrowthstackV17.tif')

intactBiomasschangesmoothmeanstack <- stack('./CodeOutputs/intactBiomasschangesmoothmeanstackV17.tif')
intactBiomasschangesmoothmaxstack <- stack('./CodeOutputs/intactBiomasschangesmoothmaxstackV17.tif')
intactBiomasschangetrendmeanstack <- stack('./CodeOutputs/intactBiomasschangetrendmeanstackV17.tif')

AGC_LVOD_amazoniamethodsmeancrop <- stack('./CodeOutputs/AGC_LVOD_amazoniamethodsmeancropV18.tif')
intactforestcrop <- stack('./CodeOutputs/intactforestcropV16.tif')
mapbiomasforestfraccrop <- stack('./CodeOutputs/mapbiomasforestfraccropV16.tif')

modelAGCchangestack <- stack()
modelAGCchange_incIFstack <- stack()

#create different stacks of modelled AGC values per year, including different processes covering:
#deforestation, SF growth, edge effects, non-edge degradation and old-growth (intact) forest change from adjacent LVOD grid cells

modelAGCstack <- AGC_LVOD_amazoniamethodsmeancrop[[1]]
modelAGC_incIFstack <- AGC_LVOD_amazoniamethodsmeancrop[[1]]
modelAGC_onlydeforestationstack <- AGC_LVOD_amazoniamethodsmeancrop[[1]]
modelAGC_deforestation_nonedgedegradationstack <- AGC_LVOD_amazoniamethodsmeancrop[[1]]
modelAGC_deforestation_alldegradationstack <- AGC_LVOD_amazoniamethodsmeancrop[[1]]

for(i in 1:8){
  modelAGCstack <- stack(modelAGCstack,modelAGCstack[[i]]-deforestedBiomassmapbiomasstack[[i]]-degradedBiomassstack[[i]]-edgeBiomassstack[[i]]+SFBiomassstack[[i]])
  
  modelAGCchangestack <- stack(modelAGCchangestack,(-1)*deforestedBiomassmapbiomasstack[[i]]-degradedBiomassstack[[i]]-edgeBiomassstack[[i]]+SFBiomassstack[[i]])
  
  modelAGC_incIFstack <- stack(modelAGC_incIFstack,modelAGC_incIFstack[[i]]-deforestedBiomassmapbiomasstack[[i]]-degradedBiomassstack[[i]]-edgeBiomassstack[[i]]+SFBiomassstack[[i]]+mean(intactBiomasschangesmoothmeanstack[[i]],intactBiomasschangesmoothmaxstack[[i]],intactBiomasschangetrendmeanstack[[i]]))
  modelAGC_onlydeforestationstack<- stack(modelAGC_onlydeforestationstack, modelAGC_onlydeforestationstack[[i]]-deforestedBiomassmapbiomasstack[[i]])
  modelAGC_deforestation_nonedgedegradationstack <-  stack(modelAGC_deforestation_nonedgedegradationstack,modelAGC_deforestation_nonedgedegradationstack[[i]]-deforestedBiomassmapbiomasstack[[i]]-degradedBiomassstack[[i]])
  modelAGC_deforestation_alldegradationstack <- stack(modelAGC_deforestation_alldegradationstack, modelAGC_deforestation_alldegradationstack[[i]]-deforestedBiomassmapbiomasstack[[i]]-degradedBiomassstack[[i]]-edgeBiomassstack[[i]])
  modelAGCchange_incIFstack <- stack(modelAGCchange_incIFstack,(-1)*deforestedBiomassmapbiomasstack[[i]]-degradedBiomassstack[[i]]-edgeBiomassstack[[i]]+SFBiomassstack[[i]]+mean(intactBiomasschangesmoothmeanstack[[i]],intactBiomasschangesmoothmaxstack[[i]],intactBiomasschangetrendmeanstack[[i]]))
  }

modelNetAGCchange <- modelAGC_incIFstack[[9]]-modelAGC_incIFstack[[1]]
modelNetAGCchangenoIF <- modelAGCstack[[9]]-modelAGCstack[[1]]
LVODNetAGCchange <- AGC_LVOD_amazoniamethodsmeancrop[[9]]-AGC_LVOD_amazoniamethodsmeancrop[[1]]

modelNetAGCchange_onlydeforestationstack <- modelAGC_onlydeforestationstack[[9]]-modelAGC_onlydeforestationstack[[1]]
modelNetAGCchange_deforestation_nonedgedegradationstack <- modelAGC_deforestation_nonedgedegradationstack[[9]]-modelAGC_deforestation_nonedgedegradationstack[[1]]
modelNetAGCchange_deforestation_alldegradationstack <- modelAGC_deforestation_alldegradationstack[[9]]-modelAGC_deforestation_alldegradationstack[[1]]

NetAGCchangediffLVODminusmodel <- LVODNetAGCchange-modelNetAGCchange
NetAGCchangediffmodelminusLVOD <- modelNetAGCchange-LVODNetAGCchange

#make mask of areas with <90% intact forest in 2018
notintactmask <- intactforestcrop[[8]]<0.9
notintactmask[notintactmask<1] <- NA 


#Correlation stats

#forest frac correlation

mapbiomasforestfracchange <-  mapbiomasforestfraccrop[[8]]-mapbiomasforestfraccrop[[1]]  

par(pty='s')
LVODchangenotintactvals <- getValues(LVODNetAGCchange*notintactmask)
modelchangenotintactvals <- getValues(mapbiomasforestfracchange*notintactmask)
plot(modelchangenotintactvals,LVODchangenotintactvals,pch=16,col=alpha('black',0.05),xlim=c(-0.5,0.5),ylim=c(-40,25),xlab=bquote("Forest frac. cover change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),ylab=bquote("LVOD AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))
LVODmodellinmod <- lm(LVODchangenotintactvals~modelchangenotintactvals)
#abline(0,1,lty=2)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(LVODmodellinmod$coefficients[1],LVODmodellinmod$coefficients[2])
summary(LVODmodellinmod)
pearson <- stats::cor(LVODchangenotintactvals,modelchangenotintactvals, use="complete.obs")
print(paste0('pearsons:',pearson))
MAE <- mean(abs(LVODchangenotintactvals-modelchangenotintactvals),na.rm=T)
print(paste0('MAE: ',MAE))
bias <- mean((modelchangenotintactvals-LVODchangenotintactvals),na.rm=T)
print(paste0('bias: ',bias))
relbias <- bias/mean(LVODchangenotintactvals,na.rm=T)
print(paste0('relbias: ',relbias))


######Absolute changes correlation: spatial correlation of AGC changes between beginning of 2011 and beginning of 2019

LVODvsModel_abs_mat <- matrix(data=NA,ncol=6,nrow=7)

LVODvsModel_abs_mat[1,] <- c('abs','Full model',	'degradation + deforestation+SF',	'degradation + deforestation',	'non-edge degr. + deforestation',	'deforestation')
LVODvsModel_abs_mat[,1] <- c('abs','R2',	'RSE',	'pearsons r',	'MAE',	'bias','relbias') 


#full model
par(pty='s')
LVODchangenotintactvals <- getValues(LVODNetAGCchange*notintactmask)
modelchangenotintactvals <- getValues(modelNetAGCchange*notintactmask)
plot(modelchangenotintactvals,LVODchangenotintactvals,pch=16,col=alpha('black',0.05),xlim=c(-40,25),ylim=c(-40,25),xlab=bquote("Modelled AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),ylab=bquote("LVOD AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))
LVODmodellinmod <- lm(LVODchangenotintactvals~modelchangenotintactvals)
abline(0,1,lty=2)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(LVODmodellinmod$coefficients[1],LVODmodellinmod$coefficients[2])
LVODmodellinmodsummary <-  summary(LVODmodellinmod)
r2 <- LVODmodellinmodsummary$r.squared
RSE <- LVODmodellinmodsummary$sigma
pearson <- stats::cor(LVODchangenotintactvals,modelchangenotintactvals, use="complete.obs")
print(paste0('pearsons:',pearson))
MAE <- mean(abs(LVODchangenotintactvals-modelchangenotintactvals),na.rm=T)
print(paste0('MAE: ',MAE))
bias <- mean((modelchangenotintactvals-LVODchangenotintactvals),na.rm=T)
print(paste0('bias: ',bias))
relbias <- bias/mean(LVODchangenotintactvals,na.rm=T)
print(paste0('relbias: ',relbias))

LVODvsModel_abs_mat[2:7,2] <- c(r2,RSE,pearson,MAE,bias,relbias)



#no IF
par(pty='s')
LVODchangenotintactvals <- getValues(LVODNetAGCchange*notintactmask)
modelchangenotintactvals <- getValues(modelNetAGCchangenoIF*notintactmask)
plot(modelchangenotintactvals,LVODchangenotintactvals,pch=16,col=alpha('black',0.05),xlim=c(-40,25),ylim=c(-40,25),xlab=bquote("Modelled AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),ylab=bquote("LVOD AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))
LVODmodellinmod <- lm(LVODchangenotintactvals~modelchangenotintactvals)
abline(0,1,lty=2)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(LVODmodellinmod$coefficients[1],LVODmodellinmod$coefficients[2])
LVODmodellinmodsummary <-  summary(LVODmodellinmod)
r2 <- LVODmodellinmodsummary$r.squared
RSE <- LVODmodellinmodsummary$sigma
pearson <- stats::cor(LVODchangenotintactvals,modelchangenotintactvals, use="complete.obs")
print(paste0('pearsons:',pearson))
MAE <- mean(abs(LVODchangenotintactvals-modelchangenotintactvals),na.rm=T)
print(paste0('MAE: ',MAE))
bias <- mean((modelchangenotintactvals-LVODchangenotintactvals),na.rm=T)
print(paste0('bias: ',bias))
relbias <- bias/mean(LVODchangenotintactvals,na.rm=T)
print(paste0('relbias: ',relbias))

LVODvsModel_abs_mat[2:7,3] <- c(r2,RSE,pearson,MAE,bias,relbias)

#no IF, no SF
par(pty='s')
LVODchangenotintactvals <- getValues(LVODNetAGCchange*notintactmask)
modelchangenotintactvals <- getValues(modelNetAGCchange_deforestation_alldegradationstack*notintactmask)
plot(modelchangenotintactvals,LVODchangenotintactvals,pch=16,col=alpha('black',0.05),xlim=c(-40,25),ylim=c(-40,25),xlab=bquote("Modelled AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),ylab=bquote("LVOD AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))
LVODmodellinmod <- lm(LVODchangenotintactvals~modelchangenotintactvals)
abline(0,1,lty=2)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(LVODmodellinmod$coefficients[1],LVODmodellinmod$coefficients[2])
LVODmodellinmodsummary <-  summary(LVODmodellinmod)
r2 <- LVODmodellinmodsummary$r.squared
RSE <- LVODmodellinmodsummary$sigma
pearson <- stats::cor(LVODchangenotintactvals,modelchangenotintactvals, use="complete.obs")
print(paste0('pearsons:',pearson))
MAE <- mean(abs(LVODchangenotintactvals-modelchangenotintactvals),na.rm=T)
print(paste0('MAE: ',MAE))
bias <- mean((modelchangenotintactvals-LVODchangenotintactvals),na.rm=T)
print(paste0('bias: ',bias))
relbias <- bias/mean(LVODchangenotintactvals,na.rm=T)
print(paste0('relbias: ',relbias))

LVODvsModel_abs_mat[2:7,4] <- c(r2,RSE,pearson,MAE,bias,relbias)

#no IF, no SF, no edge degr
par(pty='s')
LVODchangenotintactvals <- getValues(LVODNetAGCchange*notintactmask)
modelchangenotintactvals <- getValues(modelNetAGCchange_deforestation_nonedgedegradationstack*notintactmask)
plot(modelchangenotintactvals,LVODchangenotintactvals,pch=16,col=alpha('black',0.05),xlim=c(-40,25),ylim=c(-40,25),xlab=bquote("Modelled AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),ylab=bquote("LVOD AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))
LVODmodellinmod <- lm(LVODchangenotintactvals~modelchangenotintactvals)
abline(0,1,lty=2)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(LVODmodellinmod$coefficients[1],LVODmodellinmod$coefficients[2])
LVODmodellinmodsummary <-  summary(LVODmodellinmod)
r2 <- LVODmodellinmodsummary$r.squared
RSE <- LVODmodellinmodsummary$sigma
pearson <- stats::cor(LVODchangenotintactvals,modelchangenotintactvals, use="complete.obs")
print(paste0('pearsons:',pearson))
MAE <- mean(abs(LVODchangenotintactvals-modelchangenotintactvals),na.rm=T)
print(paste0('MAE: ',MAE))
bias <- mean((modelchangenotintactvals-LVODchangenotintactvals),na.rm=T)
print(paste0('bias: ',bias))
relbias <- bias/mean(LVODchangenotintactvals,na.rm=T)
print(paste0('relbias: ',relbias))

LVODvsModel_abs_mat[2:7,5] <- c(r2,RSE,pearson,MAE,bias,relbias)

#no IF, no SF, no degr
par(pty='s')
LVODchangenotintactvals <- getValues(LVODNetAGCchange*notintactmask)
modelchangenotintactvals <- getValues(modelNetAGCchange_onlydeforestationstack*notintactmask)
plot(modelchangenotintactvals,LVODchangenotintactvals,pch=16,col=alpha('black',0.05),xlim=c(-40,25),ylim=c(-40,25),xlab=bquote("Modelled AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),ylab=bquote("LVOD AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))
LVODmodellinmod <- lm(LVODchangenotintactvals~modelchangenotintactvals)
abline(0,1,lty=2)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(LVODmodellinmod$coefficients[1],LVODmodellinmod$coefficients[2])
LVODmodellinmodsummary <-  summary(LVODmodellinmod)
r2 <- LVODmodellinmodsummary$r.squared
RSE <- LVODmodellinmodsummary$sigma
pearson <- stats::cor(LVODchangenotintactvals,modelchangenotintactvals, use="complete.obs")
print(paste0('pearsons:',pearson))
MAE <- mean(abs(LVODchangenotintactvals-modelchangenotintactvals),na.rm=T)
print(paste0('MAE: ',MAE))
bias <- mean((modelchangenotintactvals-LVODchangenotintactvals),na.rm=T)
print(paste0('bias: ',bias))
relbias <- bias/mean(LVODchangenotintactvals,na.rm=T)
print(paste0('relbias: ',relbias))

LVODvsModel_abs_mat[2:7,6] <- c(r2,RSE,pearson,MAE,bias,relbias)

write.table(LVODvsModel_abs_mat,'D:/RECCAP/stats/LVODvsmodel_absolute_stats_V17_corr.csv',sep=',',row.names = F,col.names = F)


######Correlation stats trends

fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates)}}
fun2=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$p.value) }}



LVODtrendAGCchange=calc(AGC_LVOD_amazoniamethodsmeancrop[[1:9]], fun1)
#AGC_trend_intact <- raster::mask(AGC_trend,intactmask)
LVODtrendAGCchangemean <- cellStats(LVODtrendAGCchange,mean,na.rm=T)
LVODtrendAGCchangenotintactmean <- cellStats(LVODtrendAGCchange*notintactmask,mean,na.rm=T)

LVODtrendAGCchange_sig=calc(AGC_LVOD_amazoniamethodsmeancrop[[1:9]], fun2)
#AGC_trend_intact_sig <- raster::mask(AGC_trend_sig,intactmask)

LVODtrendAGCchange_sigmask <- LVODtrendAGCchange_sig
LVODtrendAGCchange_sigmask[LVODtrendAGCchange_sig<0.05] <- 1
LVODtrendAGCchange_sigmask[LVODtrendAGCchange_sig>=0.05] <- 0
LVODtrendAGCchange_sigmask[LVODtrendAGCchange_sigmask==0] <- NA

modeltrendAGCchange <- calc(modelAGC_incIFstack,fun1)
modeltrendAGCchangenoIF <- calc(modelAGCstack,fun1)

modeltrendAGCchangemean <- cellStats(modeltrendAGCchange,mean,na.rm=T)
modeltrendAGCchangenotintactmean <- cellStats(modeltrendAGCchange*notintactmask,mean,na.rm=T)

modeltrendAGCchange_onlydeforestationstack <- calc(modelAGC_onlydeforestationstack,fun1)
modeltrendAGCchange_deforestation_nonedgedegradationstack <- calc(modelAGC_deforestation_nonedgedegradationstack,fun1)
modeltrendAGCchange_deforestation_alldegradationstack <- calc(modelAGC_deforestation_alldegradationstack,fun1)

LVODvsModel_trend_mat <- matrix(data=NA,ncol=6,nrow=7)

LVODvsModel_trend_mat[1,] <- c('trend','Full model',	'degradation + deforestation+SF',	'degradation + deforestation',	'non-edge degr. + deforestation',	'deforestation')
LVODvsModel_trend_mat[,1] <- c('trend','R2',	'RSE',	'pearsons r',	'MAE',	'bias','relbias') 


#full model
par(pty='s')
LVODchangenotintactvals <- getValues(LVODtrendAGCchange*notintactmask)
modelchangenotintactvals <- getValues(modeltrendAGCchange*notintactmask)
plot(modelchangenotintactvals,LVODchangenotintactvals,pch=16,col=alpha('black',0.05),xlim=c(-6,5),ylim=c(-6,5),xlab=bquote("Modelled AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),ylab=bquote("LVOD AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))
LVODmodellinmod <- lm(LVODchangenotintactvals~modelchangenotintactvals)
abline(0,1,lty=2)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(LVODmodellinmod$coefficients[1],LVODmodellinmod$coefficients[2])
LVODmodellinmodsummary <-  summary(LVODmodellinmod)
r2 <- LVODmodellinmodsummary$r.squared
RSE <- LVODmodellinmodsummary$sigma
pearson <- stats::cor(LVODchangenotintactvals,modelchangenotintactvals, use="complete.obs")
print(paste0('pearsons:',pearson))
MAE <- mean(abs(LVODchangenotintactvals-modelchangenotintactvals),na.rm=T)
print(paste0('MAE: ',MAE))
bias <- mean((modelchangenotintactvals-LVODchangenotintactvals),na.rm=T)
print(paste0('bias: ',bias))
relbias <- bias/mean(LVODchangenotintactvals,na.rm=T)
print(paste0('relbias: ',relbias))

LVODvsModel_trend_mat[2:7,2] <- c(r2,RSE,pearson,MAE,bias,relbias)

#no IF
par(pty='s')
LVODchangenotintactvals <- getValues(LVODtrendAGCchange*notintactmask)
modelchangenotintactvals <- getValues(modeltrendAGCchangenoIF*notintactmask)
plot(modelchangenotintactvals,LVODchangenotintactvals,pch=16,col=alpha('black',0.05),xlim=c(-10,5),ylim=c(-10,5),xlab=bquote("Modelled AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),ylab=bquote("LVOD AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))
LVODmodellinmod <- lm(LVODchangenotintactvals~modelchangenotintactvals)
abline(0,1,lty=2)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(LVODmodellinmod$coefficients[1],LVODmodellinmod$coefficients[2])
LVODmodellinmodsummary <-  summary(LVODmodellinmod)
r2 <- LVODmodellinmodsummary$r.squared
RSE <- LVODmodellinmodsummary$sigma
pearson <- stats::cor(LVODchangenotintactvals,modelchangenotintactvals, use="complete.obs")
print(paste0('pearsons:',pearson))
MAE <- mean(abs(LVODchangenotintactvals-modelchangenotintactvals),na.rm=T)
print(paste0('MAE: ',MAE))
bias <- mean((modelchangenotintactvals-LVODchangenotintactvals),na.rm=T)
print(paste0('bias: ',bias))
relbias <- bias/mean(LVODchangenotintactvals,na.rm=T)
print(paste0('relbias: ',relbias))

LVODvsModel_trend_mat[2:7,3] <- c(r2,RSE,pearson,MAE,bias,relbias)

#no IF, no SF
par(pty='s')
LVODchangenotintactvals <- getValues(LVODtrendAGCchange*notintactmask)
modelchangenotintactvals <- getValues(modeltrendAGCchange_deforestation_alldegradationstack*notintactmask)
plot(modelchangenotintactvals,LVODchangenotintactvals,pch=16,col=alpha('black',0.05),xlim=c(-40,25),ylim=c(-40,25),xlab=bquote("Modelled AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),ylab=bquote("LVOD AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))
LVODmodellinmod <- lm(LVODchangenotintactvals~modelchangenotintactvals)
abline(0,1,lty=2)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(LVODmodellinmod$coefficients[1],LVODmodellinmod$coefficients[2])
LVODmodellinmodsummary <-  summary(LVODmodellinmod)
r2 <- LVODmodellinmodsummary$r.squared
RSE <- LVODmodellinmodsummary$sigma
pearson <- stats::cor(LVODchangenotintactvals,modelchangenotintactvals, use="complete.obs")
print(paste0('pearsons:',pearson))
MAE <- mean(abs(LVODchangenotintactvals-modelchangenotintactvals),na.rm=T)
print(paste0('MAE: ',MAE))
bias <- mean((modelchangenotintactvals-LVODchangenotintactvals),na.rm=T)
print(paste0('bias: ',bias))
relbias <- bias/mean(LVODchangenotintactvals,na.rm=T)
print(paste0('relbias: ',relbias))

LVODvsModel_trend_mat[2:7,4] <- c(r2,RSE,pearson,MAE,bias,relbias)


#no IF, no SF, no edge degr
par(pty='s')
LVODchangenotintactvals <- getValues(LVODtrendAGCchange*notintactmask)
modelchangenotintactvals <- getValues(modeltrendAGCchange_deforestation_nonedgedegradationstack*notintactmask)
plot(modelchangenotintactvals,LVODchangenotintactvals,pch=16,col=alpha('black',0.05),xlim=c(-40,25),ylim=c(-40,25),xlab=bquote("Modelled AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),ylab=bquote("LVOD AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))
LVODmodellinmod <- lm(LVODchangenotintactvals~modelchangenotintactvals)
abline(0,1,lty=2)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(LVODmodellinmod$coefficients[1],LVODmodellinmod$coefficients[2])
LVODmodellinmodsummary <-  summary(LVODmodellinmod)
r2 <- LVODmodellinmodsummary$r.squared
RSE <- LVODmodellinmodsummary$sigma
pearson <- stats::cor(LVODchangenotintactvals,modelchangenotintactvals, use="complete.obs")
print(paste0('pearsons:',pearson))
MAE <- mean(abs(LVODchangenotintactvals-modelchangenotintactvals),na.rm=T)
print(paste0('MAE: ',MAE))
bias <- mean((modelchangenotintactvals-LVODchangenotintactvals),na.rm=T)
print(paste0('bias: ',bias))
relbias <- bias/mean(LVODchangenotintactvals,na.rm=T)
print(paste0('relbias: ',relbias))

LVODvsModel_trend_mat[2:7,5] <- c(r2,RSE,pearson,MAE,bias,relbias)


#no IF, no SF, no degr
par(pty='s')
LVODchangenotintactvals <- getValues(LVODtrendAGCchange*notintactmask)
modelchangenotintactvals <- getValues(modeltrendAGCchange_onlydeforestationstack*notintactmask)
plot(modelchangenotintactvals,LVODchangenotintactvals,pch=16,col=alpha('black',0.05),xlim=c(-6,5),ylim=c(-6,5),xlab=bquote("Modelled AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),ylab=bquote("LVOD AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))
LVODmodellinmod <- lm(LVODchangenotintactvals~modelchangenotintactvals)
abline(0,1,lty=2) 
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(LVODmodellinmod$coefficients[1],LVODmodellinmod$coefficients[2])
LVODmodellinmodsummary <-  summary(LVODmodellinmod)
r2 <- LVODmodellinmodsummary$r.squared
RSE <- LVODmodellinmodsummary$sigma
pearson <- stats::cor(LVODchangenotintactvals,modelchangenotintactvals, use="complete.obs")
print(paste0('pearsons:',pearson))
MAE <- mean(abs(LVODchangenotintactvals-modelchangenotintactvals),na.rm=T)
print(paste0('MAE: ',MAE))
bias <- mean((modelchangenotintactvals-LVODchangenotintactvals),na.rm=T)
print(paste0('bias: ',bias))
relbias <- bias/mean(LVODchangenotintactvals,na.rm=T)
print(paste0('relbias: ',relbias))

LVODvsModel_trend_mat[2:7,6] <- c(r2,RSE,pearson,MAE,bias,relbias)

write.table(LVODvsModel_trend_mat,'D:/RECCAP/stats/LVODvsmodel_trend_stats_V17_corr.csv',sep=',',row.names = F,col.names = F)

