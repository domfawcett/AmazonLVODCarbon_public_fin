#Author: D. Fawcett
#
#Description:
#This code generates the figures 2 and 3 in the paper

#Outputs:
#-Figure 2 a-e and 3 a and b
#
#Date: 25/10/2022


library(raster)
library(trend)
library(scales)
library(gridExtra)
library(rasterVis)
library(maptools)
library(tls)
library(ggplot2)
library(reshape2)
library(rgdal) 
library(RColorBrewer)

#adjust directory

setwd("D:/RECCAP")

amazonBasin <- readOGR(dsn = './Datasets', layer = "amazonBasinCelso")


##############
#Visualise and analyse outputs

amazonBasin <- readOGR(dsn = './Datasets', layer = "amazonBasinCelso")

#deforested biomass
deforestedBiomassmapbiomasstack <- stack('./CodeOutputs/deforestedbiomassstackv17.tif')
deforestedBiomassmapbiomasPlusSDstack <- stack('./CodeOutputs/deforestedbiomassPlusSDstackv17.tif')
deforestedBiomassmapbiomasMinusSDstack <- stack('./CodeOutputs/deforestedbiomassMinusSDstackv17.tif')

#degraded Biomass
degradedBiomassstack <- stack('./CodeOutputs/degradedbiomassstackv17.tif')
degradedBiomassPlusSDstack <- stack('./CodeOutputs/degradedbiomassPlusSDstackv17.tif')
degradedBiomassMinusSDstack <- stack('./CodeOutputs/degradedbiomassMinusSDstackv17.tif')


#edge Biomass change
edgeBiomassstack <- stack('./CodeOutputs/edgebiomasschangestackv17.tif')
edgeBiomassPlusSDstack <- stack('./CodeOutputs/edgebiomassPlusSDchangestackv17.tif')
edgeBiomassMinusSDstack <- stack('./CodeOutputs/edgebiomassMinusSDchangestackv17.tif')

#secondary forest growth change
SFBiomassstack <- stack('./CodeOutputs/SFBiomassgrowthstackv17.tif')

#intact forest changes
intactBiomasschangemethodsmeanstack <- stack('./CodeOutputs/intactBiomasschangemethodsmeanstackV18.tif')
intactBiomasschangemethodsmeanPlusSDstack <- stack('./CodeOutputs/intactBiomasschangemethodsmeanPlusSDstackV18.tif')
intactBiomasschangemethodsmeanMinusSDstack <- stack('./CodeOutputs/intactBiomasschangemethodsmeanMinusSDstackV18.tif')

intactBiomasschangesmoothmeanstack <- stack('./CodeOutputs/intactBiomasschangesmoothmeanstackv17.tif')
intactBiomasschangesmoothmaxstack <- stack('./CodeOutputs/intactBiomasschangesmoothmaxstackv17.tif')
intactBiomasschangetrendmeanstack <- stack('./CodeOutputs/intactBiomasschangetrendmeanstackv17.tif')

#AGC LVOD changes
AGC_LVOD_amazoniasmoothmeancrop <- stack('./CodeOutputs/AGC_LVOD_amazoniasmoothmeancropV18.tif')
AGC_LVOD_amazoniasmoothmaxcrop <- stack('./CodeOutputs/AGC_LVOD_amazoniasmoothmaxcropV18.tif')
AGC_LVOD_amazoniatrendmeancrop <- stack('./CodeOutputs/AGC_LVOD_amazoniatrendmeancropV18.tif')

AGC_LVOD_amazoniamethodsmeancrop <- stack('./CodeOutputs/AGC_LVOD_amazoniamethodsmeancropV18.tif')
AGC_LVOD_amazoniamethodsmeancropPlusSD <- stack('./CodeOutputs/AGC_LVOD_amazoniamethodsmeancropPlusSDV18.tif')

intactforestcrop <- stack('./CodeOutputs/intactforestcropV16.tif')

totalMask <- raster('./CodeOutputs/totalMaskV18.tif')


#time series of summed values

subsetmask <- totalMask


intactBiomasschangemethodsmeanTS <- cellStats(raster::mask(intactBiomasschangemethodsmeanstack,subsetmask)*area(intactBiomasschangemethodsmeanstack[[1]])*100,sum,na.rm=T)/10^9
intactBiomasschangemethodsmeanPlusSDTS <- cellStats(raster::mask(intactBiomasschangemethodsmeanPlusSDstack,subsetmask)*area(intactBiomasschangemethodsmeanPlusSDstack[[1]])*100,sum,na.rm=T)/10^9
intactBiomasschangemethodsmeanMinusSDTS <- cellStats(raster::mask(intactBiomasschangemethodsmeanMinusSDstack ,subsetmask)*area(intactBiomasschangemethodsmeanMinusSDstack[[1]])*100,sum,na.rm=T)/10^9

intactBiomasschangetrendmeanTS <- cellStats(raster::mask(intactBiomasschangetrendmeanstack,subsetmask)*area(intactBiomasschangetrendmeanstack[[1]])*100,sum,na.rm=T)/10^9
intactBiomasschangesmoothmaxTS <- cellStats(raster::mask(intactBiomasschangesmoothmaxstack,subsetmask)*area(intactBiomasschangesmoothmaxstack[[1]])*100,sum,na.rm=T)/10^9
intactBiomasschangesmoothmeanTS <- cellStats(raster::mask(intactBiomasschangesmoothmeanstack,subsetmask)*area(intactBiomasschangesmoothmeanstack[[1]])*100,sum,na.rm=T)/10^9


SFBiomassTS <- cellStats(raster::mask(SFBiomassstack,subsetmask)*area(SFBiomassstack[[1]])*100,sum,na.rm=T)/10^9

deforestedBiomassmapbiomasTS <- cellStats(raster::mask(deforestedBiomassmapbiomasstack,subsetmask)*area(deforestedBiomassmapbiomasstack[[1]])*100,sum,na.rm=T)/10^9
deforestedBiomassmapbiomasPlusSDTS <- cellStats(raster::mask(deforestedBiomassmapbiomasPlusSDstack,subsetmask)*area(deforestedBiomassmapbiomasPlusSDstack[[1]])*100,sum,na.rm=T)/10^9
deforestedBiomassmapbiomasMinusSDTS <- cellStats(raster::mask(deforestedBiomassmapbiomasMinusSDstack,subsetmask)*area(deforestedBiomassmapbiomasMinusSDstack[[1]])*100,sum,na.rm=T)/10^9

degradedBiomassTS <- cellStats(raster::mask(degradedBiomassstack,subsetmask)*area(degradedBiomassstack[[1]])*100,sum,na.rm=T)/10^9
degradedBiomassPlusSDTS<- cellStats(raster::mask(degradedBiomassPlusSDstack,subsetmask)*area(degradedBiomassPlusSDstack[[1]])*100,sum,na.rm=T)/10^9
degradedBiomassMinusSDTS <- cellStats(raster::mask(degradedBiomassMinusSDstack,subsetmask)*area(degradedBiomassMinusSDstack[[1]])*100,sum,na.rm=T)/10^9

edgeBiomassTS <- cellStats(raster::mask(edgeBiomassstack,subsetmask)*area(edgeBiomassstack[[1]])*100,sum,na.rm=T)/10^9
edgeBiomassPlusSDTS <- cellStats(raster::mask(edgeBiomassPlusSDstack,subsetmask)*area(edgeBiomassPlusSDstack[[1]])*100,sum,na.rm=T)/10^9
edgeBiomassMinusSDTS <- cellStats(raster::mask(edgeBiomassMinusSDstack,subsetmask)*area(edgeBiomassMinusSDstack[[1]])*100,sum,na.rm=T)/10^9


#get max annual mean deforestation
deforestedBiomassmapbiomasmean <- calc(deforestedBiomassmapbiomasstack,mean)

maxdeforestedBiomassmapbiomasmean <- cellStats(deforestedBiomassmapbiomasmean*area(deforestedBiomassmapbiomasmean)*100,max)

maxdeforestedTgCyr <- maxdeforestedBiomassmapbiomasmean/10^6

#secondary forest standard deviation, using 95% CI reported for average growth rate in Heinrich et al. 2021, assuming normal distribution

SFrelSD = 1/3/1.96

#combine uncertainty (SD) from biomass CCI map and from use of different LVOD indices
intactBiomasschangemethodmeans <- intactBiomasschangemethodsmeanTS
intactBiomasschangemethodSDs <- apply(cbind(intactBiomasschangetrendmeanTS,intactBiomasschangesmoothmeanTS,intactBiomasschangesmoothmaxTS),1,sd) 
intactBiomasschangemethodCCImapSDs <- abs(intactBiomasschangemethodsmeanPlusSDTS-intactBiomasschangemethodsmeanTS) 
intactBiomasschangeuncertainty <- sqrt(intactBiomasschangemethodCCImapSDs^2+intactBiomasschangemethodSDs^2)


AGCTSdf <- data.frame(years=seq(2011.5,2018.5,1),deforestedMapbiomas=-deforestedBiomassmapbiomasTS,
                      deforestedMapbiomasMinusSD=-deforestedBiomassmapbiomasMinusSDTS,
                      deforestedMapbiomasPlusSD=-deforestedBiomassmapbiomasPlusSDTS,
                      degradedBiomass=-degradedBiomassTS,
                      degradedBiomassMinusSD=-degradedBiomassMinusSDTS,
                      degradedBiomassPlusSD=-degradedBiomassPlusSDTS,
                      degradedEdges=-edgeBiomassTS,
                      degradedEdgesMinusSD=-edgeBiomassMinusSDTS,
                      degradedEdgesPlusSD=-edgeBiomassPlusSDTS,
                      SFgrowth=SFBiomassTS,
                      intactChangeBiomass=intactBiomasschangemethodmeans,
                      intactChangeBiomassSD=intactBiomasschangeuncertainty)#,LVODAGC=LVODAGCTS)


#make Figure 3 a

pdf(file = "./Figures/FinFigs/Fig3a_V18.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 5) # The height of the plot in inches

p2new <- ggplot(AGCTSdf)+ 
  theme_classic() +
  ylab(bquote("modelled AGC change " ~ "["~ Pg ~ C  ~yr^{-1}~ "]")) +
  xlab('year') +
  scale_x_continuous(breaks=seq(2011,2019,1),limits=c(2011, 2019))+
  theme(text = element_text(size=18),axis.text.x=element_text(angle=45,hjust=1))+
  geom_hline(yintercept=0,lty=2)+
  geom_line(aes(years, deforestedMapbiomas),color='#EE6677')+
  geom_line(aes(years, SFgrowth),color='#66CCEE')+
  

  geom_ribbon(aes(years, SFgrowth,ymin = SFgrowth*(1-SFrelSD),
                  ymax = SFgrowth*(1+SFrelSD)), alpha = 0.2,fill='#66CCEE')+
  geom_line(aes(years, intactChangeBiomass),color='#4477AA')+
  geom_ribbon(aes(years, intactChangeBiomass,ymin = intactChangeBiomass - intactChangeBiomassSD,
                  ymax = intactChangeBiomass + intactChangeBiomassSD), alpha = 0.2,fill='#4477AA')+
  geom_ribbon(aes(years, deforestedMapbiomas,ymin = deforestedMapbiomasMinusSD,
                 ymax =deforestedMapbiomasPlusSD), alpha = 0.2,fill='#EE6677')+
  geom_line(aes(years, degradedBiomass),color='#AA3377')+
  geom_ribbon(aes(years, degradedBiomass,ymin = degradedBiomassMinusSD,
                  ymax =degradedBiomassPlusSD), alpha = 0.2,fill='#AA3377')+
  geom_line(aes(years, degradedEdges),color='#CCBB44')+
  geom_ribbon(aes(years, degradedEdges,ymin = degradedEdgesMinusSD,
                  ymax =degradedEdgesPlusSD), alpha = 0.2,fill='#CCBB44')

plot(p2new)

dev.off()

#seperate figure legend
plot(AGCTSdf$years,AGCTSdf$deforestedMapbiomas,type='n')
legend(2012,-0.14, legend=c("deforestation", "degradation (edge)","degradation (non-edge)","secondary forest growth","old-growth forest change"),
       col=c("#EE6677", "#CCBB44","#AA3377","#66CCEE","#4477AA"), lty=1,lwd=1.5,box.lty=0)




#secondary forest trend

SFgrowth <- lm(formula = AGCTSdf$SFgrowth ~ c(1:8))
summary(SFgrowth)

#get L-VOD AGC time series and uncertainty
LVODAGCTSsmoothmean <- cellStats(raster::mask(AGC_LVOD_amazoniasmoothmeancrop,subsetmask)*area(AGC_LVOD_amazoniasmoothmeancrop)*100,sum,na.rm=T)/10^9

LVODAGCTSsmoothmax <- cellStats(raster::mask(AGC_LVOD_amazoniasmoothmaxcrop,subsetmask)*area(AGC_LVOD_amazoniasmoothmaxcrop)*100,sum,na.rm=T)/10^9

LVODAGCTStrendmean <- cellStats(raster::mask(AGC_LVOD_amazoniatrendmeancrop,subsetmask)*area(AGC_LVOD_amazoniatrendmeancrop)*100,sum,na.rm=T)/10^9

LVODAGCTSmethodsmean <- cellStats(AGC_LVOD_amazoniamethodsmeancrop*area(AGC_LVOD_amazoniamethodsmeancrop)*100,sum,na.rm=T)/10^9
LVODAGCTSmethodsmeanPlusSD <- cellStats(AGC_LVOD_amazoniamethodsmeancropPlusSD*area(AGC_LVOD_amazoniamethodsmeancropPlusSD)*100,sum,na.rm=T)/10^9


#calculate means and uncertainty of LVOD and modelled AGC
methodmeans <- LVODAGCTSmethodsmean

methodSDs <- apply(cbind(LVODAGCTSsmoothmean,LVODAGCTSsmoothmax,LVODAGCTStrendmean),1,sd) 

methodCCISDs <- abs(LVODAGCTSmethodsmeanPlusSD-LVODAGCTSmethodsmean ) 

methodtotSDs <- sqrt(methodSDs^2+methodCCISDs^2)


changemethodmeans <- methodmeans[2:9]-methodmeans[1:8]


LVODAGCTSsmoothmeanchange <- LVODAGCTSsmoothmean[2:9]-LVODAGCTSsmoothmean [1:8]
LVODAGCTSsmoothmaxchange <- LVODAGCTSsmoothmax[2:9]-LVODAGCTSsmoothmax[1:8]
LVODAGCTStrendmeanchange <- LVODAGCTStrendmean[2:9]-LVODAGCTStrendmean[1:8]

changemethodSDs <- apply(cbind(LVODAGCTSsmoothmeanchange,LVODAGCTSsmoothmaxchange,LVODAGCTStrendmeanchange),1,sd) #sqrt(methodSDs[2:9]^2+methodSDs[1:8]^2) 

LVODAGCTSmethodsmeanPlusSDchange <- LVODAGCTSmethodsmeanPlusSD[2:9]-LVODAGCTSmethodsmeanPlusSD[1:8]
LVODAGCTSmethodsmeanchange <- LVODAGCTSmethodsmean[2:9]-LVODAGCTSmethodsmean[1:8]

changemethodCCISDs <- abs(LVODAGCTSmethodsmeanPlusSDchange-LVODAGCTSmethodsmeanchange)

changemethodtotSDs <- sqrt(changemethodCCISDs^2+changemethodSDs^2)

#calculate total AGC stocks per year from modelled changes and 2011 L-VOD AGC as baseline


modelAGCchange <- vector(length=8)
modelAGCchange_combSD <- vector(length=8)

modelAGC_inclIF <- vector(length=9)
modelAGC_inclIF[1] <- methodmeans[1]

modelAGC_inclIFmin <- vector(length=9)
modelAGC_inclIFmin[1] <- methodmeans[1]

modelAGC_inclIFmax <- vector(length=9)
modelAGC_inclIFmax[1] <- methodmeans[1]

modelAGC_inclIFcumulativeErrors <- vector(length=9)
modelAGC_inclIFcumulativeErrors[1] <- 0

modelAGCchange_inclIF <- vector(length=8)
modelAGCchange_inclIF_combSD <- vector(length=8)

modelAGCstack <- AGC_LVOD_amazoniamethodsmeancrop[[1]]
modelAGC_incIFstack <- AGC_LVOD_amazoniamethodsmeancrop[[1]]
modelAGCchangestack <- stack()
modelAGCchange_incIFstack <- stack()


for(i in 1:8){
 modelAGC_inclIF[i+1] <- modelAGC_inclIF[i]-(deforestedBiomassmapbiomasTS[i])-(degradedBiomassTS[i])-(edgeBiomassTS[i])+SFBiomassTS[i]+intactBiomasschangemethodmeans[i]
 deforestedSD <- deforestedBiomassmapbiomasPlusSDTS[i]-deforestedBiomassmapbiomasTS[i]
 degradedSD <- degradedBiomassPlusSDTS[i]-degradedBiomassTS[i]
 edgeSD <- edgeBiomassPlusSDTS[i]-edgeBiomassTS[i]
 SFSD <- SFBiomassTS[i]*SFrelSD
 intactSD <- intactBiomasschangeuncertainty[i]
 
 #cumulative errors up to current year (i+1)
 modelAGC_inclIFcumulativeErrors[i+1] <- sqrt(modelAGC_inclIFcumulativeErrors[i]^2+(sqrt(deforestedSD^2+degradedSD^2+edgeSD^2+SFSD^2+intactSD^2))^2)
 modelAGCchange_inclIF_combSD[i] <- sqrt(deforestedSD^2+degradedSD^2+edgeSD^2+SFSD^2+intactSD^2)
 
 
 modelAGCchange[i] <- ((-1)*deforestedBiomassmapbiomasTS[i])-(degradedBiomassTS[i])-(edgeBiomassTS[i])+SFBiomassTS[i]#+intactBiomasschangeTS[i]
 modelAGCchange_inclIF[i] <- (-1)*(deforestedBiomassmapbiomasTS[i])-(degradedBiomassTS[i])-(edgeBiomassTS[i])+SFBiomassTS[i]+intactBiomasschangemethodmeans[i]
 
 modelAGCchangestack <- stack(modelAGCchangestack,(-1)*deforestedBiomassmapbiomasstack[[i]]-degradedBiomassstack[[i]]-edgeBiomassstack[[i]]+SFBiomassstack[[i]])
 
modelAGC_incIFstack <- stack(modelAGC_incIFstack,modelAGC_incIFstack[[i]]-deforestedBiomassmapbiomasstack[[i]]-degradedBiomassstack[[i]]-edgeBiomassstack[[i]]+SFBiomassstack[[i]]+mean(intactBiomasschangesmoothmeanstack[[i]],intactBiomasschangesmoothmaxstack[[i]],intactBiomasschangetrendmeanstack[[i]]))
modelAGCchange_incIFstack <- stack(modelAGCchange_incIFstack,(-1)*deforestedBiomassmapbiomasstack[[i]]-degradedBiomassstack[[i]]-edgeBiomassstack[[i]]+SFBiomassstack[[i]]+mean(intactBiomasschangesmoothmeanstack[[i]],intactBiomasschangesmoothmaxstack[[i]],intactBiomasschangetrendmeanstack[[i]]))

}

#get minimum and maximum for error bar ribbon display
modelAGC_inclIFmin <- modelAGC_inclIF-modelAGC_inclIFcumulativeErrors
modelAGC_inclIFmax <- modelAGC_inclIF+modelAGC_inclIFcumulativeErrors


modelAGCchangeTScomb <- modelAGCchange_inclIF#(-1)*deforestedBiomassmapbiomasTS-degradedBiomassTS-edgeBiomassTS+SFBiomassTS+intactBiomasschangemethodmeans
modelAGCchangeTScombmax <-  modelAGCchange_inclIF+modelAGCchange_inclIF_combSD#(-1)*deforestedBiomassmapbiomasPlusSDTS-degradedBiomassPlusSDTS-edgeBiomassPlusSDTS+SFBiomassTS*0.66666+(intactBiomasschangemethodmeans-intactBiomasschangemethodSDs)
modelAGCchangeTScombmin <- modelAGCchange_inclIF-modelAGCchange_inclIF_combSD#(-1)*deforestedBiomassmapbiomasMinusSDTS-degradedBiomassMinusSDTS-edgeBiomassMinusSDTS+SFBiomassTS*1.33333+(intactBiomasschangemethodmeans+intactBiomasschangemethodSDs)

changemethodcompdf <- data.frame(years=seq(2011.5,2018.5,1),changemethodmeans=changemethodmeans[1:8],changemethodtotSDs=changemethodtotSDs[1:8],modelAGCchangeTScomb,modelAGCchangeTScombmin,modelAGCchangeTScombmax)

#make Figure S4

pdf(file = "./Figures/FinFigs/Supplementary/AnnualAGCchangeModelvsLVOD_V18.pdf",   
    width = 6, 
    height = 5) 


p <- ggplot(changemethodcompdf, aes(x = years, y = changemethodmeans)) + 
  theme_classic() +
  ylab(bquote("AGC change " ~ "["~ Pg ~ C  ~yr^{-1}~ "]")) +
  xlab('year') +
  scale_x_continuous(breaks=seq(2011,2019,1),limits=c(2011, 2019))+
  theme(text = element_text(size=18),axis.text.x=element_text(angle=45,hjust=1))+
  
  geom_hline(yintercept=mean(changemethodmeans),lty=2,col=  '#009988')+
  geom_hline(yintercept=mean(modelAGCchangeTScomb),lty=2,col= '#EE7733')+
  geom_hline(yintercept=0,lty=2)+
  geom_line(aes(years, changemethodmeans),color= '#009988')+
  geom_line(aes(years, modelAGCchangeTScomb),color= '#EE7733')+
  geom_ribbon(aes(ymin = changemethodmeans - changemethodtotSDs,
                  ymax = changemethodmeans + changemethodtotSDs), alpha = 0.2,fill= '#009988')+
geom_ribbon(aes(ymin = modelAGCchangeTScombmin,
                ymax = modelAGCchangeTScombmax), alpha = 0.2,fill= '#EE7733')

plot(p)

dev.off()

#make figure 3 b

pdf(file = "./Figures/FinFigs/Fig3b_V18_nolvoderror.pdf",   
    width = 6,
    height = 5) 

methodcompdf <- data.frame(years=seq(2011,2019,1),methodmeans=methodmeans[1:9],methodSDs=methodSDs[1:9],modelAGC_inclIF)

p <- ggplot(methodcompdf, aes(x = years, y = methodmeans)) + #vegetation crown plot
  #ggtitle(paste(classnames[i],', n=',size,sep=''))+
  theme_classic() +
  ylab(bquote("AGC " ~ "["~ Pg ~ C ~"]")) +
  xlab('year') +
  scale_x_continuous(breaks=seq(2011,2019,1),limits=c(2011, 2019))+
  theme(text = element_text(size=18),axis.text.x=element_text(angle=45,hjust=1))+
  #ylim(0.9,1.3)+
  geom_line(col='#009988') +
  geom_line(aes(years, modelAGC_inclIF),color= '#EE7733',lty=2)+

  geom_ribbon(aes(ymin = modelAGC_inclIFmin,
                ymax = modelAGC_inclIFmax), alpha = 0.2,fill= '#EE7733')
plot(p)

dev.off()

#make seperate legend
plot(methodcompdf$years,methodcompdf$methodmeans)
legend(2015,55.5, legend=c("L-VOD AGC", "modelled AGC"),
       col=c('#009988', '#EE7733'), lty=c(1,2), cex=0.8,box.lty=0)
 

#calculate AGC change from L-VOD AGC between 2011, 2012 and 2018
maxchange=methodmeans[2]-methodmeans[9]
maxchangeSD=sqrt(methodSDs[2]^2+methodSDs[9]^2)

maxincrease=methodmeans[2]-methodmeans[1]
maxincreaseSD=sqrt(methodSDs[2]^2+methodSDs[1]^2)

#calculate AGC change maps from model and L-VOD as well as differences

modelNetAGCchange <- modelAGC_incIFstack[[9]]-modelAGC_incIFstack[[1]]
LVODNetAGCchange <- AGC_LVOD_amazoniamethodsmeancrop[[9]]-AGC_LVOD_amazoniamethodsmeancrop[[1]]


NetAGCchangediffLVODminusmodel <- LVODNetAGCchange-modelNetAGCchange
NetAGCchangediffmodelminusLVOD <- modelNetAGCchange-LVODNetAGCchange

modelNetAGCchangecutoff <- (modelNetAGCchange<(-1))+(modelNetAGCchange>1) #cutoff value 1 
LVODNetAGCchangecutoff <- (LVODNetAGCchange<(-1))+(LVODNetAGCchange>1) #cutoff value 1 

modelNetAGCchangecutoff[modelNetAGCchangecutoff<1] <- NA
LVODNetAGCchangecutoff[LVODNetAGCchangecutoff<1] <- NA

#remove too small AGC changes
relNetAGCchangediffLVOD <- (NetAGCchangediffLVODminusmodel/abs(modelNetAGCchange))*modelNetAGCchangecutoff*100 #percentage LVOD deviates from model
relNetAGCchangediffmodel <- (NetAGCchangediffmodelminusLVOD/abs(LVODNetAGCchange))*LVODNetAGCchangecutoff*100 #percentage model deviates from LVOD


notintactmask <- intactforestcrop[[8]]<=0.9
notintactmask[notintactmask<1] <- NA

intactmask <- intactforestcrop[[8]]>0.9
intactmask[intactmask<1] <- NA


par(pty='s')
LVODchangenotintactvals <- getValues(LVODNetAGCchange*notintactmask)
modelchangenotintactvals <- getValues(modelNetAGCchange*notintactmask)
plot(modelchangenotintactvals,LVODchangenotintactvals,pch=16,col=alpha('black',0.05),xlim=c(-40,25),ylim=c(-40,25),xlab=bquote("Modelled AGC change 2010-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),ylab=bquote("LVOD AGC change 2010-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))
LVODmodellinmod <- lm(LVODchangenotintactvals~modelchangenotintactvals)
abline(0,1,lty=2)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(LVODmodellinmod$coefficients[1],LVODmodellinmod$coefficients[2])

modelNetAGCchangelog <- log(abs(modelNetAGCchange),10)

source("./Code/divergePlotFuns.R")

modelNetAGCchangeminmax <- getMinMax(modelNetAGCchange)

#build color scale with steps to visualise small changes
breakvals <- seq(0,modelNetAGCchangeminmax[2])
breaklogvals <- breakvals^2#log10(breakvals)
breaklogsteps <- breaklogvals
breaklogsteps <- breaklogsteps[-1]
breaklogsteps <- breaklogsteps*(modelNetAGCchangeminmax[2]/breaklogvals[42])
breaklogstepsfin <- c(rev(breaklogsteps[-1]*(-1)),breaklogsteps)

labelvals <- seq(0,40,10)
labelvalssqrd <- labelvals^2
labellogvals <- labelvalssqrd*(modelNetAGCchangeminmax[2]/labelvalssqrd[5])

myTheme <- rasterTheme(region = brewer.pal(n=11,'RdBu'),panel.background=list(col="lightgrey"))

writeRaster(modelNetAGCchange,"./CodeOutputs/modelNetAGCchangeTotal.tif")

pdf(file = "./Figures/FinFigs/Fig2aV17RdBu.pdf",   
    width = 6,
    height = 5) 

#plot log spaced color bar for combined model net change
levelplot(modelNetAGCchange,main=bquote("Modelled AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),at=breaklogstepsfin,colorkey = list(at = breaklogstepsfin),margin=FALSE,maxpixels = 2e10,par.settings=myTheme)+latticeExtra::layer(sp.polygons(amazonBasin,col='darkgrey'))#list(panel.background=list(col="lightgrey")))

dev.off()

#plot maps of biomass change per process

#Plot deforestation biomass change

deforestedBiomasstotchange <- calc(deforestedBiomassmapbiomasstack,sum)

writeRaster(deforestedBiomasstotchange*(-1),"./CodeOutputs/deforestedAGCchangeTotal.tif")

pdf(file = "./Figures/FinFigs/Fig2bV17RdBu.pdf",   
    width = 6, 
    height = 5) 
levelplot(-deforestedBiomasstotchange ,main=bquote("Modelled AGC change deforestation 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),at=breaklogstepsfin,colorkey = list(at = breaklogstepsfin),margin=FALSE,maxpixels = 2e10,par.settings=myTheme)+latticeExtra::layer(sp.polygons(amazonBasin,col='darkgrey'))#list(panel.background=list(col="lightgrey")))
dev.off()
#get maximum deforested
maxchange=cellStats(deforestedBiomasstotchange,max)
maxchangeperyear=maxchange/8

#Plot degraded biomass change


degradedBiomasstotchange <- calc(degradedBiomassstack,sum)+calc(edgeBiomassstack,sum)

writeRaster(degradedBiomasstotchange*(-1),"./CodeOutputs/degradedAGCchangeTotal.tif")

pdf(file = "./Figures/FinFigs/Fig2cV17RdBu.pdf",   
    width = 6, 
    height = 5) 
levelplot(-degradedBiomasstotchange,main=bquote("Modelled AGC change degradation 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),at=breaklogstepsfin,colorkey = list(at = breaklogstepsfin),margin=FALSE,maxpixels = 2e10,par.settings=myTheme)+latticeExtra::layer(sp.polygons(amazonBasin,col='darkgrey'))#list(panel.background=list(col="lightgrey")))
dev.off()

#Plot SF biomass change

SFBiomasstotchange <- calc(SFBiomassstack,sum)

writeRaster(SFBiomasstotchange,"./CodeOutputs/SFAGCchangeTotal.tif")

pdf(file = "./Figures/FinFigs/Fig2dV17RdBu.pdf",  
    width = 6, 
    height = 5) 
levelplot(SFBiomasstotchange,main=bquote("Modelled AGC change SF growth 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),at=breaklogstepsfin,colorkey = list(at = breaklogstepsfin),margin=FALSE,maxpixels = 2e10,par.settings=myTheme)+latticeExtra::layer(sp.polygons(amazonBasin,col='darkgrey'))#list(panel.background=list(col="lightgrey")))
dev.off()
#Plot intact biomass change


intactBiomasstotchange <- (calc(intactBiomasschangesmoothmeanstack,sum)+calc(intactBiomasschangesmoothmaxstack,sum)+calc(intactBiomasschangetrendmeanstack,sum))/3

writeRaster(intactBiomasstotchange,"./CodeOutputs/intactBiomasstotchangeTotal.tif")

pdf(file = "./Figures/FinFigs/Fig2eV17RdBu.pdf",  
    width = 6, 
    height = 5) 
levelplot(intactBiomasstotchange,main=bquote("Modelled AGC change intact 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),at=breaklogstepsfin,colorkey = list(at = breaklogstepsfin),margin=FALSE,maxpixels = 2e10,par.settings=myTheme)+latticeExtra::layer(sp.polygons(amazonBasin,col='darkgrey'))#list(panel.background=list(col="lightgrey")))
dev.off()

