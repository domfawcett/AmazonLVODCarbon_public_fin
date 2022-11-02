#Author: D. Fawcett
#
#Description:
#This code generates Amazon country specific changes in AGC attributed to different processes using outputs from analyseLVOD_ccgfilt_diff_Amazon_Main_V18_fin.R
#Excluding 2011 which was a la Nina year and saw a large increase in AGC

#Outputs:
#-statistics of gross and net biomass changes and trends for each country and the entire Amazon (SI Table S4)
#Supplementary Information Figures S10
#
#Date: 16/08/2022

library(mblm)
library(viridis)
library(rgdal)
library(Kendall)
library(trend)
library(scales)
library(ncdf4)
library(raster)
library(gridExtra)
library(rasterVis)
library(maptools)
library(raster)
library(tls)
library(epiR)
library(ggplot2)
library(rgdal)
  

setwd("D:/RECCAP")

#load previously generated datasets

totalMask <- raster('./CodeOutputs/totalMaskV18.tif')

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

SFBiomassstack <- stack('./CodeOutputs/SFBiomassgrowthstackv17.tif')


intactBiomasschangemethodsmeanstack <- stack('./CodeOutputs/intactBiomasschangemethodsmeanstackV18.tif')
intactBiomasschangemethodsmeanPlusSDstack <- stack('./CodeOutputs/intactBiomasschangemethodsmeanPlusSDstackV18.tif')
intactBiomasschangemethodsmeanMinusSDstack <- stack('./CodeOutputs/intactBiomasschangemethodsmeanMinusSDstackV18.tif')

intactBiomasschangesmoothmeanstack <- stack('./CodeOutputs/intactBiomasschangesmoothmeanstackv17.tif')
intactBiomasschangesmoothmaxstack <- stack('./CodeOutputs/intactBiomasschangesmoothmaxstackv17.tif')
intactBiomasschangetrendmeanstack <- stack('./CodeOutputs/intactBiomasschangetrendmeanstackv17.tif')

AGC_LVOD_amazoniamethodsmeancrop <- stack('./CodeOutputs/AGC_LVOD_amazoniamethodsmeancropV18.tif')
AGC_LVOD_amazoniasmoothmeancrop <- stack('./CodeOutputs/AGC_LVOD_amazoniasmoothmeancropV18.tif')
AGC_LVOD_amazoniasmoothmaxcrop <- stack('./CodeOutputs/AGC_LVOD_amazoniasmoothmaxcropV18.tif')
AGC_LVOD_amazoniatrendmeancrop <- stack('./CodeOutputs/AGC_LVOD_amazoniatrendmeancropV18.tif')

AGC_LVOD_amazoniamethodsmeancrop <- stack('./CodeOutputs/AGC_LVOD_amazoniamethodsmeancropV18.tif')
AGC_LVOD_amazoniamethodsmeancropPlusSD <- stack('./CodeOutputs/AGC_LVOD_amazoniamethodsmeancropPlusSDV18.tif')


intactforestcrop <- stack('./CodeOutputs/intactforestcropV16.tif')
mapbiomasforestfraccrop <- stack('./CodeOutputs/mapbiomasforestfraccropV16.tif')


#lists of generated plots
changeprocessplotlist <- list()
changeprocessplotlistnoIF <- list()
changeprocessplotlistwithIF <- list()
changeprocessbarplotlistwithIF <- list()
changeprocessbarplotannuallistwithIF <- list()

changetotalplotlist <- list()
totalplotlist <- list()

#get amazon country shapes
amazonCountries <- readOGR(dsn = './Datasets/South_America', layer = "South_America")
selectCountries <- c(1,2,4,6,8,9,11,14)
selectCountrynames <- c('Bolivia','Brazil','Colombia','Ecuador','French Guiana (France)','Guyana','Suriname','Peru','Venezuela','Amazon')
 
#prepare output tables
countrystatstable <- matrix(data=NA,nrow=10,ncol=22)
countrystatstrendtable <- matrix(data=NA,nrow=10,ncol=18)

colnames <- c('name',
              'ncells',
                           'LVODfirstAGC',
                           'LVODfirstAGCSD',
              'modelnetDeforestation',
              'modelnetDeforestationSD',
              'modelnetDegradedEdge',
              'modelnetDegradedEdgeSD',
              'modelnetDegradedNonEdge',
              'modelnetDegradedNonEdgeSD',
              'modelnetSF',
              'modelnetSFSD',
              'modelnetIF',
              'modelnetIFSD',
              'modelnetAGCchange',
              'modelnetAGCchangeSD',
              'modelrelnetAGCchange',
              'modelrelnetAGCchangeSD',
                           'LVODnetAGCchange',
                            'LVODnetAGCchangeSD',
              'LVODrelnetAGCchange',
              'LVODrelnetAGCchangeSD')


colnamestrends <- c('name',
              'ncells',
              'LVODfirstAGC',
              'LVODfirstAGCSD',
'LVODtrend',
'LVODtrendCI1',
'LVODtrendCI2',
'LVODreltrend',
'LVODreltrendCI1',
'LVODreltrendCI2',
'LVODtrendsig',
'modeltrend',
'modeltrendCI1',
'modeltrendCI2',
'modelreltrend',
'modelreltrendCI1',
'modelreltrendCI2',
'modeltrendsig'
)

countrystatsdf <- data.frame(countrystatstable)
names(countrystatsdf) <- colnames
countrystatsdf$name <- selectCountrynames

countrystatstrenddf <- data.frame(countrystatstrendtable)
names(countrystatstrenddf) <- colnamestrends
countrystatstrenddf$name <- selectCountrynames

#iterate through countries and total Amazon
for(i in 1:10){

if(i<10){
 countryshape <- amazonCountries[amazonCountries@data$COUNTRY == selectCountrynames[i],]  
 countryshapecrop <-  crop(countryshape,amazonBasin)
}else{
  countryshape <- amazonBasin
  countryshapecrop <- amazonBasin
}

  
subsetmask <- raster::mask(totalMask,countryshapecrop)



#time series of AGC change (2012 to 2019)


deforestedBiomassmapbiomasTS <- (cellStats(raster::mask(deforestedBiomassmapbiomasstack,subsetmask)*area(deforestedBiomassmapbiomasstack)*100,sum,na.rm=T)/10^9)[2:8]
deforestedBiomassmapbiomasPlusSDTS <- (cellStats(raster::mask(deforestedBiomassmapbiomasPlusSDstack,subsetmask)*area(deforestedBiomassmapbiomasstack)*100,sum,na.rm=T)/10^9)[2:8]
deforestedBiomassmapbiomasMinusSDTS <- (cellStats(raster::mask(deforestedBiomassmapbiomasMinusSDstack,subsetmask)*area(deforestedBiomassmapbiomasstack)*100,sum,na.rm=T)/10^9)[2:8]

degradedBiomassTS <- (cellStats(raster::mask(degradedBiomassstack,subsetmask)*area(degradedBiomassstack)*100,sum,na.rm=T)/10^9)[2:8]
degradedBiomassPlusSDTS <- (cellStats(raster::mask(degradedBiomassPlusSDstack,subsetmask)*area(degradedBiomassPlusSDstack)*100,sum,na.rm=T)/10^9)[2:8]
degradedBiomassMinusSDTS <- (cellStats(raster::mask(degradedBiomassMinusSDstack,subsetmask)*area(degradedBiomassMinusSDstack)*100,sum,na.rm=T)/10^9)[2:8]

edgeBiomassTS <- (cellStats(raster::mask(edgeBiomassstack,subsetmask)*area(edgeBiomassstack)*100,sum,na.rm=T)/10^9)[2:8]
edgeBiomassPlusSDTS <- (cellStats(raster::mask(edgeBiomassPlusSDstack,subsetmask)*area(edgeBiomassPlusSDstack)*100,sum,na.rm=T)/10^9)[2:8]
edgeBiomassMinusSDTS <- (cellStats(raster::mask(edgeBiomassMinusSDstack,subsetmask)*area(edgeBiomassMinusSDstack)*100,sum,na.rm=T)/10^9)[2:8]

SFBiomassTS <- (cellStats(raster::mask(SFBiomassstack,subsetmask)*area(SFBiomassstack)*100,sum,na.rm=T)/10^9)[2:8]

intactBiomasschangesmoothmeanTS <- (cellStats(raster::mask(intactBiomasschangesmoothmeanstack,subsetmask)*area(intactBiomasschangesmoothmeanstack)*100,sum,na.rm=T)/10^9)[2:8]
intactBiomasschangesmoothmaxTS <- (cellStats(raster::mask(intactBiomasschangesmoothmaxstack,subsetmask)*area(intactBiomasschangesmoothmaxstack)*100,sum,na.rm=T)/10^9)[2:8]
intactBiomasschangetrendmeanTS <- (cellStats(raster::mask(intactBiomasschangetrendmeanstack,subsetmask)*area(intactBiomasschangetrendmeanstack)*100,sum,na.rm=T)/10^9)[2:8]

intactBiomasschangemethodsmeanTS <- (cellStats(raster::mask(intactBiomasschangemethodsmeanstack,subsetmask)*area(intactBiomasschangemethodsmeanstack[[1]])*100,sum,na.rm=T)/10^9)[2:8]
intactBiomasschangemethodsmeanPlusSDTS <- (cellStats(raster::mask(intactBiomasschangemethodsmeanPlusSDstack,subsetmask)*area(intactBiomasschangemethodsmeanPlusSDstack[[1]])*100,sum,na.rm=T)/10^9)[2:8]
intactBiomasschangemethodsmeanMinusSDTS <- (cellStats(raster::mask(intactBiomasschangemethodsmeanMinusSDstack ,subsetmask)*area(intactBiomasschangemethodsmeanMinusSDstack[[1]])*100,sum,na.rm=T)/10^9)[2:8]

#combine the uncertainties related to biomass map (CCISD) and LVOD index used (SD)

#for intact change
intactBiomasschangemethodSDs <- apply(cbind(intactBiomasschangetrendmeanTS,intactBiomasschangesmoothmeanTS,intactBiomasschangesmoothmaxTS),1,sd) 
intactBiomasschangemethodCCImapSDs <- abs(intactBiomasschangemethodsmeanPlusSDTS-intactBiomasschangemethodsmeanTS) 
intactBiomasschangeuncertainty <- sqrt(intactBiomasschangemethodCCImapSDs^2+intactBiomasschangemethodSDs^2)

#for total LVOD AGC change
LVODAGCTSmethodsmean <- (cellStats(raster::mask(AGC_LVOD_amazoniamethodsmeancrop,subsetmask)*area(AGC_LVOD_amazoniamethodsmeancrop)*100,sum,na.rm=T)/10^9)[2:9]
LVODAGCTSmethodsmeanPlusSD <- (cellStats(raster::mask(AGC_LVOD_amazoniamethodsmeancropPlusSD,subsetmask)*area(AGC_LVOD_amazoniamethodsmeancropPlusSD)*100,sum,na.rm=T)/10^9)[2:9]

LVODAGCTSsmoothmean <- (cellStats(raster::mask(AGC_LVOD_amazoniasmoothmeancrop,subsetmask)*area(AGC_LVOD_amazoniasmoothmeancrop)*100,sum,na.rm=T)/10^9)[2:9]
#LVODAGCchangeTSsmoothmean <- LVODAGCTSsmoothmean[2:10]-LVODAGCTSsmoothmean[1:9]

LVODAGCTSsmoothmax <- (cellStats(raster::mask(AGC_LVOD_amazoniasmoothmaxcrop,subsetmask)*area(AGC_LVOD_amazoniasmoothmaxcrop)*100,sum,na.rm=T)/10^9)[2:9]
#LVODAGCchangeTSsmoothmax <- LVODAGCTSsmoothmax[2:10]-LVODAGCTSsmoothmax[1:9]

LVODAGCTStrendmean <- (cellStats(raster::mask(AGC_LVOD_amazoniatrendmeancrop,subsetmask)*area(AGC_LVOD_amazoniatrendmeancrop)*100,sum,na.rm=T)/10^9)[2:9]
#LVODAGCchangeTStrendmean <- LVODAGCTStrendmean[2:10]-LVODAGCTStrendmean[1:9]


#for total LVOD AGC
LVODAGCTSCCISDs <- abs(LVODAGCTSmethodsmeanPlusSD-LVODAGCTSmethodsmean)

methodmeans <- apply(cbind(LVODAGCTSsmoothmean,LVODAGCTSsmoothmax,LVODAGCTStrendmean),1,mean) 
methodSDs <- apply(cbind(LVODAGCTSsmoothmean,LVODAGCTSsmoothmax,LVODAGCTStrendmean),1,sd) 
methodtotSDs <- sqrt(methodSDs^2+LVODAGCTSCCISDs^2)


#secondary forest standard deviation, using 95% CI reported for average growth rate in Heinrich et al. 2021, assuming normal distribution

SFrelSD = 1/3/1.96


#vary axis depending on country
if(i>3 & i<8){
  plotymin <- (-0.03)
  plotymax <- (0.03)
  plotymin2 <- (-0.5)
  plotymax2 <- (0.7)
  plotymin3 <- (-0.01)
  plotymax3 <- (0.005)
  plotymin4 <- (-30)
  plotymax4 <- (30)
  plotymin5 <- (-400)
  plotymax5 <- (160)  
  
}else if(i!=2){
  plotymin <- (-0.08)
  plotymax <- (0.075)
  plotymin2 <- (-0.1)
  plotymax2 <- (0.1)
  plotymin3 <- (-0.02)
  plotymax3 <- (0.01)
  plotymin4 <- (-130)
  plotymax4 <- (90)
  plotymin5 <- (-40)
  plotymax5 <- (40)
}else{
  plotymin <- (-0.3)
  plotymax <- (0.5)
  plotymin2 <- (-0.5)
  plotymax2 <- (0.7)
  plotymin3 <- (-0.3)
  plotymax3 <- (0.05)
  plotymin4 <- (-2000)
  plotymax4 <- (500)
  plotymin5 <- (-400)
  plotymax5 <- (160)
}


AGCTSdf <- data.frame(years=seq(2012.5,2018.5,1),deforestedMapbiomas=-deforestedBiomassmapbiomasTS,
                      deforestedMapbiomasMinusSD=-deforestedBiomassmapbiomasMinusSDTS,
                      deforestedMapbiomasPlusSD=-deforestedBiomassmapbiomasPlusSDTS,
                      degradedBiomass=-degradedBiomassTS,
                      degradedBiomassMinusSD=-degradedBiomassMinusSDTS,
                      degradedBiomassPlusSD=-degradedBiomassPlusSDTS,
                      degradedEdges=-edgeBiomassTS,
                      degradedEdgesMinusSD=-edgeBiomassMinusSDTS,
                      degradedEdgesPlusSD=-edgeBiomassPlusSDTS,
                      SFgrowth=SFBiomassTS,
                      intactChangeBiomass=intactBiomasschangemethodsmeanTS,
                      intactChangeBiomassSD=intactBiomasschangeuncertainty)

AGCTSSDdf <- data.frame(years=seq(2012.5,2018.5,1),
                        deforestedMapbiomasSD=abs(AGCTSdf$deforestedMapbiomas-AGCTSdf$deforestedMapbiomasPlusSD),
                        degradedBiomassSD=abs(AGCTSdf$degradedBiomass-AGCTSdf$degradedBiomassPlusSD),
                        degradedEdgesSD=abs(AGCTSdf$degradedEdges-AGCTSdf$degradedEdgesPlusSD),
                        SFgrowthSD=abs(AGCTSdf$SFgrowth*SFrelSD),
                        intactChangeBiomassSD=intactBiomasschangeuncertainty
                        )

#add stats to data frame

#functions to estimate trends and significance
fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates)}}
fun2=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$p.value) }}

#calculate trend
LVODAGCchangetrend <- fun1(methodmeans)
LVODAGCchangetrendsig <- fun2(methodmeans)

#confidence interval of trend
LVODAGCchangeCI1 <- sens.slope(methodmeans)$conf.int[1]
LVODAGCchangeCI2 <- sens.slope(methodmeans)$conf.int[2]


#net LVOD AGC changes and uncertainties: 2012 to 2019
LVODnetAGCchange <- methodmeans[8]-methodmeans[1]

LVODnetAGCTSmethodsmeanPlusSDchange <- LVODAGCTSmethodsmeanPlusSD[8]-LVODAGCTSmethodsmeanPlusSD[1]
LVODnetAGCTSmethodsmeanchange <- LVODAGCTSmethodsmean[8]-LVODAGCTSmethodsmean[1]

LVODnetAGCchangeCCISD <- abs(LVODnetAGCTSmethodsmeanPlusSDchange-LVODnetAGCTSmethodsmeanchange)

LVODnetAGCchangeTSsmoothmean <- LVODAGCTSsmoothmean[8]-LVODAGCTSsmoothmean[1]
LVODnetAGCchangeTSsmoothmax <- LVODAGCTSsmoothmax[8]-LVODAGCTSsmoothmax[1]
LVODnetAGCchangeTStrendmean <- LVODAGCTStrendmean[8]-LVODAGCTStrendmean[1]

LVODnetAGCchangeSD <- sd(c(LVODnetAGCchangeTSsmoothmean,LVODnetAGCchangeTSsmoothmax,LVODnetAGCchangeTStrendmean))

LVODnetAGCchangetotSD <- sqrt(LVODnetAGCchangeCCISD^2+LVODnetAGCchangeSD^2)#sqrt(methodSDs[1]^2+methodSDs[9]^2)


#calculate total gross changes of processes
modeldeforestationtotchange <- sum(AGCTSdf$deforestedMapbiomas)
modeldegradationtotchange <- sum(AGCTSdf$degradedBiomass)
modeledgestotchange <- sum(AGCTSdf$degradedEdges)
modelSFtotchange <- sum(AGCTSdf$SFgrowth)

#calculate total net change for intact forest
modelintacttotchange <- sum(AGCTSdf$intactChangeBiomass)

#get combined uncertainties of gross change
modeldeforestationtotSD <- sqrt(sum(AGCTSSDdf$deforestedMapbiomasSD^2))
modeldegradationtotSD <- sqrt(sum(AGCTSSDdf$degradedBiomassSD^2))
modeledgestotSD <- sqrt(sum(AGCTSSDdf$degradedEdgesSD^2))
modelSFtotSD <- sqrt(sum(AGCTSSDdf$SFgrowthSD^2))
modelintacttotSD <- sqrt(sum(AGCTSSDdf$intactChangeBiomassSD^2))

#get total modelled change by combining processes
modelAGCannualtotchange <- AGCTSdf$deforestedMapbiomas+AGCTSdf$degradedBiomass+AGCTSdf$degradedEdges+AGCTSdf$SFgrowth+AGCTSdf$intactChangeBiomass

modelAGCtotchange <- modeldeforestationtotchange+modeldegradationtotchange+modeledgestotchange+modelSFtotchange+modelintacttotchange

modelAGCtotSD <- sqrt(modeldeforestationtotSD^2+modeldegradationtotSD^2+modeledgestotSD^2+modelSFtotSD^2+modelintacttotSD^2)

modelAGCannualSD <- mean(sqrt(AGCTSSDdf$deforestedMapbiomasSD^2+AGCTSSDdf$degradedBiomassSD^2+AGCTSSDdf$degradedEdgesSD^2+AGCTSSDdf$SFgrowthSD^2+AGCTSSDdf$intactChangeBiomassSD^2),na.rm=T)

modelAGCtot <- methodmeans
for(j in 1:7){
modelAGCtot[j+1] <- modelAGCtot[j]+modelAGCannualtotchange[j]  
}

#trends and CI of modelled change
modelAGCchangetrend <- fun1(modelAGCtot)
modelAGCchangetrendsig <- fun2(modelAGCtot)

modelAGCchangeCI1 <- sens.slope(modelAGCtot)$conf.int[1]
modelAGCchangeCI2 <- sens.slope(modelAGCtot)$conf.int[2]

#total AGC change dataframe for table
countrystatsdf$LVODfirstAGC[i] <- methodmeans[1]
countrystatsdf$LVODfirstAGCSD[i] <- methodtotSDs[1]
countrystatsdf$LVODnetAGCchange[i] <- LVODnetAGCchange
countrystatsdf$LVODnetAGCchangeSD[i] <- LVODnetAGCchangetotSD
countrystatsdf$LVODrelnetAGCchange[i] <- LVODnetAGCchange/methodmeans[1]*100
countrystatsdf$LVODrelnetAGCchangeSD[i] <- LVODnetAGCchangetotSD/methodmeans[1]*100
countrystatsdf$modelnetAGCchange[i] <- modelAGCtotchange
countrystatsdf$modelnetAGCchangeSD[i] <- modelAGCtotSD
countrystatsdf$modelrelnetAGCchange[i] <- modelAGCtotchange/methodmeans[1]*100
countrystatsdf$modelrelnetAGCchangeSD[i] <- modelAGCtotSD/methodmeans[1]*100
countrystatsdf$modelnetDeforestation[i] <- modeldeforestationtotchange
countrystatsdf$modelnetDeforestationSD[i] <- modeldeforestationtotSD
countrystatsdf$modelnetDegradedEdge[i] <- modeledgestotchange
countrystatsdf$modelnetDegradedEdgeSD[i] <- modeledgestotSD
countrystatsdf$modelnetDegradedNonEdge[i] <- modeldegradationtotchange
countrystatsdf$modelnetDegradedNonEdgeSD[i] <- modeldegradationtotSD
countrystatsdf$modelnetSF[i] <- modelSFtotchange
countrystatsdf$modelnetSFSD[i] <- modelSFtotSD
countrystatsdf$modelnetIF[i] <- modelintacttotchange
countrystatsdf$modelnetIFSD[i] <- modelintacttotSD
countrystatsdf$ncells[i] <- sum(getValues(subsetmask),na.rm=T)

#AGC trends dataframe for table

countrystatstrenddf$ncells[i] <- sum(getValues(subsetmask),na.rm=T)
countrystatstrenddf$LVODfirstAGC[i] <- methodmeans[1]
countrystatstrenddf$LVODfirstAGCSD[i] <- methodtotSDs[1]
countrystatstrenddf$LVODtrend[i] <- LVODAGCchangetrend
countrystatstrenddf$LVODtrendCI1[i] <- LVODAGCchangeCI1
countrystatstrenddf$LVODtrendCI2[i] <- LVODAGCchangeCI2
countrystatstrenddf$LVODreltrend[i] <- LVODAGCchangetrend/methodmeans[1]*100
countrystatstrenddf$LVODreltrendCI1[i] <- LVODAGCchangeCI1/methodmeans[1]*100
countrystatstrenddf$LVODreltrendCI2[i] <- LVODAGCchangeCI2/methodmeans[1]*100
countrystatstrenddf$LVODtrendsig[i] <- LVODAGCchangetrendsig
countrystatstrenddf$modeltrend[i] <- modelAGCchangetrend
countrystatstrenddf$modeltrendCI1[i] <- modelAGCchangeCI1
countrystatstrenddf$modeltrendCI2[i] <- modelAGCchangeCI2
countrystatstrenddf$modelreltrend[i] <- modelAGCchangetrend/methodmeans[1]*100
countrystatstrenddf$modelreltrendCI1[i] <- modelAGCchangeCI1/methodmeans[1]*100
countrystatstrenddf$modelreltrendCI2[i] <- modelAGCchangeCI2/methodmeans[1]*100
countrystatstrenddf$modeltrendsig[i] <- modelAGCchangetrendsig


#bar plot total changes

barname <- c('total','deforestation','degradation','SF growth','IF change','L-VOD')
colvec <- c('grey30','red','brown','green','darkgreen','blue')

modelchangedf <- data.frame(id=c(5,1,2,3,4,6),name=barname,cols=colvec,
                            vals=1000*c(modelAGCtotchange,modeldeforestationtotchange,modeledgestotchange+modeldegradationtotchange,modelSFtotchange,modelintacttotchange,LVODnetAGCchange),
                            errors=1000*c(modelAGCtotSD,modeldeforestationtotSD,sqrt(modeledgestotSD^2+modeldegradationtotSD^2),modelSFtotSD,modelintacttotSD,LVODnetAGCchangetotSD))
          


modelchangedfordered  <- modelchangedf[order(modelchangedf$id,decreasing=T),]

barP<-ggplot(data=modelchangedf,aes(x=id, y=vals,fill=name))+
  geom_bar(stat="identity", position=position_dodge())+
  ylab(bquote("AGC change " ~ "["~ Tg ~ C  ~ "]")) +
  xlab("") +
  scale_fill_manual("legend", values = c("deforestation" = "#EE6677", "degradation" = "#aa7633", "SF growth" = '#66CCEE',"IF change"='#4477AA',"total"="#ababab","L-VOD"='#009988'))+
  geom_errorbar(aes(ymin=vals-errors, ymax=vals+errors), width=.2,
                position=position_dodge(.9))+#+ 
  ggtitle(selectCountrynames[i])+
  ylim(plotymin4,plotymax4)+
  geom_hline(yintercept=0,lty=2)+
  theme_bw()+
theme(legend.position='none',axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank())

  

changeprocessbarplotlistwithIF[[i]] <- barP



p2new <- ggplot(AGCTSdf)+#, aes(x = years, y = changemethodmeans)) + 
  #ggtitle(paste(classnames[i],', n=',size,sep=''))+
  theme_classic() +
  ggtitle(selectCountrynames[i])+
  ylab(bquote("modelled AGC change " ~ "["~ Pg ~ C  ~yr^{-1}~ "]")) +
  xlab('year') +
  # ylim(-0.15,0.15)+
  scale_x_continuous(breaks=seq(2012,2019,1),limits=c(2012, 2019))+
  theme(text = element_text(size=18),axis.text.x=element_text(angle=45,hjust=1))+
  ylim(plotymin3,plotymax3)+
  
  geom_hline(yintercept=0,lty=2)+
 
  geom_line(aes(years, deforestedMapbiomas),color='#EE6677')+
  geom_line(aes(years, SFgrowth),color='#66CCEE')+
  
  
  geom_ribbon(aes(years, SFgrowth,ymin = SFgrowth*(1-SFrelSD),
                  ymax = SFgrowth*(1+SFrelSD)), alpha = 0.2,fill='#66CCEE')+
  geom_ribbon(aes(years, deforestedMapbiomas,ymin = deforestedMapbiomasMinusSD,
                  ymax =deforestedMapbiomasPlusSD), alpha = 0.2,fill='#EE6677')+
  geom_line(aes(years, degradedBiomass),color='#AA3377')+
  geom_ribbon(aes(years, degradedBiomass,ymin = degradedBiomassMinusSD,
                  ymax =degradedBiomassPlusSD), alpha = 0.2,fill='#AA3377')+
  geom_line(aes(years, degradedEdges),color='#CCBB44')+
  geom_ribbon(aes(years, degradedEdges,ymin = degradedEdgesMinusSD,
                  ymax =degradedEdgesPlusSD), alpha = 0.2,fill='#CCBB44')

changeprocessplotlistnoIF[[i]] <- p2new

p2newinclIF <- ggplot(AGCTSdf)+
  theme_classic() +
  ggtitle(selectCountrynames[i])+
  ylab(bquote("modelled AGC change " ~ "["~ Pg ~ C  ~yr^{-1}~ "]")) +
  xlab('year') +
  # ylim(-0.15,0.15)+
  scale_x_continuous(breaks=seq(2012,2019,1),limits=c(2012, 2019))+
  theme(text = element_text(size=18),axis.text.x=element_text(angle=45,hjust=1))+
  ylim(plotymin,plotymax)+
  
  geom_hline(yintercept=0,lty=2)+

  
  geom_line(aes(years, intactChangeBiomass),color='#4477AA')+
  geom_ribbon(aes(years, intactChangeBiomass,ymin = intactChangeBiomass - intactChangeBiomassSD,
                  ymax = intactChangeBiomass + intactChangeBiomassSD), alpha = 0.2,fill='#4477AA')+
  geom_line(aes(years, deforestedMapbiomas),color='#EE6677')+
  geom_ribbon(aes(years, SFgrowth,ymin = SFgrowth*(1-SFrelSD),
                  ymax = SFgrowth*(1+SFrelSD)), alpha = 0.2,fill='#66CCEE')+
  geom_line(aes(years, SFgrowth),color='#66CCEE')+
  geom_ribbon(aes(years, deforestedMapbiomas,ymin = deforestedMapbiomasMinusSD,
                  ymax =deforestedMapbiomasPlusSD), alpha = 0.2,fill='#EE6677')+
  geom_line(aes(years, degradedBiomass),color='#AA3377')+
  geom_ribbon(aes(years, degradedBiomass,ymin = degradedBiomassMinusSD,
                  ymax =degradedBiomassPlusSD), alpha = 0.2,fill='#AA3377')+
  geom_line(aes(years, degradedEdges),color='#CCBB44')+
  geom_ribbon(aes(years, degradedEdges,ymin = degradedEdgesMinusSD,
                  ymax =degradedEdgesPlusSD), alpha = 0.2,fill='#CCBB44')
   
changeprocessplotlistwithIF[[i]] <- p2newinclIF

}


#convert to teragram C
countrystatsdfTgC <-countrystatsdf
countrystatsdfTgC[,c(-1,-2,-17,-18,-21,-22)] <- countrystatsdfTgC[,c(-1,-2,-17,-18,-21,-22)]*1000 

countrystatstrenddfTgC <-countrystatstrenddf
countrystatstrenddfTgC[,c(-1,-2,-8,-9,-10,-11,-15,-16,-17,-18)] <- countrystatstrenddfTgC[,c(-1,-2,-8,-9,-10,-11,-15,-16,-17,-18)]*1000 

#make Figure 4

pdf(file = "./Figures/FinFigs/Fig4_excl2011_V18.pdf",   
    width = 10, 
    height = 9) 

grid.arrange(grobs=changeprocessbarplotlistwithIF[1:9],nrow=3,ncol=3)

dev.off()

#export tables of data to make Extended Data Tables 2 and 3
write.table(countrystatsdfTgC,'D:/RECCAP/stats/countryAGCchangestatsuncertaintiesJRC_excl2011_V18.csv',sep=',',row.names=F)
write.table(countrystatstrenddfTgC,'D:/RECCAP/stats/countryAGCchangestatsTrenduncertaintiesJRC_excl2011_V18.csv',sep=',',row.names=F)

#make Figures S5 and S6

pdf(file = "./Figures/FinFigs/Supplementary/FigS5_excl2011_V18.pdf",   
    width = 16, 
    height = 15) 
grid.arrange(grobs=changeprocessplotlistnoIF[1:9],nrow=3,ncol=3)
dev.off()

pdf(file = "./Figures/FinFigs/Supplementary/FigS6_excl2011_V18.pdf",   # The directory you want to save the file in
    width = 16, # The width of the plot in inches
    height = 15) # The height of the plot in inches
grid.arrange(grobs=changeprocessplotlistwithIF[1:9],nrow=3,ncol=3)
dev.off()

