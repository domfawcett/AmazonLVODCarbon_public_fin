#Author: D. Fawcett
#
#Description:
#This code 1. generates modelled values of AGC stocks and change and 2. compares them to annual L-VOD AGC

#Outputs:
#-Figure 2 a-e and 3 a and b, Supplementary Information Figure S5, S8, S15, S22
#
#Date: 19/10/2022


library(mblm)
library(viridis)
library(rgdal)
library(raster)
library(Kendall)
library(trend)
library(scales)
library(ncdf4)
library(gridExtra)
library(rasterVis)
library(maptools)
library(tls)
library(epiR)
library(ggplot2)
library(reshape2)
library(rgdal) 
library(RColorBrewer)


setwd("D:/RECCAP")

amazonBasin <- readOGR(dsn = './Datasets', layer = "amazonBasinCelso")

########
#read L-VOD AGC calculated using different L-VOD indices (smoothed data mean, smoothed data max and trend data mean), as well as +- SD of the ESA CCI biomass map

#LVOD based AGC ccgfilt smooth mean
AGC_LVOD_amazoniasmoothmean <- stack("./Outputs/AGC_LVOD_amazoniaBB_LVODccgfiltannual4wetmonthsmoothmean_EXTV3_Hui_inclDec2010_2ndpoly_Jan_Apr.tif")
AGC_LVOD_amazoniasmoothmeanPlusSD <- stack("./Outputs/AGC_LVOD_amazoniaBB_LVODccgfiltannual4wetmonthsmoothmean_EXTV3_Hui_inclDec2010_2ndpoly_Jan_Apr_plusSD.tif")
AGC_LVOD_amazoniasmoothmeanMinusSD <- stack("./Outputs/AGC_LVOD_amazoniaBB_LVODccgfiltannual4wetmonthsmoothmean_EXTV3_Hui_inclDec2010_2ndpoly_Jan_Apr_minusSD.tif")

AGC_LVOD_amazoniasmoothmeancrop <- raster::mask(crop(AGC_LVOD_amazoniasmoothmean,amazonBasin),amazonBasin)
AGC_LVOD_amazoniasmoothmeancropPlusSD <- raster::mask(crop(AGC_LVOD_amazoniasmoothmeanPlusSD,amazonBasin),amazonBasin)
AGC_LVOD_amazoniasmoothmeancropMinusSD <- raster::mask(crop(AGC_LVOD_amazoniasmoothmeanMinusSD,amazonBasin),amazonBasin)


#LVOD based AGC ccgfilt smooth max
AGC_LVOD_amazoniasmoothmax <- stack("./Outputs/AGC_LVOD_amazoniaBB_LVODccgfiltannual4wetmonthsmoothmax_EXTV3_Hui_inclDec2010_2ndpoly_Jan_Apr.tif")
AGC_LVOD_amazoniasmoothmaxPlusSD <- stack("./Outputs/AGC_LVOD_amazoniaBB_LVODccgfiltannual4wetmonthsmoothmax_EXTV3_Hui_inclDec2010_2ndpoly_Jan_Apr_plusSD.tif")
AGC_LVOD_amazoniasmoothmaxMinusSD <- stack("./Outputs/AGC_LVOD_amazoniaBB_LVODccgfiltannual4wetmonthsmoothmax_EXTV3_Hui_inclDec2010_2ndpoly_Jan_Apr_minusSD.tif")


AGC_LVOD_amazoniasmoothmaxcrop <- raster::mask(crop(AGC_LVOD_amazoniasmoothmax,amazonBasin),amazonBasin)
AGC_LVOD_amazoniasmoothmaxcropPlusSD <- raster::mask(crop(AGC_LVOD_amazoniasmoothmaxPlusSD,amazonBasin),amazonBasin)
AGC_LVOD_amazoniasmoothmaxcropMinusSD <- raster::mask(crop(AGC_LVOD_amazoniasmoothmaxMinusSD,amazonBasin),amazonBasin)

#LVOD based AGC ccgfilt trend mean
AGC_LVOD_amazoniatrendmean <- stack("./Outputs/AGC_LVOD_amazoniaBB_LVODccgfiltannual4wetmonthtrendmean_EXTV3_Hui_inclDec2010_2ndpoly_Jan_Apr.tif")
AGC_LVOD_amazoniatrendmeanPlusSD <- stack("./Outputs/AGC_LVOD_amazoniaBB_LVODccgfiltannual4wetmonthtrendmean_EXTV3_Hui_inclDec2010_2ndpoly_Jan_Apr_plusSD.tif")
AGC_LVOD_amazoniatrendmeanMinusSD <- stack("./Outputs/AGC_LVOD_amazoniaBB_LVODccgfiltannual4wetmonthtrendmean_EXTV3_Hui_inclDec2010_2ndpoly_Jan_Apr_minusSD.tif")


AGC_LVOD_amazoniatrendmeancrop <- raster::mask(crop(AGC_LVOD_amazoniatrendmean,amazonBasin),amazonBasin)
AGC_LVOD_amazoniatrendmeancropPlusSD <- raster::mask(crop(AGC_LVOD_amazoniatrendmeanPlusSD,amazonBasin),amazonBasin)
AGC_LVOD_amazoniatrendmeancropMinusSD <- raster::mask(crop(AGC_LVOD_amazoniatrendmeanMinusSD,amazonBasin),amazonBasin)

#mean of the three methods
AGC_LVOD_amazoniamethodsmeancrop <- (AGC_LVOD_amazoniasmoothmeancrop+AGC_LVOD_amazoniasmoothmaxcrop+AGC_LVOD_amazoniatrendmeancrop)/3
AGC_LVOD_amazoniamethodsmeancropPlusSD <- (AGC_LVOD_amazoniasmoothmeancropPlusSD+AGC_LVOD_amazoniasmoothmaxcropPlusSD+AGC_LVOD_amazoniatrendmeancropPlusSD)/3
AGC_LVOD_amazoniamethodsmeancropMinusSD <- (AGC_LVOD_amazoniasmoothmeancropMinusSD+AGC_LVOD_amazoniasmoothmaxcropMinusSD+AGC_LVOD_amazoniatrendmeancropMinusSD)/3


######
#make masks

#wetland masking
nonwetmaskresamp <- raster("./Datasets/AmazonNeverWetMaskResamp.tif")

wetmaskreference <- nonwetmaskresamp
wetmaskreference[wetmaskreference>0] <- 1
wetmaskreference[wetmaskreference==0] <- NA

#75% nonwet parameter subject to change, tradeoff between excluding flooded areas and masking too many pixels
wetmask25perclessflood <- nonwetmaskresamp>0.75
wetmask25perclessflood[is.na(wetmask25perclessflood)] <- 1
wetmask25perclessflood[wetmask25perclessflood<1] <- NA

#LVOD data available in every year

hasdata <- !is.na(calc(AGC_LVOD_amazoniamethodsmeancrop[[1:9]],sum))
hasdata[hasdata==0] <- NA

annualchangesformask <- AGC_LVOD_amazoniamethodsmeancrop[[2:10]]-AGC_LVOD_amazoniamethodsmeancrop[[1:9]]

test <- getValues(annualchangesformask)

#mask regions with changes over 20 Mg ha yr threshold
anomalousIncreaseMask <- calc(annualchangesformask,max,na.rm=T)<20
anomalousIncreaseDisp <- calc(annualchangesformask,max,na.rm=T)>=20
anomalousIncreaseDisp[anomalousIncreaseDisp<1] <- NA

anomalousIncreaseMask[anomalousIncreaseMask<1] <- NA

partialMask <- 
totalMask <- anomalousIncreaseMask*wetmask25perclessflood*hasdata
noMask <- AGC_LVOD_amazoniamethodsmeancrop[[1]]
noMask[!is.na(noMask)] <- 1

#write combination of masks for reference in other codes
#writeRaster(totalMask,'./CodeOutputs/totalMaskV18.tif',overwrite=T)

#display mask extent
levelplot(totalMask,col.regions=magma(100),margin=FALSE,maxpixels = 2e10,colorkey=F,par.settings=list(panel.background=list(col="lightgrey")))+latticeExtra::layer(sp.polygons(amazonBasin,col='black'))

#display threshold excluded areas
levelplot(anomalousIncreaseDisp,col.regions=magma(100),margin=FALSE,maxpixels = 2e10,colorkey=F,par.settings=list(panel.background=list(col="lightgrey")))+latticeExtra::layer(sp.polygons(amazonBasin,col='black'))


########
#mask L-VOD datasets
AGC_LVOD_amazoniamethodsmeancropPlusSD  <- raster::mask(AGC_LVOD_amazoniamethodsmeancropPlusSD ,totalMask)

AGC_LVOD_amazoniamethodsmeancropMinusSD <- raster::mask(AGC_LVOD_amazoniamethodsmeancropMinusSD,totalMask)

AGC_LVOD_amazoniamethodsmeancrop <- raster::mask(AGC_LVOD_amazoniamethodsmeancrop,totalMask)

AGC_LVOD_amazoniasmoothmeancrop <- raster::mask(AGC_LVOD_amazoniasmoothmeancrop,totalMask)

AGC_LVOD_amazoniasmoothmaxcrop <- raster::mask(AGC_LVOD_amazoniasmoothmaxcrop,totalMask)

AGC_LVOD_amazoniatrendmeancrop <- raster::mask(AGC_LVOD_amazoniatrendmeancrop,totalMask)


#write key datasets for use in other code

writeRaster(AGC_LVOD_amazoniasmoothmeancrop[[1:9]],'./CodeOutputs/AGC_LVOD_amazoniasmoothmeancropV18.tif',overwrite=T)
writeRaster(AGC_LVOD_amazoniasmoothmaxcrop[[1:9]],'./CodeOutputs/AGC_LVOD_amazoniasmoothmaxcropV18.tif',overwrite=T)
writeRaster(AGC_LVOD_amazoniatrendmeancrop[[1:9]],'./CodeOutputs/AGC_LVOD_amazoniatrendmeancropV18.tif',overwrite=T)
writeRaster(AGC_LVOD_amazoniamethodsmeancrop[[1:9]],'./CodeOutputs/AGC_LVOD_amazoniamethodsmeancropV18.tif',overwrite=T)
writeRaster(AGC_LVOD_amazoniamethodsmeancropPlusSD[[1:9]],'./CodeOutputs/AGC_LVOD_amazoniamethodsmeancropPlusSDV18.tif',overwrite=T)
writeRaster(AGC_LVOD_amazoniamethodsmeancropMinusSD[[1:9]],'./CodeOutputs/AGC_LVOD_amazoniamethodsmeancropMinusSDV18.tif',overwrite=T)




#Time series of mean L-VOD AGC

LVODAGCTSmethodsmean <- cellStats(AGC_LVOD_amazoniamethodsmeancrop*area(AGC_LVOD_amazoniamethodsmeancrop)*100,sum,na.rm=T)/10^9
LVODAGCTSmethodsmeanPlusSD <- cellStats(AGC_LVOD_amazoniamethodsmeancropPlusSD*area(AGC_LVOD_amazoniamethodsmeancropPlusSD)*100,sum,na.rm=T)/10^9
LVODAGCTSmethodsmeanMinusSD <- cellStats(AGC_LVOD_amazoniamethodsmeancropMinusSD*area(AGC_LVOD_amazoniamethodsmeancropMinusSD)*100,sum,na.rm=T)/10^9


methodmeans <- apply(cbind(LVODAGCTSsmoothmean,LVODAGCTSsmoothmax,LVODAGCTStrendmean),1,mean) 
methodSDs <- apply(cbind(LVODAGCTSsmoothmean,LVODAGCTSsmoothmax,LVODAGCTStrendmean),1,sd) 


###########
#load land cover datasets and resample to common resolution and extents (~1km and ~25km)

#mapbiomas forest fractional cover
mapbiomasforestfrac <- stack("D:/RECCAP/Datasets/mapbiomas/mapbiomasC2_1km_2009_2018_forestFractions.tif")[[3:10]]

mapbiomasforestfracresamp <- raster::resample(mapbiomasforestfrac,AGC_LVOD_amazoniamethodsmeancrop[[1]])
mapbiomasforestfracresamp <- mapbiomasforestfracresamp

mapbiomasforestfraccrop <- raster::mask(raster::mask(crop(mapbiomasforestfracresamp,amazonBasin),amazonBasin),totalMask)

#write raster for use in other code
writeRaster(mapbiomasforestfraccrop,'./CodeOutputs/mapbiomasforestfraccropV16.tif',overwrite=T)


#mapbiomas SF extent

mapbiomasSFfrac <- stack("D:/RECCAP/Datasets/mapbiomas/forestSFC2Frac_1km_2009_2018.tif")[[2:10]]

mapbiomasSFfracresamp <- raster::resample(mapbiomasSFfrac,AGC_LVOD_amazoniamethodsmeancrop[[1]])

mapbiomasSFfraccrop <- raster::mask(raster::mask(crop(mapbiomasSFfracresamp,amazonBasin),amazonBasin),totalMask)


#mapbiomas SF biomass, increase calculated in GEE

mapbiomasSFbiomassincrease <- stack("D:/RECCAP/Datasets/mapbiomas/SFbiomassgrowthAmazon_1km_2010_2018_V5.tif")[[2:9]]

#mapbiomas SF age

mapbiomasSFage <- stack("D:/RECCAP/Datasets/mapbiomas/forestSFC2Ages_1km_2009_2018.tif")[[2:10]]


#deforestation of intact forest
mapbiomasdeforestation <- stack("D:/RECCAP/Datasets/mapbiomas/mapbiomasC2_deforestationintactfrac_1km_2010_2018_NoSF_V4JRC.tif")[[2:9]]

mapbiomasdeforestationresamp <- raster::resample(mapbiomasdeforestation,AGC_LVOD_amazoniamethodsmeancrop[[1]])

mapbiomasdeforestationcrop <- raster::mask(raster::mask(crop(mapbiomasdeforestationresamp,amazonBasin),amazonBasin),totalMask)

#deforestation of degraded forest
mapbiomasdeforestationdegraded <- stack("D:/RECCAP/Datasets/mapbiomas/mapbiomasC2_deforestationdegradedfrac_1km_2010_2018_NoSF_V4JRC.tif")[[2:9]]#2011-2018

mapbiomasdeforestationdegradedresamp <- raster::resample(mapbiomasdeforestationdegraded,AGC_LVOD_amazoniamethodsmeancrop[[1]])

mapbiomasdeforestationdegradedcrop <- raster::mask(raster::mask(crop(mapbiomasdeforestationdegradedresamp,amazonBasin),amazonBasin),totalMask)

#deforestation of SF
mapbiomasdeforestationSF <- stack("D:/RECCAP/Datasets/mapbiomas/mapbiomasC2_deforestationSFfrac_1km_2010_2018_V4.tif")[[2:9]]#2011-2018

mapbiomasdeforestationSFresamp <- raster::resample(mapbiomasdeforestationSF,AGC_LVOD_amazoniamethodsmeancrop[[1]])

mapbiomasdeforestationSFcrop <- raster::mask(raster::mask(crop(mapbiomasdeforestationSFresamp,amazonBasin),amazonBasin),totalMask)

#deforestation of SF biomass, calculated in GEE
mapbiomasdeforestationSFbiomass <- stack("D:/RECCAP/Datasets/mapbiomas/SFbiomassdeforestedAmazon_1km_2010_2018_V5.tif")[[2:9]]#2011-2018


#make plot of different deforested areas (Extended Data Fig. 5)

mapbiomasdeforestationtotTS <- cellStats(mapbiomasdeforestationcrop*area(mapbiomasdeforestationcrop),sum,na.rm=T)

mapbiomasdeforestationdegradedtotTS <- cellStats(mapbiomasdeforestationdegradedcrop*area(mapbiomasdeforestationdegradedcrop),sum,na.rm=T)

mapbiomasdeforestationSFtotTS <- cellStats(mapbiomasdeforestationSFcrop*area(mapbiomasdeforestationSFcrop),sum,na.rm=T)

deforestationdf <- data.frame(years=seq(2011,2018,1),oldgrowth=mapbiomasdeforestationtotTS,degraded=mapbiomasdeforestationdegradedtotTS,secondary=mapbiomasdeforestationSFtotTS)

deforestationdfmelt <- melt(deforestationdf,id='years')

ggplot(data=deforestationdfmelt, aes(x=years,y=value, fill=variable)) +
  ylab(bquote('deforested area ['~ ~km^2~']'))+
  scale_fill_manual(labels = c("old-growth", "degraded","secondary"),values=c("#8dd3c7", "#ffffb3","#bebada"))+
  geom_bar(stat="identity")



#cumulative loss factor for edges

edgescumulativelossfactor <- stack("D:/RECCAP/Datasets/celsoEdge/edgeAgeLossFactorCumulativeAmazon_1km_2009_2017_V4.tif")[[2:9]]#2010-2017

#cumulative loss factor for edges that were deforested
#note: years indicative of years of detected deforestation, prior year of age is used
deforestededgescumulativelossfactor <- stack("D:/RECCAP/Datasets/celsoEdge/deforestedEdgeAgeLossFactorCumulativeAmazon_1km_2010_2018_V4.tif")[[2:9]]#2011-2018

#age of SF that was deforested
#note: years indicative of years of detected deforestation, prior year of age is used
deforestedSFage <- stack("D:/RECCAP/Datasets/mapbiomas/deforestedSFageAmazon_1km_2010_2018_V4.tif")[[2:9]]


#JRC new degraded areas, edges excluded 

JRCdegradationnotedge<- stack("D:/RECCAP/Datasets/JRCdata/degradedforestnotedge_1km_JRCcelso_2010_2018.tif")[[2:9]]

JRCdegradationnotedge1kmresamp <- raster::resample(JRCdegradationnotedge,mapbiomasforestfrac[[1]])

JRCdegradationnotedgeresamp <- raster::resample(JRCdegradationnotedge1kmresamp,AGC_LVOD_amazoniamethodsmeancrop[[1]])


JRCdegradationonlyedge<- stack("D:/RECCAP/Datasets/JRCdata/degradedforestonlyedge_1km_JRCcelso_2010_2018.tif")[[2:9]]

JRCdegradationonlyedge1kmresamp <- raster::resample(JRCdegradationonlyedge,mapbiomasforestfrac[[1]])

JRCdegradationonlyedgeresamp <- raster::resample(JRCdegradationonlyedge1kmresamp,AGC_LVOD_amazoniamethodsmeancrop[[1]])


#JRC degraded baseline for years 2010-2018, currently not used

# JRCdegradationbase<- stack("D:/RECCAP/Datasets/JRCdata/degradedforestonly_1km_JRC_baseline_2010_2018.tif")[[1:8]]
# 
# JRCdegradationbaseresamp <- resample(JRCdegradationbase,AGC_LVOD_amazoniamethodsmeancrop[[1]])
# JRCdegradationbasecrop <- mask(mask(crop(JRCdegradationbaseresamp,amazonBasin),amazonBasin),totalMask)

#JRC degraded baseline for years 2010-2018, excluding SF

JRCdegradationbasenoSFnoEdge<- stack("D:/RECCAP/Datasets/JRCdata/degradedforestonlynoSFnoEdge_1km_JRC_baseline_2010_2018.tif")[[2:9]]
# 
#JRCdegradationbasenoSFnoEdge1kmresamp <- resample(JRCdegradationbasenoSFnoEdge,mapbiomasforestfrac[[1]])
# 
 JRCdegradationbaseresampnoSFnoEdge <- raster::resample(JRCdegradationbasenoSFnoEdge,AGC_LVOD_amazoniamethodsmeancrop[[1]])
 JRCdegradationbasenoSFnoEdgecrop <- raster::mask(raster::mask(crop(JRCdegradationbaseresampnoSFnoEdge,amazonBasin),amazonBasin),totalMask)


#degraded forest edges

forestedgesfrac <- stack("D:/RECCAP/Datasets/celsoEdge/forestAnthroEdgeAmazonFrac_ExclRecovered_1km_2010_2018.tif")[[2:9]]
forestedgesfracresamp <- raster::resample(forestedgesfrac,AGC_LVOD_amazoniamethodsmeancrop[[1]])
forestedgesfraccrop <- raster::mask(raster::mask(crop(forestedgesfracresamp,amazonBasin),amazonBasin),totalMask)

forestedgeslossfactor <- stack("D:/RECCAP/Datasets/celsoEdge/forestEdgeAgeLossFactorAmazon_ExclRecovered_1km_2009_2018.tif")[[3:10]]

forestedgescumulativelossfactor <- stack("D:/RECCAP/Datasets/celsoEdge/forestEdgeAgeLossFactorCumulativeAmazon_1km_2009_2018.tif")[[3:10]]


#intact forest based on mapbiomas forest cover removing SF, edges and disturbed forest pixels
intactforest1km <- stack("D:/RECCAP/Datasets/mapbiomas/intactforest_1km_mapbiomas_2010_2018.tif")[[2:9]]

intactforest <- raster::resample(intactforest1km,AGC_LVOD_amazoniamethodsmeancrop[[1]])
intactforestcrop <- raster::mask(raster::mask(crop(intactforest,amazonBasin),amazonBasin),totalMask)

#write raster for use in other code
writeRaster(intactforestcrop,'./CodeOutputs/intactforestcropV16.tif',overwrite=T)

##########
#read intact AGC reference calculated from ESA CCI biomass map adjacent cells with >90% old growth forest cover
#generated in calculateIntactForestReferenceV16.R

intactAGCrefadj <- raster("D:/RECCAP/Datasets/ESACCIv2_2017amazonTot1km_JRCintactAdjRef.tif")
intactAGCrefadjplusSD <- raster("D:/RECCAP/Datasets/ESACCIv2_2017amazonTot1km_JRCintactAdjRefPlusSD.tif")
intactAGCrefadjminusSD <- raster("D:/RECCAP/Datasets/ESACCIv2_2017amazonTot1km_JRCintactAdjRefMinusSD.tif")

#########

#functions to derive modelled AGC change for each process

getDeforestedBiomass <- function(AGCmap,deforestedmap,inmask=NULL){
  if(is.null(inmask)){
    inmask <- noMask
  }

  neighbourhoodintactAGC <- AGCmap

  deforestedAGC <- neighbourhoodintactAGC*deforestedmap
  
  deforestedAGCmaskedresamp <- raster::mask(raster::resample(deforestedAGC,AGC_LVOD_amazoniamethodsmeancrop[[1]]),totalMask)
  
  return(raster::mask(deforestedAGCmaskedresamp,inmask))
}


getDeforestedDegradedBiomass <- function(AGCmap,deforesteddegradedmap,cumulativeedgelossfactorsmap,inmask=NULL){
  if(is.null(inmask)){
    inmask <- noMask
  }

  neighbourhoodintactAGC <- AGCmap

  edgedegradationref <-  neighbourhoodintactAGC*cumulativeedgelossfactorsmap
  deforesteddegradedAGC <- (neighbourhoodintactAGC-edgedegradationref)*deforesteddegradedmap
  
  deforesteddegradedAGCmaskedresamp <- raster::mask(raster::resample(deforesteddegradedAGC,AGC_LVOD_amazoniamethodsmeancrop[[1]]),totalMask)
  return(raster::mask(deforesteddegradedAGCmaskedresamp,inmask))
}


getDeforestedSFBiomass <- function(deforestedSFbiomass,inmask=NULL){
  if(is.null(inmask)){
    inmask <- noMask
  }
  
  SFdeforestedAGC <-  deforestedSFbiomass
  
  deforestedSFAGCmaskedresamp <- raster::mask(raster::resample(SFdeforestedAGC,AGC_LVOD_amazoniamethodsmeancrop[[1]]),totalMask)
  return(raster::mask(deforestedSFAGCmaskedresamp,inmask))
}

getDegradedBiomass <- function(AGCmap,degradednotedgemap,degradedonlyedgemap,cumulativeedgelossfactorsmap,inmask=NULL){
  if(is.null(inmask)){
    inmask <- noMask
  }

  neighbourhoodintactAGC <- AGCmap
  
  #get reference value for local edge biomass
  edgedegradationref <-  neighbourhoodintactAGC*cumulativeedgelossfactorsmap
  edgebiomassref <- (neighbourhoodintactAGC-edgedegradationref)
  
  degradededgediff <- edgebiomassref-edgebiomassref*0.64684 #factor relating degraded forest to intact
  degradededgeAGC <- degradededgediff*degradedonlyedgemap
  
  degradedintactdiff <- neighbourhoodintactAGC-neighbourhoodintactAGC*0.64684
  degradedintactAGC <- degradedintactdiff*degradednotedgemap
  
  degradedAGC <- degradedintactAGC+degradededgeAGC
  
  degradedAGCmaskedresamp <- raster::mask(raster::resample(degradedAGC,AGC_LVOD_amazoniamethodsmeancrop[[1]]),totalMask)
  
  return(raster::mask(degradedAGCmaskedresamp,inmask))
}


getEdgeBiomass <- function(AGCmap,edgefracmap,edgelossfactormap,inmask=NULL){
  if(is.null(inmask)){
    inmask <- noMask
  }

  neighbourhoodintactAGC <- AGCmap

  edgedegradationref <-  neighbourhoodintactAGC*edgelossfactormap
  
  #calculate biomass loss by edges per pixel
  edgedegradationtot <- edgedegradationref*edgefracmap
  
  edgedegradationtotresampmasked <- raster::mask(raster::resample(edgedegradationtot,AGC_LVOD_amazoniamethodsmeancrop[[1]]),totalMask)
  
  return(raster::mask(edgedegradationtotresampmasked,inmask))
}

#derive change in intact forest per pixel from AGC differences
getIntactChange <- function(AGCmapold,AGCmap,intactmap,nonintactchange,inmask=NULL){

  if(is.null(inmask)){
    inmask <- noMask
  }

  fullyintact <- intactmap>0.9 
  fullyintact[fullyintact==0] <- NA
  AGCintact <- AGCmap*fullyintact
  AGCintact[!is.finite(AGCintact)] <- NA
  AGCintactold <- AGCmapold*fullyintact
  AGCintactold[!is.finite(AGCintactold)] <- NA
  
  AGCintactchange <- AGCintact-AGCintactold-nonintactchange
  
  #focal operation to get IF reference values for mixed cells
  fwModel <- focalWeight(AGCintactchange, 2.5, type='circle')#circle of values to consider, needs to be large enough to contain reference values
  fwModel[fwModel>0] <- 1
  fwModel[fwModel==0] <- NA
  neighbourhoodintactAGCchange <- focal(AGCintactchange,fwModel,na.rm=T,fun=median,NAonly=T,pad=T,padValue=NA)

  #second pass focal operation with bigger window size to include areas in greater distance from intact forest gridcells
  fwModel <- focalWeight(AGCintactchange, 5, type='circle')
  fwModel[fwModel>0] <- 1
  fwModel[fwModel==0] <- NA
  neighbourhoodintactAGC2ndpasschange <- focal(neighbourhoodintactAGCchange,fwModel,na.rm=T,fun=median,NAonly=T,pad=T,padValue=NA)

  intactmap[intactmap>0.9] <- 1
  
  intactdiff <- neighbourhoodintactAGC2ndpasschange
  intactdiffAGC <- raster::mask(intactdiff*intactmap,totalMask)
  return(raster::mask(intactdiffAGC,inmask))
}


getSFBiomassIncrease <- function(SFbiomassincrease,inmask=NULL){
  if(is.null(inmask)){
    inmask <- noMask
  }
  
  SFgrowth <- SFbiomassincrease
  
  SFgrowthresamp <- raster::resample(SFgrowth,AGC_LVOD_amazoniamethodsmeancrop[[1]])
  
  SFgrowthresampcrop <- raster::mask(raster::mask(crop(SFgrowthresamp,amazonBasin),amazonBasin),totalMask)
  
  return(raster::mask(SFgrowthresampcrop,inmask))
}

###########
#compute raster maps and time series of modelled AGC loss

#time series of AGC loss
deforestedBiomassmapbiomasTS <- vector(length=8)
deforestedBiomassmapbiomasPlusSDTS <- vector(length=8)
deforestedBiomassmapbiomasMinusSDTS <- vector(length=8)

deforestedSFBiomassTS <- vector(length=8)
degradedBiomassTS <- vector(length=8)
degradedBiomassPlusSDTS <- vector(length=8)
degradedBiomassMinusSDTS <- vector(length=8)

intactBiomasschangemethodsmeanTS <- vector(length=8)
intactBiomasschangemethodsmeanPlusSDTS  <- vector(length=8)
intactBiomasschangemethodsmeanMinusSDTS  <- vector(length=8)

intactBiomasschangesmoothmeanTS <- vector(length=8)
intactBiomasschangesmoothmaxTS <- vector(length=8)
intactBiomasschangetrendmeanTS <- vector(length=8)

edgeBiomassTS <- vector(length=8)
edgeBiomassPlusSDTS <- vector(length=8)
edgeBiomassMinusSDTS <- vector(length=8)

SFBiomassTS <- vector(length=8)

LVODAGCTS <- cellStats(AGC_LVOD_amazoniamethodsmeancrop*area(AGC_LVOD_amazoniamethodsmeancrop)*100,sum,na.rm=T)/10^9

#raster time series of AGC loss
deforestedBiomassmapbiomasstack <- stack()
deforestedBiomassmapbiomasPlusSDstack <- stack()
deforestedBiomassmapbiomasMinusSDstack <- stack()

deforestedSFBiomassstack <- stack()
degradedBiomassstack <-  stack()
degradedBiomassPlusSDstack <-  stack()
degradedBiomassMinusSDstack <-  stack()

nonintactBiomasschange <- stack()

intactBiomasschangemethodsmeanstack <- stack()
intactBiomasschangemethodsmeanPlusSDstack <- stack()
intactBiomasschangemethodsmeanMinusSDstack <- stack()

intactBiomasschangesmoothmeanstack <- stack()
intactBiomasschangesmoothmaxstack <- stack()
intactBiomasschangetrendmeanstack <- stack()

edgeBiomassstack <-  stack()
edgeBiomassPlusSDstack <-  stack()
edgeBiomassMinusSDstack <-  stack()


SFBiomassstack <-  stack()

subsetmask <- totalMask

for(i in 1:8){ #(2011-2018)


deforestedBiomassmapbiomasstack <- stack(deforestedBiomassmapbiomasstack,getDeforestedBiomass(intactAGCrefadj,mapbiomasdeforestation[[i]])+
                                           getDeforestedDegradedBiomass(intactAGCrefadj,mapbiomasdeforestationdegraded[[i]],deforestededgescumulativelossfactor[[i]])+
                                           getDeforestedSFBiomass(mapbiomasdeforestationSFbiomass[[i]]))

deforestedBiomassmapbiomasPlusSDstack<- stack(deforestedBiomassmapbiomasPlusSDstack,getDeforestedBiomass(intactAGCrefadjplusSD,mapbiomasdeforestation[[i]])+
                                                getDeforestedDegradedBiomass(intactAGCrefadjplusSD,mapbiomasdeforestationdegraded[[i]],deforestededgescumulativelossfactor[[i]])+
                                                getDeforestedSFBiomass(mapbiomasdeforestationSFbiomass[[i]]))

deforestedBiomassmapbiomasMinusSDstack<- stack(deforestedBiomassmapbiomasMinusSDstack,getDeforestedBiomass(intactAGCrefadjminusSD,mapbiomasdeforestation[[i]])+
                                                getDeforestedDegradedBiomass(intactAGCrefadjminusSD,mapbiomasdeforestationdegraded[[i]],deforestededgescumulativelossfactor[[i]])+
                                                getDeforestedSFBiomass(mapbiomasdeforestationSFbiomass[[i]]))


deforestedBiomassmapbiomasTS[i] <- cellStats(raster::mask(deforestedBiomassmapbiomasstack[[i]],subsetmask)*area(deforestedBiomassmapbiomasstack[[i]])*100,sum,na.rm=T)/10^9
deforestedBiomassmapbiomasPlusSDTS[i] <- cellStats(raster::mask(deforestedBiomassmapbiomasPlusSDstack[[i]],subsetmask)*area(deforestedBiomassmapbiomasPlusSDstack[[i]])*100,sum,na.rm=T)/10^9
deforestedBiomassmapbiomasMinusSDTS[i] <- cellStats(raster::mask(deforestedBiomassmapbiomasMinusSDstack[[i]],subsetmask)*area(deforestedBiomassmapbiomasMinusSDstack[[i]])*100,sum,na.rm=T)/10^9

degradedBiomassstack <-  stack(degradedBiomassstack,getDegradedBiomass(intactAGCrefadj,JRCdegradationnotedge1kmresamp[[i]],JRCdegradationonlyedge1kmresamp[[i]],edgescumulativelossfactor[[i]]))
degradedBiomassPlusSDstack <-  stack(degradedBiomassPlusSDstack,getDegradedBiomass(intactAGCrefadjplusSD,JRCdegradationnotedge1kmresamp[[i]],JRCdegradationonlyedge1kmresamp[[i]],edgescumulativelossfactor[[i]]))
degradedBiomassMinusSDstack <-  stack(degradedBiomassMinusSDstack,getDegradedBiomass(intactAGCrefadjminusSD,JRCdegradationnotedge1kmresamp[[i]],JRCdegradationonlyedge1kmresamp[[i]],edgescumulativelossfactor[[i]]))

degradedBiomassTS[i] <- cellStats(raster::mask(degradedBiomassstack[[i]],subsetmask)*area(degradedBiomassstack[[i]])*100,sum,na.rm=T)/10^9
degradedBiomassPlusSDTS[i] <- cellStats(raster::mask(degradedBiomassPlusSDstack[[i]],subsetmask)*area(degradedBiomassPlusSDstack[[i]])*100,sum,na.rm=T)/10^9
degradedBiomassMinusSDTS[i] <- cellStats(raster::mask(degradedBiomassMinusSDstack[[i]],subsetmask)*area(degradedBiomassMinusSDstack[[i]])*100,sum,na.rm=T)/10^9

edgeBiomassstack <- stack(edgeBiomassstack,getEdgeBiomass(intactAGCrefadj,forestedgesfrac[[i]],forestedgeslossfactor[[i]]))
edgeBiomassPlusSDstack <- stack(edgeBiomassPlusSDstack,getEdgeBiomass(intactAGCrefadjplusSD,forestedgesfrac[[i]],forestedgeslossfactor[[i]]))
edgeBiomassMinusSDstack <- stack(edgeBiomassMinusSDstack,getEdgeBiomass(intactAGCrefadjminusSD,forestedgesfrac[[i]],forestedgeslossfactor[[i]]))

edgeBiomassTS[i] <- cellStats(raster::mask(edgeBiomassstack[[i]],subsetmask)*area(edgeBiomassstack[[i]])*100,sum,na.rm=T)/10^9
edgeBiomassPlusSDTS[i] <- cellStats(raster::mask(edgeBiomassPlusSDstack[[i]],subsetmask)*area(edgeBiomassPlusSDstack[[i]])*100,sum,na.rm=T)/10^9
edgeBiomassMinusSDTS[i] <- cellStats(raster::mask(edgeBiomassMinusSDstack[[i]],subsetmask)*area(edgeBiomassMinusSDstack[[i]])*100,sum,na.rm=T)/10^9

SFBiomassstack <- stack(SFBiomassstack,getSFBiomassIncrease(mapbiomasSFbiomassincrease[[i]]))
SFBiomassTS[i] <- cellStats(raster::mask(SFBiomassstack[[i]],subsetmask)*area(SFBiomassstack[[i]])*100,sum,na.rm=T)/10^9

#sum of all non-intact forest modelled change
nonintactBiomasschange <- stack(nonintactBiomasschange,SFBiomassstack[[i]]-deforestedBiomassmapbiomasstack[[i]]-degradedBiomassstack[[i]]-edgeBiomassstack[[i]])

intactBiomasschangemethodsmeanstack <- stack(intactBiomasschangemethodsmeanstack ,getIntactChange(AGC_LVOD_amazoniamethodsmeancrop[[i]],AGC_LVOD_amazoniamethodsmeancrop[[1+i]],intactforestcrop[[i]],nonintactBiomasschange[[i]]))
intactBiomasschangemethodsmeanTS[i] <- cellStats(raster::mask(intactBiomasschangemethodsmeanstack[[i]],subsetmask)*area(intactBiomasschangemethodsmeanstack[[i]])*100,sum,na.rm=T)/10^9

intactBiomasschangemethodsmeanPlusSDstack <- stack(intactBiomasschangemethodsmeanPlusSDstack,getIntactChange(AGC_LVOD_amazoniamethodsmeancropPlusSD[[i]],AGC_LVOD_amazoniamethodsmeancropPlusSD[[1+i]],intactforestcrop[[i]],nonintactBiomasschange[[i]]))
intactBiomasschangemethodsmeanPlusSDTS[i] <- cellStats(raster::mask(intactBiomasschangemethodsmeanPlusSDstack[[i]],subsetmask)*area(intactBiomasschangemethodsmeanPlusSDstack[[i]])*100,sum,na.rm=T)/10^9

intactBiomasschangemethodsmeanMinusSDstack <- stack(intactBiomasschangemethodsmeanMinusSDstack ,getIntactChange(AGC_LVOD_amazoniamethodsmeancropMinusSD[[i]],AGC_LVOD_amazoniamethodsmeancropMinusSD[[1+i]],intactforestcrop[[i]],nonintactBiomasschange[[i]]))
intactBiomasschangemethodsmeanMinusSDTS[i] <- cellStats(raster::mask(intactBiomasschangemethodsmeanMinusSDstack[[i]],subsetmask)*area(intactBiomasschangemethodsmeanMinusSDstack[[i]])*100,sum,na.rm=T)/10^9

intactBiomasschangesmoothmeanstack <- stack(intactBiomasschangesmoothmeanstack,getIntactChange(AGC_LVOD_amazoniasmoothmeancrop[[i]],AGC_LVOD_amazoniasmoothmeancrop[[1+i]],intactforestcrop[[i]],nonintactBiomasschange[[i]]))
intactBiomasschangesmoothmeanTS[i] <- cellStats(raster::mask(intactBiomasschangesmoothmeanstack[[i]],subsetmask)*area(intactBiomasschangesmoothmeanstack[[i]])*100,sum,na.rm=T)/10^9

intactBiomasschangesmoothmaxstack <- stack(intactBiomasschangesmoothmaxstack,getIntactChange(AGC_LVOD_amazoniasmoothmaxcrop[[i]],AGC_LVOD_amazoniasmoothmaxcrop[[1+i]],intactforestcrop[[i]],nonintactBiomasschange[[i]]))
intactBiomasschangesmoothmaxTS[i] <- cellStats(raster::mask(intactBiomasschangesmoothmaxstack[[i]],subsetmask)*area(intactBiomasschangesmoothmaxstack[[i]])*100,sum,na.rm=T)/10^9

intactBiomasschangetrendmeanstack <- stack(intactBiomasschangetrendmeanstack,getIntactChange(AGC_LVOD_amazoniatrendmeancrop[[i]],AGC_LVOD_amazoniatrendmeancrop[[1+i]],intactforestcrop[[i]],nonintactBiomasschange[[i]]))
intactBiomasschangetrendmeanTS[i] <- cellStats(raster::mask(intactBiomasschangetrendmeanstack[[i]],subsetmask)*area(intactBiomasschangetrendmeanstack[[i]])*100,sum,na.rm=T)/10^9


}


#write intermediate outputs for use in other code

#deforested biomass
writeRaster(deforestedBiomassmapbiomasstack,'./CodeOutputs/deforestedbiomassstackV17.tif',overwrite=T)
writeRaster(deforestedBiomassmapbiomasPlusSDstack,'./CodeOutputs/deforestedbiomassPlusSDstackV17.tif',overwrite=T)
writeRaster(deforestedBiomassmapbiomasMinusSDstack,'./CodeOutputs/deforestedbiomassMinusSDstackV17.tif',overwrite=T)


#degraded Biomass
writeRaster(degradedBiomassstack,'./CodeOutputs/degradedbiomassstackV17.tif',overwrite=T)
writeRaster(degradedBiomassPlusSDstack,'./CodeOutputs/degradedbiomassPlusSDstackV17.tif',overwrite=T)
writeRaster(degradedBiomassMinusSDstack,'./CodeOutputs/degradedbiomassMinusSDstackV17.tif',overwrite=T)


#edge Biomass change
writeRaster(edgeBiomassstack,'./CodeOutputs/edgebiomasschangestackV17.tif',overwrite=T)
writeRaster(edgeBiomassPlusSDstack,'./CodeOutputs/edgebiomassPlusSDchangestackV17.tif',overwrite=T)
writeRaster(edgeBiomassMinusSDstack,'./CodeOutputs/edgebiomassMinusSDchangestackV17.tif',overwrite=T)


#SF growth
writeRaster(SFBiomassstack,'./CodeOutputs/SFBiomassgrowthstackV17.tif',overwrite=T)

#intact forest changes
writeRaster(intactBiomasschangemethodsmeanstack,'./CodeOutputs/intactBiomasschangemethodsmeanstackV18.tif',overwrite=T)
writeRaster(intactBiomasschangemethodsmeanPlusSDstack ,'./CodeOutputs/intactBiomasschangemethodsmeanPlusSDstackV18.tif',overwrite=T)
writeRaster(intactBiomasschangemethodsmeanMinusSDstack ,'./CodeOutputs/intactBiomasschangemethodsmeanMinusSDstackV18.tif',overwrite=T)

writeRaster(intactBiomasschangesmoothmeanstack,'./CodeOutputs/intactBiomasschangesmoothmeanstackv17.tif')
writeRaster(intactBiomasschangesmoothmaxstack,'./CodeOutputs/intactBiomasschangesmoothmaxstackv17.tif')
writeRaster(intactBiomasschangetrendmeanstack,'./CodeOutputs/intactBiomasschangetrendmeanstackv17.tif')


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

SFBiomassstack <- stack('./CodeOutputs/SFBiomassgrowthstackv17.tif')


intactBiomasschangemethodsmeanstack <- stack('./CodeOutputs/intactBiomasschangemethodsmeanstackV18.tif')
intactBiomasschangemethodsmeanPlusSDstack <- stack('./CodeOutputs/intactBiomasschangemethodsmeanPlusSDstackV18.tif')
intactBiomasschangemethodsmeanMinusSDstack <- stack('./CodeOutputs/intactBiomasschangemethodsmeanMinusSDstackV18.tif')

intactBiomasschangesmoothmeanstack <- stack('./CodeOutputs/intactBiomasschangesmoothmeanstackv17.tif')
intactBiomasschangesmoothmaxstack <- stack('./CodeOutputs/intactBiomasschangesmoothmaxstackv17.tif')
intactBiomasschangetrendmeanstack <- stack('./CodeOutputs/intactBiomasschangetrendmeanstackv17.tif')

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


modelAGCchangeTScomb <- modelAGCchange_inclIF
modelAGCchangeTScombmax <-  modelAGCchange_inclIF+modelAGCchange_inclIF_combSD
modelAGCchangeTScombmin <- modelAGCchange_inclIF-modelAGCchange_inclIF_combSD

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

pdf(file = "./Figures/FinFigs/Fig3b_V18_nolvoderror.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 5) # The height of the plot in inches

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
  #geom_ribbon(aes(ymin = methodmeans - methodSDs,
  #                ymax = methodmeans + methodSDs), alpha = 0.2,fill='#009988')+
  geom_ribbon(aes(ymin = modelAGC_inclIFmin,
                ymax = modelAGC_inclIFmax), alpha = 0.2,fill= '#EE7733')
plot(p)

dev.off()

#make seperate legend
plot(methodcompdf$years,methodcompdf$methodmeans)
legend(2015,55.5, legend=c("L-VOD AGC", "modelled AGC"),
       col=c('#009988', '#EE7733'), lty=c(1,2), cex=0.8,box.lty=0)
 

#

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

levelplot(intactmask,col.regions=magma(100),margin=FALSE,maxpixels = 2e10,colorkey=F,par.settings=list(panel.background=list(col="lightgrey")))+latticeExtra::layer(sp.polygons(amazonBasin,col='black'))


par(pty='s')
LVODchangenotintactvals <- getValues(LVODNetAGCchange*notintactmask)
modelchangenotintactvals <- getValues(modelNetAGCchange*notintactmask)
plot(modelchangenotintactvals,LVODchangenotintactvals,pch=16,col=alpha('black',0.05),xlim=c(-40,25),ylim=c(-40,25),xlab=bquote("Modelled AGC change 2010-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),ylab=bquote("LVOD AGC change 2010-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"))
LVODmodellinmod <- lm(LVODchangenotintactvals~modelchangenotintactvals)
abline(0,1,lty=2)
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(LVODmodellinmod$coefficients[1],LVODmodellinmod$coefficients[2])

source("D:/RECCAP/Code/divergePlotFuns.R")


modelNetAGCchangelog <- log(abs(modelNetAGCchange),10)

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

pdf(file = "./Figures/FinFigs/Fig2aV17RdBu.pdf",   
    width = 6,
    height = 5) 

#plot log spaced color bar for combined model net change
levelplot(modelNetAGCchange,main=bquote("Modelled AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),at=breaklogstepsfin,colorkey = list(at = breaklogstepsfin),margin=FALSE,maxpixels = 2e10,par.settings=myTheme)+latticeExtra::layer(sp.polygons(amazonBasin,col='darkgrey'))#list(panel.background=list(col="lightgrey")))

dev.off()

#plot maps of biomass change per process

#Plot deforestation biomass change

deforestedBiomasstotchange <- calc(deforestedBiomassmapbiomasstack,sum)

pdf(file = "./Figures/FinFigs/Fig2bV17RdBu.pdf",   
    width = 6, 
    height = 5) 
#diverge0(levelplot(-deforestedBiomasstotchange ,main=bquote("Modelled AGC change deforestation 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),at=breaklogstepsfin,margin=FALSE,maxpixels = 2e10,par.settings=list(panel.background=list(col="lightgrey"))),colorRampPalette(c('red','white','blue')))+latticeExtra::layer(sp.polygons(amazonBasin,col='darkgrey'))
levelplot(-deforestedBiomasstotchange ,main=bquote("Modelled AGC change deforestation 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),at=breaklogstepsfin,colorkey = list(at = breaklogstepsfin),margin=FALSE,maxpixels = 2e10,par.settings=myTheme)+latticeExtra::layer(sp.polygons(amazonBasin,col='darkgrey'))#list(panel.background=list(col="lightgrey")))
dev.off()
#get maximum deforested
maxchange=cellStats(deforestedBiomasstotchange,max)
maxchangeperyear=maxchange/8

#Plot degraded biomass change


degradedBiomasstotchange <- calc(degradedBiomassstack,sum)+calc(edgeBiomassstack,sum)

pdf(file = "./Figures/FinFigs/Fig2cV17RdBu.pdf",   
    width = 6, 
    height = 5) 
#diverge0(levelplot(-degradedBiomasstotchange,main=bquote("Modelled AGC change degradation 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),at=breaklogstepsfin,margin=FALSE,maxpixels = 2e10,par.settings=list(panel.background=list(col="lightgrey"))),colorRampPalette(c('red','white','blue')))+latticeExtra::layer(sp.polygons(amazonBasin,col='darkgrey'))
levelplot(-degradedBiomasstotchange,main=bquote("Modelled AGC change degradation 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),at=breaklogstepsfin,colorkey = list(at = breaklogstepsfin),margin=FALSE,maxpixels = 2e10,par.settings=myTheme)+latticeExtra::layer(sp.polygons(amazonBasin,col='darkgrey'))#list(panel.background=list(col="lightgrey")))
dev.off()

#Plot SF biomass change

SFBiomasstotchange <- calc(SFBiomassstack,sum)

pdf(file = "./Figures/FinFigs/Fig2dV17RdBu.pdf",  
    width = 6, 
    height = 5) 
#diverge0(levelplot(SFBiomasstotchange,main=bquote("Modelled AGC change SF growth 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),at=breaklogstepsfin,margin=FALSE,maxpixels = 2e10,par.settings=list(panel.background=list(col="lightgrey"))),colorRampPalette(c('red','white','blue')))+latticeExtra::layer(sp.polygons(amazonBasin,col='darkgrey'))
levelplot(SFBiomasstotchange,main=bquote("Modelled AGC change SF growth 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),at=breaklogstepsfin,colorkey = list(at = breaklogstepsfin),margin=FALSE,maxpixels = 2e10,par.settings=myTheme)+latticeExtra::layer(sp.polygons(amazonBasin,col='darkgrey'))#list(panel.background=list(col="lightgrey")))
dev.off()
#Plot intact biomass change


intactBiomasstotchange <- (calc(intactBiomasschangesmoothmeanstack,sum)+calc(intactBiomasschangesmoothmaxstack,sum)+calc(intactBiomasschangetrendmeanstack,sum))/3

pdf(file = "./Figures/FinFigs/Fig2eV17RdBu.pdf",  
    width = 6, 
    height = 5) 
#diverge0(levelplot(intactBiomasstotchange,main=bquote("Modelled AGC change intact 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),at=breaklogstepsfin,margin=FALSE,maxpixels = 2e10,par.settings=list(panel.background=list(col="lightgrey"))),colorRampPalette(c('red','white','blue')))+latticeExtra::layer(sp.polygons(amazonBasin,col='darkgrey'))
levelplot(intactBiomasstotchange,main=bquote("Modelled AGC change intact 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),at=breaklogstepsfin,colorkey = list(at = breaklogstepsfin),margin=FALSE,maxpixels = 2e10,par.settings=myTheme)+latticeExtra::layer(sp.polygons(amazonBasin,col='darkgrey'))#list(panel.background=list(col="lightgrey")))
dev.off()


#make figure of difference between modelled and L-VOD AGC changes
#Fig. S5
 
#masking >90% intact areas

LVODNetAGCchangeminmax <- getMinMax(raster::mask((LVODNetAGCchange-(modelNetAGCchange)),notintactmask))
# 
pdf(file = "./Figures/FinFigs/ExtDataFig1corr.pdf",  
    width = 6, 
    height = 5) 
diverge0(levelplot(raster::mask((LVODNetAGCchange-(modelNetAGCchange)),notintactmask),main=bquote("L-VOD minus modelled net AGC change 2011-2019" ~ "["~ Mg ~ C ~ ha^{-1}  ~"]"),at=seq(LVODNetAGCchangeminmax[1],LVODNetAGCchangeminmax[2], len = 100),margin=FALSE,maxpixels = 2e10,par.settings=list(panel.background=list(col="lightgrey"))),'RdBu')+latticeExtra::layer(sp.polygons(amazonBasin,col='black'))
dev.off()

#relative difference, capped at 200%

relNetAGCchangediffLVOD[relNetAGCchangediffLVOD>200] <- 200
relNetAGCchangediffLVOD[relNetAGCchangediffLVOD<(-200)] <- (-200)

relNetAGCchangediffLVODnotintact <-  raster::mask(relNetAGCchangediffLVOD,notintactmask)

relNetAGCchangediffminmax <- getMinMax(relNetAGCchangediffLVODnotintact)

diverge0(levelplot(relNetAGCchangediffLVODnotintact,main=bquote("rel. diff. LVOD minus modelled net AGC change 2011-2019" ~ "[ %"~ Mg ~ C ~ ha^{-1}  ~"]"),at=seq(relNetAGCchangediffminmax[1],relNetAGCchangediffminmax[2], len = 100),margin=FALSE,maxpixels = 2e10,par.settings=list(panel.background=list(col="lightgrey"))),colorRampPalette(c('red','white','blue')))+latticeExtra::layer(sp.polygons(amazonBasin,col='black'))

