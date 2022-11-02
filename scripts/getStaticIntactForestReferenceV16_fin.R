
#Author: D. Fawcett
#
#Description:
#This code uses the ESA CCI biomass map of 2017 and precomputed values of intact forest cover to generate a reference biomass map of intact forest AGC for each 1 km cell used in further processing

#Outputs:
#-intact forest reference AGC maps, plus and minus SD of the original maps


#Date: 16/08/2022

library(mblm)
library(viridis)
library(rgdal)
library(Kendall)
library(trend)
library(raster)
library(scales)
library(ncdf4)
library(gridExtra)
library(rasterVis)
library(maptools)
library(tls)
library(epiR)
library(ggplot2)
library(rgdal) 


###########
#ESA CCI Biomass products, calculate intact forest reference values using local median filtering

ESACCIbiomass100m <- raster("D:/RECCAP/Datasets/ESACCIv2_biomass_AmazonTot_fix.tif")
ESACCIbiomass1km <- raster::resample(ESACCIbiomass100m,mapbiomasforestfrac[[1]],filename="D:/RECCAP/Datasets/ESACCIv2_biomass_AmazonTot1km.tif",overwrite=T)
ESACCIbiomass100mplusSD <- raster("D:/RECCAP/Datasets/ESACCIv2_biomass_AmazonTot_plusSD_fix.tif")
ESACCIbiomass1kmplusSD <- raster::resample(ESACCIbiomass100mplusSD,mapbiomasforestfrac[[1]],filename="D:/RECCAP/Datasets/ESACCIv2_biomass_AmazonTot1km_plusSD.tif",overwrite=T)
ESACCIbiomass100mminusSD <- raster("D:/RECCAP/Datasets/ESACCIv2_biomass_AmazonTot_minusSD_fix.tif")
ESACCIbiomass1kmminusSD <- raster::resample(ESACCIbiomass100mminusSD,mapbiomasforestfrac[[1]],filename="D:/RECCAP/Datasets/ESACCIv2_biomass_AmazonTot1km_minusSD.tif",overwrite=T)

# 
ESACCIbiomass1km <- raster("D:/RECCAP/Datasets/ESACCIv2_biomass_AmazonTot1km.tif")
ESACCIbiomass1km <- ESACCIbiomass1km*0.47 #convert AGB to AGC

ESACCIbiomass1kmplusSD <- raster("D:/RECCAP/Datasets/ESACCIv2_biomass_AmazonTot1km_plusSD.tif")
ESACCIbiomass1kmplusSD  <- ESACCIbiomass1kmplusSD*0.47 #convert AGB to AGC

ESACCIbiomass1kmminusSD <- raster("D:/RECCAP/Datasets/ESACCIv2_biomass_AmazonTot1km_minusSD.tif")
ESACCIbiomass1kmminusSD <- ESACCIbiomass1kmminusSD*0.47 #convert AGB to AGC
# 


# #make intact focal reference biomass map
# 

#intact forest based on mapbiomas forest cover removing SF, edges and disturbed forest pixels
intactforest1km <- stack("D:/RECCAP/Datasets/mapbiomas/intactforest_1km_mapbiomas_2010_2018.tif")[[2:9]]


fullyintact <- intactforest1km[[8]]>0.9 #90% threshold of intact forest cells

fullyintact[fullyintact==0] <- NA
# 
AGCmap <- ESACCIbiomass1km
# 

#perform filtering to generate reference biomass
AGCintact <- AGCmap*fullyintact
AGCintact[!is.finite(AGCintact)] <- NA
fwModel <- focalWeight(AGCintact, 1, type='circle')#circle of values to consider
fwModel[fwModel>0] <- 1
fwModel[fwModel==0] <- NA
neighbourhoodintactAGC <- focal(AGCintact,fwModel,na.rm=T,fun=median,NAonly=T,pad=T,padValue=NA,filename="D:/RECCAP/Datasets/ESACCIv2_2017amazonTot1km_JRCintactAdjRef.tif",overwrite=T)

#perform filtering to generate reference biomass plus SD
AGCmap <- ESACCIbiomass1kmplusSD
AGCintact <- AGCmap*fullyintact
AGCintact[!is.finite(AGCintact)] <- NA
fwModel <- focalWeight(AGCintact, 1, type='circle')#circle of values to consider
fwModel[fwModel>0] <- 1
fwModel[fwModel==0] <- NA
neighbourhoodintactAGC <- focal(AGCintact,fwModel,na.rm=T,fun=median,NAonly=T,pad=T,padValue=NA,filename="D:/RECCAP/Datasets/ESACCIv2_2017amazonTot1km_JRCintactAdjRefPlusSD.tif",overwrite=T)

#perform filtering to generate reference biomass minus SD

AGCmap <- ESACCIbiomass1kmminusSD
AGCintact <- AGCmap*fullyintact
AGCintact[!is.finite(AGCintact)] <- NA
fwModel <- focalWeight(AGCintact, 1, type='circle')#circle of values to consider
fwModel[fwModel>0] <- 1
fwModel[fwModel==0] <- NA
neighbourhoodintactAGC <- focal(AGCintact,fwModel,na.rm=T,fun=median,NAonly=T,pad=T,padValue=NA,filename="D:/RECCAP/Datasets/ESACCIv2_2017amazonTot1km_JRCintactAdjRefMinusSD.tif",overwrite=T)
