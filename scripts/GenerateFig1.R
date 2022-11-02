#Author: D. Fawcett
#
#Description:
#This code calculates trends in L-VOD AGC, compares land-cover fractions in 2018 and for areas of increasing and decreasing AGC and displays land cover fraction maps and breakdowns per country


#Outputs:
#-Figure 1
#-Figure S10


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


setwd("D:/RECCAP")

amazonBasin <- readOGR(dsn = './Datasets', layer = "amazonBasinCelso")

#read L-VOD time series
AGC_LVOD_amazoniasmoothmeancrop <- stack('./CodeOutputs/AGC_LVOD_amazoniasmoothmeancropV18.tif')
AGC_LVOD_amazoniasmoothmaxcrop <- stack('./CodeOutputs/AGC_LVOD_amazoniasmoothmaxcropV18.tif')
AGC_LVOD_amazoniatrendmeancrop <- stack('./CodeOutputs/AGC_LVOD_amazoniatrendmeancropV18.tif')

AGC_LVOD_amazoniamethodsmeancrop <- stack('./CodeOutputs/AGC_LVOD_amazoniamethodsmeancropV18.tif')

totalMask <- raster('./CodeOutputs/totalMaskV18.tif')

#Time series of mean L-VOD AGC

LVODAGCTSmethodsmean <- cellStats(AGC_LVOD_amazoniamethodsmeancrop*area(AGC_LVOD_amazoniamethodsmeancrop)*100,sum,na.rm=T)/10^9


###########
#load land cover datasets and resample to common resolution and extents (~1km and ~25km)

#mapbiomas forest fractional cover

mapbiomasforestfraccrop <- stack('./CodeOutputs/mapbiomasforestfraccropV16.tif')

#mapbiomas SF extent

mapbiomasSFfrac <- stack("D:/RECCAP/Datasets/mapbiomas/forestSFC2Frac_1km_2009_2018.tif")[[2:10]]

mapbiomasSFfracresamp <- raster::resample(mapbiomasSFfrac,AGC_LVOD_amazoniamethodsmeancrop[[1]])

mapbiomasSFfraccrop <- raster::mask(raster::mask(crop(mapbiomasSFfracresamp,amazonBasin),amazonBasin),totalMask)


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


#intact forest based on mapbiomas forest cover removing SF, edges and disturbed forest pixels
intactforest1km <- stack("D:/RECCAP/Datasets/mapbiomas/intactforest_1km_mapbiomas_2010_2018.tif")[[2:9]]

intactforest <- raster::resample(intactforest1km,AGC_LVOD_amazoniamethodsmeancrop[[1]])
intactforestcrop <- raster::mask(raster::mask(crop(intactforest,amazonBasin),amazonBasin),totalMask)


#############

#get 2018 state of landcover and the net change


mapbiomasdeforestationnet <- mapbiomasforestfraccrop[[8]]-mapbiomasforestfraccrop[[1]]

JRCdegradationcropsum <- JRCdegradationbasenoSFnoEdgecrop[[8]]#calc(bullockdegradationnotedgeresamp,sum)

mapbiomasSFfraccropmean <- mapbiomasSFfraccrop[[9]]#calc(mapbiomasSFfraccrop,mean)

IFfraccropmean <- intactforestcrop[[8]]#calc(intactforest,mean)

mapbiomasnonforestfraccropinitial <- abs(mapbiomasforestfraccrop[[1]]-1) #areas that were initially not forest

edgesdegradationcropmean <- forestedgesfraccrop[[8]]#mean(forestedgesfraccrop)#[[9]]#calc(forestedgesfraccrop,mean)


##############
#make figure showing different forest covers per cell
landcoverdisplaystack <- stack(JRCdegradationcropsum+edgesdegradationcropmean,mapbiomasSFfraccropmean,IFfraccropmean,1-(JRCdegradationcropsum+edgesdegradationcropmean+mapbiomasSFfraccropmean+IFfraccropmean))


names(landcoverdisplaystack) <- c("degraded forest cover","secondary forest cover","intact forest cover","non-forest cover")
raster.names <-  c("degraded forest cover","secondary forest cover","intact forest cover","non-forest cover")
levelplot(landcoverdisplaystack,names.attr=raster.names)

landcoverstack <- stack(JRCdegradationcropsum+edgesdegradationcropmean,mapbiomasSFfraccropmean,IFfraccropmean)


#make figure showing surface cover per country

#get amazon country shapes
amazonCountries <- readOGR(dsn = './Datasets/South_America', layer = "South_America")
#selectCountries <- c(1,2,4,6,8,9,11,14)
selectCountrynames <- c('Bolivia','Brazil','Colombia','Ecuador','French Guiana (France)','Guyana','Suriname','Peru','Venezuela')

countrymat <- matrix(data=NA, nrow=9, ncol=5)

countrycoverdf <- data.frame(name=character(0),degraded=numeric(0),sf=numeric(0),oldgrowth=numeric(0),nonforest=numeric(0))

for(i in 1:9){
countryshape <- amazonCountries[amazonCountries@data$COUNTRY == selectCountrynames[i],]  
             
extractedvalues <- extract(landcoverdisplaystack,countryshape,df=T, fun=sum,na.rm=T)
extractedrefpixels <- extract(landcoverdisplaystack,countryshape,df=T)
pixelnrs <- nrow(extractedrefpixels[complete.cases(extractedrefpixels),])

coverfractions <- extractedvalues[2:5]/pixelnrs
forestcover <- sum(coverfractions)
dftoadd <- data.frame(name=selectCountrynames[i],degraded=coverfractions[1],sf=coverfractions[2],oldgrowth=coverfractions[3],nonforest=coverfractions[4])
countrycoverdf <- rbind(countrycoverdf,dftoadd)

}

names(countrycoverdf) <-  c("country","degraded forest","secondary forest","old-growth forest","non-forest")

countrycoverdfmelted <- melt(countrycoverdf,id.vars='country',measure.vars=c("degraded forest","secondary forest","old-growth forest","non-forest"))
countrycoverdfmelted$variable<-  factor(countrycoverdfmelted$variable,levels = c('non-forest', "degraded forest",'secondary forest','old-growth forest'))

p1 <- ggplot(countrycoverdfmelted,aes(fill=variable,y=value,x=country)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = c('grey', "#997732",'#66CCEE','#4477AA'))+#scale_fill_viridis(discrete = T, option = "E") +
  ylim(c(0,1))+
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size=18))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab("grid-cells covered fraction 2018")+
  xlab("")+
  labs(fill='land cover type') 


plot(p1)



colnames(countrymat) <- c('name','degraded','secondary','old-growth','non-forest')
           
#############
#calculate trends



intactmask <- intactforestcrop[[8]]
intactmask[intactmask>0.9] <- 1
intactmask[intactmask<1] <- NA

source("./Code/divergePlotFuns.R")


fun1=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$estimates)}}
fun2=function(t) { if (!is.finite(sum(t))){ return(NA) } else { m = sens.slope(t); return(m$p.value) }}



AGC_trend=calc(AGC_LVOD_amazoniamethodsmeancrop[[1:9]], fun1)
AGC_trend_intact <- raster::mask(AGC_trend,intactmask)

AGC_trend_sig=calc(AGC_LVOD_amazoniamethodsmeancrop[[1:9]], fun2)
AGC_trend_intact_sig <- raster::mask(AGC_trend_sig,intactmask)

AGC_trend_intact_sig_mask <- AGC_trend_intact_sig
AGC_trend_intact_sig_mask[ AGC_trend_intact_sig<0.05] <- 1
AGC_trend_intact_sig_mask[ AGC_trend_intact_sig>=0.05] <- 0



############
#trends vs landcover types

plot(getValues(AGC_trend),getValues(mapbiomasdeforestationnet),pch=16,col = alpha('black', 0.1))

plot(getValues(AGC_trend),getValues(JRCdegradationcropsum),pch=16,col = alpha('black', 0.1))

plot(getValues(AGC_trend),getValues(mapbiomasSFfraccropmean),pch=16,col = alpha('black', 0.1))


sigloss <- (AGC_trend*(AGC_trend_sig<0.05))<(0)
siglossmask <- sigloss
siglossmask[siglossmask<1] <- NA

siglossN <- cellStats(sigloss,sum)

siggain <- (AGC_trend*(AGC_trend_sig<0.05))>(0)
siggainmask <- siggain
siggainmask[siggainmask<1] <- NA

sigainN <- cellStats(siggain,sum)

siglossgainraster <- siggain-sigloss
siglossgainrasterfordisp <- (totalMask*0)+siglossgainraster
siglossgainraster[siglossgainraster==0] <- NA

siglossgainmask <- siglossgainraster
siglossgainmask[!is.na(siglossgainraster)] <- 1
#siglossgainmaskbackground

diverge0(levelplot(siglossgainrasterfordisp,at=seq(-1, 1, len = 10),margin=FALSE,maxpixels = 2e10,par.settings=list(panel.background=list(col="lightgrey"))),colorRampPalette(c('red','white','blue')))+latticeExtra::layer(sp.polygons(amazonBasin,col='black'))#+layer(sp.polygons(brazilianAmazon,col='black'))

#display significant trends in LVOD AGC
 
trendminmax <- getMinMax(AGC_trend*siglossgainmask)
diverge0(levelplot(totalMask*0+AGC_trend*siglossgainmask,at=seq(trendminmax[1], trendminmax[2], len = 100),main='AGC trends 2011-2019',margin=FALSE,maxpixels = 2e10,par.settings=list(panel.background=list(col="lightgrey"))),colorRampPalette(c('red','white','blue')))+latticeExtra::layer(sp.polygons(amazonBasin,col='black'))#+layer(sp.polygons(brazilianAmazon,col='black'))
diverge0(levelplot(AGC_trend,at=seq(trendminmax[1], trendminmax[2], len = 100),main='AGC trends 2011-2019',margin=FALSE,maxpixels = 2e10,par.settings=list(panel.background=list(col="lightgrey"))),colorRampPalette(c('red','white','blue')))+latticeExtra::layer(sp.polygons(amazonBasin,col='black'))#+layer(sp.polygons(brazilianAmazon,col='black'))

sigtrendsfordispout <- totalMask*0+AGC_trend*siglossgainmask

#write all trends and significant trend rasters
writeRaster(sigtrendsfordispout,'./CodeOutputs/LVODAGCsigtrends2011_2019_V18.tif')
writeRaster(AGC_trend,'./CodeOutputs/LVODAGCalltrends2011_2019_V18.tif')

aresigtrends <- !is.na(sigtrendsfordispout)
aresigtrends[aresigtrends==0] <- NA
writeRaster(aresigtrends,'./CodeOutputs/LVODAGCsigtrendsmask2011_2019_V18.tif')

##########
#calculate and plot land cover in 2018 for increasing and decreasing cells

nonforestforuncertainty <- mapbiomasnonforestfraccropinitial#+SFincreasecompensation 

landcovergainlossstack <- stack(mapbiomasdeforestationnet,JRCdegradationcropsum,edgesdegradationcropmean,mapbiomasSFfraccropmean,IFfraccropmean)
siglossfracs <- cellStats(raster::mask(landcovergainlossstack,siglossmask),mean,na.rm=T)
siggainfracs <- cellStats(raster::mask(landcovergainlossstack,siggainmask),mean,na.rm=T)

#round very small values to 0 for visualisation
siggainfracs[abs(siggainfracs)<0.001] <- 0
siglossfracs[abs(siglossfracs)<0.001] <- 0


#calculate non-forest area (remainder)
nonforestsigloss <- 1-sum(siglossfracs[2:5])
nonforestsiggain <- 1-sum(siggainfracs[2:5])

#combined bars
lctype <- c('net forest change (2011-2018)','degraded forest (non-edge)','degraded forest (edge)','secondary forest','old-growth forest','non-forest')#,'secondary forest') 


comblossfracs <- data.frame(lctype,cbind(c(siglossfracs[1:5],nonforestsigloss),c(siggainfracs[1:5],nonforestsiggain)))



  comblossfracscurrent <- data.frame(lctype=rep(comblossfracs$lctype,2),subsetname=c(rep("loss" , 6) , rep("gain" , 6)),value=unlist(c(comblossfracs[,2:3])))
  #comblossfracscurrent$lctype <-  factor( comblossfracscurrent$lctype,levels = c("deforestation (2011-2018)", "degraded (non-edge)",'degraded (edge)','secondary forest','intact forest','uncertainty'))
  comblossfracscurrent$lctype <-  factor( comblossfracscurrent$lctype,levels = c('non-forest', "degraded forest (non-edge)",'degraded forest (edge)','secondary forest','old-growth forest',"net forest change (2011-2018)"))
  
  p1 <- ggplot(comblossfracscurrent,aes(fill=lctype,y=value,x=subsetname)) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(values = c('grey', "#AA3377",'#CCBB44','#66CCEE','#4477AA',"#EE6677"))+#scale_fill_viridis(discrete = T, option = "E") +
    ylim(c(-0.05,1))+
    geom_hline(yintercept=0,lty=2)+
    theme_bw() +
    theme_classic() +
    theme(text = element_text(size=18))+
    ylab("land cover fraction 2018")+
    xlab("")
  
  
plot(p1)


#calculate fraction of disturbed forest in 2018

disturbedfrac2018=sum(comblossfracscurrent[2:4,3])/sum(comblossfracscurrent[2:5,3])

#calculate reduction of forest area since 2011 in loss cells
percreduction <- abs(comblossfracs$X1[1])/sum(abs(comblossfracs$X1))*100


