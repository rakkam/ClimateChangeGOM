##--------------------------------------------------------------------------##--------------------------------------------------------------------------
## Climate change effects on deep-water coral communities in the NW Atlantic
## Script 5: Calculates delta values of all indices of change
## Author: Maria Rakka
##--------------------------------------------------------------------------

## Libraries
library(terra)
library(tidyverse)

## Set directories
setwd("...")

#Get a file that has x, y

load("./Data/joined_env_data_bioOracle.RData")

# Get a random raster to copy crs

slayer<-rast("./Data/misc/romsNWA_deltas2.tif")

##-----------------------------------------------------------------------------
# Taxonomic richness extract median and 25,75 percentiles for each dataset
# and estimate deltas
##----------------------------------------------------------------------------

## Create function and run it to estimate median, 25, 75 percentiles
calc_ric<-function(myname){
  myfile1=readRDS(file=file.path(paste("./Outputs/",myname,
                                       "_preds_full1000.RData",sep="")))
 rich_trial=lapply(myfile1,rowSums)
 rich_matrix=do.call(rbind, rich_trial)
 rich_med=apply(rich_matrix, 2, median)
 rich_25=apply(rich_matrix,2,quantile,probs=0.25,na.rm=TRUE)
 rich_75=apply(rich_matrix,2,quantile,probs=0.75,na.rm=TRUE)
return(cbind(rich_med,rich_25, rich_75))}


rich_pres<-calc_ric("present")
rich_bio<-calc_ric("BioOracle")
rich_ROMS<-calc_ric("ROMS_NWA2")

## Calculate delta values and turn result into raster

bOracle_delta_srich<-(rich_bio-rich_pres)*100/rich_pres
ROMS_delta_srich<-(rich_ROMS-rich_pres)*100/rich_pres

bOracle_delta_srich_tab<-data.frame(bOracle_delta_srich) %>% 
  mutate(x=joined_env_dat_bio$x,y=joined_env_dat_bio$y) %>% 
  relocate(x,y)
bOracle_delta_srich_rast<-rast(bOracle_delta_srich_tab,
                               extent=ext(slayer),type="xyz")

crs(bOracle_delta_srich_rast) <- as.character(crs(slayer))

ROMS_delta_srich_tab<-data.frame(ROMS_delta_srich) %>% 
  mutate(x=joined_env_dat_bio$x,y=joined_env_dat_bio$y) %>% 
  relocate(x,y)
ROMS_delta_srich_rast<-rast(ROMS_delta_srich_tab,
                               extent=ext(slayer),type="xyz")

crs(ROMS_delta_srich_rast) <- as.character(crs(slayer))

writeRaster(bOracle_delta_srich_rast,file="./Outputs/spatial_files/bOracle_Sperc.tif",overwrite=TRUE)
writeRaster(ROMS_delta_srich_rast,file="./Outputs/spatial_files/ROMS_Sperc.tif",overwrite=TRUE)

## Check against depth

# Get depth raster
depth_ras<-rast("./Data/misc/depth_raster.tif")

mydepth<-extract(depth_ras,cbind(joined_env_dat_bio$x,joined_env_dat_bio$y))

bOracle_delta_srich_tab %>% 
  mutate(depth=mydepth$depth) %>% 
  ggplot(aes(y=rich_med,x=abs(depth)))+
  geom_point()+
  theme_light()+
  ggtitle("BioOracle")+ylab("Srich change (%)")

ROMS_delta_srich_tab %>% 
  mutate(depth=mydepth$depth) %>% 
  ggplot(aes(y=rich_med,x=abs(depth)))+
  geom_point()+
  theme_light()+
  ggtitle("ROMS")+ylab("Srich change (%)")

##--------------------------------------------------------------
# Calculate bray-curtis dissimilarity index
##--------------------------------------------------------------

## Libraries
library(vegan)
library(proxy)

## Calculate bray-curtis and perform cell-wise comparison (present vs climate
## scenarios)

freq_oc_pres<-readRDS("./Outputs/freq_occ_present.RData")
freq_oc_BioOracle<-readRDS("./Outputs/freq_occ_BioOracle.RData")
freq_oc_ROMS<-readRDS("./Outputs/freq_occ_ROMS_NWA2.RData")

CellWise_comp<-function(n,x){
  trial=rbind(freq_oc_pres[n,],x[n,])
  return(vegdist(trial,method="bray"))
}

diss_pres_BioOracle<-sapply(1:nrow(freq_oc_pres),CellWise_comp,x=freq_oc_BioOracle)
diss_pres_ROMS<-sapply(1:nrow(freq_oc_pres),CellWise_comp,x=freq_oc_ROMS)

# Convert results to spatial raster files
diss_pres_tab<-data.frame(x=joined_env_dat_bio$x, 
                          y=joined_env_dat_bio$y,
                          diss_pres_BioOracle=diss_pres_BioOracle,
                          diss_pres_ROMS=diss_pres_ROMS)

diss_pres_rast <- rast(diss_pres_tab,extent=ext(slayer),type="xyz")
crs(diss_pres_rast) <- as.character(crs(slayer))
plot(diss_pres_rast)

writeRaster(diss_pres_rast,
            file="./Outputs/spatial_files/dissimilarity_indexFreqOcc.tif",
            overwrite=TRUE)


##--------------------------------------------------------------
# Calculate delta values for trait diversity indices and CWM
##--------------------------------------------------------------

# Trait diversity indices

Bio_metrics<-rast("BioOracle_trait_div_metrics_raster.tif")
ROMS_metrics<-rast("ROMS_NWA2_trait_div_metrics_raster.tif")
Pres_metrics<-rast("present_trait_div_metrics_raster.tif")

Bio_delta_ras<-Bio_metrics-Pres_metrics
writeRaster(Bio_delta_ras,file="Bio_Trait_div_deltas.tif",overwrite=TRUE)
ROMS_delta_ras<-ROMS_metrics-Pres_metrics
writeRaster(ROMS_delta_ras,file="ROMS2_Trait_div_deltas.tif",overwrite=TRUE)
plot(ROMS_delta_ras)

# Community weighted mean

present_cwm<-rast("present_cwm_raster.tif")
bio_cwm<-rast("BioOracle_cwm_raster.tif")
ROMS_cwm<-rast("ROMS_NWA2_cwm_raster.tif")

Bio_delta_ras<-bio_cwm-present_cwm
writeRaster(Bio_delta_ras,file="Bio_CWM_deltas.tif",overwrite=TRUE)
ROMS_delta_ras<-ROMS_cwm-present_cwm
writeRaster(ROMS_delta_ras,file="ROMS2_CWM_deltas.tif",overwrite=TRUE)
