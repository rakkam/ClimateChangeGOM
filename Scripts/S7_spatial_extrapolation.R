##--------------------------------------------------------------------------
## Climate change effects on deep-water coral communities in the NW Atlantic
## Script 7: Evaluates the extent of spatial extrapolation of the model
## Author: Maria Rakka
##--------------------------------------------------------------------------

## Libraries

library(tidyverse)
library(terra)
library(raster)
library(dsmextra)
library(magrittr)


## Get extrapolation data

load("./Data/joined_env_data_bioOracle.RData")
load("./Data/joined_env_data_ROMS2.RData")

dat_bio<-joined_env_dat_bio %>% 
  dplyr::select(temp_mean=temp_RCP85,
         sal_mean=sal_RCP85,
         vel_mean=vel_RCP85,mud,aspect2000,
         tpi20km,x,y,cell_nb,sam_ef) %>% 
  drop_na()

dat_ROMS<-joined_env_dat_ROMS %>% 
  dplyr::select(temp_mean=temp_RCP85,
                sal_mean=sal_RCP85,
                vel_mean=vel_RCP85,mud,aspect2000,
                tpi20km,x,y,cell_nb,sam_ef) %>% 
  drop_na()

## Get training data

load("./Data/models_thin_100_samples_250_chains_4.RData")

traindat<-cbind(model_envdat,model_expdes) %>% 
  dplyr::select(sal_mean,temp_mean,vel_mean,tpi20km,aspect2000,mud,x=Lon,y=Lat)

#Get a file for CRS

env_layer<-rast("./Data/misc/romsNWA_deltas2.tif")


## Extrapolation for Bio Oracle

# Set tibble options
options(tibble.width = Inf) # All tibble columns shown
options(pillar.neg = FALSE) # No colouring negative numbers
options(pillar.subtle = TRUE)
options(pillar.sigfig = 4)


# Set knitr options
knitr::opts_chunk$set(echo = TRUE)

# Define projected coordinate system
aftt_crs <- crs(env_layer)

# Define environmental covariates of interest
covariates.spermwhale <- c("sal_mean","temp_mean","vel_mean","tpi20km","aspect2000","mud")

testextrapolation <- compute_extrapolation(samples = traindat,
                                           covariate.names = covariates.spermwhale,
                                           prediction.grid = dat_bio,
                                           coordinate.system = aftt_crs)
summary(testextrapolation)

# Make raster map

# Add land on map
load("./Data/misc/shapefiles_for_basemap.RData")


dataformap<-testextrapolation$data$all %>% 
  mutate(mycat=case_when(ExDet<0~"Univariate",
                         ExDet>1 ~"Combinatorial",
                         ExDet>0 & ExDet<1 ~ "Analogous",
                         TRUE~NA)) %>% 
  dplyr::select(x,y,mycat)

ggplot()+
geom_raster(aes(x=x,y=y,fill=mycat),data=dataformap)+
  scale_fill_manual(values=c("#FED439FF", "#FD7446FF","#709AE1FF","#D2AF81FF" ))+
  geom_polygon(data =bs_land , aes(x = long, 
                                   y = lat,group=group),fill = "white",col="grey30")+
  coord_equal()+
  theme_minimal()+
  xlab("Lon")+ylab("Lat")+
  theme(legend.title=element_blank())+
  ggtitle("SP5.8.5")


## Extrapolation ROMS

testextrapolation <- compute_extrapolation(samples = traindat,
                                           covariate.names = covariates.spermwhale,
                                           prediction.grid = dat_ROMS,
                                           coordinate.system = aftt_crs)
summary(testextrapolation)

# Make raster map

dataformap<-testextrapolation$data$all %>% 
  mutate(mycat=case_when(ExDet<0~"Univariate",
                         ExDet>1 ~"Combinatorial",
                         ExDet>0 & ExDet<1 ~ "Analogous",
                         TRUE~NA)) %>% 
  dplyr::select(x,y,mycat)

ggplot()+
  geom_raster(aes(x=x,y=y,fill=mycat),data=dataformap)+
  scale_fill_manual(values=c("#FED439FF", "#FD7446FF","#709AE1FF","#D2AF81FF" ))+
  geom_polygon(data =bs_land , aes(x = long, 
                                   y = lat,group=group),fill = "white",col="grey30")+
  coord_equal()+
  theme_minimal()+
  xlab("Lon")+ylab("Lat")+
  theme(legend.title=element_blank())+
  ggtitle("RCP8.5")
