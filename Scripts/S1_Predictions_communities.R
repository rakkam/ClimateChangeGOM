##--------------------------------------------------------------------------
## Climate change effects on deep-water coral communities in the NW Atlantic
## Script 1: Uses HMSC model to predict community composition under different
##           climate scenarios
## Author: Maria Rakka
##--------------------------------------------------------------------------


##----------------------------------------------------------------------------
## Part 1: Predict species communities under 2 climate scenarios 
##----------------------------------------------------------------------------


## Libraries

library(Hmsc)
library(tidyverse)
library(terra)

## Set directories
setwd("...")

## Get and prepare environmental data

load("./Data/joined_env_data_bioOracle.RData")
load("./Data/joined_env_data_ROMS2.RData")
load("./joined_env_data_present.RData")

joined_env_dat_bio<-joined_env_dat_bio %>% 
  dplyr::select(contains("RCP85"),mud,aspect2000,tpi20km,x,y,sam_ef,cell_nb)

joined_env_dat_ROMS<-joined_env_dat_ROMS %>% 
  dplyr::select(contains("RCP85"),mud,aspect2000,tpi20km,x,y,sam_ef,cell_nb)

## Get a random raster to copy crs

env_layer<-rast("./Data/misc/romsNWA_deltas2.tif")

## Load model and make predictions

load("./Data/models_thin_100_samples_250_chains_4.RData")

set.seed(1)

# Make function to predict species composition
trial_fun<-function(scen,x,...){
  trial_env1<-x %>% 
    dplyr::select(...,aspect2000,mud, tpi20km,sam_ef) %>% 
    rename(temp_mean=contains("temp"),
           sal_mean=contains("sal"),
           vel_mean=contains("vel")) %>% 
    as.data.frame()
  
  xydat<-x %>% 
    dplyr::select(x,y) %>%
    as.matrix()
  
  gradient_trial=prepareGradient(models$pres.abs.quad,
                                 XDataNew=trial_env1,
                                 sDataNew=list("ID"=xydat))
  predict_com=predict(models$pres.abs.quad,Gradient=gradient_trial,predictEtaMean = TRUE)
  
  print("finished predictions")
  saveRDS(predict_com,file=paste("./Outputs/",scen,"_preds_full1000.RData",sep=""))
  print("saved predictions")
  EpredY=Reduce("+",predict_com)/length(predict_com)
  getS = function(p){return(rowSums(p))}
  aS = simplify2array(lapply(X = predict_com, FUN = getS))
  
  Srich = apply(aS, 1, mean)
  sdS = sqrt(apply(aS, 1, var))
  S5 = apply(aS > 5, 1, mean)
  S10 = apply(aS > 10, 1, mean)
  sp_ric_dat<-data.frame(x=x$x, y=x$y, Srich,sdS,S5,S10) 
  sp_ric_ras <- rast(sp_ric_dat,extent=ext(env_layer),type="xyz")
  crs(sp_ric_ras) <- as.character(crs(env_layer))
  saveRDS(sp_ric_ras,file=paste("./Outputs/spatial_files/",scen,"_metrics_raster.RData",sep=""))
  print("saved spatial metrics")
  my_metrics<-data.frame(Srich,sdS,S5,S10,cell_nb=x$cell_nb)
  saveRDS(my_metrics,file=paste("./Outputs/",scen,"_metrics_table.RData",sep=""))
}

# Use function to predict species communities (takes long time to run and needs
# a lot of memory)

trial_fun(scen="bioOracle",x=joined_env_dat_bio,temp_RCP85,sal_RCP85,vel_RCP85)
trial_fun(scen="ROMS_NWA2",x=joined_env_dat_ROMS,temp_RCP85,sal_RCP85,vel_RCP85)
trial_fun(scen="present",x=analogue_dat,temp_mean,sal_mean,vel_mean)

# Read/Check predictions
trial_ras<-readRDS("./Outputs/spatial_files/ROMS_NWA2_metrics_raster.RData")
plot(trial_ras$Srich)

##----------------------------------------------------------------------------
## Part 2: Calculate frequency of occurrence for each scenario and prepare maps
## of delta values
##----------------------------------------------------------------------------

##Libraries

library(paletteer)
library(patchwork)
library(purrr) 

## Read lists of predictions and calculate frequency of occurrence
## Note that full lists are large files and heavy to process

my_giant_list_ROMS<-readRDS("./Outputs/ROMS_NWA2_preds_full1000.RData")
EpredY_ROMS=Reduce("+",my_giant_list_ROMS)/length(my_giant_list_ROMS)
saveRDS(EpredY_ROMS,file="./Outputs/freq_occ_ROMS_NWA2.RData")
rm(my_giant_list_ROMS)

my_giant_list_Bio<-readRDS("./Outputs/BioOracle_preds_full1000.RData")
EpredY_BioOracle=Reduce("+",my_giant_list_Bio)/length(my_giant_list_Bio)
saveRDS(EpredY_BioOracle,file="./Outputs/freq_occ_BioOracle.RData")
rm(my_giant_list_Bio)

my_giant_list_pres<-readRDS("./Outputs/present_preds_full1000.RData")
EpredY_pres=Reduce("+",my_giant_list_pres)/length(my_giant_list_pres)
saveRDS(EpredY_pres,file="./Outputs/freq_occ_present.RData")

## Make maps for each species

# Write function to make maps
make_freq_plot<-function(species_name,mydat){
  ggplot(data = mydat)+ 
    geom_raster(aes(x = x, y = y, fill = get(species_name))) +
    scale_fill_paletteer_c("pals::coolwarm",direction=-1,limits = c(-0.7,0.7))+
    theme_void() +
    theme(legend.position = "bottom",
          legend.title=element_blank(),
          panel.border = element_rect(colour = "black",fill=NA),
          plot.margin = unit(c(0,30,0,0), "pt"))+
    coord_equal()+
    ggtitle(species_name)
}

# Make Maps

deltafreq_Bio<-EpredY_BioOracle-EpredY_pres

deltafreq_Bio_tab<-data.frame(deltafreq_Bio) %>% 
  mutate(x=joined_env_dat_bio$x,y=joined_env_dat_bio$y)

myspecies<-names(deltafreq_Bio_tab)[1:30]
allplots<-lapply(myspecies,make_freq_plot,mydat=deltafreq_Bio_tab)

pdf("./Outputs/graphs/SupMat-BioOracle_SpeciesDelta_FreqOc.pdf")
reduce(allplots[1:15], `+`)+
  plot_layout(ncol = 3,guides = 'collect')&
  theme(legend.position = 'bottom')
reduce(allplots[16:30], `+`)+
  plot_layout(ncol = 3,guides = 'collect')&
  theme(legend.position = 'bottom')
dev.off()


deltafreq_ROMS<-EpredY_ROMS-EpredY_present

deltafreq_ROMS_tab<-data.frame(deltafreq_ROMS) %>% 
  mutate(x=joined_env_dat_bio$x,y=joined_env_dat_bio$y)

myspecies<-names(deltafreq_ROMS_tab)[1:30]
ROMSplots<-lapply(myspecies,make_freq_plot,mydat=deltafreq_ROMS_tab)

pdf("./Outputs/graphs/SupMat-ROMS2_SpeciesDelta_FreqOc.pdf")
reduce(ROMSplots[1:15], `+`)+
  plot_layout(ncol = 3,guides = 'collect')&
  theme(legend.position = 'bottom')
reduce(ROMSplots[16:30], `+`)+
  plot_layout(ncol = 3,guides = 'collect')&
  theme(legend.position = 'bottom')
dev.off()