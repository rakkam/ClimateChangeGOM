##--------------------------------------------------------------------------##--------------------------------------------------------------------------
## Climate change effects on deep-water coral communities in the NW Atlantic
## Script 5: Packs the outputs of script 4 into spatial files
## Author: Maria Rakka
##--------------------------------------------------------------------------


## Libraries

library(purrr)
library(tidyverse)
library(terra)
library(abind)


## Set directories
setwd("...")

#Get a file that has x, y

load("./Data/joined_env_data_bioOracle.RData")

# Get a random raster to copy crs

env_layer<-rast("./Data/misc/romsNWA_deltas2.tif")

# function to load, merge metrics and create spatial files from them

make_spat_files<-function(myfile_part1,myfile_part2,myname){
  load(file.path(paste("./Outputs/Beluga_outputs/Future_indices/",myfile_part2,sep="")))
  final_list2<-final_list
  load(file.path(paste("./Outputs/Beluga_outputs/Future_indices/",myfile_part1,sep="")))

  #Indices
  trait_div_indices<-as.data.frame(abind(final_list[[1]],final_list2[[1]],along=1))
  trait_div_indices$x<-joined_env_dat_bio$x
  trait_div_indices$y<-joined_env_dat_bio$y
  trait_div_indices<-trait_div_indices %>% 
  relocate(x,y)

 trait_div_ras <- rast(trait_div_indices,extent=ext(env_layer),type="xyz") 
 crs(trait_div_ras) <- as.character(crs(env_layer))
 plot(trait_div_ras$feve_median)

 saveRDS(trait_div_ras,file=file.path(paste("./Outputs/spatial_files/",
                                  myname,"_trait_div_metrics_raster.RData",sep="")))
 writeRaster(trait_div_ras,
             file=file.path(paste("./Outputs/spatial_files/metrics/",
                                     myname,"_trait_div_metrics_raster.tif",sep="")),overwrite=TRUE)
 print("estimated and saved trait indices")

 #CWM
 trait_cwm<-as.data.frame(abind(final_list[[2]],final_list2[[2]],along=1))
 trait_cwm$x<-joined_env_dat_bio$x
 trait_cwm$y<-joined_env_dat_bio$y
 trait_cwm<-trait_cwm %>% 
  relocate(x,y)
 
 trait_cwm_ras <- rast(trait_cwm,extent=ext(env_layer),type="xyz")
 crs(trait_cwm_ras) <- as.character(crs(env_layer))

 #setwd("C:/Users/maria/Projects/OFI/Climate_change_v2/Data/spatial_files")
 saveRDS(trait_cwm_ras,file=file.path(paste("./Outputs/spatial_files/",
                                            myname,"_cwm_raster.RData",sep="")))
 writeRaster(trait_cwm_ras,file=file.path(paste("./Outputs/spatial_files/metrics/",
                                     myname,"_cwm_raster.tif",sep="")),overwrite=TRUE)
 print("estimated and saved trait cwm")          
}

make_spat_files(myfile_part1="BioOracle_final_trait_indices_part1.RData",
                myfile_part2="BioOracle_final_trait_indices_part2.RData",myname="BioOracle")

make_spat_files(myfile_part1="ROMS_NWA2_final_trait_indices_part1.RData",
                myfile_part2="ROMS_NWA2_final_trait_indices_part2.RData",myname="ROMS_NWA2")


## Slightly changing function to work with present indices


make_spat_files<-function(myfile_part1,myfile_part2,myname){
  load(file.path(paste("./Outputs/Beluga_outputs/Present_indices/",myfile_part2,sep="")))
  final_list2<-final_list
  load(file.path(paste("./Outputs/Beluga_outputs/Present_indices/",myfile_part1,sep="")))
  
  #Indices
  trait_div_indices<-as.data.frame(abind(final_list[[1]],final_list2[[1]],along=1))
  trait_div_indices$x<-joined_env_dat_bio$x
  trait_div_indices$y<-joined_env_dat_bio$y
  trait_div_indices<-trait_div_indices %>% 
    relocate(x,y)
  
  trait_div_ras <- rast(trait_div_indices,extent=ext(env_layer),type="xyz") 
  crs(trait_div_ras) <- as.character(crs(env_layer))
  plot(trait_div_ras$feve_median)
  
  #setwd("C:/Users/maria/Projects/OFI/Climate_change_v2/Data/spatial_files")
  saveRDS(trait_div_ras,file=file.path(paste("./Outputs/spatial_files/",
                                             myname,"_trait_div_metrics_raster.RData",sep="")))
  writeRaster(trait_div_ras,
              file=file.path(paste("./Outputs/spatial_files/metrics/",
                                   myname,"_trait_div_metrics_raster.tif",sep="")),overwrite=TRUE)
  print("estimated and saved trait indices")
  
  #CWM
  trait_cwm<-as.data.frame(abind(final_list[[2]],final_list2[[2]],along=1))
  trait_cwm$x<-joined_env_dat_bio$x
  trait_cwm$y<-joined_env_dat_bio$y
  trait_cwm<-trait_cwm %>% 
    relocate(x,y)
  
  trait_cwm_ras <- rast(trait_cwm,extent=ext(env_layer),type="xyz")
  crs(trait_cwm_ras) <- as.character(crs(env_layer))
  
  #setwd("C:/Users/maria/Projects/OFI/Climate_change_v2/Data/spatial_files")
  saveRDS(trait_cwm_ras,file=file.path(paste("./Outputs/spatial_files/",
                                             myname,"_cwm_raster.RData",sep="")))
  writeRaster(trait_cwm_ras,file=file.path(paste("./Outputs/spatial_files/metrics/",
                                                 myname,"_cwm_raster.tif",sep="")),overwrite=TRUE)
  print("estimated and saved trait cwm")          
}

make_spat_files(myfile_part1="present_final_trait_indices_part1.RData",
                myfile_part2="present_final_trait_indices_part2.RData",myname="present")


