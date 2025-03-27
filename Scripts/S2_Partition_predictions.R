##--------------------------------------------------------------------------
## Climate change effects on deep-water coral communities in the NW Atlantic
## Script 2: Partitions every list of predictions into smaller chunks to 
##           facilitate the estimate of indices of change
## Author: Maria Rakka
##--------------------------------------------------------------------------


## Libraries
library(tidyverse)

## Make result reproducible
set.seed(1)

## Set working directories
setwd("...")

## Partition predictions in two files with 300 samples each

# Make function to split lists
split_my_list<-function(huge_list,myname){
  
  sel_elems<-sample(1:length(huge_list),300)
  huge_list_100<-huge_list[sel_elems]
  rm(huge_list)
  magic_num<-(dim(huge_list_100[[1]])[1])/2
  
  magic_vec1<-1:floor(magic_num)
  magic_vec2<-(max(magic_vec1)+1):dim(huge_list_100[[1]])[1]
  
  my_huge_list_part1<-map(huge_list_100,~.[magic_vec1,])
  my_huge_list_part2<-map(huge_list_100,~.[magic_vec2,])

  saveRDS(my_huge_list_part1,
          file=file.path(paste("./Outputs/partitions/",myname,
                               "_preds_full300_part1.RData",sep="")))
  saveRDS(my_huge_list_part2,
          file=file.path(paste("./Outputs/partitions/",myname,
                               "_preds_full300_part2.RData",sep="")))
  }


# Read and split list for present conditions 
huge_list_present<-readRDS("./Outputs/present_preds_full1000.RData")

split_my_list(huge_list_present, myname="Present")
rm(huge_list_present)

# Read and split list for climate scenario BioOracle

my_giant_list_Bio<-readRDS("bioOracle_preds_full1000.RData")

split_my_list(my_giant_list_Bio, myname="BioOracle")
rm(my_giant_list_Bio)

# Read and split list for climate scenario ROMS-NWA
my_giant_list_NWA<-readRDS("ROMS_NWA2_preds_full1000.RData")

split_my_list(my_giant_list_NWA, myname="ROMS_NWA2")
rm(my_giant_list_NWA)