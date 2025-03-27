##--------------------------------------------------------------------------
## Climate change effects on deep-water coral communities in the NW Atlantic
## Script 4: Uses predictions from HMSC model to estimate metrics of change
##           under different climate scenarios
## Attention: This script requires a lot of memory and can only run in 
##            facilities of High Performance Computing (takes approximately
##            4h for 8 clusters with 65G memory each)
## Author: Maria Rakka
##--------------------------------------------------------------------------

## Libraries
  
library(fundiversity)
library(tidyverse)
library(doParallel)
library(future)
library(future.apply)
library(abind)

## Prepare directories
setwd("..")
localDir="."
outputDir=file.path(localDir,"Outputs")

## Read in traits from model
load(file="traits_trial.RData")

my_giant_list<-readRDS("./Data/ROMS_NWA_preds_full1000.RData")


## Trait indices are computed in parallel
mytrait_fun<-function(x){
  fric<-fd_fric(sample_traits,x, stand = TRUE) 
  feve<-fd_feve(sample_traits,x)
  fdis<-fd_fdis(sample_traits,x)
  #return(as.matrix(cbind(sric=sp_rich,
  #                       fric=fric$FRic,
  #                       feve=feve$FEve,
  #                       fdis=fdis$FDis)))
  return(cbind(fric=fric$FRic,feve=feve$FEve,fdis=fdis$FDis))
}

## Using 8 clusters
plan(multisession, workers=8)

trial_indexes<-future_lapply(my_giant_list,mytrait_fun)
save(trial_indexes,file=file.path(outputDir,"all_trait_div_indices.RData"))

a <- do.call(abind, c(trial_indexes,list(along=3)))
trait_median<-apply(a, 1:2, median,na.rm=TRUE)
trait_25<-apply(a,1:2,quantile,probs=0.25,na.rm=TRUE)
trait_75<-apply(a,1:2,quantile,probs=0.75,na.rm=TRUE)

final_res<-cbind(trait_median,trait_25,trait_75)

colnames(final_res)<-c("fric_median",
                       "feve_median",
                       "fdis_median",
                       "fric_25",
                       "feve_25",
                       "fdis_25",
                       "fric_75",
                       "feve_75",
                       "fdis_75")

## CWM needs matrix multiplication and is not parallelized

myCWM<-function(x){
  sp_rich<-rowSums(x)
  return((x %*% as.matrix(sample_traits))/
           matrix(rep(sp_rich, dim(sample_traits)[2]), ncol = dim(sample_traits)[2]))
}

plan(multisession, workers=8)
CWM_indexes<-future_lapply(my_giant_list,myCWM)

cwm_array <- do.call(abind, c(CWM_indexes,list(along=3)))
cwm_median<-apply(cwm_array, 1:2, median,na.rm=TRUE)
cwm_25<-apply(cwm_array,1:2,quantile,probs=0.25,na.rm=TRUE)
cwm_75<-apply(cwm_array,1:2,quantile,probs=0.75,na.rm=TRUE)

final_cwm<-cbind(cwm_median,cwm_25,cwm_75)

colnames(final_cwm)<-c("PC1_median",
                       "PC2_median",
                       "PC3_median",
                       "PC1_25",
                       "PC2_25",
                       "PC3_25",
                       "PC1_75",
                       "PC2_75",
                       "PC3_75")
final_list<-list(final_res,final_cwm)

save(final_list,
     file=file.path(outputDir,"final_trait_indices.RData"))