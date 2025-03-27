##--------------------------------------------------------------------------
## Climate change effects on deep-water coral communities in the NW Atlantic
## Script 6: Calculates interquantile range to evaluate model uncertainty
## Author: Maria Rakka
##--------------------------------------------------------------------------

## Libraries
library(terra)
library(tidyverse)
library(patchwork)


## Set directories 
setwd("..")

## Calculate interquantile range

all.files<-list.files(path="./Outputs/Spatial_files")

cwm_files<-all.files[grep("CWM_deltas",all.files)]
trait_div_files<-all.files[grep("div_deltas",all.files)]
sperc<-all.files[grep("Sperc",all.files)]

interquant_range<-function(myfile){
  myraster=rast(myfile)
  lower_quant=myraster[[grep("25",names(myraster))]]
  higher_quant=myraster[[grep("75",names(myraster))]]
  #med_quant=myraster[[grep("median",names(myraster))]]
  int_range=abs(higher_quant-lower_quant)
  mod_suf=case_when(str_detect(myfile,"ROMS")~"ROMS",
                      str_detect(myfile,"Bio|bO")~"BioOracle",
                      TRUE~"present")
  ind_suf=case_when(str_detect(myfile,"CWM")~"_cwm",
                    str_detect(myfile,"Trait")~"_trait_div",
                    TRUE~"_Sric")
  writeRaster(int_range,file=file.path(paste("Uncertainty/",mod_suf,
                                             ind_suf,
                                             "_iqrange.tif",sep="")),
              overwrite=TRUE)
  print("raster stored")
  #return((higher_quant-lower_quant))
}

interquant_range("BioOracle_trait_div_metrics_raster.tif")
lapply(cwm_files,interquant_range)
lapply(trait_div_files,interquant_range)
lapply(sperc,interquant_range)

## Evaluate relationship between IQR and depth

# Get depth raster
depth_ras<-rast("./Data/misc/depth_raster.tif")

unc_files<-list.files(path="Outputs/spatial_files/Uncertainty")
#[-grep("graphs",list.files())]
#cwm_files<-unc_files[grep("cwm",unc_files)]
#trait_files<-unc_files[grep("trait_div",unc_files)]

test1<-rast(x=file.path(paste("./Outputs/spatial_files/Uncertainty/",unc_files[3],sep="")))
depth_ras_masked<-mask(depth_ras,test1)

get_depth_graph<-function(x){
  modelname=case_when(str_detect(x,"ROMS")~"RCP8.5",
                      str_detect(x,"racle")~"SPS5.8.5",
                      TRUE~"Present")
  myras<-rast(file.path(paste("./Outputs/spatial_files/Uncertainty/",x,sep="")))
  if(any(str_detect(names(myras),"fric|feve|fdis"))){
    names(myras)<-c("Fric_75", "Feve_75", "Fdis_75", "frao_75")}else{
      names(myras)<-names(myras)
    }
  depth_vals=values(depth_ras_masked[[1]])
  test1_vals=values(myras)
  mydf<-as.data.frame(cbind(depth_vals,test1_vals)) %>% 
    drop_na()
  mygraph<-mydf%>%
    rename_at(vars(contains("_")),~word(., 1, sep = "_")) %>%
    dplyr::select(!any_of(contains("frao"))) %>%
    #relocate(any_of(c("fric","feve"))) %>% 
    pivot_longer(cols=!lyr1,names_to="index",values_to="val") %>%
    ggplot(aes(x=abs(lyr1),y=val))+
    geom_point()+
    facet_wrap(~index,ncol=1,scales="free_y")+
    theme_light()+xlab("Depth (m)")+
    theme(strip.background =element_rect(colour="gray",fill="white"),
          strip.text=element_text(colour="black"))+
    ggtitle(modelname)+ylab("Interquantile range")
  ggsave(plot=mygraph,filename=file.path(paste("./Outputs/graphs/",
                                               modelname,
                                               str_remove(x,".tif"),".jpg",sep="")),
         width = 5, height = 8)
  print("graph done")
}

lapply(unc_files[-c(2,5)],get_depth_graph)


## Slightly different graphs for Taxonomic richness (to adjust size of output
## graph)

get_depth_graph<-function(x){
  modelname=case_when(str_detect(x,"ROMS")~"RCP8.5",
                      str_detect(x,"racle")~"SPS5.8.5",
                      TRUE~"Present")
  myras<-rast(file.path(paste("./Outputs/spatial_files/Uncertainty/",x,sep="")))
  depth_vals=values(depth_ras_masked[[1]])
  test1_vals=values(myras)
  mydf<-as.data.frame(cbind(depth_vals,test1_vals)) %>% 
    drop_na()
  mygraph<-mydf%>%
    rename_at(vars(contains("_")),~word(., 1, sep = "_")) %>%
    dplyr::select(!any_of(contains("frao"))) %>%
    #relocate(any_of(c("fric","feve"))) %>% 
    pivot_longer(cols=!lyr1,names_to="index",values_to="val") %>% 
    ggplot(aes(x=abs(lyr1),y=val))+
    geom_point()+
    facet_wrap(~index,ncol=1,scales="free_y")+
    theme_light()+xlab("Depth (m)")+
    theme(strip.background =element_rect(colour="gray",fill="white"),
          strip.text=element_text(colour="black"))+
    ggtitle(modelname)+ylab("Interquantile range")
  ggsave(plot=mygraph,filename=file.path(paste("Ind_graphs/",
                                               modelname,
                                               str_remove(x,".tif"),".jpg",sep="")),
         width = 5, height = 5)
  print("graph done")
}

lapply(unc_files[c(2,5)],get_depth_graph)