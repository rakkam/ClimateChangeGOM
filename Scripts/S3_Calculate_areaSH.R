##--------------------------------------------------------------------------
## Climate change effects on deep-water coral communities in the NW Atlantic
## Script 3: Calculates metrics of changes in the suitable habitat for each
##           coral genus
## Author: Maria Rakka
##--------------------------------------------------------------------------

## libraries
library(tidyverse)
library(terra)
library(ggpubr)


## Set directories
setwd("...")

## Calculate the total area of suitable habitat for each coral genus

# Function to calculate the median and quartiles
calc_area_fun<-function(myfile){
  myfile1<-readRDS(file.path(paste("./Outputs/partitions/",myfile,"_part1.RData",sep="")))
  myfile2<-readRDS(file.path(paste("./Outputs/partitions/",myfile,"_part2.RData",sep="")))
  
  mymatrix1<-do.call(rbind,myfile1 %>% 
                       map(.,colSums)) 
  mymatrix2<-do.call(rbind,myfile2 %>% 
                       map(.,colSums)) 
  
  mymat<-mymatrix1+mymatrix2
  
  sumdat<-mymat%>%
    as.data.frame() %>% 
    summarize_all(list(mymed=median,
                       q25=~quantile(.x,0.25),
                       q75=~quantile(.x,0.75))) %>% 
    pivot_longer(cols=everything(),names_sep="_",names_to=c("species","var")) %>% 
    pivot_wider(id_cols="species",names_from="var",values_from="value")
  
  return(list(mymat,sumdat))
}

present<-calc_area_fun("present_preds_full300")
NWA_ROMS<-calc_area_fun("ROMS_NWA2_preds_full300")
BioOracle<-calc_area_fun("BioOracle_preds_full300")

figA<-present[[2]]%>% 
  ggplot(aes(y=reorder(species,mymed),x=mymed))+
  geom_point(size=0.9)+
  geom_errorbarh(aes(xmin=q25,xmax=q75),alpha=0.7,linewidth=0.4)+
  geom_point(data=BioOracle[[2]],aes(y=species,x=mymed),colour="darkorange",
             position=position_nudge(x = 0, y = 0.1),size=0.9)+
  geom_errorbarh(data=BioOracle[[2]],
                 aes(xmin=q25,xmax=q75),colour="darkorange",
                 position=position_nudge(x = 0, y = 0.1),alpha=0.7,linewidth=0.4)+
  theme_light()+
  theme(axis.title  = element_blank(),
        axis.text.y=element_text(face="italic"))+
  geom_point(y=10.5,x=57000,colour="darkorange",size=1)+
  annotate(geom="text",y=10.5,x=67000,label="SPS5.8.5",colour="gray43",size=3)+
  geom_point(y=12,x=57000,colour="black",size=1)+
  annotate(geom="text",y=12,x=67000,label="Present",colour="gray43",size=3)

figB<-present[[2]]%>% 
  ggplot(aes(y=reorder(species,mymed),x=mymed))+
  geom_point(size=0.9)+
  geom_errorbarh(aes(xmin=q25,xmax=q75),alpha=0.7,linewidth=0.4)+
  geom_point(data=NWA_ROMS[[2]],aes(y=species,x=mymed),colour="darkorange",
             position=position_nudge(x = 0, y = 0.1),size=0.9)+
  geom_errorbarh(data=NWA_ROMS[[2]],
                 aes(xmin=q25,xmax=q75),colour="darkorange",
                 position=position_nudge(x = 0, y = 0.1),alpha=0.7,linewidth=0.4)+
  theme_light()+
  theme(axis.title.y = element_blank(),
        axis.text.y=element_text(face="italic"))+
  labs(x=bquote('Area of suitable habitat' ~(km^2)))+
  geom_point(y=10.5,x=57000,colour="darkorange",size=1)+
  annotate(geom="text",y=10.5,x=67000,label="RCP8.5",colour="gray43",size=3)+
  geom_point(y=12,x=57000,colour="black",size=1)+
  annotate(geom="text",y=12,x=67000,label="Present",colour="gray43",size=3)


pdf("./Outputs/graphs/Fig5.pdf",width=5,height=8)
ggarrange(figA,figB, ncol=1)
dev.off()

## Calculate the average depth

# Get environmental variables for cell numbers
load("./Data/joined_env_data_present.RData")

# Get depth raster
depth_ras<-rast("Data/misc/depth_raster.tif")

# Get a vector with cells that have presences in each dataset and each species
# rownames have been maintained in the two datasets


myfile="present_preds_full300"
myfile1<-readRDS(file.path(paste(myfile,"_part1.RData",sep="")))
myfile2<-readRDS(file.path(paste(myfile,"_part2.RData",sep="")))


give_me_rownames<-function(my_matrix,mycol){
return(rownames(my_matrix)[my_matrix[, mycol] == 1])
}

give_me_info<-function(myfile,species_nb){
  myfile1<-readRDS(file.path(paste(myfile,"_part1.RData",sep="")))
  myfile2<-readRDS(file.path(paste(myfile,"_part2.RData",sep="")))
  
  firstbit<-lapply(myfile1,FUN=give_me_rownames,mycol=species_nb)
  secondbit<-lapply(myfile2,FUN=give_me_rownames,mycol=species_nb)
  allbits<-map2(firstbit,secondbit,c)%>% 
    map(.,as.numeric)
  
  mytable<-lapply(1:100,function(x){
    mypres<-analogue_dat[allbits[[x]],]
    return(mypres$cell_nb)})
  
  mydepths<-lapply(mytable,terra::extract,x=depth_ras,xy=TRUE)
  names(mydepths)<-paste("dataset",1:100,sep="_")
  return(map_df(mydepths, ~ .x |> select(x,y,depth),.id="mydataset"))
  
}

# Get general graph with median and quantiles among 1000 datasets for each
# scenario

present_latlon_depth<-lapply(1:30,give_me_info,myfile="present_preds_full300")
names(present_latlon_depth)<-colnames(myfile1[[1]])

present_sums<-map_df(present_latlon_depth, ~.x %>% 
  mutate(depth=abs(depth)) %>%
  group_by(mydataset) %>% 
  summarize_at(vars(x,y,depth),mean) %>% 
  ungroup() %>% 
  #summarize_at(vars(x,y,depth),list(mymean=mean,mysd=sd)),.id="Genus")
    summarize_at(vars(x,y,depth),list(mymed=median,
                                      q25=~quantile(.x,0.25),
                                      q75=~quantile(.x,0.75))),.id="Genus")


BioOracle_latlon_depth<-lapply(1:30,give_me_info,myfile="BioOracle_preds_full300")
names(BioOracle_latlon_depth)<-colnames(myfile1[[1]])

BioOracle_sums<-map_df(BioOracle_latlon_depth, ~.x%>% 
                       mutate(depth=abs(depth)) %>%
                       group_by(mydataset) %>% 
                       summarize_at(vars(x,y,depth),mean) %>% 
                       ungroup() %>% 
                       summarize_at(vars(x,y,depth),list(mymed=median,
                                                         q25=~quantile(.x,0.25),
                                                         q75=~quantile(.x,0.75))),.id="Genus")

figD1<-present_sums%>% 
  ggplot(aes(y=reorder(Genus,depth_mymed,decreasing=TRUE),x=depth_mymed))+
  geom_point(size=0.9)+
  geom_errorbarh(aes(xmin=depth_q25,xmax=depth_q75),alpha=0.7,linewidth=0.4)+
  geom_point(data=BioOracle_sums,aes(y=Genus,x=depth_mymed),colour="darkorange",
             position=position_nudge(x = 0, y = 0.1),size=0.9)+
  geom_errorbarh(data=BioOracle_sums,
                 aes(xmin=depth_q25,xmax=depth_q75),colour="darkorange",
                 position=position_nudge(x = 0, y = 0.1),alpha=0.7,linewidth=0.4)+
  theme_light()+
  theme(axis.title = element_blank(),
        axis.text.y=element_text(face="italic"))+
  geom_point(y=10.5,x=500,colour="darkorange",size=1)+
  annotate(geom="text",y=10.5,x=800,label="SPS5.8.5",colour="gray43",size=3)+
  geom_point(y=12,x=500,colour="black",size=1)+
  annotate(geom="text",y=12,x=800,label="Present",colour="gray43",size=3)+
  xlim(c(0,3100))

ROMS_latlon_depth<-lapply(1:30,give_me_info,myfile="ROMS_NWA2_preds_full300")
names(ROMS_latlon_depth)<-colnames(myfile1[[1]])

ROMS_sums<-map_df(ROMS_latlon_depth, ~.x%>% 
                         mutate(depth=abs(depth)) %>%
                         group_by(mydataset) %>% 
                         summarize_at(vars(x,y,depth),mean) %>% 
                         ungroup() %>% 
                         summarize_at(vars(x,y,depth),list(mymed=median,
                                                           q25=~quantile(.x,0.25),
                                                           q75=~quantile(.x,0.75))),.id="Genus")
figD2<-present_sums%>% 
  ggplot(aes(y=reorder(Genus,depth_mymed,decreasing=TRUE),x=depth_mymed))+
  geom_point(size=0.9)+
  geom_errorbarh(aes(xmin=depth_q25,xmax=depth_q75),alpha=0.7,linewidth=0.4)+
  geom_point(data=ROMS_sums,aes(y=Genus,x=depth_mymed),colour="darkorange",
             position=position_nudge(x = 0, y = 0.1),size=0.9)+
  geom_errorbarh(data=ROMS_sums,
                 aes(xmin=depth_q25,xmax=depth_q75),colour="darkorange",
                 position=position_nudge(x = 0, y = 0.1),alpha=0.7,linewidth=0.4)+
  theme_light()+
  theme(axis.title.y = element_blank(),
        axis.text.y=element_text(face="italic"))+xlab("Depth (m)")+
  geom_point(y=10.5,x=500,colour="darkorange",size=1)+
  annotate(geom="text",y=10.5,x=800,label="RCP8.5",colour="gray43",size=3)+
  geom_point(y=12,x=500,colour="black",size=1)+
  annotate(geom="text",y=12,x=800,label="Present",colour="gray43",size=3)+
  xlim(c(0,3100))

pdf("./Outputs/graphs/Fig6.pdf",width=5.2,height=8)
ggarrange(figD1,figD2, ncol=1)
#ggsave("Suitable_habitat.pdf", figA, width=5,height=9,device=cairo_pdf)
dev.off()

save(BioOracle_latlon_depth,ROMS_latlon_depth,present_latlon_depth,
        file="./Outputs/depth_lat_lon.RData")

## Calculate centroid of distribution (lat)

# Transform y and x from UTM to lat lon
transform_xy<-function(x){
v <- vect(cbind(x$x_mymed,x$y_mymed), 
          crs="+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs")
v_proj<-project(v, "+proj=longlat +datum=WGS84")
lonlat <- geom(v_proj)[, c("x", "y")]
x<-x %>% 
  mutate(x_mymed=lonlat[,1],y_mymed=lonlat[,2])

v <- vect(cbind(x$x_q25,x$y_q25), 
          crs="+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs")
v_proj<-project(v, "+proj=longlat +datum=WGS84")
lonlat <- geom(v_proj)[, c("x", "y")]
x<-x %>% 
  mutate(x_q25=lonlat[,1],y_q25=lonlat[,2])

v <- vect(cbind(x$x_q75,x$y_q75), 
          crs="+proj=utm +zone=19 +datum=WGS84 +units=m +no_defs")
v_proj<-project(v, "+proj=longlat +datum=WGS84")
lonlat <- geom(v_proj)[, c("x", "y")]
x<-x %>% 
  mutate(x_q75=lonlat[,1],y_q75=lonlat[,2])

return(x)
}

present_sums_tsf<-transform_xy(present_sums)
BioOracle_sums_tsf<-transform_xy(BioOracle_sums)
ROMS_sums_tsf<-transform_xy(ROMS_sums)


lat1<-present_sums_tsf%>% 
  ggplot(aes(y=reorder(Genus,y_mymed,decreasing=TRUE),x=y_mymed))+
  geom_point(size=0.9)+
  geom_errorbarh(aes(xmin=y_q25,xmax=y_q75),alpha=0.7,linewidth=0.4)+
  geom_point(data=BioOracle_sums_tsf,aes(y=Genus,x=y_mymed),colour="darkorange",
             position=position_nudge(x = 0, y = 0.1),size=0.9)+
  geom_errorbarh(data=BioOracle_sums_tsf,
                 aes(xmin=y_q25,xmax=y_q75),colour="darkorange",
                 position=position_nudge(x = 0, y = 0.1),alpha=0.7,linewidth=0.4)+
  theme_light()+
  #ggtitle("Present vs NW ROMS")+
  theme(axis.title = element_blank())+
  geom_point(y=25,x=42,colour="darkorange",size=1)+
  annotate(geom="text",y=25,x=42.3,label="SP5.8.5",colour="gray43",size=3)+
  geom_point(y=27,x=42,colour="black",size=1)+
  annotate(geom="text",y=27,x=42.3,label="Present",colour="gray43",size=3)+
  xlim(39.8,42.75)

lat2<-present_sums_tsf%>% 
  ggplot(aes(y=reorder(Genus,y_mymed,decreasing=TRUE),x=y_mymed))+
  geom_point(size=0.9)+
  geom_errorbarh(aes(xmin=y_q25,xmax=y_q75),alpha=0.7,linewidth=0.4)+
  geom_point(data=ROMS_sums_tsf,aes(y=Genus,x=y_mymed),colour="darkorange",
             position=position_nudge(x = 0, y = 0.1),size=0.9)+
  geom_errorbarh(data=ROMS_sums_tsf,
                 aes(xmin=y_q25,xmax=y_q75),colour="darkorange",
                 position=position_nudge(x = 0, y = 0.1),alpha=0.7,linewidth=0.4)+
  theme_light()+
  #ggtitle("Present vs NW ROMS")+
  theme(axis.title.y = element_blank())+xlab("Latitude")+
  geom_point(y=25,x=42,colour="darkorange",size=1)+
  annotate(geom="text",y=25,x=42.3,label="RCP8.5",colour="gray43",size=3)+
  geom_point(y=27,x=42,colour="black",size=1)+
  annotate(geom="text",y=27,x=42.3,label="Present",colour="gray43",size=3)+
  xlim(39.8,42.75)

pdf("./Outputs/graphs/Appendix6-Species-Lat.pdf",width=5.2,height=8)
ggarrange(lat1,lat2, ncol=1)
dev.off()