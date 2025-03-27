# ClimateChangeGOM
This page includes the data and code for the manuscript: Climate change drives bathymetric shifts in species and trait diversity of deep-sea benthic communities 


## Scripts

### S1_Predictions_communities ###
Here, we use the developed HMSC model (*models_thin_100_samples_250_chains_4.RData*) to predict coral communities under three climate scenarios. For each climate scenario, the script produces a list with 1000 predictions that include presence/absence of each of the 30 coral genera for each grid cell in our study area and stores it in the file *[model_name]_preds_full1000.RData*. It also produces a file with the following metrics for each cell of the study area: *Srich* which is the total number of coral genera, *sdS* which is the standard deviation of the total number of coral genera, as well as *S5* and *S1*0 which represent the probability of having a total of 5 or 10 species respectively. The metrics are stored as table (*[model_name]_metrics_table.RData*) and as raster files (*[model_name]_metrics_raster.RData*). Subsequently, we use the lists of predictions to calculate the probability of occurrence of each coral genus in each grid cell, and store it in a matrix (*freq_occ_[model_name].RData*). Finally, the script produces a map with the delta values of the frequency of occurrence of each coral genus between present conditions and each climate scenario. 

### S2_Partition_predictions ###
In this script, we partition the lists of the full predictions *[model_name]_preds_full1000.RData* into two smaller files to facilitate further processing and calculations of other indices of change (*./partitions/[model]_preds_full300_part1.RData* and *./partitions/[model]_preds_full300_part2.RData*). Due to their large size, partitions are stored in an external repository, and can be accessed through [this link](https://www.dropbox.com/scl/fo/5x03nla4e67zy5b962b0s/APyVOpYiOzsRInqrwXJCSEo?rlkey=rqg8gnov6bh1r5bdkjnj1i60h&st=iz8n2uua&dl=0)

### S3_Calculate_areaSH ###
In this script we take predictions (*[model]/preds_full300.RData*) as inputs and compute metrics of change of the spatial distributions of the 30 coral genera within the study area. Firstly, we compute the total area of suitable habitat for each of the 1000 produced datasets of each climate scenario (present, RCP5.8, SPS5.8.5), and subsequently estimate the median, as well as the 25th and 75th percentiles to provide a measure of uncertainty. We proceed by repeating the same process to estimate the median, 25th and 75th percentile of the average depth and centroid for each genus and climate scenario. The outputs are figures of the three metrics, as well as a dataset with the results (*depth_lat_lon.RData*)

### S4_Calculate_metrics ###
In this script, we use the partitions of the predicted data (*./partitions/[model]_preds_full300_part1.RData*) to estimate metrics of change, including trait diversity indices (Functional richness, functional evenness and functional dispersion), as well as community weighted mean (CWM) of the coral traits in each cell of our study area and for each climate scenario. The trait diversity indices are stored in the file *[model]_all_trait_div_indices_partx.RData*, while CWMs are stored in the file *[model]_final_trait_indices_partx.RData*, both stored in an external repository, due to their large size (these files can be accessed through [this link](https://www.dropbox.com/scl/fo/5x03nla4e67zy5b962b0s/APyVOpYiOzsRInqrwXJCSEo?rlkey=rqg8gnov6bh1r5bdkjnj1i60h&st=iz8n2uua&dl=0)) . This script requires a lot of memory and can only run in facilities of High Performance Computing. The script was run for each of the partitions separately, to optimize processing time

### S4-2_Pack_beluga_outputs ###
This is a continuation of script 4. Here we obtain the outputs from script 4, merge the partitions in one file and transform the outputs into spatial files that can be later used to calculate delta values.

### S5_compare_metrics ###
Here we estimate delta values for all metrics, i.e. the differences of these metrics between present conditions and future climate scenarios and save them as spatial raster files. First we focus on taxonomic richness. We estimate the median, 25th, and 75th percentiles and then the delta values, and  we store the output delta values in two raster files (*bOracle_Sperc.tif*, and *ROMS_Sperc.tif*). We subsequently calculate the Bray-Curtis dissimilarity index and perform a cell-wise comparison of coral communities between present conditions and climate change scenarios. The output is stored in a raster file (*dissimilarity_indexFreqOcc.tif*).Lastly, we calculate delta values for all trait diversity indices and the community weighted mean. The outputs are stored in raster files: *[model]_ROMS_Trait_div_deltas.tif* for trait diversity indices and *[model]_CWM_deltas.tif* for community weighted mean

### S6_evaluate_uncertainty ###
In this script we estimate the interquantile range (IQR) for each metric, which is used to evaluate uncertainty. The IQRs are stored in spatial raster files that can be found on the folder *Outputs/Spatial_files/Uncertainty*. Subsequently, we make figures of the IQR range versus depth, and store then in the *Outputs/graphs folder*.

### S7_spatial_extrapolation ###
Uses the package dsmextra to evaluate the number of cells in the study area that have analogous conditions to those that were used to train the Hmsc model. Takes as input the environmental data for all the study area (*Data/joined_env_data_bioOracle.RData*) and environmental data that were used to train the model (*Data/joined_env_data_ROMS2.RData*). Produces a map with the locations of analogous and non-analogous areas, and a table with their percentages.

