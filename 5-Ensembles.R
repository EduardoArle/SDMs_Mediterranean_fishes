################################################################################
#                                Make ensembles                                #
################################################################################

#load packages
library(sf); library(terra); library(pROC); library(PresenceAbsence)

#list WDs
wd_projections <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Model_projections'
wd_eval <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Evaluation'
wd_ens <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Ensembles'
wd_pr_pa <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/pr_PA'
wd_plots <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Plots_ensemble'
  
#determine CRS we are using
CRS <- "GEOGCRS[\"WGS 84 (CRS84)\",\n    DATUM[\"World Geodetic System 1984\",\n        ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n            LENGTHUNIT[\"metre\",1]]],\n    PRIMEM[\"Greenwich\",0,\n        ANGLEUNIT[\"degree\",0.0174532925199433]],\n    CS[ellipsoidal,2],\n        AXIS[\"geodetic longitude (Lon)\",east,\n            ORDER[1],\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        AXIS[\"geodetic latitude (Lat)\",north,\n            ORDER[2],\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n    USAGE[\n        SCOPE[\"unknown\"],\n        AREA[\"World\"],\n        BBOX[-90,-180,90,180]],\n    ID[\"OGC\",\"CRS84\"]]"

#list species (for now I have only minT options)
setwd(wd_projections)

sps_list <- unique(gsub("_minT_\\d+\\.grd$", "",
                 gsub('Pred_', '', list.files(pattern = '.grd$'))))

# #create colour ramp to represent the values
colramp <- colorRampPalette(c("#9e0142", "#d53e4f", "#f46d43",
                              "#fdae61", "#fee08b", "#ffffbf",
                              "#e6f598", "#abdda4", "#66c2a5",
                              "#3288bd", "#08306b"))

### Install adapted function myGradientLegend
### (script modified_function_gradientLegend)

#create vectors to store AUC and TSS results
AUC_minT <- numeric()
TSS_minT <- numeric()

#loop through species ensembling the resulting maps

for(i in 1:length(sps_list))
{
  #load all results per species (for now I have only minT options)
  setwd(wd_projections)
  projs_list <- list.files(pattern = paste0(sps_list[i], ".*\\.grd$"))
  projs_rasts <- lapply(projs_list, rast)
  
  #load evaluation metrics (for now I have only minT options)
  setwd(wd_eval)
  eval <- read.csv(paste0('eval_model_', sps_list[i], '_minT.csv'))
  
  #binarise each map according to its evaluation metric 
  projs_bin <- projs_rasts
  for(j in 1:length(projs_rasts))
  {
    #get model ID of each successful prediction (AUC and TSS satisfactory)
    rast_name <- names(projs_rasts[[j]])
    model_ID <- sub("id_(\\d+)__sp.*", "\\1", rast_name)
    
    #find threshold for the specific raster
    th <- eval$threshold[as.character(eval$modelID) == model_ID]
    
    #binarise
    projs_bin[[j]] <- ifel(projs_bin[[j]] < th, 0, 1)
  }
  
  #sum all layers and calculate rate of agreement per cell
  ens_minT <- sum(rast(projs_bin)) / length(projs_bin) 
  ##ens_meanT <- sum(rast(projs_bin)) / length(projs_bin)  #change accordingly
  ##ens_maxT <- sum(rast(projs_bin)) / length(projs_bin)   #change accordingly

  #save ensemble
  setwd(wd_ens)
  writeRaster(ens_minT, 
              paste0('Ensemble_', sps_list[i], '_minT.tif'), filetype = 'GTiff')
  
  ## Calculate metrics of the ensemble ##
  
  #load species presences and PAs
  setwd(wd_pr_pa)
  occ <- read.csv(paste0(sps_list[i], ' .csv'))
  
  #make spatial object
  occ_sf <- st_as_sf(occ, coords = c('lon', 'lat'), crs = CRS)
  
  #extract precictions in locations of occurrences
  pred_minT <- extract(ens_minT, occ_sf)[,2]
  
  ##  AUC  ##
  
  roc_obj <- roc(occ$occurrence, pred_minT)
  auc_val_minT <- roc_obj$auc

  #store result
  AUC_minT[i] <- auc_val_minT
  
  ##  TSS  ##
  
  #use PresenceAbsence package to calculate across thresholds
  df_minT <- data.frame(ID = 1:nrow(occ), obs = occ$occurrence,
                            pred = pred_minT)
  
  #calculate accuracy measures for many thresholds
  acc_minT <- presence.absence.accuracy(df_minT, threshold = seq(0, 1, 0.01))
  
  #find threshold with maximum TSS
  best_row_minT <- acc_minT[which.max(acc_minT$sensitivity +
                                            acc_minT$specificity - 1), ]
  tss_val_minT <- best_row_minT$sensitivity + best_row_minT$specificity - 1
  
  TSS_minT[i] <- tss_val_minT
  
  #save map plots with values
  setwd(wd_plots)
  pdf(file = paste0(sps_list[i], "_minT.pdf"))
  
  #plot background med
  plot(ens_minT, col = colramp(200),
       main = sps_list[i], font.main = 3, cex.main = 2,
       box = NA, legend = NA, axes = NA)
  
  #add evaluation metrics
  text(1, 33, paste0('AUC = ', round(auc_val_minT, 2)), cex = 1.2)
  text(1, 31.5, paste0('TSS = ', round(tss_val_minT, 2)), cex = 1.2)
  
  #plot legend
  myGradientLegend(valRange = c(0, 1),
                   pos=c(0.3,0.1,0.8,0.12),
                   color = colramp(200),
                   side = 1,
                   n.seg = 0,
                   values = c(0, 1),
                   cex = 2)

  dev.off()
}
