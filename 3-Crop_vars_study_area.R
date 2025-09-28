################################################################################
#                      Preparing BioOracle variables                           #
################################################################################


#load packages
library(sf); library(terra)

#list WDs
wd_raw_vars_pr <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Variables/BioOracle/Raw_present'
wd_vars_pr <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Variables/BioOracle/Processed_present'
wd_raw_vars_2090 <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Variables/BioOracle/Raw_2090_ssp585'
wd_vars_2090 <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Variables/BioOracle/Processed_2090_ssp585'

#determine CRS we are using
CRS <- "GEOGCRS[\"WGS 84 (CRS84)\",\n    DATUM[\"World Geodetic System 1984\",\n        ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n            LENGTHUNIT[\"metre\",1]]],\n    PRIMEM[\"Greenwich\",0,\n        ANGLEUNIT[\"degree\",0.0174532925199433]],\n    CS[ellipsoidal,2],\n        AXIS[\"geodetic longitude (Lon)\",east,\n            ORDER[1],\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        AXIS[\"geodetic latitude (Lat)\",north,\n            ORDER[2],\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n    USAGE[\n        SCOPE[\"unknown\"],\n        AREA[\"World\"],\n        BBOX[-90,-180,90,180]],\n    ID[\"OGC\",\"CRS84\"]]"

#make a polygon with mediterranean boundaries
coords_med <- matrix(c(-6, 30, 38, 30, 38, 46, 2, 46, -6, 37, -6, 30),  
                     ncol = 2, byrow = TRUE)

poly_med <- st_polygon(list(coords_med))
study_area <- st_sf(geometry = st_sfc(poly_med, crs = CRS))


###################
##### PRESENT #####
###################

#load variables
setwd(wd_raw_vars_pr)

vars <- lapply(list.files(), rast)
names(vars) <- lapply(vars, names)

#mask, crop and save variables by study area
setwd(wd_vars_pr)
crop_vars_pr <- list() #to check
for(i in 1:length(vars))
{
  mask_vars <- mask(vars[[i]], study_area)
  crop_vars_pr[[i]] <- crop(mask_vars, study_area)
  writeRaster(crop_vars_pr[[i]],
              paste0(names(crop_vars_pr[[i]]),'.tif'))
  print(i)
}



####################
####### 2090 #######
####################

#load variables
setwd(wd_raw_vars_2090)

vars_2090 <- lapply(list.files(), rast)
names(vars_2090) <- lapply(vars_2090, names)

#mask, crop and save variables by study area
setwd(wd_vars_2090)
crop_vars_2090 <- list() #to check
for(i in 1:length(vars_2090))
{
  mask_vars_2090 <- mask(vars_2090[[i]], study_area)
  crop_vars_2090[[i]] <- crop(mask_vars_2090, study_area)
  writeRaster(crop_vars_2090[[i]],
              paste0(names(crop_vars_2090[[i]]),'.tif'))
  print(i)
}




