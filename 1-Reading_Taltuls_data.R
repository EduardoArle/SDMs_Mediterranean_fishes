#load packages
library(data.table); library(sf); library(terra)

#list WDs
wd_med <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med'
wd_data <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Data'
wd_proc <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Processes_data'
wd_vars <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Variables_present'
wd_maps <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Maps'

#laod one variable for CRS and ploting
setwd(wd_vars)
var <- rast(list.files()[2])

#manually input CRS. Not coming automatically. Perhaps .nc issue?
crs(var) <- "GEOGCRS[\"WGS 84 (CRS84)\",\n    DATUM[\"World Geodetic System 1984\",\n        ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n            LENGTHUNIT[\"metre\",1]]],\n    PRIMEM[\"Greenwich\",0,\n        ANGLEUNIT[\"degree\",0.0174532925199433]],\n    CS[ellipsoidal,2],\n        AXIS[\"geodetic longitude (Lon)\",east,\n            ORDER[1],\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        AXIS[\"geodetic latitude (Lat)\",north,\n            ORDER[2],\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n    USAGE[\n        SCOPE[\"unknown\"],\n        AREA[\"World\"],\n        BBOX[-90,-180,90,180]],\n    ID[\"OGC\",\"CRS84\"]]"

#make ID raster
ID_rast <- var
non_na_cells <- which(!is.na(values(ID_rast)))
new_vals <- rep(NA, ncell(ID_rast))
new_vals[non_na_cells] <- seq_along(non_na_cells)
ID_rast[] <- new_vals
names(ID_rast) <- 'ID'

#make a polygon with mediterranean boundaries
coords_med <- matrix(c(-6, 30, 38, 30, 38, 46, 2, 46, -6, 37, -6, 30),  
                     ncol = 2, byrow = TRUE)

poly_med <- st_polygon(list(coords_med))
study_area_med <- st_sf(geometry = st_sfc(poly_med, crs = crs(var)))

#crop ID raster to mediterranean resolution
cropped <- crop(ID_rast, study_area_med)   
ID_med <- mask(cropped, study_area_med) 

#load data
setwd(wd_data)
load('Clean_Medata_lis.Rdata')

#understand the structure of the data
n_entries <- nrow(medata)
n_sps <- length(unique(medata$species))
sps_list <- unique(medata$species)

#split data by species
sps <- split(medata, medata$species)

#### no species has at least 10 unique absences and 10 unique presences

#select only species with minimum data adequacy for ENMs (30 pr, unique cellID)
sps_sel <- list()
for(i in 1:length(all_sps))
{
  #make a vector o presences and absences
  sps[[i]]$occurrence <- ifelse(sps[[i]]$sp.n == 0, 0, 1)
  
  #create a spatial object
  sps_sf <- st_as_sf(sps[[i]], coords = c('lon', 'lat'), crs = crs(var))
  
  #extract cell ID 
  ID <- extract(ID_rast, sps_sf)
  sps[[i]]$ID <- ID[,2]
  
  #split data into presence and absence
  sps_pr <- sps[[i]][sps[[i]]$occurrence == 1,]
  sps_abs <- sps[[i]][sps[[i]]$occurrence == 0,]
  
  #elimintate occurrences on land (NAs in the ID)
  sps_pr <- sps_pr[complete.cases(sps_pr$ID),]
  sps_abs <- sps_abs[complete.cases(sps_abs$ID),]
  
  #delete absences in same cell as presences
  if(nrow(sps_pr) > 0){
    sps_abs <- sps_abs[-which(sps_abs$cellID %in% sps_pr$cellID),]
  }
  
  #delete duplicates
  sps_pr <- unique(as.data.table(sps_pr), by = 'ID')
  sps_abs <- unique(as.data.table(sps_abs), by = 'ID')
  
  #join pr and abs
  sps_all <- rbind(sps_pr, sps_abs)
  
  #check whether there are enough presences and absences to run ENMs
  if(length(which(sps_all$occurrence == "1")) >= 30){
    sps_sel[[length(sps_sel) + 1]] <- sps_all
    
    #write table
    setwd(wd_proc)
    write.csv(sps_all, paste0(unique(sps_all$species), '.csv'), row.names = F)
    
    #make a spatial object
    sps_all <- st_as_sf(sps_all, coords = c('lon', 'lat'), crs = crs(var))
    
    #check if all points are in the med region
    test_med <- extract(ID_med, sps_all)
    
    if(NA %in% test_med){stop('Points out of the region!')}
    
    #save images of the maps
    setwd(wd_maps)
    pdf(file=paste0(unique(sps_all$species),".pdf"))
    
    #plot background med
    plot(ID_med, col = 'grey80',
         main = gsub('\\.', ' ',unique(sps_all$species)))
    
    #plot occ
    plot(sps_all[which(sps_all$occurrence == 1),], add = T,
         pch = 19, col = 'darkgreen')
    plot(sps_all[which(sps_all$occurrence == 0),], add = T,
         pch = 19, col = 'orange')
    
    dev.off()
    
  }
  print(i)
}

#make list of species with at lest 30 presences in unique cells (var res)
species <- character()
n_pr <- numeric()
n_abs <- numeric()

for(i in 1:length(sps_sel))
{
  species[i] <- unique(sps_sel[[i]]$species)
  n_pr[i] <- length(which(sps_sel[[i]]$occurrence == 1))
  n_abs[i] <- length(which(sps_sel[[i]]$occurrence == 0))
}

#make table
sps_occ <- data.frame(species = species, n_presences = n_pr, n_absences = n_abs)

#save
setwd(wd_med)
write.csv(sps_occ, 'Species_occurrence.csv', row.names = F)

