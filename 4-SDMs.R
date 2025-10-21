################################################################################
#                                Run models                                    #
################################################################################

#load packages
library(sf); library(terra); library(data.table); library(sdm); library(usdm)
library(sp); library(raster) #for the sdm pakcage

#list WDs
wd_pr_pa <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/pr_PA'
wd_vars_pr <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Variables/BioOracle/Processed_present'
wd_vars_2090 <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Variables/BioOracle/Processed_2090_ssp585'
wd_SDM_data <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/SDM_data'
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Models'
wd_eval <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Evaluation'
wd_projections <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Model_projections'
wd_proj_future <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Model_projections_2090'


#load present variable layers
setwd(wd_vars_pr)
vars <- lapply(list.files(pattern = '.tif$'), raster)

#name variables
names(vars) <- lapply(vars, names)

#stack variables 
vars_minT <- stack(vars$so_min, vars$phyc_mean, vars$sws_mean, vars$thetao_min)
vars_meanT <- stack(vars$so_min, vars$phyc_mean, vars$sws_mean, vars$thetao_mean)
vars_maxT <- stack(vars$so_min, vars$phyc_mean, vars$sws_mean, vars$thetao_max)

#load 2090 variable layers
setwd(wd_vars_2090)
vars_2090 <- lapply(list.files(pattern = '.tif$'), raster)

#name variables
names(vars_2090) <- lapply(vars_2090, names)

#stack variables 
vars_minT_fut <- stack(vars_2090$so_min, vars_2090$phyc_mean,
                       vars_2090$sws_mean, vars_2090$thetao_min)
vars_meanT_fut <- stack(vars_2090$so_min, vars_2090$phyc_mean,
                        vars_2090$sws_mean, vars_2090$thetao_mean)
vars_maxT_fut <- stack(vars_2090$so_min, vars_2090$phyc_mean,
                       vars_2090$sws_mean, vars_2090$thetao_max)

#load occurrence data
setwd(wd_pr_pa)
occ_data <- lapply(list.files(), read.csv)

#name objects
names(occ_data) <- sub(' .csv', '', list.files())

#harmonise the occurrence table with the requirements of the package
occ_data_2 <- lapply(occ_data, function(x){
  data.frame(occurrence = x$occurrence,
             lon = x$lon,
             lat = x$lat)})

#create spatial points data frames and inform the geographic system
occ_data_sp <- occ_data_2 

for(i in 1:length(occ_data_sp))
{
  coordinates(occ_data_sp[[i]]) <- ~ lon + lat
  proj4string(occ_data_sp[[i]]) <- crs(vars[[1]])
  
  print(i)
}

#extract values from all points
vals_minT <- lapply(occ_data_sp, function(x){
  as.data.frame(extract(vars_minT, x))
})

vals_meanT <- lapply(occ_data_sp, function(x){
  as.data.frame(extract(vars_meanT, x))
})

vals_maxT <- lapply(occ_data_sp, function(x){
  as.data.frame(extract(vars_maxT, x))
})


#prepare data frames (only species with enough data and no corel issues)
vals_minT_2 <- list()
for(i in 1:length(vals_minT))
{
  #check whether there are at least half PA as pr
  if(length(which(occ_data_2[[i]]$occurrence == 0)) >= 
     length(which(occ_data_2[[i]]$occurrence == 1)) / 2){
    
    #check whether the variables are not too correlate
    cor <- vifcor(vals_minT[[i]], th = 0.7)
    
    if(length(cor@excluded) == 0){
      
      vals_minT_2[[length(vals_minT_2) + 1]] <- cbind(occ_data_2[[i]],
                                                      vals_minT[[i]])
      
      names(vals_minT_2)[length(vals_minT_2)] <- names(vals_minT)[i]
    }
  }
  print(i)
}

vals_meanT_2 <- list()
for(i in 1:length(vals_meanT))
{
  #check whether there are at least half PA as pr
  if(length(which(occ_data_2[[i]]$occurrence == 0)) >= 
     length(which(occ_data_2[[i]]$occurrence == 1)) / 2){
    
    #check whether the variables are not too correlate
    cor <- vifcor(vals_meanT[[i]], th = 0.7)
    
    if(length(cor@excluded) == 0){
      
      vals_meanT_2[[length(vals_meanT_2) + 1]] <- cbind(occ_data_2[[i]],
                                                        vals_meanT[[i]])
      
      names(vals_meanT_2)[length(vals_meanT_2)] <- names(vals_meanT)[i]
    }
  }
  print(i)
}

vals_maxT_2 <- list()
for(i in 1:length(vals_maxT))
{
  #check whether there are at least half PA as pr
  if(length(which(occ_data_2[[i]]$occurrence == 0)) >= 
     length(which(occ_data_2[[i]]$occurrence == 1)) / 2){
    
    #check whether the variables are not too correlate
    cor <- vifcor(vals_maxT[[i]], th = 0.7)
    
    if(length(cor@excluded) == 0){
      
      vals_maxT_2[[length(vals_maxT_2) + 1]] <- cbind(occ_data_2[[i]],
                                                      vals_maxT[[i]])
      
      names(vals_maxT_2)[length(vals_maxT_2)] <- names(vals_maxT)[i]
    }
  }
  print(i)
}

#prepare data object
data_obj_minT <- lapply(vals_minT_2, function(x){
  sdmData(formmula = occurrence ~ . + coords(lon+lat), train = x)
})
                        
data_obj_meanT <- lapply(vals_meanT_2, function(x){
  sdmData(formmula = occurrence ~ . + coords(lon+lat), train = x)
})

data_obj_maxT <- lapply(vals_maxT_2, function(x){
  sdmData(formmula = occurrence ~ . + coords(lon+lat), train = x)
})

#save sdmData objects
setwd(wd_SDM_data)

mapply(function(x, y) {
  write.sdm(x, paste0('sdmData_', y, '_minT'))},
  data_obj_minT, names(data_obj_minT))

mapply(function(x, y) {
  write.sdm(x, paste0('sdmData_', y, '_meanT'))},
  data_obj_meanT, names(data_obj_meanT))

mapply(function(x, y) {
  write.sdm(x, paste0('sdmData_', y, '_maxT'))},
  data_obj_maxT, names(data_obj_maxT))

#run models
sdm_minT <- lapply(data_obj_minT, function(x){
  sdm(occurrence ~ ., data = x,
      methods = c('brt', 'cart', 'fda', 'glm', 'mars', 'maxlike',
                  'mda', 'gam', 'rf', 'svm'), 
      replication = 'cv', cv.folds = 5, n = 5)
})

sdm_meanT <- lapply(data_obj_meanT, function(x){
  sdm(occurrence ~ ., data = x,
      methods = c('brt', 'cart', 'fda', 'glm', 'mars', 'maxlike',
                  'mda', 'gam', 'rf', 'svm'), 
      replication = 'cv', cv.folds = 5, n = 5)
})

sdm_maxT <- lapply(data_obj_maxT, function(x){
  sdm(occurrence ~ ., data = x,
      methods = c('brt', 'cart', 'fda', 'glm', 'mars', 'maxlike',
                  'mda', 'gam', 'rf', 'svm'), 
      replication = 'cv', cv.folds = 5, n = 5)
})


#save model objects
setwd(wd_models)

mapply(function(x, y) {
  write.sdm(x, paste0('models_', y, '_minT'))},
  sdm_minT, names(sdm_minT))

mapply(function(x, y) {
  write.sdm(x, paste0('models_', y, '_meanT'))},
  sdm_meanT, names(sdm_meanT))

mapply(function(x, y) {
  write.sdm(x, paste0('models_', y, '_maxT'))},
  sdm_maxT, names(sdm_maxT))


#get model evaluation
eval_minT <- lapply(sdm_minT, function(x){
  getEvaluation(x, w = 1:250,
                wtest='test.dep', 
                stat=c('AUC','TSS','th'), opt = 2)
})

eval_meanT <- lapply(sdm_meanT, function(x){
  getEvaluation(x, w = 1:250,
                wtest='test.dep', 
                stat=c('AUC','TSS','th'), opt = 2)
})


eval_maxT <- lapply(sdm_maxT, function(x){
  getEvaluation(x, w = 1:250,
                wtest='test.dep', 
                stat=c('AUC','TSS','th'), opt = 2)
})


#save evaluation metrics
setwd(wd_eval)

mapply(function(x, y) {
  write.csv(x, paste0('eval_model_', y, '_minT.csv'), row.names = F)},
  eval_minT, names(eval_minT))

mapply(function(x, y) {
  write.csv(x, paste0('eval_model_', y, '_meanT.csv'), row.names = F)},
  eval_meanT, names(eval_meanT))

mapply(function(x, y) {
  write.csv(x, paste0('eval_model_', y, '_maxT.csv'), row.names = F)},
  eval_maxT, names(eval_maxT))


#select models with TSS higher than 0.5 and AUC higher than 0.7
sel_minT <- lapply(eval_minT, function(x){
  x[which(x$TSS >= 0.5 & x$AUC >= 0.7),]
})

sel_meanT <- lapply(eval_meanT, function(x){
  x[which(x$TSS >= 0.5 & x$AUC >= 0.7),]
})

sel_maxT <- lapply(eval_maxT, function(x){
  x[which(x$TSS >= 0.5 & x$AUC >= 0.7),]
})


###################
##### PRESENT #####
###################


#project all selected models
setwd(wd_projections)

for(i in 1:length(sdm_minT))
{
  for(j in 1:nrow(sel_minT[[i]]))
  {
    predict(sdm_minT[[i]], id = sel_minT[[i]]$modelID[j], newdata = vars_minT, 
            filename = paste0('Pred_', names(sdm_minT[i]), '_minT_',
                              sel_minT[[i]]$modelID[j], '.grd'))
  }
}


for(i in 1:length(sdm_meanT))
{
  for(j in 1:nrow(sel_meanT[[i]]))
  {
    predict(sdm_meanT[[i]], id = sel_meanT[[i]]$modelID[j], newdata = vars_meanT, 
            filename = paste0('Pred_', names(sdm_meanT[i]), '_meanT_',
                              sel_meanT[[i]]$modelID[j], '.grd'))
  }
}


for(i in 1:length(sdm_maxT))
{
  for(j in 1:nrow(sel_maxT[[i]]))
  {
    predict(sdm_maxT[[i]], id = sel_maxT[[i]]$modelID[j], newdata = vars_maxT, 
            filename = paste0('Pred_', names(sdm_maxT[i]), '_maxT_',
                              sel_maxT[[i]]$modelID[j], '.grd'))
    print(j)
  }
}










###################
#####   2090  #####
###################


#project models
setwd(wd_proj_future)
pred_2090 <- list()
for(i in 1:200)
{
  pred_2090[[i]] <- predict(sdm, id = i, newdata = vars_2090, 
                          filename = paste0('Epinephelus_aeneus_', i, '_',
                                            eval$Algorithm[i], '.grd'))
}







##### load and process projections


###################
##### PRESENT #####
###################

#load present 
setwd(wd_projections)
pred_pr <- lapply(list.files(pattern = ".grd$"), raster)
names(pred_pr) <- gsub("[^0-9]", "", list.files(pattern = ".grd$"))

#binarise projections according to threshold
pred_pr_bin <- list()
for(i in 1:length(pred_pr))
{
  pred_pr_bin[[i]] <- pred_pr[[i]]
  a <- which(eval$modelID == as.numeric(names(pred_pr)[i]))
  th <- eval$threshold[a]
  pred_pr_bin[[i]][] <- ifelse(pred_pr_bin[[i]][] >= th, 1, 0)
  print(i)
}

names(pred_pr_bin) <- names(pred_pr)

#select models with TSS >= 0.4 and AUC >= 0.7
pred_pr_sel <- pred_pr_bin[names(pred_pr_bin) %in% 
                                 as.character(eval$modelID[which(
                                   eval$TSS >= 0.4 & eval$AUC >= 0.7)])]

#stack all selected projections projections
pred_pr_stack <- stack(pred_pr_sel)

#sum all layers and calculate percentage of agreement
pred_pr_ens <- sum(pred_pr_stack) / length(pred_pr_sel) * 100




###################
#####   2090  #####
###################


  
#load 2090
setwd(wd_proj_future)
pred_2090 <- lapply(list.files(pattern = ".grd$"), raster)
names(pred_2090) <- gsub("[^0-9]", "", list.files(pattern = ".grd$"))

#binarise projections according to threshold
pred_2090_bin <- list()
for(i in 1:length(pred_2090))
{
  pred_2090_bin[[i]] <- pred_2090[[i]]
  a <- which(eval$modelID == as.numeric(names(pred_2090)[i]))
  th <- eval$threshold[a]
  pred_2090_bin[[i]][] <- ifelse(pred_2090_bin[[i]][] >= th, 1, 0)
  print(i)
}

names(pred_2090_bin) <- names(pred_2090)

#select models with TSS >= 0.4 and AUC >= 0.7
pred_2090_sel <- pred_2090_bin[names(pred_2090_bin) %in% 
                                 as.character(eval$modelID[which(
                                   eval$TSS >= 0.4 & eval$AUC >= 0.7)])]
  

#stack all selected projections projections
pred_2090_stack <- stack(pred_2090_sel)

#sum all layers and calculate percentage of agreement
pred_2090_ens <- sum(pred_2090_stack) / length(pred_2090_sel) * 100


###### PLOTS ######

#create colour ramp to represent the values
colramp <- colorRampPalette(c("#d53e4f", "#f46d43",
                              "#fdae61", "#fee08b", "#ffffbf",
                              "#e6f598", "#abdda4", "#66c2a5",
                              "#238b45"))


plot(pred_pr_ens, main = 'Present', col = colramp(200))

plot(pred_2090_ens, main = '2090', col = colramp(200), zlim = c(0, 100))

t <- pred_pr_ens / pred_2090_ens

i=51
plot(pred_pr[[i]])
plot(pred_2090[[i]])


i=42
t <- pred_pr[[i]] - pred_2090[[i]]
t

i=3
vars[[i]]
vars_2090[[i]]

t <- vars[[i]] - vars_2090[[2]]
t
plot(t)


###### test projecting models to fake layers with absurd values and see what happens

#create absurd variable values
vars_absurd <- vars
vars_absurd[[1]] <- vars_absurd[[1]] - 5
vars_absurd[[2]] <- vars_absurd[[2]] + 7
vars_absurd[[3]] <- vars_absurd[[3]] - 11

setwd('/Users/carloseduardoaribeiro/Documents/Collaborations/Uri Roll/Tests_projection')

test_pr <- predict(sdm, id = 1, newdata = vars, filename = 'Test_pr.grd')

test_2090 <- predict(sdm, id = 1, newdata = vars_2090,
                              filename = 'Test_2090.grd')

test_absurd_values <- predict(sdm, id = 1, newdata = vars_absurd,
                              filename = 'Test_absurd_values.grd')

plot(test_absurd_values)
