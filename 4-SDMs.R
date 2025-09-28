################################################################################
#                                Run models                                    #
################################################################################

#load packages
library(sf); library(terra); library(data.table); library(sdm)
library(sp); library(raster) #for the sdm pakcage

#list WDs
wd_pr_pa <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/pr_PA'
wd_vars_pr <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Variables/BioOracle/Processed_present'
wd_vars_2090 <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Variables/BioOracle/Processed_2090_ssp585'
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Models'
wd_evaluation <- '/Users/carloseduardoaribeiro/Documents/Post-doc/ENMs_med/Model_evaluation'
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
  proj4string(occ_data_sp[[i]]) <- crs(vars)
  
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
for(i in length(vals_minT))
{
  #check whether there are at least half PA as pr
  if(length(which(occ_data_2[[i]]$occurrence == 0)) >= 
     length(which(occ_data_2[[i]]$occurrence == 1)) / 2){
    
    #check whether the variables are not too correlate
  }
}

cbind(data_model2$occurrence, vals)

#fix col names
names(vals2)[1] <- 'occurrence'

#prepare data object
data_obj <- sdmData(formmula = occurrence ~ . + coords(lon+lat),
                    train = vals2)

#run models
sdm <- sdm(occurrence ~ ., data = data_obj,
           methods = c('fda', 'glm', 'mars', 'maxlike',
                       'mda', 'gam', 'rf', 'svm'), 
               replication = 'cv', cv.folds = 5, n = 5)

#run models
sdm <- sdm(occurrence ~ ., data = data_obj,
           methods = c('glm', 'maxlike',
                       'gam', 'rf'), 
           replication = 'cv', cv.folds = 5, n = 5)

gui(sdm)

#save model objects
setwd(wd_models)
write.sdm(sdm, 'Epinephelus_aeneus')


#get model evaluation
eval <- getEvaluation(sdm, w = 1:100,
                          wtest='test.dep', 
                          stat=c('AUC','TSS','th'), opt = 2)


#include column informing algorithm
eval$Algorithm <- NA
eval$Algorithm[1:25] <- 'fda'
eval$Algorithm[26:50] <- 'glm'
eval$Algorithm[51:75] <- 'mars'
eval$Algorithm[76:100] <- 'maxlike'
eval$Algorithm[101:125] <- 'mda'
eval$Algorithm[126:150] <- 'gam'
eval$Algorithm[151:175] <- 'rf'
eval$Algorithm[176:200] <- 'svm'

#save evaluation metrics
setwd(wd_evaluation)
write.csv(eval, 'Eval_Epinephelus_aeneus.csv')


#list models with TSS higher than 0.5
sel_TSS <- which(eval$TSS >= 0.4)



###################
##### PRESENT #####
###################


#project models
setwd(wd_projections)

pred_pr <- list()
for(i in 1:100)
{
  pred_pr[[i]] <- predict(sdm, id = i, newdata = vars, 
              filename = paste0('Epinephelus_aeneus_', i, '_',
                                eval$Algorithm[i], '.grd'))
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
