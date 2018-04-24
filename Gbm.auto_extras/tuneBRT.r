# tune BRT code to find the optimal parameter values (i.e, lowest mean deviance with n.trees > 1000).
# Zack Oyafuso & Erik Franklin.

##############################
## Tune BRT settings
## Create all combination of
## (1) learning rate (lrt)
## (2) tree complexity (tc)
## (3) bagging rate (bag_rate)
##############################
tune_settings = expand.grid(lrt = c(0.001, 0.005, 0.01), tc = 1:3, bag_rate = c(0.50, 0.75))

###########################################
## Run BRT model under each setting in the tune_settings df
###########################################

############################################
## Register Cores for Parallel Processing
#############################################
cl = makeCluster(detectCores())
registerDoParallel(cl, cores = detectCores())
clusterExport(cl = cl, varlist = 'fishCompDF')

############################################
## Partition BRT work to workers
############################################
x = foreach(i = 1:nrow(tune_settings), .combine = rbind) %dopar% {
  library(dismo); library(gbm)
  test = gbm.step(data = fishCompDF,
                  gbm.x = predictorVars,
                  gbm.y = responseVar, #response
                  family = "gaussian",
                  tree.complexity = tune_settings$tc[i],
                  learning.rate = tune_settings$lrt[i],
                  bag.fraction = tune_settings$bag_rate[i])

  print( c(test$n.trees, test$cv.statistics$deviance.mean) )
}
stopCluster(cl)

#Append number of trees and mean deviance to tune_settings df
tune_settings[,c('n_trees', 'dev_mean')] = x

##############################################
## The "optimal" settings is the combination
## of bagging rate, tree complexity, and learning rate
## that (1) produces > 1000 trees (as reccommended by Elith et al. 2008)
## and (2) has the lowest mean deviance
##############################################

accept_mods = subset(tune_settings, n_trees > 1000)
best_settings = accept_mods[which.min(accept_mods$dev_mean),]

##############################################
## Fit BRT model with the "optimal" settings
##############################################
brt_opt =  gbm.step(data = fishCompDF,
                    gbm.x = predictorVars, #covariates
                    gbm.y = responseVar,
                    family = "gaussian",
                    tree.complexity = best_settings$tc,
                    learning.rate = best_settings$lrt,
                    bag.fraction = best_settings$bag_rate)
