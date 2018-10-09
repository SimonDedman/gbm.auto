#' Find the optimal parameter values for BRTs
#' (lowest mean deviance with n.trees > 1000)
#' Zack Oyafuso & Erik Franklin <=2017
#' Simon Dedman, 2016, simondedman@gmail.com,
#' github.com/SimonDedman/gbm.auto
#'
#' @param samples Explanatory and response variables to predict from. Keep col
#' names short, no odd characters, starting numerals or terminal periods. Spaces
#' may be converted to periods in directory names, underscores won't. Can be a
#' subset
#' @param expvar List of names or column numbers of explanatory variables in
#' 'samples': c(1,3,6) or c("Temp","Sal"). No default
#' @param resvar Name or column number of response variable in samples: 12,
#' "Rockfish". No default. Column name is ideally species name
#' @param tc Permutations of tree complexity to attempt, single or vector,
#' biggest number <= number of explanatory variables e.g. c(2,7)
#' @param lr Permutations of learning rate to attmpt. Single or vector.
#' @param bf Permutations of bag fraction to attempt. Single vector or "CHECK" to test
#' @param family Probability distribution family: bernoulli poisson laplace gaussian
#' @param runbest Option to run BRT with best param combo, default false
#'
#' @return Optionally runs best BRT and prints best parameter combination
#' @export
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @examples None
#' @import gbm
#' @import dismo
#' @import parallel
#' @import foreach
#' @import doParallel

gbm.tune <- function(
  lr = c(0.001, 0.005, 0.01), # Permutations of learning rate to attmpt. Single or vector.
  tc = 1:3, # Permutations of tree complexity to attempt, single or vector,
  # biggest number <= number of explanatory variables e.g. c(2,7)
  bf = "CHECK", # Permutations of bag fraction to attempt. Single vector or "CHECK" to test
  samples, # Explanatory and response variables to predict from. Keep col
  # names short, no odd characters, starting numerals or terminal periods. Spaces
  # may be converted to periods in directory names, underscores won't. Can be a
  # subset
  expvar, # List of names or column numbers of explanatory variables in samples
  # e.g. c(1,3,6) or c("Temp","Sal"). No default
  resvar, # Name or column number of response variable in samples: 12,
  # "Rockfish". No default. Column name is ideally species name
  family = "gaussian", #Probability distribution family: bernoulli poisson laplace gaussian
  runbest = FALSE){ #run BRT of best parameter combo? Default FALSE

  library(dismo); library(gbm); library(parallel); library(foreach); library(doParallel)

  #SD TODO####
  # currently tries preset permutations of tc lr bf Which is better than manually
  # tc could be based on n of vars?
  # what happens if a BRT doesn't run? If code can learn from crash & continue then
  # lr & bf could be set large enough to probably crash and iteratively shrink
  # change name of this function

  #Get minimum bf with gbm.bfcheck if requested####
  if(bf == "CHECK") if(family == "gaussian"){
    minbf <- gbm.bfcheck(samples = samples, resvar = resvar, ZI = "CHECK")[2]
  } else {
    minbf <- gbm.bfcheck(samples = samples, resvar = resvar, ZI = "CHECK")[1]
  }

  #set min bf based on bfcheck & untouchable ceiling of 1####
  if(minbf < 0.5) { #open 1st IF
    minbf <- 0.5 # default to 0.5 if acceptable bf is <0.5
  } else { #close 1st IF open 1st ELSE
    if(minbf > 0.9 && minbf < 1){ #if minbf is between 0.9 & 1 #open 2nd IF
      minbf <- min(ceiling(minbf*100) / 100, 0.99) # round up 0.1 or 0.99 never 1. Still probably won't work!
    } else { #open 2nd ELSE within 2nd IF
      minbf <- ceiling(minbf*10) / 10 # else round up to nearest 0.1
    } #close 2nd ELSE
  } # close 1st IF. minbf is single min acceptable bf else 0.5.

  #create vector of bf permutations####
  headroom <- 1-minbf #space between minbf and 1
  headincrement <- headroom/3 # thirds of that space
  bf <- c(minbf, # create a vector of 3 values, min bf + 1/3 & 2/3
          minbf + (headincrement * 1),
          minbf + (headincrement * 2))


  #Tune BRT settings####
  ## Create all combination of
  ## (1) learning rate (lr)
  ## (2) tree complexity (tc)
  ## (3) bagging rate (bf)
  tune_settings = expand.grid(lr, tc, bf)

  #Run BRT model per setting in tune_settings####
  #Register Cores for Parallel Processing####
  cl = makeCluster(detectCores())
  registerDoParallel(cl, cores = detectCores())
  #samples <- samples #delete if unrequired
  clusterExport(cl = cl, varlist = 'samples')

  #Partition BRTs to workers####
  x = foreach(i = 1:nrow(tune_settings), .combine = rbind) %dopar% {
    library(dismo); library(gbm)
    test = gbm.step(data = samples,
                    gbm.x = expvar,
                    gbm.y = resvar, #response
                    family = family,
                    tree.complexity = tune_settings$tc[i],
                    learning.rate = tune_settings$lr[i],
                    bag.fraction = tune_settings$bf[i])

    print( c(test$n.trees, test$cv.statistics$deviance.mean) )
  }
  stopCluster(cl)

  #Append n.trees mean dev to tune_settings####
  tune_settings[,c('n.trees', 'dev_mean')] = x

  #Select best settings####
  # The "optimal" settings is the combination of bf tc & lr that
  # 1: produces > 1000 trees (as recommended by Elith et al. 2008), and
  # 2: has the lowest mean deviance
  if(min(tune_settings$n.trees < 1000)) { #if no treees>1000 produces
    best_settings = tune_settings[which.min(tune_settings$dev_mean),] #make best settings anyway
    print("No BRT combinations produced >= 1000 trees") #let the user know
  } else { #else subset for only trees > 1000
    accept_mods = subset(tune_settings, n.trees > 1000)
    best_settings = accept_mods[which.min(accept_mods$dev_mean),]
  }

  #Fit BRT with best settings####
  if(runbest){ #if runbest BRT params selected, do so
    brt_opt =  gbm.step(data = samples,
                        gbm.x = expvar, #covariates
                        gbm.y = resvar,
                        family = "gaussian",
                        tree.complexity = best_settings$tc,
                        learning.rate = best_settings$lr,
                        bag.fraction = best_settings$bf)
  }
  return(best_settings) #return best settings combo to user
} # close function
