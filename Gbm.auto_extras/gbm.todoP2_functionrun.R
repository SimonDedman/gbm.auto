## Run gbm.cons
library("devtools")
install_github("SimonDedman/gbm.auto") # update gbm.auto to latest
library("gbm.auto")
mygrids <- gbm.auto::grids
Juveniles <- gbm.auto::Juveniles
Adult_Females <- gbm.auto::Adult_Females

setwd("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Testbench")
# source('~/Dropbox/Galway/Analysis/R/gbm.auto/gbm.cons.R')
# Juveniles <- read.csv("Hauls&J&Preds&Enviros_Trimmed_ISonly_newdata_oldbkuporder&enviros&rays.csv")
# Adult_Females <- read.csv("F_Mat_plus_LPUE_plus_Enviro_IS_AllSp.csv") # single or vector of samples data csv files corresponding to subsets

gbm.cons(mygrids = mygrids, # csv file (/ & location) of gridded lat long data to predict to
         subsets = c("Juveniles","Adult_Females"), # Subset names
         alerts = TRUE,                # play sounds to mark progress steps
         map = TRUE,                   # produce maps
         BnW = TRUE,
         resvars = c(44:47,11:14),
         gbmautos = TRUE,
         expvars = list(c(4:11,15,17,21,25,29,37),
                        c(4:11,15,18,22,26,30,38),
                        c(4:11,15,19,23,27,31),
                        c(4:11,15,20,24,28,32,39),
                        4:10,
                        4:10,
                        4:10,
                        4:10),
         tcs = list(c(2,14),
                    c(2,14),
                    13,
                    c(2,14),
                    c(2,6),
                    c(2,6),
                    6,
                    c(2,6)),
         lrs = list(c(0.01,0.005),
                    c(0.01,0.005),
                    0.005,
                    c(0.01,0.005),
                    0.005,
                    0.005,
                    0.001,
                    0.005),
         ZIs = rep(TRUE,8),
         savegbms = rep(FALSE, 8),
         varints = rep(FALSE, 8),
         RSBs = rep(FALSE, 8),
         BnWs = rep(FALSE, 8),
         zeroes = rep(FALSE,8))


# without doing gmb.auto runs
gbm.cons(mygrids = subsets = c("Juveniles","Adult_Females"), # Subset names
         subsets = c("Juveniles","Adult Females"), # Subset names
         gbmautos = FALSE, resvars = c(44:47,11:14))
