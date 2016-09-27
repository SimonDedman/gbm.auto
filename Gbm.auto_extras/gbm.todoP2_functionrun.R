## Run gbm.cons
setwd("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Testbench")
source('~/Dropbox/Galway/Analysis/R/gbm.auto/gbm.cons.R')
gbm.cons(grids = "grids_Enviro_HansLPUE_MI&MMOlog_MIscallopVMS_MMOWhelk_MMOScal_Dist2Srvy_Preds_IS_NA_HansE.csv", # csv file (/ & location) of gridded lat long data to predict to
         subsets = c("Juveniles","Adult Females"), # Subset names
         conssamples = c("Hauls&J&Preds&Enviros_Trimmed_ISonly_newdata_oldbkuporder&enviros&rays.csv","F_Mat_plus_LPUE_plus_Enviro_IS_AllSp.csv"), # single or vector of samples data csv files corresponding to subsets
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


# without going gmb.auto runs
gbm.cons(grids = "grids_Enviro_HansLPUE_MI&MMOlog_MIscallopVMS_MMOWhelk_MMOScal_Dist2Srvy_Preds_IS_NA_HansE.csv", # csv file (/ & location) of gridded lat long data to predict to
         subsets = c("Juveniles","Adult Females"), # Subset names
         conssamples = c("Hauls&J&Preds&Enviros_Trimmed_ISonly_newdata_oldbkuporder&enviros&rays.csv","F_Mat_plus_LPUE_plus_Enviro_IS_AllSp.csv"), # single or vector of samples data csv files corresponding to subsets
         alerts = TRUE,                # play sounds to mark progress steps
         gbmautos = FALSE,
         resvars = c(44:47,11:14))