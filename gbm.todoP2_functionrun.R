## Run gbm.cons
setwd("/home/simon/Dropbox/Galway/Project Sections/2. Spatial subsets inc fishery data/Data/Maps/Testbench")
source('~/Dropbox/Galway/Analysis/R/gbm.auto/gbm.cons.R')
gbm.cons(grids = "grids_Enviro_HansLPUE_MI&MMOlog_MIscallopVMS_MMOWhelk_MMOScal_Dist2Srvy_Preds_IS_NA_HansE.csv", # csv file (/ & location) of gridded lat long data to predict to
         subsets = c("Juveniles","Adult Females"), # Subset names
         conssamples = c("Hauls&J&Preds&Enviros_Trimmed_ISonly_newdata_oldbkuporder&enviros&rays.csv","F_Mat_plus_LPUE_plus_Enviro_IS_AllSp.csv"), # single or vector of samples data csv files corresponding to subsets
         alerts = TRUE,                # play sounds to mark progress steps
         map = TRUE,                   # produce maps
         BnW = TRUE,
         expvars = list(c(4:11,15,17,21,25,29,37),
                        c(4:11,15,18,22,26,30,38),
                        c(4:11,15,19,23,27,31),
                        c(4:11,15,20,24,28,32,39),
                        4:10,
                        4:10,
                        4:10,
                        4:10),
         resvars = c(44:47,11:14),
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
                    0.005,
                    0.005),
         ZIs = rep(TRUE,8),
         savegbms = rep(FALSE, 8),
         varints = rep(FALSE, 8),
         RSBs = rep(FALSE, 8),
         BnWs = rep(FALSE, 8),
         zeroes = rep(FALSE,8))





####TESTBENCH####
# lrs: NEED IF(ISNULL) solution for individual entries? Or have user manually enter gbm.auto defaults?

# subsets: Juveniles & ADult Females
(resvars)[1:(length(resvars)/length(subsets))] # 44:47 correct
length(resvars)/length(subsets) # =4 = currently "species"
names(samples)[11:14] # correct. remember samples = samples[2] = adult females.
names(samples)[(resvars)[1:(length(resvars)/length(subsets))]] # that works
# but only does 1st 4. What if resvar order of species is different in different subsets?
# resvars: 8
# subsets: 2
# group size: 4: resvars/subsets
# 1:GS = 1:4
# GS+1:GS*2
#
# 1+(GS*(n-1)):(GS*n) = 1+(4*0):4*1 = 1:4 (n=1)
# n=2: 1+(4*(1)):4*2 = 5:8
# So, loop through 1+(GS*(n-1)):(GS*n) for length(subsets) n's
# creates a list (?) of ranges for resvars
# "resvarrange", resvarrange[[listnumber]]
names(samples)[(resvars)[resvarrange[[i]]]]
names(samples)[(resvars)[resvarrange[[2]]]][2]

assign(paste("All_", k, sep = ""),rep(0,dim(mygrids)[1]))
assign(paste("All_", "Cuckoo", sep = ""),rep(0,dim(mygrids)[1]))
juves <- read.csv(paste("./Juveniles/Individual Predators/Cuckoo/Abundance_Preds_only.csv", sep = ""), header = TRUE)
matfs <- read.csv(paste("./Mature Females plus Hans' F/Cuckoo/Abundance_Preds_only.csv", sep = ""), header = TRUE)
assign(paste("All_", "Cuckoo", sep = ""),get(paste("All_", "Cuckoo", sep = "")) + juves[,3])
getwd()
length(mygrids)
dim(mygrids)[1]

# mysamples is i loop dependent.... right?
for (k in names(mysamples)[(resvars)[resvarrange[[length(subsets)]]]]) {  # name list from last mysamples & last subset resvarnames list,
  # Cuckoo Thornback Blonde Spotted
# before subsets loop, in k species loop, make grids-length blank object
assign(paste("All_", k, sep = ""),rep(0,dim(mygrids)[1])) # e.g. All_Cuckoo
# loop subsets
  for (i in 1:length(subsets)) {
    # replace All_Cuckoo w/ All_Cuckoo + Juveniles_Cuckoo
assign(paste("All_", k, sep = ""),get(paste("All_", k, sep = "")) + get(paste(subsets[i], "_", k, sep = ""))[,3])
 } # then reloop i & add Adult Females_Cuckoo
} # then reloop

tcs = list()
for (g in 1:length(resvars)) {tcs[[g]] <- c(2,length(expvars[[g]]))}
