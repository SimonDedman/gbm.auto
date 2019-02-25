# parallel computing####
# https://www.r-bloggers.com/parallel-r-loops-for-windows-and-linux/
# https://github.com/tobigithub/R-parallel/wiki/R-parallel-Errors
# http://www.vikparuchuri.com/blog/monitoring-progress-inside-foreach-loop/ 
#windows####
library(foreach)
library(doSNOW)
cl<-makeCluster(4) #number of CPU cores. win=4? Should be 8?
registerDoSNOW(cl)

foreach(i = c("F.Atlantic.sharpnose.shark",
              "M.Atlantic.sharpnose.shark",
              "F.blacknose.shark",
              "M.blacknose.shark",
              "F.blacktip.shark",
              "M.blacktip.shark",
              "F.bull.shark",
              "M.bull.shark",
              "F.sandbar.shark",
              "M.sandbar.shark",
              "M.scalloped.hammerhead",
              "F.smoothhound.sharks",
              "F.Southern.stingray",
              "F.spinner.shark",
              "M.spinner.shark")) %dopar% {
                
                gbm.auto(expvar = expvars,
                         resvar = i,
                         samples = sharks,
                         grids = gridsspring,
                         shape = mobilebay,
                         lr = 0.0005,
                         bf = 0.8,
                         gaus = TRUE,
                         savegbm = FALSE,
                         BnW = FALSE,
                         simp = TRUE)
              }
stopCluster(cl)

#linux####
library(foreach)
library(doMC)
registerDoMC(16)  #number of CPU cores. linux=16
writeLines(c(""), "log.txt") #sends r output to a file. "" blanks file first

foreach(i = c("F.Atlantic.sharpnose.shark",
              "M.Atlantic.sharpnose.shark",
              "F.blacknose.shark",
              "M.blacknose.shark",
              "F.blacktip.shark",
              "M.blacktip.shark",
              "F.bull.shark",
              "M.bull.shark",
              "F.sandbar.shark",
              "M.sandbar.shark",
              "M.scalloped.hammerhead",
              "F.smoothhound.sharks",
              "F.Southern.stingray",
              "F.spinner.shark",
              "M.spinner.shark")) %dopar% {
                
                sink("log.txt", append = TRUE) #sends otherwise hidden messages to log file which can be opened at any time
                
                gbm.auto(expvar = expvars,
                         resvar = i,
                         samples = sharks,
                         grids = gridsspring,
                         shape = mobilebay,
                         lr = 0.0005,
                         bf = 0.8,
                         gaus = TRUE,
                         savegbm = FALSE,
                         BnW = FALSE,
                         simp = TRUE)
              }