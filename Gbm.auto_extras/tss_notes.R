#Run 5 no mangrove####
library(gbm.auto)
setwd("/media/Seagate/Work/PostDoc Work/Kroetz & Dedman Sawfish BRT/")
setwd("/home/simon/Documents/Si Work/PostDoc Work/Kroetz & Dedman Sawfish BRT/")
setwd("X:/PostDoc Work/Kroetz & Dedman Sawfish BRT")
samples <- read.csv("AndreaData2019.07.27.csv")
samples$Against.Shoreline <- as.factor(samples$Against.Shoreline)
samples$Shell <- as.factor(samples$Shell)
samples$Oyster <- as.factor(samples$Oyster)
samples$Seagrass <- as.factor(samples$Seagrass)
samples$Mesh.Inch <- as.factor(samples$Mesh.Inch)
samples$Net.Length.Ft <- as.factor(samples$Net.Length.Ft)
expvars5 <- c("Superregion", "Tidal.State", "Water.Temp.C", "Salinity",
              "Depth.M", "DO.MgL", "Grain.Size.MM.Log",
              "Month","Year")
gbm.auto(grids = NULL,
         samples = samples,
         expvar = expvars5,
         resvar = "P.Pectinata.CPUE",
         tc = length(expvars5),
         lr = 0.01,
         bf = 0.5,
         n.trees = 50,
         ZI = "CHECK",
         fam1 = "bernoulli",
         fam2 = "gaussian",
         simp = F,
         gridslat = 5,
         gridslon = 4,
         multiplot = F,
         cols = grey.colors(1, 1, 1),
         linesfiles = F,
         smooth = F, #changed
         savegbm = F,
         loadgbm = NULL,
         varint = F,
         map = TRUE,
         shape = NULL,
         RSB = F,
         BnW = F,
         alerts = TRUE,
         pngtype = "cairo-png",
         gaus = F)

# has to be done on the results of a prediction run or can be done on a model run only?
# evaluate() help allows for a model & predictor variables...

'e' gives you a number of things including AUC,
the true positive rate (TPR; known as sensitivity),
and the true negative rate (TNR; known as specificity).
TSS can be calculated from the TNR and TPR (last line of code)


DataInput.kfolds <- gbm.step(data=DataInput_train,
                             gbm.x= gbm.x,
                             gbm.y = gbm.y,
                             family="bernoulli",
                             tree.complexity=tc,
                             learning.rate = lr,
                             bag.fraction = 0.6)
preds <- predict.gbm(DataInput.kfolds, # same as gbm.predict.grids, model,
                     DataInput_test, # grids, new data
                     n.trees=DataInput.kfolds$gbm.call$best.trees, #same
                     type="response") #same
dev <- calc.deviance(obs=DataInput_test$PresAbs,
                     pred=preds,
                     calc.mean=TRUE)
d <- cbind(DataInput_test$PresAbs,
           preds)
pres <- d[d[,1] == 1, 2]
abs <- d[d[,1] == 0, 2]
e <- evaluate(p=pres,
              a=abs)
e
auc <- e@auc
tss <- max(e@TPR + e@TNR-1)



e <- evaluate(p = ,
              a = ,
              model = get(Bin_Best_Model))


# see also https://rdrr.io/cran/dismo/man/plotEval.html
plot(e, 'ROC')
plot(e, 'kappa')
plot(e, 'FPR')
plot(e, 'prevalence')


ModelEvaluation {dismo}	R Documentation
Class "ModelEvaluation"
Description
Class to store results of model cross-validation with presence/absence (0/1) data

Slots
presence:
  presence data used

absence:
  absence data used

np:
  number of presence points

na:
  number of absence points

auc:
  Area under the receiver operator (ROC) curve

pauc:
  p-value for the AUC (for the Wilcoxon test W statistic)

cor:
  Correlation coefficient

pcor:
  p-value for correlation coefficient

t:
  vector of thresholds used to compute confusion matrices

confusion:
  confusion matrices

prevalence:
  Prevalence

ODP:
  Overall diagnostic power

CCR:
  Correct classification rate

TPR:
  True positive rate

TNR:
  True negative rate

FPR:
  False positive rate

FNR:
  False negative rate

PPP:
  Positive predictive power

NPP:
  Negative predictive power

MCR:
  Misclassification rate

OR:
  Odds-ratio

kappa:
  Cohen's kappa
