#' gbm.utils: roc, calibration & gbm.predict.grids functions bundle
#'
#' [Note: Man page named after 1st function in script i.e. roc but script
#' contains 3. Should be no need for users to interact with these directly]
#'
#' roc: adapted from Ferrier, Pearce and Watson's code, by J.Elith, see: Hanley,
#' J.A. & McNeil, B.J. (1982) The meaning and use of the area under a Receiver
#' Operating Characteristic (ROC) curve. Radiology, 143, 29-36. Also Pearce, J.
#' & Ferrier, S. (2000) Evaluating the predictive performance of habitat models
#' developed using logistic regression. Ecological Modelling, 133, 225-245. This
#' is the non-parametric calculation for area under the ROC curve, using the
#' fact that a MannWhitney U statistic is closely related to the area. In dismo,
#' this is used in the gbm routines, but not elsewhere (see evaluate).
#'
#' calibration: j elith/j leathwick 17th March 2005. Calculates calibration
#' statistics for either binomial or count data but the family argument must be
#' specified for the latter a conditional test for the latter will catch most
#' failures to specify the family
#'
#' Gbm.predict.grids: J.Elith / J.Leathwick, March 07. To make predictions to
#' sites or grids. If to sites, the predictions are written to the R workspace.
#' If to grid, the grids are written to a nominated directory and optionally
#' also plotted in R. New data (new.dat) must be a data frame with column names
#' identical to names for all variables in the model used for prediction.
#' pred.vec is a vector of -9999's, the length of the scanned full grid (i.e.
#' without nodata values excluded).filepath must specify the whole path as a
#' character vector,but without the final file name - eg "c:/gbm/"
#'
#' @param obsdat Observed data
#' @param preddat Predicted data
#'
#' @return roc & calibration stats internally within gbm runs e.g. in gbm.auto; gbm.predict.grids powers the predictive mapping element of gbm.map
#' @importFrom gbm predict.gbm
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @examples None
roc <- function(obsdat, preddat) {
# code adapted from Ferrier, Pearce and Watson's code, by J.Elith
# see: Hanley, J.A. & McNeil, B.J. (1982) The meaning and use of the area
# under a Receiver Operating Characteristic (ROC) curve. Radiology, 143, 29-36
#
# Pearce, J. & Ferrier, S. (2000) Evaluating the predictive performance of habitat models developed using logistic regression.
# Ecological Modelling, 133, 225-245.
# this is the non-parametric calculation for area under the ROC curve, using the fact that a MannWhitney U statistic is closely related to
# the area. In dismo, this is used in the gbm routines, but not elsewhere (see evaluate).

    if (length(obsdat) != length(preddat))
        stop("obs and preds must be equal lengths")
    n.x <- length(obsdat[obsdat == 0])
    n.y <- length(obsdat[obsdat == 1])
    xy <- c(preddat[obsdat == 0], preddat[obsdat == 1])
    rnk <- rank(xy)
    wilc <- ((n.x * n.y) + ((n.x * (n.x + 1))/2) - sum(rnk[1:n.x]))/(n.x * n.y)
    return(round(wilc, 4))
}


calibration <- function(obs, preds, family = "binomial")  {
# j elith/j leathwick 17th March 2005
# calculates calibration statistics for either binomial or count data but the family argument must be specified for the latter
# a conditional test for the latter will catch most failures to specify the family

if (family == "bernoulli") family <- "binomial"
pred.range <- max(preds) - min(preds)
if (pred.range > 1.2 & family == "binomial") {
print(paste0("range of response variable is ", round(pred.range, 2)), quote = F)
print("check family specification", quote = F)
return()
}
if (family == "binomial") {
pred <- preds + 1e-005
pred[pred >= 1] <- 0.99999
mod <- glm(obs ~ log((pred)/(1 - (pred))), family = binomial)
lp <- log((pred)/(1 - (pred)))
a0b1 <- glm(obs ~ offset(lp) - 1, family = binomial)
miller1 <- 1 - pchisq(a0b1$deviance - mod$deviance, 2)
ab1 <- glm(obs ~ offset(lp), family = binomial)
miller2 <- 1 - pchisq(a0b1$deviance - ab1$deviance, 1)
miller3 <- 1 - pchisq(ab1$deviance - mod$deviance, 1)
}
if (family == "poisson") {
mod <- glm(obs ~ log(preds), family = poisson)
lp <- log(preds)
a0b1 <- glm(obs ~ offset(lp) - 1, family = poisson)
miller1 <- 1 - pchisq(a0b1$deviance - mod$deviance, 2)
ab1 <- glm(obs ~ offset(lp), family = poisson)
miller2 <- 1 - pchisq(a0b1$deviance - ab1$deviance, 1)
miller3 <- 1 - pchisq(ab1$deviance - mod$deviance, 1)
}
calibration.result <- c(mod$coef, miller1, miller2, miller3)
names(calibration.result) <- c("intercept", "slope", "testa0b1", "testa0|b1", "testb1|a")
return(calibration.result)
}


gbm.predict.grids <- function(model, new.dat, want.grids = F, preds2R = T, sp.name = "preds", pred.vec = NULL, filepath = NULL,
                               num.col = NULL, num.row = NULL, xll = NULL, yll = NULL, cell.size = NULL, no.data = NULL, plot = F,
                               full.grid = T, part.number = NULL, part.row = NULL, header = T)
{
# J.Elith / J.Leathwick, March 07
# to make predictions to sites or grids. If to sites, the predictions are written to the R workspace. If to grid,
# the grids are written to a nominated directory and optionally also plotted in R
#
# new data (new.dat) must be a data frame with column names identical to names for all variables in the model used for prediction.
# pred.vec is a vector of -9999's, the length of the scanned full grid (i.e. without nodata values excluded).
# filepath must specify the whole path as a character vector,but without the final file name - eg "c:/gbm/"

temp <- predict.gbm(model, new.dat, n.trees = model$gbm.call$best.trees, type = "response")

if (want.grids)
{
newname <- paste0(filepath, sp.name,".asc")
full.pred <- pred.vec
full.pred[as.numeric(row.names(new.dat))] <- temp
if (header) {
write(paste0("ncols          ",num.col),newname)
write(paste0("nrows          ",num.row),newname,append = T)
write(paste0("xllcorner      ",xll),newname,append = T)
write(paste0("yllcorner      ",yll),newname,append = T)
write(paste0("cellsize       ",cell.size),newname,append = T)
write(paste0("NODATA_value ",no.data),newname,append = T)
}
  if (full.grid) {
         full.pred.mat <- matrix(full.pred, nrow = num.row, ncol = num.col, byrow = T)
	if (plot)
	{
	image(z = t(full.pred.mat)[, nrow(full.pred.mat):1], zlim =  c(0,1), col = rev(topo.colors(12)))
	}
	write.table(full.pred.mat, newname, sep = " ", append = T, row.names = F, col.names = F)
	#also write to R directory, if required:
	if (preds2R) {assign(sp.name, temp, pos = 1)}
	}
	else{
         full.pred.mat <- matrix(full.pred, nrow = part.row, ncol = num.col, byrow = T)
	write.table(full.pred.mat, newname, sep = " ", append = T, row.names = F, col.names = F)
	if (preds2R) {assign(paste0(sp.name, part.number), temp, pos = 1)}
	}
}
else{
assign(sp.name, temp, pos = 1)
}
}
