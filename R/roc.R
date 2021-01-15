#' roc
#'
#' Internal use only. Adapted from Ferrier, Pearce and Watson's code, by J.Elith
#' , see: Hanley, J.A. & McNeil, B.J. (1982) The meaning and use of the area
#' under a Receiver Operating Characteristic (ROC) curve. Radiology, 143, 29-36.
#' Also Pearce, J. & Ferrier, S. (2000) Evaluating the predictive performance of
#'  habitat models developed using logistic regression. Ecological Modelling,
#'  133, 225-245. This is the non-parametric calculation for area under the ROC
#'  curve, using the fact that a MannWhitney U statistic is closely related to
#'  the area. In dismo, this is used in the gbm routines, but not elsewhere (see
#'  evaluate).
#'
#' @param obsdat Observed data.
#' @param preddat Predicted data.
#'
#' @return roc & calibration stats internally within gbm runs e.g. in gbm.auto.
#' @importFrom gbm predict.gbm
#' @importFrom grDevices topo.colors
#' @importFrom graphics image
#' @importFrom stats binomial glm pchisq poisson
#' @importFrom utils write.table
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @author Jane Elith
#' @author John Leathwick
#'
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
