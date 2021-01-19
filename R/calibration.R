#' calibration
#'
#' Internal use only. Jane Elith/John Leathwick 17th March 2005. Calculates
#' calibration statistics for either binomial or count data but the family
#' argument must be specified for the latter a conditional test for the latter
#' will catch most failures to specify the family.
#'
#' @param obs Observed data.
#' @param preds Predicted data.
#' @param family Statistical distribution family.
#'
#' @return roc & calibration stats internally within gbm runs e.g. in gbm.auto.
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#'
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
