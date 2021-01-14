#' Calculates minimum Bag Fraction size for gbm.auto
#'
#' Provides minimum bag fractions for gbm.auto, preventing failure
#' due to bf & samples rows limit. Simon Dedman, 2016, simondedman@gmail.com,
#' github.com/SimonDedman/gbm.auto
#'
#' @param samples Samples dataset, same as gbm.auto.
#' @param resvar Response variable column in samples.
#' @param ZI Are samples zero-inflated? TRUE/FALSE/"CHECK".
#' @param grv Dummy param for package testing for CRAN, ignore.
#'
#' @return Prints minimum Bag Fraction size for gbm.auto.
#' @export
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @examples
#' data(samples)
#' gbm.bfcheck(samples = samples, resvar = "Cuckoo")
gbm.bfcheck <- function(
  samples, # samples dataset, same as gbm.auto
  resvar, # response variable column in samples
  ZI = "CHECK") # are samples zero-inflated? TRUE/FALSE/"CHECK"
  # grv = NULL) # addresses devtools::check's no visible binding for global variable https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html#globals
{

# gbm.bfcheck: provides minimum bag fractions for gbm.auto,
# preventing failure due to bf & samples rows limit
# Simon Dedman, 2016, simondedman@gmail.com, github.com/SimonDedman/gbm.auto

utils::globalVariables("grv") # addresses devtools::check's no visible binding for global variable https://cran.r-project.org/web/packages/data.table/vignettes/datatable-importing.html#globals

# Minimum bag fraction for binary
minbfbin <- 21/nrow(samples)
print(
  paste0("  binary bag fraction must be at least ",
         round(minbfbin, 3),
         ". n = ",
         nrow(samples)
         )
      )

# if user has asked code to check for ZI, check it & set new ZI status
if (ZI == "CHECK") if (sum(samples[,resvar] == 0, na.rm = TRUE) / length(samples[,resvar]) >= 0.5) ZI = TRUE else ZI = FALSE
logem <- log1p(samples[,resvar])
dont  <- samples[,resvar]
if (ZI) {samples$grv <- logem} else {samples$grv <- dont}
grv_yes <- subset(samples, grv > 0) # nonzero subset for gaussian BRTs

# Minimum bag fraction for Gaussian
minbfgaus <- 21/nrow(grv_yes)
print(
  paste0("Gaussian bag fraction must be at least ",
         round(minbfgaus,3),
         ". n = ",
         nrow(grv_yes)
         )
      )

# return bin & gaus min bfs
return(c(minbfbin,minbfgaus))
}
