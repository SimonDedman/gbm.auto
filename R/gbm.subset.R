#' Subset data into two groups
#'
#' Returns the variable value corresponding to the 0 value on the lineplots,
#' which should be the optimal place to split the dataset into 2 subsets.
#'
#'  loop varnames are BinLineLoop_VAR.csv & GausLineLoop_VAR.csv
#'  normal varnames are Bin_Best_line_VAR.csv & Gaus_Best_line_VAR.csv
#'
#'  just use average between the last negative & first positive point
#'  unless any points fall on zero

#' @param x Vector of variable names.
#' @param fams statistical data distribution family names to be modelled by gbm.
#' @param loop Is the folder a gbm.loop output?
#'
#' @export
#' @importFrom utils read.csv
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @examples None
#'
gbm.subset <- function(x, #Vector of variable names
                       fams = c("Bin", "Gaus"), #family names modelled by gbm
                       loop = FALSE){ #is the folder a gbm.loop output?
  subsetsplits <- list() #create blank list object
  for (j in fams) { #loop through families
    for (i in x) { #loop through variable names' files
      if (loop) {if (file.exists(paste0(j, "LineLoop_", i, ".csv"))) { #if file exists
        tmp <- read.csv(paste0(j, "LineLoop_", i, ".csv")) #read in csv file
        tmp <- tmp[,c(1, length(tmp) - 2)]} #keep only X & averageY
      } else {#if not loop
        if (file.exists(paste0(j, "_Best_line_", i, ".csv"))) {#if file exists
          tmp <- read.csv(paste0(j, "_Best_line_", i, ".csv"))}} #read in csv file

      if (exists("tmp")) { #if csv file was read (x names used in gbm may not have generated files)
        if (!is.na(match(0,sign(tmp[,2])))) { # if there's an exact 0 value in the Y column,
          subsetsplits[[i]] <- tmp[match(0,sign(tmp[,2])),1] #set the corresponding X value as a list item named i
        } else {# if there isn't (normal)
          row1 <- which(diff(sign(tmp[,2])) != 0) #gives the last row before crossing the Y=0 intercept
          row2 <- row1 + 1 #first row after
          subsetsplits[[i]] <- mean(c(tmp[row1,1], #set list value i as average of the points
                                      tmp[row2,1])) #before & after intercept crossing
        } #close 'point on 0 line' ifelse
        rm(tmp) #remove tmp if it was created
      } #close if exists tmp
    } #close var.names loop
  } #close fam.names loop
  return(subsetsplits) #return object
} #close function
