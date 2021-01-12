#' gbm.predict.grids
#'
#' Internal use only. J.Elith / J.Leathwick, March 07. To make predictions to
#' sites or grids. If to sites, the predictions are written to the R workspace.
#' If to grid, the grids are written to a nominated directory and optionally
#' also plotted in R. New data (new.dat) must be a data frame with column names
#' identical to names for all variables in the model used for prediction.
#' pred.vec is a vector of -9999's, the length of the scanned full grid (i.e.
#' without nodata values excluded).filepath must specify the whole path as a
#' character vector,but without the final file name - eg "c:/gbm/"
#' @param model Gbm model object.
#' @param new.dat New data to predict to. Must be a data frame with column names
#'  identical to names for all variables in the model used for prediction.
#' @param want.grids If used, grids are written to a nominated directory and
#' optionally also plotted in R.
#' @param preds2R Write predictions to R directory.
#' @param sp.name Species prediction name
#' @param pred.vec A vector of -9999's, the length of the scanned full grid
#' (i.e. without nodata values excluded).
#' @param filepath Must specify the whole path as a character vector,but without
#'  the final file name - eg "c:/gbm/".
#' @param num.col Number of columns for grids.
#' @param num.row Number of rows for grids.
#' @param xll X latitude longitude corner.
#' @param yll Y latitude longitude corner.
#' @param cell.size Cell size.
#' @param no.data No data value.
#' @param plot Plot output.
#' @param full.grid Predict to full grid else to part grid.
#' @param part.number Number of grid subset part.
#' @param part.row Number of rows of grid subset part.
#' @param header Want.grids header included.
#'
#' @return Powers the predictive mapping element of gbm.map.
#' @author Simon Dedman, \email{simondedman@@gmail.com}
#' @author Jane Elith
#' @author John Leathwick
#' @examples None
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

  if (want.grids) {
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
      if (plot) {image(z = t(full.pred.mat)[, nrow(full.pred.mat):1], zlim =  c(0,1), col = rev(topo.colors(12)))}
      write.table(full.pred.mat, newname, sep = " ", append = T, row.names = F, col.names = F)
      #also write to R directory, if required:
      if (preds2R) {assign(sp.name, temp, pos = 1)}
    } else {
      full.pred.mat <- matrix(full.pred, nrow = part.row, ncol = num.col, byrow = T)
      write.table(full.pred.mat, newname, sep = " ", append = T, row.names = F, col.names = F)
      if (preds2R) {assign(paste0(sp.name, part.number), temp, pos = 1)}
    }
  } else {
    assign(sp.name, temp, pos = 1)
  }
}
