#' plotMolecules
#' 
#' Plotting functions for spatially resolved transcriptomics data.
#' 
#' Function to plot molecule-based spatially resolved transcriptomics data
#' stored in a \code{SpatialExperiment} object.
#' 
#' This function generates a plot in spatial coordinates (e.g. x-y coordinates
#' on a tissue slide), for a selected molecule.
#' 
#' 
#' @param spe (SpatialExperiment) Input data, assumed to be a
#'   \code{SpatialExperiment} object.
#' 
#' @param molecule (character) Name of mRNA molecule to plot (assumed to match
#'   one of the row names of \code{rowData}).
#' 
#' @param x_coord (character) Name of column in \code{spatialCoords} containing
#'   x-coordinates of the cell centroids. Default = NULL, which selects the
#'   first column.
#' 
#' @param y_coord (character) Name of column in \code{spatialCoords} containing
#'   y-coordinates of the cell centroids. Default = NULL, which selects the
#'   second column.
#' 
#' @param palette (character) Color palette, provided as a vector of length 2
#'   for the low and high range. Default = \code{c("gray90", "navy")}.
#' 
#' @param size (numeric) Point size for \code{geom_point()}. Default = 0.3.
#' 
#' 
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, labels, formatting, etc).
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment counts
#' @importFrom ggplot2 ggplot aes_string geom_point scale_color_gradient
#'   coord_fixed ggtitle theme_void
#' 
#' @export
#' 
#' @examples
#' library(STexampleData)
#' spe <- seqFISH_mouseEmbryo()
#' plotMolecules(spe, molecule = "Sox2")
#' 
plotMolecules <- function(spe, 
                          molecule =  NULL, 
                          x_coord = NULL, y_coord = NULL, 
                          palette = c("gray90", "navy"), 
                          size = 0.3) {
  
  if (!is.null(molecule)) stopifnot(is.character(molecule))
  
  if (is.null(x_coord)) x_coord <- colnames(spatialCoords(spe))[1]
  if (is.null(y_coord)) y_coord <- colnames(spatialCoords(spe))[2]
  
  mRNA_counts <- as.numeric(counts(spe)[molecule, ])
  stopifnot(length(mRNA_counts) == ncol(spe))
  
  # providing a single value e.g. "navy" will create a vector c("gray95", "navy")
  palette <- .get_pal(palette)
  
  df <- cbind.data.frame(spatialCoords(spe), sum = mRNA_counts)
  
  p <- ggplot(df, aes_string(x = x_coord, y = y_coord, color = "sum")) + 
    geom_point(size = size) + 
    scale_color_gradient(low = palette[1], high = palette[2], trans = "sqrt") + 
    coord_fixed() + 
    ggtitle(molecule) + 
    theme_void()
  
  p
}

