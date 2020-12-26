#' plotMolecules
#' 
#' Plots for spatially resolved transcriptomics datasets
#' 
#' Function to plot molecule-based spatially resolved transcriptomics data in
#' spatial (x-y) coordinates.
#' 
#' This function generates a plot showing counts for a given molecule in the x-y
#' coordinates of the tissue slide.
#' 
#' 
#' @param spe (SpatialExperiment) Input data object.
#' 
#' @param molecule (character) Name of mRNA molecule to plot (matching to one of
#'   the row names of 'rowData').
#' 
#' @param x_coord (character) Name of column in 'spatialData' containing
#'   x-coordinates of the cell centroids. Default = 'x_coord'.
#' 
#' @param y_coord (character) Name of column in 'spatialData' containing
#'   y-coordinates of the cell centroids. Default = 'y_coord'.
#' 
#' @param palette (character) Color palette for points. Options are a single
#'   color name (e.g. 'red', 'navy', etc), or a vector of length two containing
#'   color names for each end of the scale. Default = 'navy'.
#' 
#' 
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, customized formatting, etc).
#' 
#' 
#' @importFrom SpatialExperiment spatialData
#' @importFrom SingleCellExperiment colData counts
#' @importFrom ggplot2 ggplot aes_string geom_point scale_color_gradient
#'   coord_fixed theme_void ggtitle
#' @importFrom methods as
#' 
#' @export
#' 
#' @examples
#' # library(ggspavis)
#' # library(STdata)
#' # spe <- load_data("seqFISH_mouseEmbryo")
#' # plotMolecules(spe, molecule = "Sox2")
#' 
plotMolecules <- function(spe, 
                          molecule =  NULL, 
                          x_coord = "x_coord", y_coord = "y_coord", 
                          palette = NULL) {
  
  # set up color palette
  if (is.null(palette)) {
    palette <- "navy"
  }
  if (!is.null(palette) && length(palette) == 1) {
    # if providing a single color name (e.g. 'navy' or 'red'), combine with gray
    # for color scale; else if length(palette) > 1, use palette as provided
    # (i.e. expecting a vector of length 2)
    palette <- c("gray90", palette)
  }
  
  df_plot <- spatialData(spe)[, c(x_coord, y_coord), drop = FALSE]
  mRNA_counts <- as.numeric(counts(spe)[molecule, , drop = FALSE])
  stopifnot(length(mRNA_counts) == nrow(df_plot))
  df_plot <- cbind(df_plot, sum = mRNA_counts)
  
  df_plot <- as.data.frame(df_plot)
  
  p <- ggplot(df_plot, aes_string(x = "x_coord", y = "y_coord", color = "sum")) + 
    geom_point(size = 0.5) + 
    scale_color_gradient(low = palette[1], high = palette[2], trans = "sqrt") + 
    coord_fixed() + 
    ggtitle(molecule) + 
    theme_void()
  
  p
}

