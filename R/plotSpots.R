#' plotSpots
#' 
#' Plots for spatially resolved transcriptomics datasets
#' 
#' Function to plot spot-based spatially resolved transcriptomics data in
#' spatial (x-y) coordinates.
#' 
#' This function generates a plot showing spatial coordinates (spots) in the x-y
#' coordinates of the tissue slide, with optional colors for cluster labels,
#' ground truth labels, or other quantities.
#' 
#' 
#' @param spe (SpatialExperiment) Input data object.
#' 
#' @param x_coord (character) Name of column in spatialData containing
#'   x-coordinates. Default = "x".
#' 
#' @param y_coord (character) Name of column in spatialData containing
#'   y-coordinates. Default = "y".
#' 
#' @param in_tissue (logical) Whether to plot only spots over tissue (TRUE) or
#'   all spots (FALSE). Default = TRUE.
#' 
#' @param discrete (character) Name of column in colData containing discrete
#'   labels (e.g. cluster labels or ground truth labels) to show with colors.
#'   Default = NULL.
#' 
#' @param continuous (character) Name of column in colData containing continuous
#'   values (e.g. total UMI counts) to show with colors. Default = NULL.
#' 
#' @param palette (character) Color palette. Options for discrete labels are
#'   "libd_layer_colors", "Okabe-Ito", or a vector of hex codes for a custom
#'   palette. Default = "libd_layer_colors". Options for continuous values are
#'   "navy", or a vector of length two containing custom colors. Default =
#'   "navy".
#' 
#' @param flip_xy_Visium (logical) Whether to switch x and y coordinates and
#'   reverse y coordinates (sometimes required for 10x Visium platform if using
#'   the default `pxl_row_in_fullres` and `pxl_col_in_fullres` columns as
#'   coordinates). Default = FALSE.
#' 
#' 
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, customized formatting, etc).
#' 
#' 
#' @importFrom rlang sym "!!"
#' @importFrom SpatialExperiment spatialData
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot aes geom_point coord_fixed ggtitle
#'   scale_color_manual scale_color_gradient theme_bw theme element_blank
#' 
#' @export
#' 
#' @examples
#' # library(ggspavis)
#' # library(STdata)
#' # spe <- load_data("Visium_humanDLPFC")
#' # plotSpots(spe, discrete = "ground_truth", palette = "libd_layer_colors")
#' 
plotSpots <- function(spe, 
                      x_coord = "x", y_coord = "y", in_tissue = TRUE, 
                      discrete = NULL, continuous = NULL, 
                      palette = NULL, flip_xy_Visium = FALSE) {
  
  # note: using quasiquotation to allow custom variable names in ggplot ("sym" and "!!")
  
  x_coord_sym <- sym(x_coord)
  y_coord_sym <- sym(y_coord)
  if (!is.null(discrete)) discrete_sym <- sym(discrete)
  if (!is.null(continuous)) continuous_sym <- sym(continuous)
  
  if (!is.null(palette) && palette == "libd_layer_colors") {
    palette <- c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", "#FFD700", "#FF7F00", "#1A1A1A", "#666666")
  } else if (!is.null(palette) && palette == "Okabe-Ito") {
    palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  } else if (!is.null(palette) && palette == "navy") {
    palette <- c("gray95", "navy")
  }
  
  if (in_tissue) {
    spe <- spe[, spatialData(spe)$in_tissue == 1]
  }
  
  if (flip_xy_Visium) {
    x_coord_tmp <- spatialData(spe)[, y_coord]
    y_coord_tmp <- -1 * spatialData(spe)[, x_coord] + max(spatialData(spe)[, x_coord]) + 1
    df_plot <- data.frame(x_coord = x_coord_tmp, y_coord = y_coord_tmp)
  }
  
  df_plot <- spatialData(spe)[, c(x_coord, y_coord), drop = FALSE]
  
  if (!is.null(discrete)) {
    df_plot <- cbind(df_plot, colData(spe)[, discrete, drop = FALSE])
  }
  if (!is.null(continuous)) {
    df_plot <- cbind(df_plot, colData(spe)[, continuous, drop = FALSE])
  }
  
  df_plot <- as.data.frame(df_plot)
  
  p <- ggplot(df_plot, aes(x = !!x_coord_sym, y = !!y_coord_sym)) + 
    geom_point(size = 0.5) + 
    coord_fixed() + 
    ggtitle("Spatial coordinates") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  if (!is.null(discrete)) {
    p <- p + aes(color = !!discrete_sym) + 
      scale_color_manual(values = palette)
  }
  
  if (!is.null(continuous)) {
    p <- p + aes(color = !!continuous_sym) + 
      scale_color_gradient(low = palette[1], high = palette[2])
  }
  
  p
  
}

