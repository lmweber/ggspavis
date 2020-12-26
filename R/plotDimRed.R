#' plotDimRed
#' 
#' Plots for spatially resolved transcriptomics datasets
#' 
#' Function to plot spatially resolved transcriptomics data in reduced dimension
#' coordinates.
#' 
#' This function generates a plot showing spatial coordinates (spots) in reduced
#' dimension coordinates (PCA or UMAP), with optional colors for cluster labels,
#' ground truth labels, or other quantities.
#' 
#' 
#' @param spe (SpatialExperiment) Input data object.
#' 
#' @param type (character) Type of reduced dimension plot. Options are "PCA" or
#'   "UMAP". Default = "UMAP.
#' 
#' @param x_axis (character) Name of column in reducedDim slot containing
#'   x-coordinates. Default = "UMAP1".
#' 
#' @param y_axis (character) Name of column in reducedDim slot containing
#'   y-coordinates. Default = "UMAP2".
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
#' 
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, customized formatting, etc).
#' 
#' 
#' @importFrom rlang sym "!!"
#' @importFrom SingleCellExperiment reducedDim colData
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab ggtitle
#'   scale_color_manual scale_color_gradient theme_bw theme element_blank
#' 
#' @export
#' 
#' @examples
#' # to do
#' 
plotDimRed <- function(spe, type = "UMAP", 
                       x_axis = "UMAP1", y_axis = "UMAP2", 
                       discrete = NULL, continuous = NULL, 
                       palette = NULL) {
  
  # note: using quasiquotation to allow custom variable names in ggplot ("sym" and "!!")
  
  x_axis_sym <- sym(x_axis)
  y_axis_sym <- sym(y_axis)
  if (!is.null(discrete)) discrete_sym <- sym(discrete)
  if (!is.null(continuous)) continuous_sym <- sym(continuous)
  
  if (!is.null(palette) && palette == "libd_layer_colors") {
    palette <- c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", "#FFD700", "#FF7F00", "#1A1A1A", "#666666")
  } else if (!is.null(palette) && palette == "Okabe-Ito") {
    palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  } else if (!is.null(palette) && palette == "navy") {
    palette <- c("gray95", "navy")
  }
  
  df_plot <- reducedDim(spe, type)
  
  if (!is.null(discrete)) {
    df_plot <- cbind(df_plot, colData(spe)[, discrete, drop = FALSE])
  }
  if (!is.null(continuous)) {
    df_plot <- cbind(df_plot, colData(spe)[, continuous, drop = FALSE])
  }
  
  df_plot <- as.data.frame(df_plot)
  
  p <- ggplot(df_plot, aes(x = !!x_axis_sym, y = !!y_axis_sym)) + 
    geom_point(size = 0.5) + 
    xlab(paste0(type, "1")) + 
    ylab(paste0(type, "2")) + 
    ggtitle("Reduced dimensions") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
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

