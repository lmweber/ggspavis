#' plotSpots
#' 
#' Plotting functions for spatially resolved transcriptomics data.
#' 
#' Function to plot spot-based spatially resolved transcriptomics data stored in
#' a \code{SpatialExperiment} object.
#' 
#' This function generates a plot in spatial coordinates (e.g. x-y coordinates
#' on a tissue slide), along with annotation such as cluster labels or total UMI
#' counts.
#' 
#' 
#' @param spe (SpatialExperiment) Input data, assumed to be a
#'   \code{SpatialExperiment} object.
#' 
#' @param x_coord (character) Name of column in \code{spatialCoords} containing
#'   x-coordinates. Default = "x".
#' 
#' @param y_coord (character) Name of column in \code{spatialCoords} containing
#'   y-coordinates. Default = "y".
#' 
#' @param in_tissue (logical) Whether to show only spots over tissue, or all
#'   spots. Options are TRUE (show spots over tissue; requires a column labelled
#'   "in_tissue" in \code{spatialData} identifying spots over tissue, which is
#'   the standard format for 10x Genomics Visium data), FALSE (show all spots),
#'   or a character value with the name of a column in \code{spatialData}
#'   identifying the spots to show.
#' 
#' @param annotate (character) Name of column in \code{colData} containing
#'   values to annotate spots with colors, e.g. cluster labels (discrete values)
#'   or total UMI counts (continuous values). For discrete values such as
#'   cluster labels, the column in \code{colData} should be formatted as a
#'   factor.
#' 
#' @param palette (character) Color palette for annotation. Options for discrete
#'   labels (e.g. cluster labels) are "libd_layer_colors", "Okabe-Ito", or a
#'   vector of color names or hex values. For continuous values (e.g. total UMI
#'   counts), provide a vector of length 2 for the low and high range, e.g.
#'   \code{c("gray90", "navy")}. Default = \code{"libd_layer_colors"}.
#' 
#' @param y_reverse (logical) Whether to reverse y coordinates, which is often
#'   required for 10x Genomics Visium data. Default = TRUE.
#' 
#' @param size (numeric) Point size for \code{geom_point()}. Default = 0.3.
#' 
#' 
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, labels, formatting, etc).
#' 
#' 
#' @importFrom SpatialExperiment spatialData spatialCoords
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot aes_string geom_point coord_fixed ggtitle theme_bw
#'   theme element_blank scale_y_reverse scale_color_manual scale_color_gradient
#' 
#' @export
#' 
#' @examples
#' # library(ggspavis)
#' # library(STexampleData)
#' # spe <- Visium_humanDLPFC()
#' # plotSpots(spe, annotate = "ground_truth")
#' 
plotSpots <- function(spe, 
                      x_coord = "x", y_coord = "y", 
                      in_tissue = TRUE, 
                      annotate = NULL, palette = "libd_layer_colors", 
                      y_reverse = TRUE, size = 0.3) {
  
  if (length(palette) == 1 && palette == "libd_layer_colors") {
    palette <- c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", "#FFD700", "#FF7F00", "#1A1A1A", "#666666")
  } else if (length(palette) == 1 && palette == "Okabe-Ito") {
    palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  }
  
  df <- as.data.frame(cbind(colData(spe), spatialData(spe), spatialCoords(spe)))
  
  if (in_tissue) {
    df <- df[df$in_tissue == 1, ]
  } else if (is.character(in_tissue)) {
    df <- df[df[, in_tissue] == 1, ]
  }
  
  p <- ggplot(df, aes_string(x = x_coord, y = y_coord, color = annotate)) + 
    geom_point(size = size) + 
    coord_fixed() + 
    ggtitle("Spatial coordinates") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  if (y_reverse) {
    p <- p + scale_y_reverse()
  }
  
  if (is.factor(df[, annotate])) {
    p <- p + scale_color_manual(values = palette)
  }
  
  if (is.numeric(df[, annotate])) {
    p <- p + scale_color_gradient(low = palette[1], high = palette[2])
  }
  
  p
}

