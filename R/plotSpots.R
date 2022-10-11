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
#'   x-coordinates. Default = NULL, which selects the first column.
#' 
#' @param y_coord (character) Name of column in \code{spatialCoords} containing
#'   y-coordinates. Default = NULL, which selects the second column.
#' 
#' @param sample_id (character) Name of column in \code{colData} containing
#'   sample IDs. For datasets with multiple samples, this is used to plot
#'   multiple panels (one per sample) using facetting.
#' 
#' @param in_tissue (character) Name of column in \code{colData} identifying
#'   spots over tissue, e.g. "in_tissue" for 10x Genomics Visium data. If this
#'   argument is provided, only spots over tissue will be shown. Alternatively,
#'   set to NULL to display all spots. Default = "in_tissue".
#' 
#' @param annotate (character) Name of column in \code{colData} containing
#'   values to annotate spots with colors, e.g. cluster labels (discrete values)
#'   or total UMI counts (continuous values).
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
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot aes_string facet_wrap geom_point coord_fixed
#'   ggtitle theme_bw theme element_blank scale_y_reverse scale_color_manual
#'   scale_color_gradient
#' 
#' @export
#' 
#' @examples
#' library(STexampleData)
#' spe <- Visium_humanDLPFC()
#' plotSpots(spe, annotate = "ground_truth")
#' 
plotSpots <- function(spe, 
                      x_coord = NULL, y_coord = NULL, 
                      sample_id = "sample_id", 
                      in_tissue = "in_tissue", 
                      annotate = NULL, palette = "libd_layer_colors", 
                      y_reverse = TRUE, size = 0.3) {
  
  if (!is.null(in_tissue)) stopifnot(is.character(in_tissue))
  
  if (is.null(x_coord)) x_coord <- colnames(spatialCoords(spe))[1]
  if (is.null(y_coord)) y_coord <- colnames(spatialCoords(spe))[2]
  
  n_samples <- length(table(colData(spe)[, sample_id]))
  
  # accepts "libd_layer_colors" and "Okabe-Ito"
  palette <- .get_pal(palette)
  
  df <- cbind.data.frame(colData(spe), spatialCoords(spe))
  
  if (!is.null(in_tissue)) {
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
  
  if (n_samples > 1) {
    p <- p + facet_wrap(~ sample_id)
  }
  
  if (y_reverse) {
    p <- p + scale_y_reverse()
  }
  
  if (is.factor(df[, annotate]) | is.character(df[, annotate])) {
    p <- p + scale_color_manual(values = palette)
  }
  
  if (is.numeric(df[, annotate])) {
    p <- p + scale_color_gradient(low = palette[1], high = palette[2])
  }
  
  p
}

