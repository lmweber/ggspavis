#' plotSpots
#' 
#' Plots for spatial transcriptomics datasets.
#' 
#' Functions to generate plots of spatial coordinates for spatial
#' transcriptomics datasets, optionally including cluster labels.
#' 
#' This function generates a plot showing spatial coordinates (spots) in the
#' physical x-y coordinates of the tissue slide. Cluster labels or ground truth
#' labels can be shown with colors.
#' 
#' 
#' @param spe Input object (SpatialExperiment).
#' 
#' @param x_coord Name of column in spatialCoords slot containing x-coordinates.
#'   Default = "pxl_row_in_fullres".
#' 
#' @param y_coord Name of column in spatialCoords slot containing x-coordinates.
#'   Default = "pxl_col_in_fullres".
#' 
#' @param cluster_id Name of column in colData containing cluster IDs. To plot
#'   without cluster labels, set "cluster_id = NULL". Default = "cluster_id".
#' 
#' @param palette Color palette for cluster labels. Options are
#'   "libd_layer_colors", "Okabe-Ito", or providing a vector of hex codes for a
#'   custom palette.
#' 
#' 
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, formatting).
#' 
#' 
#' @importFrom rlang sym "!!"
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot aes geom_point coord_fixed scale_y_reverse
#'   scale_color_manual ggtitle theme_bw theme element_blank
#' 
#' @export
#' 
#' @examples
#' # to do
#' 
plotSpots <- function(spe, 
                      x_coord = "pxl_row_in_fullres", y_coord = "pxl_col_in_fullres", 
                      cluster_id = NULL, 
                      palette = "libd_layer_colors") {
  
  # note: using quasiquotation to allow custom variable names in ggplot ("sym" and "!!")
  
  x_coord <- sym(x_coord)
  y_coord <- sym(y_coord)
  if (!is.null(cluster_id)) {
    cluster_id <- sym(cluster_id)
  }
  
  if (length(palette) == 1 && palette == "libd_layer_colors") {
    palette <- c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", "#FFD700", "#FF7F00", "#1A1A1A", "#666666")
  }
  else if (length(palette) == 1 && palette == "Okabe-Ito") {
    palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  }
  
  stopifnot(all(rownames(colData(spe)) == rownames(spatialCoords(spe))))
  
  # assuming barcode IDs are in first column of both colData and spatialCoords
  df <- cbind(as.data.frame(colData(spe)), as.data.frame(spatialCoords(spe))[, -1])
  
  p <- ggplot(df, aes(x = !!x_coord, y = !!y_coord)) + 
    geom_point(size = 0.5) + 
    coord_fixed() + 
    scale_y_reverse() + 
    ggtitle("Spatial coordinates (spots)") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  if (!is.null(cluster_id)) {
    p <- p + aes(color = !!cluster_id) + 
      scale_color_manual(values = palette)
  }
  
  p
  
}

