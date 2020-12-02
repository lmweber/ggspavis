#' plotDimRed
#'
#' Plots for spatial transcriptomics datasets.
#'
#' Functions to generate plots of spatial coordinates for spatial
#' transcriptomics datasets, optionally including cluster labels.
#'
#' This function generates a plot showing spatial coordinates (spots) in reduced
#' dimension space (either PCA or UMAP). Cluster labels or ground truth labels
#' can be shown with colors.
#'
#'
#' @param spe Input object (SpatialExperiment).
#'
#' @param type Type of dimension reduction to use. Options are "PCA" or "UMAP".
#'   Default = "UMAP".
#'
#' @param x_axis Name of column in reducedDim to use for x-axis. Default =
#'   "UMAP1".
#'
#' @param y_axis Name of column in reducedDim to use for y-axis. Default =
#'   "UMAP2".
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
#' @importFrom SingleCellExperiment colData reducedDim
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual xlab ylab
#'   ggtitle theme_bw theme element_blank
#'
#' @export
#'
#' @examples
#' # to do
#' 
plotDimRed <- function(spe, 
                       type = "UMAP", 
                       x_axis = "UMAP1", y_axis = "UMAP2", 
                       cluster_id = NULL, 
                       palette = "libd_layer_colors") {
  
  # note: using quasiquotation to allow custom variable names in ggplot ("sym" and "!!")
  
  x_axis <- sym(x_axis)
  y_axis <- sym(y_axis)
  if (!is.null(cluster_id)) {
    cluster_id <- sym(cluster_id)
  }
  
  if (length(palette) == 1 && palette == "libd_layer_colors") {
    palette <- c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", "#FFD700", "#FF7F00", "#1A1A1A", "#666666")
  }
  else if (length(palette) == 1 && palette == "Okabe-Ito") {
    palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  }
  
  stopifnot(all(rownames(colData(spe)) == rownames(reducedDim(spe, type))))
  
  df <- cbind(as.data.frame(colData(spe)), as.data.frame(reducedDim(spe, type)))
  
  p <- ggplot(df, aes(x = !!x_axis, y = !!y_axis)) + 
    geom_point(size = 0.5) + 
    xlab(paste0(type, "1")) + 
    ylab(paste0(type, "2")) + 
    ggtitle(paste0("Reduced dimensions (", type, ")")) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  if (!is.null(cluster_id)) {
    p <- p + aes(color = !!cluster_id) + 
      scale_color_manual(values = palette)
  }
  
  p
  
}

