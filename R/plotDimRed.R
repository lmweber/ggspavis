#' plotDimRed
#' 
#' Plotting functions for spatially resolved transcriptomics data.
#' 
#' Function to plot spot-based spatially resolved transcriptomics data stored in
#' a \code{SpatialExperiment} object.
#' 
#' This function generates a plot in reduced dimension coordinates (PCA or
#' UMAP), along with annotation such as cluster labels or total UMI counts.
#' 
#' 
#' @param spe (SpatialExperiment) Input data, assumed to be a
#'   \code{SpatialExperiment} object.
#' 
#' @param type (character) Type of reduced dimension plot. Options are "UMAP" or
#'   "PCA". Default = "UMAP".
#' 
#' @param x_axis (character) Name of column in \code{reducedDim} containing
#'   x-coordinates. Default = "UMAP1" or "PC1", depending on plot type.
#' 
#' @param y_axis (character) Name of column in \code{reducedDim} containing
#'   y-coordinates. Default = "UMAP2" or "PC2", depending on plot type.
#' 
#' @param annotate (character) Name of column in \code{colData} containing
#'   values to annotate spots with colors, e.g. cluster labels (discrete values)
#'   or total UMI counts (continuous values). For discrete values such as
#'   cluster labels, the column in \code{colData} should be formatted as a
#'   factor.
#' 
#' @param palette (character) Color palette for annotation. Options for discrete
#'   labels are "libd_layer_colors", "Okabe-Ito", or a vector of color names or
#'   hex values. For continuous values, provide a vector of length 2 for the low
#'   and high range, e.g. \code{c("gray90", "navy")}. Default =
#'   \code{"libd_layer_colors"}.
#' 
#' @param size (numeric) Point size for \code{geom_point()}. Default = 0.3.
#' 
#' 
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, labels, formatting, etc).
#' 
#' 
#' @importFrom SingleCellExperiment colData reducedDim
#' @importFrom ggplot2 ggplot aes_string geom_point xlab ylab ggtitle theme_bw
#'   theme element_blank scale_color_manual scale_color_gradient
#' 
#' @export
#' 
#' @examples
#' library(STexampleData)
#' spe <- Visium_humanDLPFC()
#' # add random data in reducedDims
#' dat <- matrix(ncol = 2, runif(ncol(spe) * 2))
#' colnames(dat) <- paste0("PC", 1:2)
#' reducedDims(spe, type = "PCA") <- list(PCA = dat)
#' plotDimRed(spe, type = "PCA")
#' 
plotDimRed <- function(spe, 
                       type = c("UMAP", "PCA"), 
                       x_axis = NULL, y_axis = NULL, 
                       annotate = NULL, palette = "libd_layer_colors", 
                       size = 0.3) {
  
  type <- match.arg(type)
  
  if (is.null(x_axis) & is.null(y_axis)) {
    if (type == "UMAP") {
      x_axis <- "UMAP1"
      y_axis <- "UMAP2"
    } else if (type == "PCA") {
      x_axis <- "PC1"
      y_axis <- "PC2"
    }
  }
  
  if (length(palette) == 1 && palette == "libd_layer_colors") {
    palette <- c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", "#FFD700", "#FF7F00", "#1A1A1A", "#666666")
  } else if (length(palette) == 1 && palette == "Okabe-Ito") {
    palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  }
  
  df <- as.data.frame(cbind(colData(spe), reducedDim(spe, type)))
  
  p <- ggplot(df, aes_string(x = x_axis, y = y_axis, color = annotate)) + 
    geom_point(size = size) + 
    xlab(paste0(type, "1")) + 
    ylab(paste0(type, "2")) + 
    ggtitle("Reduced dimensions") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  if (is.factor(df[, annotate])) {
    p <- p + scale_color_manual(values = palette)
  }
  
  if (is.numeric(df[, annotate])) {
    p <- p + scale_color_gradient(low = palette[1], high = palette[2])
  }
  
  p
}

