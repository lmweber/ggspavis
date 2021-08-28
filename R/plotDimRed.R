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
#' 
#' # use small subset of data for this example
#' # for longer examples see our online book OSTA
#' spe <- spe[, spatialData(spe)$in_tissue == 1]
#' set.seed(100)
#' n <- 200
#' spe <- spe[, sample(seq_len(ncol(spe)), n)]
#' 
#' # calculate log-transformed normalized counts
#' library(scran)
#' set.seed(100)
#' qclus <- quickCluster(spe)
#' spe <- computeSumFactors(spe, cluster = qclus)
#' spe <- logNormCounts(spe)
#' 
#' # identify top highly variable genes (HVGs)
#' is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
#' spe <- spe[!is_mito, ]
#' dec <- modelGeneVar(spe)
#' top_hvgs <- getTopHVGs(dec, prop = 0.1)
#' 
#' # run dimensionality reduction
#' library(scater)
#' set.seed(100)
#' spe <- runPCA(spe, subset_row = top_hvgs)
#' set.seed(100)
#' spe <- runUMAP(spe, dimred = "PCA")
#' colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)
#' 
#' # generate plot
#' plotDimRed(spe, type = "UMAP", annotate = "ground_truth")
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
  
  # accepts "libd_layer_colors" and "Okabe-Ito"
  palette <- .get_pal(palette)
  
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

